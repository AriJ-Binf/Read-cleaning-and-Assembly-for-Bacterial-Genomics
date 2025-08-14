#!/usr/bin/env bash
set -euo pipefail

# ----------------- Defaults -----------------
BASE_DIR="$(pwd)"
PHRED="33"                 # 33 or 64 (Trimmomatic only)
SLIDINGWINDOW="5:30"       # Trimmomatic: W:Q
LEADING="5"
TRAILING="5"
MINLEN="50"
THREADS="8"
TRIMMER="trimmomatic"      # trimmomatic | fastp

ADAPTERS=""                # Trimmomatic ILLUMINACLIP adapters (optional)
DETECT_ADAPTERS="false"    # fastp adapter autodetect (PE)

usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS] <raw_data_folder>

General:
  -b, --base-dir DIR            Base output directory (default: $BASE_DIR)
  -t, --threads N               Threads (default: $THREADS)
  -r, --trimmer {trimmomatic|fastp}   (default: $TRIMMER)

Trimmomatic (used only if --trimmer trimmomatic):
  -p, --phred <33|64>           (default: $PHRED)
  -w, --slidingwindow W:Q       (default: $SLIDINGWINDOW)
  -L, --leading N               (default: $LEADING)
  -T, --trailing N              (default: $TRAILING)
  -m, --minlen N                (default: $MINLEN)
  -A, --adapters FILE           ILLUMINACLIP adapters fasta (optional)

fastp (used only if --trimmer fastp):
  -D, --detect-adapters         Enable adapter auto-detection for PE

Notes:
  * If --trimmer fastp is selected, FastQC is skipped (fastp makes HTML/JSON QC).
  * Script processes ALL pairs matching *_1.fastq.gz in the input folder.
EOF
}

# -------- Long -> short mapping --------
args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-dir)         args+=(-b "${2:-}"); shift 2 ;;
    --threads)          args+=(-t "${2:-}"); shift 2 ;;
    --trimmer)          args+=(-r "${2:-}"); shift 2 ;;
    --phred)            args+=(-p "${2:-}"); shift 2 ;;
    --slidingwindow)    args+=(-w "${2:-}"); shift 2 ;;
    --leading)          args+=(-L "${2:-}"); shift 2 ;;
    --trailing)         args+=(-T "${2:-}"); shift 2 ;;
    --minlen)           args+=(-m "${2:-}"); shift 2 ;;
    --adapters)         args+=(-A "${2:-}"); shift 2 ;;
    --detect-adapters)  args+=(-D); shift ;;
    --help)             args+=(-h); shift ;;
    --)                 shift; break ;;
    *)                  args+=("$1"); shift ;;
  esac
done
set -- "${args[@]}" "$@"

# --------------------------- getopts parsing ---------------------------
RAW_DATA_FOLDER=""
while getopts ":b:t:r:p:w:L:T:m:A:Dh" opt; do
  case "$opt" in
    b) BASE_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    r) TRIMMER="$OPTARG" ;;
    p) PHRED="$OPTARG" ;;
    w) SLIDINGWINDOW="$OPTARG" ;;
    L) LEADING="$OPTARG" ;;
    T) TRAILING="$OPTARG" ;;
    m) MINLEN="$OPTARG" ;;
    A) ADAPTERS="$OPTARG" ;;
    D) DETECT_ADAPTERS="true" ;;
    h) usage; exit 0 ;;
    :)
      echo "Error: -$OPTARG requires an argument." >&2; usage; exit 1 ;;
    \?)
      echo "Error: unknown option: -$OPTARG" >&2; usage; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

# ----------------------------- Validation -----------------------------
if [[ $# -lt 1 ]]; then
  echo "Error: <raw_data_folder> is required." >&2; usage; exit 1
fi
RAW_DATA_FOLDER="$1"
[[ -d "$RAW_DATA_FOLDER" ]] || { echo "Error: folder not found: $RAW_DATA_FOLDER" >&2; exit 1; }

TRIMMER=$(echo "$TRIMMER" | tr '[:upper:]' '[:lower:]')

if [[ "$TRIMMER" != "trimmomatic" && "$TRIMMER" != "fastp" ]]; then
  echo "Error: --trimmer must be trimmomatic or fastp" >&2; exit 1
fi
if [[ "$TRIMMER" == "trimmomatic" ]]; then
  [[ "$PHRED" == "33" || "$PHRED" == "64" ]] || { echo "Error: --phred must be 33 or 64" >&2; exit 1; }
fi
if [[ -n "$ADAPTERS" && ! -f "$ADAPTERS" ]]; then
  echo "Error: adapters file not found: $ADAPTERS" >&2; exit 1
fi

# --------------------------- Pipeline ---------------------------
echo ">>> Settings"
echo "  Base dir       : $BASE_DIR"
echo "  Input folder   : $RAW_DATA_FOLDER"
echo "  Threads        : $THREADS"
echo "  Trimmer        : $TRIMMER"
if [[ "$TRIMMER" == "trimmomatic" ]]; then
  echo "  Trimmomatic    : -phred$PHRED SLIDINGWINDOW:$SLIDINGWINDOW LEADING:$LEADING TRAILING:$TRAILING MINLEN:$MINLEN"
  [[ -n "$ADAPTERS" ]] && echo "  Adapters       : $ADAPTERS"
else
  echo "  fastp          : detect_adapters=$DETECT_ADAPTERS (FastQC skipped)"
fi
echo

mkdir -p "$BASE_DIR"
shopt -s nullglob

for R1 in "$RAW_DATA_FOLDER"/*_1.fastq.gz; do
  SRA_ACC="$(basename "$R1" "_1.fastq.gz")"
  R2="$RAW_DATA_FOLDER/${SRA_ACC}_2.fastq.gz"
  [[ -f "$R2" ]] || { echo ">> Skipping $SRA_ACC: mate not found ($R2)" >&2; continue; }

  echo ">>> Processing $SRA_ACC ..."
  SRA_DIR="$BASE_DIR/$SRA_ACC"
  mkdir -p "$SRA_DIR/QC" "$SRA_DIR/trim" "$SRA_DIR/asm"

  OUT_R1_P="$SRA_DIR/trim/r1.paired.fq.gz"
  OUT_R2_P="$SRA_DIR/trim/r2.paired.fq.gz"
  OUT_R1_U="$SRA_DIR/trim/r1_unpaired.fq.gz"
  OUT_R2_U="$SRA_DIR/trim/r2_unpaired.fq.gz"
  SINGLETONS="$SRA_DIR/trim/singletons.fq.gz"

  # --------- QC (skip if fastp is the trimmer) ---------
  if [[ "$TRIMMER" == "trimmomatic" ]]; then
    fastqc --threads "$THREADS" --outdir "$SRA_DIR/QC" "$R1" "$R2"
  fi

  # ------------------------- Trimming ------------------------
  if [[ "$TRIMMER" == "trimmomatic" ]]; then
    PHRED_FLAG="-phred${PHRED}"
    TRIM_FLAGS=( "SLIDINGWINDOW:${SLIDINGWINDOW}" "LEADING:${LEADING}" "TRAILING:${TRAILING}" "MINLEN:${MINLEN}" )
    if [[ -n "$ADAPTERS" ]]; then
      TRIM_FLAGS=( "ILLUMINACLIP:${ADAPTERS}:2:30:10" "${TRIM_FLAGS[@]}" )
    fi

    trimmomatic PE -threads "$THREADS" "$PHRED_FLAG" \
      "$R1" "$R2" \
      "$OUT_R1_P" "$OUT_R1_U" \
      "$OUT_R2_P" "$OUT_R2_U" \
      "${TRIM_FLAGS[@]}" \
      1>"$SRA_DIR/trim/trimmo.stdout.log" \
      2>"$SRA_DIR/trim/trimmo.stderr.log"

    if [[ -s "$OUT_R1_U" || -s "$OUT_R2_U" ]]; then
      cat "$OUT_R1_U" "$OUT_R2_U" > "$SINGLETONS" || true
    else
      : > "$SINGLETONS"
    fi
    rm -f "$OUT_R1_U" "$OUT_R2_U"

  else
    FASTP_ARGS=( -w "$THREADS" )
    [[ "$DETECT_ADAPTERS" == "true" ]] && FASTP_ARGS+=( --detect_adapter_for_pe )
    fastp \
      -i "$R1" -I "$R2" \
      -o "$OUT_R1_P" -O "$OUT_R2_P" \
      --unpaired1 "$OUT_R1_U" --unpaired2 "$OUT_R2_U" \
      --json "$SRA_DIR/trim/fastp.json" \
      --html "$SRA_DIR/trim/fastp.html" \
      "${FASTP_ARGS[@]}" \
      1>"$SRA_DIR/trim/fastp.stdout.log" \
      2>"$SRA_DIR/trim/fastp.stderr.log"

    if [[ -s "$OUT_R1_U" || -s "$OUT_R2_U" ]]; then
      cat "$OUT_R1_U" "$OUT_R2_U" > "$SINGLETONS" || true
    else
      : > "$SINGLETONS"
    fi
    rm -f "$OUT_R1_U" "$OUT_R2_U"
  fi

  # ------------------------- Assembly (SPAdes only) ------------------------
  SPADES_ARGS=( -1 "$OUT_R1_P" -2 "$OUT_R2_P" -o "$SRA_DIR/asm/spades" --only-assembler -t "$THREADS" )
  [[ -s "$SINGLETONS" ]] && SPADES_ARGS+=( -s "$SINGLETONS" )
  spades.py "${SPADES_ARGS[@]}" \
    1>"$SRA_DIR/asm/spades.stdout.txt" \
    2>"$SRA_DIR/asm/spades.stderr.txt"

  ASM_CONTIGS="$SRA_DIR/asm/spades/contigs.fasta"

  # ----------------------- Post-assembly ----------------------
  if [[ -f "$ASM_CONTIGS" ]]; then
    echo -n "Contigs: "
    grep -c '>' "$ASM_CONTIGS" || true
  else
    echo "Warning: contigs.fasta not found for $SRA_ACC" >&2
  fi

  echo
done



