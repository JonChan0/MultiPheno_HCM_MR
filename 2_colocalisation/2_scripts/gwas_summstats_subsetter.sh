#!/bin/bash

# --- 📂 1. DEFINE YOUR FILE PATHS HERE ---
# Add the full paths to your .tsv.gz files inside the parentheses.
FILES=(
    "../../1_MR/2_scripts/henry25_hf_gwas_ssf/final/henry25_hf_gwas_ssf.h.tsv.gz"
    "../../1_MR/2_scripts/hcmr_NTproBNP_gwas_ssf/final/hcmr_NTproBNP_gwas_ssf.h.tsv.gz"
    "../../1_MR/2_scripts/hcmr_TnT_gwas_ssf/final/hcmr_TnT_gwas_ssf.h.tsv.gz"
    "../../1_MR/2_scripts/keaton24_dbp_gwas_ssf/final//keaton24_dbp_gwas_ssf.h.tsv.gz"
)

# --- 📁 2. DEFINE YOUR OUTPUT DIRECTORY ---
OUTPUT_DIR="../3_output/gwas_summstat_subsets/"

# --- ⚙️ 3. SCRIPT LOGIC (No user changes needed below) ---

# Check for correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "❌ Incorrect number of arguments."
    echo "Usage: $0 <chromosome> <base_pair_location>"
    echo "Example: $0 10 73650119"
    exit 1
fi

# Assign arguments to variables for clarity
INPUT_CHR=$1
INPUT_BP=$2
WINDOW=500000 # The range on either side of the input base pair

# Calculate the start and end positions for the range using bash arithmetic
START_POS=$((INPUT_BP - WINDOW))
END_POS=$((INPUT_BP + WINDOW))

# Create the output directory if it doesn't already exist
mkdir -p "$OUTPUT_DIR"

echo "Filtering for chromosome '$INPUT_CHR' in range $START_POS - $END_POS"
echo "Output will be saved to the '$OUTPUT_DIR' directory."
echo "---"

# Loop through each file path specified in the array
for file in "${FILES[@]}"
do
  # First, check if the input file actually exists
  if [ ! -f "$file" ]; then
      echo "⚠️ Warning: File not found at '$file'. Skipping."
      continue
  fi

  echo "Processing '$file'..."

  # --- Generate the new output filename ---
  # Get the base name of the input file (e.g., "file1.tsv.gz" -> "file1")
  base_name=$(basename "$file" .tsv.gz)
  
  # Construct the new filename with the desired suffix
  output_filename="${OUTPUT_DIR}/${base_name}_${INPUT_CHR}_${INPUT_BP}.tsv"

  # --- Process the file ---
  # 1. Write the header from the original file to the new output file (overwriting it)
  zcat < "$file" | head -n 1 > "$output_filename"

  # 2. Filter the data using awk and append (>>) the results to the new output file
  #
  # !! IMPORTANT !!
  # -> Replace '$1' with the column number for 'chromosome'.
  # -> Replace '$2' with the column number for 'base_pair_location'.
  #
  zcat < "$file" | awk -F'\t' -v chr="$INPUT_CHR" -v start="$START_POS" -v end="$END_POS" \
    'NR > 1 && $1 == chr && $2 >= start && $2 <= end' >> "$output_filename"

done

echo "---"
echo "✅ All files processed successfully."
