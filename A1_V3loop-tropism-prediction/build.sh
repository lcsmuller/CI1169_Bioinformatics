#!/bin/bash -e

results_dir='results'
tmp_dir='tmp'
newdb_dir='newdb'
patients_dir='patients'

function bioinformatics_pipeline() {
    # Check if $1 is a valid .fasta file
    if ! file -b --mime-type "$1" | grep -q '^text/'; then
        echo "Error: $1 is not a valid .fasta file."
        exit 1
    fi

    # Get the correceptor name from the file path
    correceptor_name=$(basename "$(dirname "$1")")
    dir=$results_dir/$correceptor_name

    mkdir -p $dir

    # Classify sequences
    hmmsearch --max -o "${dir}/results_r5.txt" --tblout "${dir}/sequence_hits_r5.txt" $tmp_dir/r5.hmm "$1"
    hmmsearch --max -o "${dir}/results_x4.txt" --tblout "${dir}/sequence_hits_x4.txt" $tmp_dir/x4.hmm "$1"

    # Generate plots for visualization
    python3 generate_plots.py "$dir"
}

mkdir -p $results_dir
mkdir -p $tmp_dir

# Concatenate sequences
awk '{print}' $newdb_dir/ccr5.fasta $newdb_dir/r5x4.fasta > $tmp_dir/all_r5.fasta
awk '{print}' $newdb_dir/cxcr4.fasta $newdb_dir/r5x4.fasta > $tmp_dir/all_x4.fasta

# Align sequences
clustalo -i $tmp_dir/all_r5.fasta -o $tmp_dir/aligned_r5.fasta --force
clustalo -i $tmp_dir/all_x4.fasta -o $tmp_dir/aligned_x4.fasta --force

# Build HMM profiles
hmmbuild $tmp_dir/r5.hmm $tmp_dir/aligned_r5.fasta
hmmbuild $tmp_dir/x4.hmm $tmp_dir/aligned_x4.fasta

for file in $patients_dir/*/*.fasta; do
    bioinformatics_pipeline "$file"
done
