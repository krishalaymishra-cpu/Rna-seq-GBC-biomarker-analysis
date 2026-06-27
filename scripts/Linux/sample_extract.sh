#!/bin/bash
word_dir=$(pwd)
Accession=SRP223526
batch_size=20     #number of files per batch
parallel_jobs=4  #change it to your own system capacity

#---Prepare directories--------
mkdir -p "$work_dir"
cd "$work_dir" || exit 1
echo "Working directory: $work_dir..."

# --- Query ENA for FASTQ links ---
echo "Fetching FASTQ links from ENA..."
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${Accession}&result=read_run&fields=run_accession,fastq_ftp&format=tsv" \
| awk 'NR>1 {print $2}' | tr ';' '\n' | sed 's|^|ftp://|' > fastq_links.txt

#----Count total-------
total=$(wc -l < fastq_links.txt)
echo "Found $total FASTQ files to process."

downloaded=0
skipped=0
batch_num=1

# --- Split into batches ---
split -l $batch_size fastq_links.txt batch_

for batch_file in batch_*; do
    echo "=== Processing batch $batch_num ($(wc -l < $batch_file) files) ==="

    # Use GNU parallel to download multiple files at once
    cat "$batch_file" | parallel -j $parallel_jobs '
        file=$(basename {});
        if [[ -f "$file" ]]; then
            echo "Skipping existing file: $file"
        else
            echo "Downloading: $file"
            wget -c {}
        fi
    '

    echo "=== Finished batch $batch_num ==="
    ((batch_num++))
done

echo "All tasks finished!"
