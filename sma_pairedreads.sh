# Script for metagenomic analysis of the gut microbiome from non-human primates
# This script is under development
# Author: Marcus V. C. Viana
# Date: 28/08/2026


############################################################
## SUMMARY OF SHOTGUN METAGENOME ANAYLSIS WORKFLOW (FROM SHORT-READS)
############################################################

# 1) Reads files and renaming
#     * Reads stored as local files
#     * Reads from NCBI SRA (sra-tools)
# 2) Raw reads quality assessment
#     * FastQC
#     * FastQC -> MultiQC
# 3) Raw reads trimming
#     * Fastp
# 4) Trimmed reads quality assessment
#     * FastQC
#     * FastQC -> MultiQC
# 5) Host decontamination (optional)
#     * NCBI Datasets
#     * Bwa-mem2 index
#     * Bwa-mem2 mem
#     * Bwa-mem2 reads
#     * Bwa-mem2 -> FastQC
#     * Bwa-mem2 -> FastQC -> MultiQC
# 6) Taxonomic abundance profile
#     * Kraken
#     * Kraken -> Bracken
#     * Kraken -> Bracken -> Comparison
#     * Kraken -> Bracken -> Krona
#     * MetaPhlAn
#     * MetaPhlAn -> Comparison
# 7) Metagenome assembly
#     * MEGAHIT
#     * MEGAHIT -> QUAST
# 8) Functional abundance profile and prophages
#     * Pybarrnap
#     * Aragorn
#     * Pyrodigal
#     * Pyrodigal -> AMRFinderPlus
#     * Pyrodigal -> dbCAN
#     * Pyrodigal -> eggNOG-mapper
#     * Pyrodigal -> VFDB (BLASTP)
# 9) Functional abundance profile of gene catalog
#     * MMseqs2 input (SeqKit)
#     * MMseqs2 easy-linclust
#     * MMseqs2 -> AMRFinderPlus
#     * MMseqs2 -> dbCAN
#     * MMseqs2 -> eggNOG-mapper
#     * MMseqs2 -> Pyrodigal -> VFDB (BLASTP)
#     * MMseqs2 -> Salmon index
#     * MMseqs2 -> Salmon quant
# 10) Binning - Single/Multi-sample - Input files
#     * Seqkit
#     * Seqkit -> SemiBin concatenate_fasta
#     * Seqkit -> SemiBin concatenate_fasta -> Minimap2 index
#     * Seqkit -> SemiBin concatenate_fasta -> Minimap2
# 11) Binning - Single/Multi-sample (Self-supervised mode)
#     * SemiBin
# 12) Bin quality control and taxonomy
#     * CheckM2
#     * GUNC
#     * GTDB-Tk
#     * QUAST
# 13) Bin functional abundance profile
#     * Aragorn
#     * Pybarrnap
#     * Pyrodigal
#     * Pyrodigal -> AMRFinderPlus
#     * Pyrodigal -> dbCAN
#     * Pyrodigal -> eggNOG-mapper
#     * Pyrodigal -> VFDB (BLASTP)
# 14) Bin mobile genetic elements
#     * MOB-suite 
#     * VIBRANT


############################################################
# 1) Reads files and renaming
############################################################

############################################################
## 1.1) Reads stored as local files

# 1. Standardize the paired-end file names of each sample to **samplename_1.fq.gz** and **samplename_2.fq.gz**.
# 2. In the working directory, create the directory **1_reads** and place the read files **inside it**.

# Create an output directory
mkdir -p 1_reads
# Put reads there

############################################################
## 1.2) Reads from NCBI SRA (SRA Tools)

# 1. Create a **tab-separated file** named **"1_reads_accessions.tsv"**.
# 2. This file **must contain** the NCBI SRA **accession number** in the first column and the **sample name** in the second column. Other columns will be ignored.
# 3. **Do not use** special characters in the sample names.
# 4. Place the **"1_reads_accessions.tsv"** file in the working directory.

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="1) SRA Tools"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS


# Delete previous file of not used reads
rm -f 1_reads_single-end.tsv

# Verify the presence of the file 1_reads_accessions.tsv with a list of accessions
if [ -f 1_reads_accessions.tsv ]; then
    echo "1) The file 1_reads_accessions.tsv was found. The sequencing reads will be downloaded."

    # Create output directory 
    mkdir -p 1_reads

    # SRA Tools
    # Activate Conda environment
    conda activate sra-tools

    # Loop through file lines
    tr -d '\r' < 1_reads_accessions.tsv |  awk '1'| \
    while IFS=$'\t' read -r accession sample others; do

        # Check for valid, complete paired-end or single-end output
        paired_ok=0
        single_ok=0
        if [ -f "1_reads/${sample}_1.fq.gz" ] && [ -f "1_reads/${sample}_2.fq.gz" ] \
            && gzip -t "1_reads/${sample}_1.fq.gz" 2>/dev/null && gzip -t "1_reads/${sample}_2.fq.gz" 2>/dev/null; then
            paired_ok=1
        elif [ -f "1_reads/${sample}.fq.gz" ] && gzip -t "1_reads/${sample}.fq.gz" 2>/dev/null; then
            single_ok=1
        fi

        if [ "$paired_ok" -eq 1 ]; then
            echo "1) Sample $sample paired files found and valid. Skipping download."
        elif [ "$single_ok" -eq 1 ]; then
            echo "1) Sample $sample single-end file found and valid. Skipping download."
        else
            echo "1) Downloading sample: $sample (accession: $accession)"

            # Remove any incomplete/corrupted files and leftover intermediates from a previous interrupted run
            rm -f "1_reads/${sample}_1.fq.gz" "1_reads/${sample}_2.fq.gz" "1_reads/${sample}.fq.gz"
            rm -rf "1_reads/${accession}"
            rm -f "1_reads/${accession}"*.fastq "1_reads/${accession}"*.fastq.gz

            # Run prefetch
            prefetch -p -O 1_reads "${accession}"

            # Run fasterq-dump
            cd 1_reads
            fasterq-dump \
            --threads $(nproc --ignore=1) \
            -p \
            --split-files "${accession}" \
            -O .
            # Delete temporary directories
            rm -rf "${accession}"
            # Compress files
            echo "Compressing fastq files."
            pigz -p $(nproc --ignore=1) ${accession}*.fastq

            # Check if the reads are paired-end
            if [ -f "${accession}_1.fastq.gz" ] && [ -f "${accession}_2.fastq.gz" ]  \
                && gzip -t "${accession}_1.fastq.gz" 2>/dev/null && gzip -t "${accession}_2.fastq.gz" 2>/dev/null; then
                echo "1) Sample ${sample} has paired-end reads."
                # Rename files
                echo "1) Renaming files."
                mv "${accession}_1.fastq.gz" "${sample}_1.fq.gz"
                mv "${accession}_2.fastq.gz" "${sample}_2.fq.gz"
                # Go back to main directory
                cd ..
            # Check if the reads are not paired-end
            elif [ -f "${accession}.fastq.gz" ] && gzip -t "${accession}.fastq.gz" 2>/dev/null; then
                echo "1) Sample ${sample} has single-end reads."
                # Renaming file
                echo "1) Renaming file."
                mv "${accession}.fastq.gz" "${sample}.fq.gz"
                # Warning about the requirement of paired-reads
                echo "1) Warning: this script only use paired-end reads. This file will not be used." 
                # Create list of not used read files
                echo -e "${accession}\t${sample}.fq.gz" >> ../1_reads_single-end.tsv
                # Go back to main directory
                cd ..
            else
                echo "1) Warning: no valid output produced for sample ${sample} (accession: ${accession}). It will be retried on the next run."
                rm -f "${accession}"*.fastq.gz
                cd ..
            fi
        fi
    done
    echo "1) Download process complete. Deactivating the environment."
    # Deactivate Conda environment
    conda deactivate
else
    echo "1) The file 1_reads_accessions.tsv was not found. Proceeding using local files."
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
# 2) Raw reads quality assessment
############################################################

############################################################
## 2.1) FastQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="2) FastQC"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "2_fastqc.tar.gz" ] && [ -f "2_fastqc.tar.gz.md5" ] && md5sum -c "2_fastqc.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (2_fastqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo,
    # since FastQC runs as a single call over all files and has no reliable per-file resume here.
    if [ -d "2_fastqc" ]; then
        echo "${workflow_step}: Found incomplete 2_fastqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "2_fastqc"
    fi
    rm -f "2_fastqc.tar.gz" "2_fastqc.tar.gz.md5"

    # Create an output directory
    mkdir -p 2_fastqc

    # Activate Conda environment
    conda activate fastqc
    # Run main software
    fastqc -t $(nproc --ignore=1) 1_reads/*.gz -o 2_fastqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="2_fastqc.tar.gz"
    itens_to_compress=(2_fastqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt


############################################################
## 2.2) FastQC -> MultiQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="2) MultiQC"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "2_fastqc_multiqc.tar.gz" ] && [ -f "2_fastqc_multiqc.tar.gz.md5" ] && md5sum -c "2_fastqc_multiqc.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (2_fastqc_multiqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
    # Clean up 2_fastqc too, since a completed MultiQC step means it's no longer needed
    rm -rf "2_fastqc" "2_fastqc.tar.gz.md5.bak" 2>/dev/null
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo.
    if [ -d "2_fastqc_multiqc" ]; then
        echo "${workflow_step}: Found incomplete 2_fastqc_multiqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "2_fastqc_multiqc"
    fi
    rm -f "2_fastqc_multiqc.tar.gz" "2_fastqc_multiqc.tar.gz.md5"

    # Activate Conda environment
    conda activate multiqc
    # Run main software
    multiqc 2_fastqc/*_fastqc.zip -o 2_fastqc_multiqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="2_fastqc_multiqc.tar.gz"
    itens_to_compress=(2_fastqc_multiqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

    # Delete the output directory
    rm -r 2_fastqc 2_fastqc_multiqc
fi

# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt


############################################################
# 3) Raw reads trimming
############################################################

############################################################
## 3.1) Fastp

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="3) Fastp"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the loop running time
start_time=$SECONDS

# Create an output directory
mkdir -p 3_fastp    

# Calculate sample size
i=1
r1_files=(1_reads/*_1.fq.gz)
sample_count=${#r1_files[@]}

# Activate Conda environment
conda activate fastp
# Loop through a list of sample files
for r1 in "${r1_files[@]}"; do

    # Obtain r2 path
    r2="${r1/_1.fq.gz/_2.fq.gz}"
    # Extract r1 file name
    r1filename=${r1##*/}
    # Extract sample name
    sample=${r1filename%%_*}

    # Skip sample if trimmed output already exists and is valid
    out1="3_fastp/${sample}_trimmed_1.fq.gz"
    out2="3_fastp/${sample}_trimmed_2.fq.gz"
    if [ -f "$out1" ] && [ -f "$out2" ] && gzip -t "$out1" 2>/dev/null && gzip -t "$out2" 2>/dev/null; then
        echo "${workflow_step} output files already exist and are valid for sample: $sample. Skipping sample."
        i=$((i + 1))
        continue
    elif [ -f "$out1" ] || [ -f "$out2" ] || [ -f "3_fastp/${sample}_trimmed_fastp.html" ] || [ -f "3_fastp/${sample}_trimmed_fastp.json" ]; then
        echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
        rm -f "$out1" "$out2" "3_fastp/${sample}_trimmed_fastp.html" "3_fastp/${sample}_trimmed_fastp.json"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    echo "${workflow_step} - Read 1 file: ${r1}"
    echo "${workflow_step} - Read 2 file: ${r2}"
    # Start counting the loop running time
    loop_start_time=$SECONDS

    # Run main software
    fastp \
        --thread $(nproc --ignore=1) \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --trim_poly_x \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --length_required 50 \
        --overrepresentation_analysis \
        --in1 "$r1" \
        --in2 "$r2" \
        --out1 "$out1" \
        --out2 "$out2" \
        --html "3_fastp/${sample}_trimmed_fastp.html" \
        --json "3_fastp/${sample}_trimmed_fastp.json"
    
    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="3_fastp.tar.gz"
itens_to_compress=(3_fastp/*.json 3_fastp/*.html)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Generate checksum files for the reads
(cd 3_fastp && for file in *.gz; do
    [ -f "${file}.md5" ] && md5sum -c "${file}.md5" >/dev/null 2>&1 && continue
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done) | tee -a 0_workflow_progress.txt
# Check file integrity
echo "${workflow_step}: Checking file integrity of reads" | tee -a 0_workflow_progress.txt
(cd 3_fastp && md5sum -c *.md5) | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt


############################################################
# 4) Trimmed reads quality assessment
############################################################

############################################################
## 4.1) FastQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="4) FastQC"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "4_fastqc.tar.gz" ] && [ -f "4_fastqc.tar.gz.md5" ] && md5sum -c "4_fastqc.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (4_fastqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo,
    # since FastQC runs as a single call over all files and has no reliable per-file resume here.
    if [ -d "4_fastqc" ]; then
        echo "${workflow_step}: Found incomplete 4_fastqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "4_fastqc"
    fi
    rm -f "4_fastqc.tar.gz" "4_fastqc.tar.gz.md5"

    # Create an output directory
    mkdir -p 4_fastqc

    # Activate Conda environment
    conda activate fastqc
    # Run main software
    fastqc -t $(nproc --ignore=1) 3_fastp/*.gz -o 4_fastqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="4_fastqc.tar.gz"
    itens_to_compress=(4_fastqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt


############################################################
## 4.2) FastQC -> MultiQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="4) MultiQC"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "4_fastqc_multiqc.tar.gz" ] && [ -f "4_fastqc_multiqc.tar.gz.md5" ] && md5sum -c "4_fastqc_multiqc.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (4_fastqc_multiqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo.
    if [ -d "4_fastqc_multiqc" ]; then
        echo "${workflow_step}: Found incomplete 4_fastqc_multiqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "4_fastqc_multiqc"
    fi
    rm -f "4_fastqc_multiqc.tar.gz" "4_fastqc_multiqc.tar.gz.md5"

    # Activate Conda environment
    conda activate multiqc
    # Run main software
    multiqc 4_fastqc/*_fastqc.zip -o 4_fastqc_multiqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="4_fastqc_multiqc.tar.gz"
    itens_to_compress=(4_fastqc_multiqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

    # Delete the output directory
    rm -r 4_fastqc 4_fastqc_multiqc
fi

# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt

############################################################
# 5) Host decontamination (Skippable)
############################################################

############################################################
## 5.1) NCBI Datasets (Download reference genomes)

# Create the tab-separated text file named "5_ref_ids.tsv"
# Column 1: The GenBank genome assembly ID of the reference genome. It will be used to download the genome.
# Column 2: The species name of the reference genome. It will name the Bwa-mem2 index in the next step. Do not use spaces or special characters.
# Place the file in the working directory.

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="5) NCBI Datasets"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create output directory
mkdir -p 5_ref_genomes

# Check whether every reference listed in 5_ref_ids.tsv already has a valid final file
# (either the raw .fasta from this step, or the .fasta.gz produced later by 5.2)
all_present=1
while IFS=$'\t' read -r ref_accession ref_name others; do
    [ -z "$ref_accession" ] && continue
    if [ -f "5_ref_genomes/${ref_accession}.fasta" ] && [ -s "5_ref_genomes/${ref_accession}.fasta" ]; then
        continue
    elif [ -f "5_ref_genomes/${ref_accession}.fasta.gz" ] && gzip -t "5_ref_genomes/${ref_accession}.fasta.gz" 2>/dev/null; then
        continue
    else
        all_present=0
        break
    fi
done < <(tr -d '\r' < 5_ref_ids.tsv | awk '1')

if [ "$all_present" -eq 1 ]; then
    echo "${workflow_step} all reference genomes already present and valid. Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left partial genomes or leftover temp directories;
    # this download is not reliably resumable per-accession, so wipe and redo entirely.
    echo "${workflow_step}: Missing or incomplete reference genome(s) found. Removing partial files and redownloading everything." | tee -a 0_workflow_progress.txt
    rm -rf 5_ref_genomes genome_data.zip ncbi_dataset README.md md5sum.txt
    mkdir -p 5_ref_genomes

    # Activate Conda environment
    conda activate datasets
    # Download compressed dehydrated directory
    datasets download genome accession \
        --assembly-version latest \
        --include genome \
        --dehydrated \
        --inputfile 5_ref_ids.tsv \
        --filename genome_data.zip
    # Unzip genome_data.tar.gz
    unzip genome_data.zip
    # Rehydrate directory
    datasets rehydrate --directory .
    # Deactivate Conda environment
    conda deactivate
    # Move assembly files to 5_ref_genomes
    mv ncbi_dataset/data/*/*.fna 5_ref_genomes
    # Remove sufix from file names and change file extension
    cd 5_ref_genomes
    rename 's/(.*?_.*?)_.*/$1.fasta/' *.fna
    cd ..
    # Delete temporary files and directory
    rm -rf genome_data.zip ncbi_dataset README.md md5sum.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 5.2) Bwa-mem2 index

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="5) Bwa-mem2 index"
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip the whole step if it already completed successfully (aggregate checksum present and valid)
if [ -f "5_bwa_index.md5" ] && md5sum -c "5_bwa_index.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (5_bwa_index.md5 verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create output directory
    mkdir -p 5_bwa_index

    # Activate Conda environment
    conda activate bwa-mem2
    # Loop through a list of file lines (Species as index name)
    tr -d '\r' < 5_ref_ids.tsv | awk '1'| \
    while IFS=$'\t' read -r ref_accession ref_name others; do

        # Skip blank/malformed lines
        [ -z "$ref_accession" ] && continue

        # Skip reference if its index files already exist and look complete
        idx_prefix="5_bwa_index/${ref_name}/${ref_name}"
        if [ -f "${idx_prefix}.0123" ] && [ -f "${idx_prefix}.bwt.2bit.64" ] && [ -f "${idx_prefix}.ann" ] \
            && [ -f "${idx_prefix}.amb" ] && [ -f "${idx_prefix}.pac" ]; then
            echo "${workflow_step} index already exists and looks complete for reference: ${ref_name}. Skipping."
            continue
        elif [ -d "5_bwa_index/${ref_name}" ]; then
            echo "${workflow_step} found incomplete index for reference: ${ref_name}. Removing partial files and reprocessing."
            rm -rf "5_bwa_index/${ref_name}"
        fi

        # Inform current sample
        echo "▶  ${workflow_step} - ${ref_name}"
        # Start counting the running time
        loop_start_time=$SECONDS

        # Create output directory
        mkdir -p 5_bwa_index/${ref_name}

        # Run main software
        bwa-mem2 index \
            -p "5_bwa_index/${ref_name}/${ref_name}" \
            "5_ref_genomes/${ref_accession}.fasta"

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ $workflow_step - $ref_name - ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    done
    # Deactivate Conda environment
    conda deactivate

    # Generate a checksum file for the indexes
    find 5_bwa_index -type f -exec md5sum {} + > 5_bwa_index.md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c 5_bwa_index.md5 | tee -a 0_workflow_progress.txt

    # Compress original reference genome files
    echo "Compressing reference genomes files" | tee -a 0_workflow_progress.txt
    for file in 5_ref_genomes/*.fasta; do
        [ -f "$file" ] || continue
        pigz -p $(nproc --ignore=1) "${file}"
    done
    # Generate checksum files for the reads
    (cd 5_ref_genomes && for file in *.fasta.gz; do
        [ -f "${file}.md5" ] && md5sum -c "${file}.md5" >/dev/null 2>&1 && continue
        echo "Processing checksum of file: ${file}"
        md5sum ${file} > ${file}.md5
    done) | tee -a 0_workflow_progress.txt
    # Check file integrity
    echo "${workflow_step}: Checking file integrity of reads" | tee -a 0_workflow_progress.txt
    (cd 5_ref_genomes && md5sum -c *.md5) | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 5.3) Bwa-mem2 mem

# Avoid literal glob pattern
shopt -s nullglob

# Create the tab-separated text file named "5_metagenomes.tsv":
# Column 1: The sample name.
# Column 2: The GenBank genome assembly ID of the reference genome. It will be used to register the specific version of the reference genome assembly.
# Column 3: The species name of the reference genome. The same name you used to create the Bwa-mem2 index.
# Column 4: The species name of the sample host (or isolation sorce). It will be used in the binning step.
# Place the file in the working directory.

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="5) Bwa-mem2 mem"
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create output directory
mkdir -p 5_bwa_reads

# Calculate sample size
i=1
sample_count=$(awk 'END {print NR}' 5_metagenomes.tsv)

# Activate conda environment
conda activate bwa-mem2
# Loop through a list of file lines
tr -d '\r' < 5_metagenomes.tsv |  awk '1' | \
while IFS=$'\t' read -r sample ref_accession ref_name isolation_source others; do

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} from ${ref_name} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Check if the reference index exists
    index_prefix="5_bwa_index/${ref_name}/${ref_name}"
    if [ ! -f "${index_prefix}.bwt.2bit.64" ]; then
        echo "5) Bwa-mem2 index not found for reference: ${ref_name} (${index_prefix}). Skipping sample: $sample" | tee -a 0_workflow_progress.txt
        i=$((i + 1))
        continue
    fi

    # Skip sample if output files already exist and are valid, complete gzip files
    r1_output_file="5_bwa_reads/${sample}_trimmed_nohost_1.fq.gz"
    r2_output_file="5_bwa_reads/${sample}_trimmed_nohost_2.fq.gz"
    if [ -f "$r1_output_file" ] && [ -f "$r2_output_file" ] \
        && gzip -t "$r1_output_file" 2>/dev/null && gzip -t "$r2_output_file" 2>/dev/null; then
        echo "5) Bwa-mem2 output file already exists and is valid for sample: $sample. Skipping sample."
        i=$((i + 1))
        continue
    elif [ -f "$r1_output_file" ] || [ -f "$r2_output_file" ]; then
        echo "5) Bwa-mem2 found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing." | tee -a 0_workflow_progress.txt
        rm -f "$r1_output_file" "$r2_output_file" \
            "5_bwa_reads/${sample}_alignment.log" \
            "5_bwa_reads/${sample}_nofilter_flagstat.txt" \
            "5_bwa_reads/${sample}_hostfilter_flagstat.txt" \
            "5_bwa_reads/${r1_output_file##*/}.md5" \
            "5_bwa_reads/${r2_output_file##*/}.md5"
    fi

    # Run main software (without intermediate files)
    bwa-mem2 mem \
        -t $(nproc --ignore=1) \
        "5_bwa_index/${ref_name}/${ref_name}" \
        "3_fastp/${sample}_trimmed_1.fq.gz" \
        "3_fastp/${sample}_trimmed_2.fq.gz" \
        2> "5_bwa_reads/${sample}_alignment.log" \
        | tee >(samtools flagstat - > "5_bwa_reads/${sample}_nofilter_flagstat.txt") \
        | samtools view -h -b -f 4 -@ $(nproc --ignore=1) - \
        | tee >(samtools flagstat - > "5_bwa_reads/${sample}_hostfilter_flagstat.txt") \
        | samtools sort -n -@ $(nproc --ignore=1) - \
        | samtools fastq -c 6 \
        -1 "5_bwa_reads/${sample}_trimmed_nohost_1.fq.gz" \
        -2 "5_bwa_reads/${sample}_trimmed_nohost_2.fq.gz" \
        -0 /dev/null \
        -s /dev/null \
        -@ $(nproc --ignore=1) \
        -

    # Validate reads files; if corrupted, remove so the next run reprocesses this sample
    echo "Validating R1 file"
    if ! gzip -t "$r1_output_file" 2>> 5_bwa_reads_corrupted.txt; then
        echo "5) Bwa-mem2 R1 output failed integrity check for sample: $sample" | tee -a 0_workflow_progress.txt
        rm -f "$r1_output_file" "$r2_output_file"
    fi
    echo "Validating R2 file"
    if [ -f "$r2_output_file" ] && ! gzip -t "$r2_output_file" 2>> 5_bwa_reads_corrupted.txt; then
        echo "5) Bwa-mem2 R2 output failed integrity check for sample: $sample" | tee -a 0_workflow_progress.txt
        rm -f "$r1_output_file" "$r2_output_file"
    fi

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Generate checksum files for the reads (skip files already verified against an existing .md5)
(cd 5_bwa_reads && for file in *.gz; do
    [ -f "${file}.md5" ] && md5sum -c "${file}.md5" >/dev/null 2>&1 && continue
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done) | tee -a 0_workflow_progress.txt
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
(cd 5_bwa_reads && md5sum -c *.md5) | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 5.4) Bwa-mem2 -> FastQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="5) FastQC"
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "5_bwa_reads_fastqc.tar.gz" ] && [ -f "5_bwa_reads_fastqc.tar.gz.md5" ] && md5sum -c "5_bwa_reads_fastqc.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (5_bwa_reads_fastqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo,
    # since FastQC runs as a single call over all files and has no reliable per-file resume here.
    if [ -d "5_bwa_reads_fastqc" ]; then
        echo "${workflow_step}: Found incomplete 5_bwa_reads_fastqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "5_bwa_reads_fastqc"
    fi
    rm -f "5_bwa_reads_fastqc.tar.gz" "5_bwa_reads_fastqc.tar.gz.md5"

    # Create an output directory
    mkdir -p 5_bwa_reads_fastqc

    # Activate Conda environment
    conda activate fastqc
    # Run main software
    fastqc -t $(nproc --ignore=1) 5_bwa_reads/*.gz -o 5_bwa_reads_fastqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="5_bwa_reads_fastqc.tar.gz"
    itens_to_compress=(5_bwa_reads_fastqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 5.5) Bwa-mem2 -> FastQC -> MultiQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="5) MultiQC"
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "5_bwa_reads_fastqc_multiqc.tar.gz" ] && [ -f "5_bwa_reads_fastqc_multiqc.tar.gz.md5" ] && md5sum -c "5_bwa_reads_fastqc_multiqc.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (5_bwa_reads_fastqc_multiqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo.
    if [ -d "5_bwa_reads_fastqc_multiqc" ]; then
        echo "${workflow_step}: Found incomplete 5_bwa_reads_fastqc_multiqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "5_bwa_reads_fastqc_multiqc"
    fi
    rm -f "5_bwa_reads_fastqc_multiqc.tar.gz" "5_bwa_reads_fastqc_multiqc.tar.gz.md5"

    # Activate Conda environment
    conda activate multiqc
    # Run main software
    multiqc 5_bwa_reads_fastqc/*_fastqc.zip -o 5_bwa_reads_fastqc_multiqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="5_bwa_reads_fastqc_multiqc.tar.gz"
    itens_to_compress=(5_bwa_reads_fastqc_multiqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

    # Delete the output directory
    rm -r 5_bwa_reads_fastqc 5_bwa_reads_fastqc_multiqc
fi

# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt


############################################################
# 6) Taxonomic abundance profile
############################################################

############################################################
## 6.1) Kraken

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) Kraken"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 6_kraken_report
mkdir -p 6_kraken_output

# Calculate sample size
i=1
r1_files=(5_bwa_reads/*_1.fq.gz)
sample_count=${#r1_files[@]}

# Activate Conda environment
conda activate kraken2
# Loop through a list of sample files
for r1 in "${r1_files[@]}"; do
    # Obtain r2 path
    r2=${r1/_1.fq.gz/_2.fq.gz}
    # Extract r1 file name
    filename=${r1##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Skip sample if output file already exists and is a valid, complete gzip file
    output_file="6_kraken_report/${sample}_kreport.tsv"
    gz_file="6_kraken_output/${sample}.kraken.gz"
    if [ -f "$output_file" ] && [ -f "$gz_file" ] && gzip -t "$gz_file" 2>/dev/null; then
        echo "${workflow_step} output files already exist and are valid for sample: $sample. Skipping sample."
        i=$((i + 1))
        continue
    elif [ -f "$output_file" ] || [ -f "$gz_file" ]; then
        echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
        rm -f "$output_file" "$gz_file"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS
    
    # Run main software
    kraken2 \
        --threads $(nproc --ignore=1) \
        --db /db/kraken/k2_pluspf_20260626/ \
        --memory-mapping \
        --paired "${r1}" "${r2}" \
        --report "6_kraken_report/${sample}_kreport.tsv" \
        --report-minimizer-data \
        --minimum-hit-groups 3 \
        --report-zero-counts \
        | pigz -p $(nproc --ignore=1) -c > "6_kraken_output/${sample}.kraken.gz"
    
    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Generate checksum files
(cd 6_kraken_output && for file in *.gz; do
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done) | tee -a 0_workflow_progress.txt
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
(cd 6_kraken_output && md5sum -c *.md5) | tee -a 0_workflow_progress.txt

# # Archive directory
# echo "Archiving directory: 6_kraken_output"
# tar -cvf 6_kraken_output.tar 6_kraken_output
# # Generate checksum of file
# echo "Processing checksum of compressed report file."
# md5sum 6_kraken_output.tar > 6_kraken_output.tar.md5

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 6.2) Kraken -> Bracken

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) Bracken"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 6_bracken_report

# Calculate sample size
i=1
reports=(6_kraken_report/*_kreport.tsv)
sample_count=${#reports[@]}

# Activate Conda environment
conda activate kraken2
# Loop through a list of sample files
for file in "${reports[@]}"; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Skip sample if output file already exists and is a valid, complete gzip file
    output_file="6_bracken_report/${sample}_breport.tsv"
    gz_file="6_bracken_report/${sample}.bracken.gz"
    bracken_tmp="6_bracken_report/${sample}.bracken"
    if [ -f "$output_file" ] && [ -f "$bracken_tmp" ] && [ -f "$gz_file" ] && gzip -t "$gz_file" 2>/dev/null; then

        echo "${workflow_step} output files already exist and are complete for sample: $sample. Skipping sample."
        i=$((i + 1))
        continue

    elif [ -f "$output_file" ] || [ -f "$bracken_tmp" ] || [ -f "$gz_file" ]; then

        echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
        rm -f "$output_file" "$bracken_tmp" "$gz_file"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Run main software
    bracken \
        -d /db/kraken/k2_pluspf_20260626/ \
        -i "6_kraken_report/${sample}_kreport.tsv" \
        -w "6_bracken_report/${sample}_breport.tsv" \
        -o "6_bracken_report/${sample}.bracken" \
        -r 100 \
        -l S \
        -t 0

    # Compress .bracken file
    pigz -k -p $(nproc --ignore=1) "6_bracken_report/${sample}.bracken"
    
    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="6_bracken_report.tar.gz"
itens_to_compress=(6_bracken_report)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 6.3) Kraken -> Bracken -> Comparison

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) Kraken -> Bracken -> Comparison"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "6_bracken_comparison.tar.gz" ] && [ -f "6_bracken_comparison.tar.gz.md5" ] && md5sum -c "6_bracken_comparison.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (6_bracken_comparison.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo,
    # since this step is a single call over all samples with no reliable partial resume.
    if [ -d "6_bracken_comparison" ]; then
        echo "${workflow_step}: Found incomplete 6_bracken_comparison directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "6_bracken_comparison"
    fi
    rm -f "6_bracken_comparison.tar.gz" "6_bracken_comparison.tar.gz.md5"

    # Create an output directory
    mkdir -p 6_bracken_comparison

    # Activate Conda environment
    conda activate kraken2
    # Combine Bracken outputs of all samples
    combine_bracken_outputs.py \
        --files 6_bracken_report/*.bracken \
        -o 6_bracken_comparison/bracken_allsamples_breport.tsv
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="6_bracken_comparison.tar.gz"
    itens_to_compress=(6_bracken_comparison)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 6.4) Kraken -> Bracken -> Krona

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) Krona"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "6_bracken_krona.tar.gz" ] && [ -f "6_bracken_krona.tar.gz.md5" ] && md5sum -c "6_bracken_krona.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (6_bracken_krona.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create an output directory
    mkdir -p 6_bracken_krona

    # Calculate sample size
    i=1
    reports=(6_bracken_report/*_breport.tsv)
    sample_count=${#reports[@]}

    # Activate Conda environment
    conda activate kraken2
    # Loop through a list of sample files
    for file in "${reports[@]}"; do
        # Extract file name
        filename=${file##*/}
        # Extract sample name
        sample=${filename%%_*}

        # Skip sample if output files already exist and look complete
        txt_file="6_bracken_krona/${sample}_bkrona.txt"
        html_file="6_bracken_krona/${sample}_bkrona.html"
        if [ -s "$txt_file" ] && [ -s "$html_file" ]; then
            echo "${workflow_step} output files already exist and are valid for sample: $sample. Skipping sample."
            i=$((i + 1))
            continue
        elif [ -f "$txt_file" ] || [ -f "$html_file" ]; then
            echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
            rm -f "$txt_file" "$html_file"
        fi

        # Inform current sample
        echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

        # Generate .txt file
        kreport2krona.py \
            -r "${file}" \
            -o "$txt_file" \
            --no-intermediate-ranks

        # Generate .html file
        ktImportText "$txt_file" \
            -o "$html_file"

        i=$((i + 1))

    done
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="6_bracken_krona.tar.gz"
    itens_to_compress=(6_bracken_krona)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

    # Delete output directory
    rm -r 6_bracken_krona
fi

# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 6.5) MetaPhlAn

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) MetaPhlAn"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 6_metaphlan

# Calculate sample size
i=1
r1_files=(5_bwa_reads/*_1.fq.gz)
sample_count=${#r1_files[@]}

# Activate Conda environment
conda activate metaphlan
# Loop through a list of sample files
for r1 in "${r1_files[@]}"; do

    # Obtain r2 path
    r2=${r1/_1.fq.gz/_2.fq.gz}
    # Extract r1 file name
    filename=${r1##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Skip sample if output files already exist and look complete
    output_file="6_metaphlan/${sample}_mprofile.txt"
    mapout_file="6_metaphlan/${sample}_metaphlan_bowtie2.bz2"
    if [ -s "$output_file" ] && [ -f "$mapout_file" ]; then
        echo "${workflow_step} output file already exists for sample: $sample ($output_file). Skipping sample."
        i=$((i + 1))
        continue
    elif [ -f "$output_file" ] || [ -f "$mapout_file" ]; then
        echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
        rm -f "$output_file" "$mapout_file"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Run main software
    metaphlan \
        "${r1}","${r2}" \
        --input_type fastq \
        --nproc $(nproc --ignore=1) \
        --verbose \
        --db_dir /db/metaphlan/ \
        --mapout "$mapout_file" \
        -o "$output_file"

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="6_metaphlan.tar.gz"
itens_to_compress=(6_metaphlan)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 6.6) MetaPhlAn -> Comparison

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) MetaPhlAn -> Comparison"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "6_metaphlan_comparison.tar.gz" ] && [ -f "6_metaphlan_comparison.tar.gz.md5" ] && md5sum -c "6_metaphlan_comparison.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (6_metaphlan_comparison.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo,
    # since this step is a single call over all samples with no reliable partial resume.
    if [ -d "6_metaphlan_comparison" ]; then
        echo "${workflow_step}: Found incomplete 6_metaphlan_comparison directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "6_metaphlan_comparison"
    fi
    rm -f "6_metaphlan_comparison.tar.gz" "6_metaphlan_comparison.tar.gz.md5"

    # Create an output directory
    mkdir -p 6_metaphlan_comparison

    # Merge abundance tables of all samples
    # Activate Conda environment
    conda activate metaphlan
    merge_metaphlan_tables.py 6_metaphlan/*_mprofile.txt \
        > 6_metaphlan_comparison/metaphlan_allsamples_mprofile.txt
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="6_metaphlan_comparison.tar.gz"
    itens_to_compress=(6_metaphlan_comparison)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt


############################################################
# 7) Metagenome assembly
############################################################

############################################################
## 7.1) MEGAHIT

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="7) MEGAHIT"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 7_megahit

# Calculate sample size
i=1
r1_files=(5_bwa_reads/*_1.fq.gz)
sample_count=${#r1_files[@]}

# Activate Conda environment
conda activate megahit
# Loop through a list of sample files
for r1 in "${r1_files[@]}"; do

    # Obtain r2 path
    r2=${r1/_1.fq.gz/_2.fq.gz}
    # Extract r1 file name
    filename=${r1##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Skip sample if output file already exists and looks complete (non-empty fasta)
    output_file="7_megahit/${sample}_megahit.fasta"
    log_file="7_megahit/${sample}_megahit.log"
    if [ -s "$output_file" ] && [ -f "$log_file" ]; then
        echo "7) MEGAHIT output file already exists and is valid for sample: $sample ($output_file). Skipping assembly."
        i=$((i + 1))
        continue
    elif [ -f "$output_file" ] || [ -f "$log_file" ]; then
        echo "7) MEGAHIT found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
        rm -f "$output_file" "$log_file"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Remove leftover output directory from a previous interrupted run
    [ -d "7_megahit/${sample}_megahit" ] && rm -rf "7_megahit/${sample}_megahit"

    # Run main software
    megahit \
        -t $(nproc --ignore=1) \
        -1 "${r1}" \
        -2 "${r2}" \
        -o "7_megahit/${sample}_megahit" \
        --min-contig-len 500

    # Move and rename assembly file
    mv "7_megahit/${sample}_megahit/final.contigs.fa" "$output_file"

    # Move and rename log file
    mv "7_megahit/${sample}_megahit/log" "$log_file"

    # Delete the sample directory
    rm -r "7_megahit/${sample}_megahit/"

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="7_megahit.tar.gz"
itens_to_compress=(7_megahit)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 7.2) MEGAHIT -> QUAST

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="7) QUAST"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "7_megahit_quast.tar.gz" ] && [ -f "7_megahit_quast.tar.gz.md5" ] && md5sum -c "7_megahit_quast.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (7_megahit_quast.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create an output directory
    mkdir -p 7_megahit_quast

    # QUAST for each sample
    # Activate Conda environment
    conda activate quast
    #Loop
    fasta_files=(7_megahit/*.fasta)
    for file in "${fasta_files[@]}"; do

        # Extract file name
        filename=${file##*/}
        # Extract sample name
        sample=${filename%%_*}

        # Skip sample if output already exists and looks complete
        report_file="7_megahit_quast/${sample}_quast/report.tsv"
        if [ -s "$report_file" ]; then
            echo "${workflow_step} output already exists and is valid for sample: $sample. Skipping."
            continue
        elif [ -d "7_megahit_quast/${sample}_quast" ]; then
            echo "${workflow_step} found incomplete output for sample: $sample. Removing partial directory and reprocessing."
            rm -rf "7_megahit_quast/${sample}_quast"
        fi

        # Inform current sample
        echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')"
        # Start counting the running time
        loop_start_time=$SECONDS

        # Create output directory
        mkdir -p "7_megahit_quast/${sample}_quast"

        # Run main software
        quast.py \
            -t $(nproc --ignore=1) \
            -m 0 \
            -o "7_megahit_quast/${sample}_quast" \
            "${file}"
    done
    # Deactivate Conda environment
    conda deactivate

    # QUAST for all samples
    all_report="7_megahit_quast/all_samples_quast/report.tsv"
    if [ -s "$all_report" ]; then
        echo "${workflow_step} combined all-samples output already exists and is valid. Skipping."
    else
        [ -d "7_megahit_quast/all_samples_quast" ] && rm -rf "7_megahit_quast/all_samples_quast"
        # Activate Conda environment
        conda activate quast
        quast.py \
            -t $(nproc --ignore=1) \
            -m 0 \
            -o "7_megahit_quast/all_samples_quast" \
            7_megahit/*.fasta
        # Deactivate Conda environment
        conda deactivate
    fi

    # Compress the output directory
    compressed_file="7_megahit_quast.tar.gz"
    itens_to_compress=(7_megahit_quast)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt


############################################################
# 8) Functional abundance profile (per sample)
############################################################

############################################################
## 8.1) Aragorn

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) Aragorn"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 8_aragorn

# Calculate the sample count to display loop progress
i=1
fasta_files=(7_megahit/*.fasta)
sample_count=${#fasta_files[@]}

# Activate Conda environment
conda activate aragorn
# Loop through a list of sample files
for file in "${fasta_files[@]}"; do

    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Skip sample if output file already exists and looks complete (non-empty)
    output_file="8_aragorn/${sample}_aragorn/${sample}_aragorn.txt"
    if [ -s "$output_file" ]; then
        echo "8) Aragorn output file already exists and is valid for sample: $sample ($output_file). Skipping sample."
        i=$((i + 1))
        continue
    elif [ -d "8_aragorn/${sample}_aragorn" ]; then
        echo "8) Aragorn found incomplete output for sample: $sample. Removing partial directory and reprocessing."
        rm -rf "8_aragorn/${sample}_aragorn"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create output directory
    mkdir -p "8_aragorn/${sample}_aragorn"

    # Run main software
    aragorn \
        "$file" \
        > "$output_file"
    
    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="8_aragorn.tar.gz"
itens_to_compress=(8_aragorn)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 8.2) Pybarrnap

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) Pybarrnap"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 8_pybarrnap

# Calculate the sample count to display loop progress
i=1
fasta_files=(7_megahit/*.fasta)
sample_count=${#fasta_files[@]}

# Activate Conda environment
conda activate pybarrnap
# Loop through a list of sample files
for file in "${fasta_files[@]}"; do

    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Skip sample if BOTH bac and arc output files already exist and look complete
    bac_output="8_pybarrnap/${sample}_pybarrnap/${sample}_bac_pybarrnap.fasta"
    arc_output="8_pybarrnap/${sample}_pybarrnap/${sample}_arc_pybarrnap.fasta"
    if [ -s "$bac_output" ] && [ -s "$arc_output" ]; then
        echo "${workflow_step} output files already exist and are valid for sample: $sample. Skipping sample."
        i=$((i + 1))
        continue
    elif [ -d "8_pybarrnap/${sample}_pybarrnap" ]; then
        echo "${workflow_step} found incomplete output for sample: $sample. Removing partial directory and reprocessing."
        rm -rf "8_pybarrnap/${sample}_pybarrnap"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create output directory
    mkdir -p "8_pybarrnap/${sample}_pybarrnap"

    # Run barrnap for bacteria
    pybarrnap \
        --threads $(nproc --ignore=1) \
        --quiet \
        --kingdom bac \
        --outseq "$bac_output" \
        "$file" \
        > "8_pybarrnap/${sample}_pybarrnap/${sample}_bac_pybarrnap.gff" \
        2> "8_pybarrnap/${sample}_pybarrnap/${sample}_bac_pybarrnap.log"
    
    # Run barrnap for archea
    pybarrnap \
        --threads $(nproc --ignore=1) \
        --quiet \
        --kingdom arc \
        --outseq "$arc_output" \
        "$file" \
        > "8_pybarrnap/${sample}_pybarrnap/${sample}_arc_pybarrnap.gff" \
        2> "8_pybarrnap/${sample}_pybarrnap/${sample}_arc_pybarrnap.log"
    
    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="8_pybarrnap.tar.gz"
itens_to_compress=(8_pybarrnap)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 8.3) Pyrodigal

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) Pyrodigal"

# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 8_pyrodigal

# Calculate the sample count to display loop progress
i=1
files=(7_megahit/*.fasta)
sample_count=${#files[@]}

# Activate Conda environment
conda activate pyrodigal
# Loop through a list of sample files
for file in "${files[@]}"; do

    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Skip sample if all output files already exist and look complete
    faa_output="8_pyrodigal/${sample}_pyrodigal/${sample}.faa"
    ffn_output="8_pyrodigal/${sample}_pyrodigal/${sample}.ffn"
    gff_output="8_pyrodigal/${sample}_pyrodigal/${sample}.gff"
    if [ -s "$faa_output" ] && [ -s "$ffn_output" ] && [ -s "$gff_output" ]; then
        echo "${workflow_step} output files already exist and are valid for sample: $sample ($faa_output). Skipping sample."
        i=$((i + 1))
        continue
    elif [ -d "8_pyrodigal/${sample}_pyrodigal" ]; then
        echo "${workflow_step} found incomplete output for sample: $sample. Removing partial directory and reprocessing."
        rm -rf "8_pyrodigal/${sample}_pyrodigal"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create output directory
    mkdir -p "8_pyrodigal/${sample}_pyrodigal"

    # Run main software
    pyrodigal \
        -j $(nproc --ignore=1) \
        -m \
        -p meta \
        --no-stop-codon \
        -f gff \
        -i "${file}" \
        -d "$ffn_output" \
        -a "$faa_output" \
        -o "$gff_output"

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="8_pyrodigal.tar.gz"
itens_to_compress=(8_pyrodigal)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 8.4) Pyrodigal -> AMRFinderPlus

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) AMRFinderPlus"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 8_pyrodigal_amrfinder

# Calculate the sample count to display loop progress
i=1
faa_files=(8_pyrodigal/*/*.faa)
sample_count=${#faa_files[@]}

# Activate Conda environment
conda activate amrfinder
# Loop through a list of sample files
for file in "${faa_files[@]}"; do

    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.faa}

    # Skip sample if output file already exists and looks complete (non-empty)
    output_file="8_pyrodigal_amrfinder/${sample}_amrfinder.tsv"
    format_faa="8_pyrodigal_amrfinder/${sample}_amrfinder_format.faa"
    nuc_out="8_pyrodigal_amrfinder/${sample}_amrfinder.fasta"
    prot_out="8_pyrodigal_amrfinder/${sample}_amrfinder.faa"
    if [ -s "$output_file" ] && [ ! -f "$format_faa" ]; then
        echo "8) AMRFinderPlus output file already exists and is valid for sample: $sample ($output_file). Skipping sample."
        i=$((i + 1))
        continue
    elif [ -f "$output_file" ] || [ -f "$format_faa" ] || [ -f "$nuc_out" ] || [ -f "$prot_out" ]; then
        echo "8) AMRFinderPlus found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
        rm -f "$output_file" "$format_faa" "$nuc_out" "$prot_out"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    # Start counting the running time
    loop_start_time=$SECONDS

    # Method 1 - Run AMRFinderPlus using nucleotide and proteins sequences
    # Format faa files for AMRFinderPlus
    gawk '{
      if ($0 ~ /^>/) {
        match($0, /^>([^ ]+)/, a)
        id=a[1]
        sub(/ID=[^;]+/, "ID="id)
        print
      } else print
    }' "${file}" > "$format_faa"
    # Run main software
    amrfinder \
        --threads $(nproc --ignore=1) \
        --database /db/amrfinder/latest \
        --plus \
        --annotation_format prodigal \
        --nucleotide "7_megahit/${sample}_megahit.fasta" \
        --protein "$format_faa"  \
        --gff "8_pyrodigal/${sample}_pyrodigal/${sample}.gff" \
        --name "${sample}" \
        --nucleotide_output "$nuc_out"\
        --protein_output "$prot_out"\
        --output "$output_file"
    # Remove intermediate file
    rm "$format_faa"

    # # Method 2 - Run AMRFinderPlus using only proteins sequences
    # # Run main software
    # amrfinder \
    # --threads $(nproc --ignore=1) \
    # --database /db/amrfinder/latest \
    # --plus \
    # --protein "${file}"  \
    # --name "${sample}" \
    # --protein_output "8_pyrodigal_amrfinder/${sample}_amrfinder.faa"\
    # --output "8_pyrodigal_amrfinder/${sample}_amrfinder.tsv"

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="8_pyrodigal_amrfinder.tar.gz"
itens_to_compress=(8_pyrodigal_amrfinder)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 8.5) Pyrodigal -> dbCAN

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) dbCAN"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 8_pyrodigal_dbcan

# Calculate the sample count to display loop progress
i=1
faa_files=(8_pyrodigal/*/*.faa)
sample_count=${#faa_files[@]}

# Activate Conda environment
conda activate dbcan
# Loop through a list of sample files
for file in "${faa_files[@]}"; do

    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.faa}

    # Skip sample if output file already exists and looks complete (non-empty)
    output_file="8_pyrodigal_dbcan/${sample}_dbcan/overview.tsv"
    if [ -s "$output_file" ]; then
        echo "${workflow_step} output file already exists and is valid for sample: $sample ($output_file). Skipping sample."
        i=$((i + 1))
        continue
    elif [ -d "8_pyrodigal_dbcan/${sample}_dbcan" ]; then
        echo "${workflow_step} found incomplete output for sample: $sample. Removing partial directory and reprocessing."
        rm -rf "8_pyrodigal_dbcan/${sample}_dbcan"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create output directory
    mkdir -p "8_pyrodigal_dbcan/${sample}_dbcan"

    # Run main software
    run_dbcan CAZyme_annotation \
        --threads $(nproc --ignore=1) \
        --db_dir /db/dbcan \
        --mode protein \
        --input_raw_data $file \
        --output_dir "8_pyrodigal_dbcan/${sample}_dbcan"

    # Delete unecessary file and directory
    rm -f "8_pyrodigal_dbcan/${sample}_dbcan/uniInput.faa"
    rm -rf "8_pyrodigal_dbcan/${sample}_dbcan/tmp"

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="8_pyrodigal_dbcan.tar.gz"
itens_to_compress=(8_pyrodigal_dbcan)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 8.6) Pyrodigal -> eggNOG-mapper

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) eggNOG-mapper"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 8_pyrodigal_emapper

# Calculate the sample count to display loop progress
i=1
faa_files=(8_pyrodigal/*/*.faa)
sample_count=${#faa_files[@]}

# Activate Conda environment
conda activate eggnog-mapper
# Loop through a list of sample files
for file in "${faa_files[@]}"; do

    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.faa}

    # Skip sample if output files already exist and look complete (non-empty, no leftover original)
    output_file="8_pyrodigal_emapper/${sample}_emapper/${sample}_emapper.tsv"
    comments_file="8_pyrodigal_emapper/${sample}_emapper/${sample}_emapper_comments.tsv"
    original_file="8_pyrodigal_emapper/${sample}_emapper/${sample}.emapper.annotations"
    if [ -s "$output_file" ] && [ -f "$comments_file" ] && [ ! -f "$original_file" ]; then
        echo "${workflow_step} output file already exists and is valid for sample: $sample ($output_file). Skipping sample."
        i=$((i + 1))
        continue
    elif [ -d "8_pyrodigal_emapper/${sample}_emapper" ]; then
        echo "${workflow_step} found incomplete output for sample: $sample. Removing partial directory and reprocessing."
        rm -rf "8_pyrodigal_emapper/${sample}_emapper"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create output directory
    mkdir -p "8_pyrodigal_emapper/${sample}_emapper"

    # Run main software
    emapper.py \
        --cpu $(nproc --ignore=1) \
        -i $file \
        --output "${sample}" \
        --output_dir "8_pyrodigal_emapper/${sample}_emapper"

    # Adjust output table
    grep -v '^##' "$original_file" \
    > "$output_file"

    # Create a table of results comments
    grep '^##' "$original_file" \
    > "$comments_file"

    # Delete original table
    rm "$original_file"

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="8_pyrodigal_emapper.tar.gz"
itens_to_compress=(8_pyrodigal_emapper)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

# ############################################################
# ## 8.7) Pyrodigal -> VFDB (BLASTP)

# # Avoid literal glob pattern
# shopt -s nullglob

# # Software name for tracking progress in 0_workflow_progress.txt
# workflow_step="8) VFDB (BLASTP)"
# # Update the file 0_workflow_progress.txt
# echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# # Start counting the running time
# start_time=$SECONDS

# # Create an output directory
# mkdir -p 8_pyrodigal_vfdb

# # Calculate the sample count to display loop progress
# i=1
# faa_files=(8_pyrodigal/*/*.faa)
# sample_count=${#faa_files[@]}

# # Activate Conda environment
# conda activate eggnog-mapper
# # Loop through a list of sample files
# for file in "${faa_files[@]}"; do

#     # Extract file name
#     filename=${file##*/}
#     # Extract sample name
#     sample=${filename%%.faa}

#     # Skip sample if output files already exist and look complete (non-empty, no leftover original)

#     # Inform current sample
#     echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
#     # Start counting the running time
#     loop_start_time=$SECONDS

#     # Create output directory
#     mkdir -p "8_pyrodigal_vfdb/${sample}_vfdb"

#     # Run main software
#     vfdb.py \
#         --cpu $(nproc --ignore=1) \
#         -i $file \
#         --output "${sample}" \
#         --output_dir "8_pyrodigal_vfdb/${sample}_vfdb"

#     # Stop counting the running time
#     loop_elapsed_time=$((SECONDS - $loop_start_time))
#     # Calculate the running time
#     loop_hours=$((loop_elapsed_time / 3600))
#     loop_minutes=$(((loop_elapsed_time % 3600) / 60))
#     loop_seconds=$((loop_elapsed_time % 60))
#     loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
#     # Show the running time
#     echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

#     i=$((i + 1))

# done
# # Deactivate Conda environment
# conda deactivate

# # Compress the output directory
# compressed_file="8_pyrodigal_vfdb.tar.gz"
# itens_to_compress=(8_pyrodigal_vfdb)
# echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
# tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# # Generate checksum file of compressed directory file
# md5sum "${compressed_file}" > "${compressed_file}".md5
# # Check file integrity
# echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
# md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# # Stop counting the running time
# elapsed_time=$((SECONDS - $start_time))
# # Calculate the running time
# hours=$((elapsed_time / 3600))
# minutes=$(((elapsed_time % 3600) / 60))
# seconds=$((elapsed_time % 60))
# running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# # Show the running time
# echo "${workflow_step} running time: ${running_time}" | tee -a 0_workflow_progress.txt
# # Update the file 0_workflow_progress.txt
# echo "${workflow_step} step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 9) Functional abundance profile (clustering)
############################################################

############################################################
## 9.1) MMseqs2 input (SeqKit)

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) MMseqs2 input"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "9_mmseqs_input.tar.gz" ] && [ -f "9_mmseqs_input.tar.gz.md5" ] && md5sum -c "9_mmseqs_input.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (9_mmseqs_input.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left partial/duplicated concatenated files
    # (this step uses append, which is not safely resumable per-sample); wipe and redo entirely.
    if [ -d "9_mmseqs_input" ]; then
        echo "${workflow_step}: Found incomplete 9_mmseqs_input directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "9_mmseqs_input"
    fi
    rm -f "9_mmseqs_input.tar.gz" "9_mmseqs_input.tar.gz.md5"

    # Create an output directory
    mkdir -p 9_mmseqs_input

    # Pattern of samples to exclude (NHP120 up to NHP130) (Pool samples)
    exclude_pattern="NHP(12[0-9]|130)_pyrodigal"

    # Build the filtered list of files once (nucleotide + corresponding protein)
    ffn_files=()
    faa_files=()
    for f in 8_pyrodigal/*/*.ffn; do
        [[ "$f" =~ $exclude_pattern ]] && continue
        ffn_files+=("$f")
        faa_files+=("${f%.ffn}.faa")
    done

    # Calculate the sample count to display loop progress (excluding filtered samples)
    i=1
    sample_count=${#ffn_files[@]}

    # Remove existing contatenated files
    rm -f 9_mmseqs_input/all_samples.ffn
    rm -f 9_mmseqs_input/all_samples.faa

    # Create contatenated files
    # Loop through a list of sample files
    for file in "${ffn_files[@]}"; do

        # Skip samples specified in the exclude pattern
        if [[ "$file" =~ $exclude_pattern ]]; then
            continue
        fi

        # Extract file name
        filename=${file##*/}
        # Extract sample name
        sample=${filename%%.ffn}

        # Inform current sample
        echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Concatenate nucleotide sequences and clean headers
        awk -v sample_awk="$sample" '/^>/ {print ">" sample_awk "_" substr($1,2); next} {print}' "$file" >> 9_mmseqs_input/all_samples.ffn
        # Corresponding protein file
        faa_file="${file%.ffn}.faa"
        # Concatenate protein sequences and clean headers
        awk -v sample_awk="$sample" '/^>/ {print ">" sample_awk "_" substr($1,2); next} {print}' "$faa_file" >> 9_mmseqs_input/all_samples.faa

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

        i=$((i + 1))

    done

    # Verification of the number of sequences
    echo "9) MMseqs2 input: Sequences on all nucleotide files: $(grep -h '^>' "${ffn_files[@]}" | wc -l)" | tee -a 0_workflow_progress.txt
    echo "9) MMseqs2 input: Sequences on concatenated nucleotide file: $(grep -h '^>' 9_mmseqs_input/all_samples.ffn | wc -l)" | tee -a 0_workflow_progress.txt
    echo "9) MMseqs2 input: Sequences on all protein files: $(grep -h '^>' "${faa_files[@]}" | wc -l)" | tee -a 0_workflow_progress.txt
    echo "9) MMseqs2 input: Sequences on concatenated protein file: $(grep -h '^>' 9_mmseqs_input/all_samples.faa | wc -l)" | tee -a 0_workflow_progress.txt

    # Extract and clean the ids from protein files
    grep "^>" 9_mmseqs_input/all_samples.faa | awk '{print $1}' | sed 's/>//' | sort > 9_mmseqs_input/all_samples_faa_ids.txt
    # Extract and clean the ids from nucleotide files
    grep "^>" 9_mmseqs_input/all_samples.ffn | awk '{print $1}' | sed 's/>//' | sort > 9_mmseqs_input/all_samples_ffn_ids.txt
    # Protein ids absent in nucleotide ids
    comm -23 9_mmseqs_input/all_samples_faa_ids.txt 9_mmseqs_input/all_samples_ffn_ids.txt > 9_mmseqs_input/all_samples_faa_ids_exclusive.txt
    # Nucleotide ids absent in protein ids
    comm -23 9_mmseqs_input/all_samples_ffn_ids.txt 9_mmseqs_input/all_samples_faa_ids.txt > 9_mmseqs_input/all_samples_ffn_ids_exclusive.txt
    # Delete id files
    rm 9_mmseqs_input/all_samples_faa_ids.txt 9_mmseqs_input/all_samples_ffn_ids.txt

    # Warn if protein and nucleotide catalogs are not perfectly matched
    n_faa_exclusive=$(wc -l < 9_mmseqs_input/all_samples_faa_ids_exclusive.txt)
    n_ffn_exclusive=$(wc -l < 9_mmseqs_input/all_samples_ffn_ids_exclusive.txt)
    if [ "$n_faa_exclusive" -gt 0 ]; then
        echo "WARNING: $n_faa_exclusive protein IDs have no matching nucleotide sequence (see all_samples_faa_ids_exclusive.txt)" | tee -a 0_workflow_progress.txt
    fi
    if [ "$n_ffn_exclusive" -gt 0 ]; then
        echo "WARNING: $n_ffn_exclusive nucleotide IDs have no matching protein sequence (see all_samples_ffn_ids_exclusive.txt)" | tee -a 0_workflow_progress.txt
    fi

    # Filter 1: Sequence size (>= 100 nt)
    echo "9) MMseqs2 input: Filtering sequences by size"
    # Activate Conda environment
    conda activate seqkit
    # Extract the corresponding nucleotide sequences
    seqkit seq \
        -m 100 9_mmseqs_input/all_samples.ffn \
        > 9_mmseqs_input/all_samples_f1.ffn
    # Extract representative sequence IDs from the Filter 1 catalog
    seqkit seq \
        -n -i 9_mmseqs_input/all_samples_f1.ffn \
        > 9_mmseqs_input/all_samples_f1_ids.txt
    # Extract the corresponding protein sequences
    seqkit grep \
        --id-regexp "^(\S+)" \
        -f 9_mmseqs_input/all_samples_f1_ids.txt \
        9_mmseqs_input/all_samples.faa \
        > 9_mmseqs_input/all_samples_f1.faa
    # Deactivate Conda environment
    conda deactivate
    # Remove the temporary file
    rm 9_mmseqs_input/all_samples_f1_ids.txt

    # Compare the number of sequences before and after filtering
    n_before=$(grep -c '^>' 9_mmseqs_input/all_samples.ffn)
    n_after=$(grep -c '^>' 9_mmseqs_input/all_samples_f1.ffn)
    n_removed=$((n_before - n_after))
    pct_removed=$(awk -v b="$n_before" -v a="$n_after" 'BEGIN{printf "%.2f", (b-a)/b*100}')
    echo "Sequences before filtering by size: $n_before" | tee -a 0_workflow_progress.txt
    echo "Sequences after filtering by size: $n_after" | tee -a 0_workflow_progress.txt
    echo "Removed sequences: $n_removed ($pct_removed%)" | tee -a 0_workflow_progress.txt

    # Compress original files
    echo "Compressing original concatenated files"
    pigz 9_mmseqs_input/all_samples.faa
    pigz 9_mmseqs_input/all_samples.ffn

    # Compress the output directory
    compressed_file="9_mmseqs_input.tar.gz"
    itens_to_compress=(9_mmseqs_input)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 9.2) MMseqs2 easy-linclust

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) MMseqs2 easy-linclust"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "9_mmseqs.tar.gz" ] && [ -f "9_mmseqs.tar.gz.md5" ] && md5sum -c "9_mmseqs.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (9_mmseqs.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left partial clustering output or temp files;
    # this step is a single call over the whole catalog with no reliable partial resume.
    if [ -d "9_mmseqs" ]; then
        echo "${workflow_step}: Found incomplete 9_mmseqs directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "9_mmseqs"
    fi
    rm -f "9_mmseqs.tar.gz" "9_mmseqs.tar.gz.md5"

    # MMseqs2
    echo "9) MMseqs2 is running on the concatenated database" | tee -a 0_workflow_progress.txt

    # Create an output directory
    mkdir -p 9_mmseqs/tmp

    # Activate Conda environment
    conda activate MMseqs2
        # Run main software
        mmseqs easy-linclust \
            9_mmseqs_input/all_samples_f1.ffn \
            9_mmseqs/all_samples_f1_nr \
            9_mmseqs/tmp \
            --threads $(nproc --ignore=1) \
            --min-seq-id 0.95 \
            -c 0.9 --cov-mode 1 \
            --split-memory-limit 90G \
            --kmer-per-seq-scale 0.3
    # Deactivate Conda environment
    conda deactivate

    # Delete intermediary file
    rm -f 9_mmseqs/*_all_seqs.fasta

    # Rename nr result file
    mv 9_mmseqs/all_samples_f1_nr_rep_seq.fasta 9_mmseqs/all_samples_f1_nr.ffn

    echo "${workflow_step}: Extracting corresponding protein sequences..."
    # Activate Conda environment
    conda activate seqkit 
    # Extract representative sequence IDs from the non-redundant sequences catalog
    seqkit seq \
        -n -i 9_mmseqs/all_samples_f1_nr.ffn \
        > 9_mmseqs/all_samples_f1_nr_ids.txt
    # Extract the corresponding protein sequences
    seqkit grep \
        --id-regexp "^(\S+)" \
        -f 9_mmseqs/all_samples_f1_nr_ids.txt \
        9_mmseqs_input/all_samples_f1.faa \
        > 9_mmseqs/all_samples_f1_nr.faa
    # Deactivate Conda environment
    conda deactivate

    echo "${workflow_step}: Number and percentage of redundant sequences"
    # Compare the number of sequences before and after clustering
    n_before=$(grep -c '^>' 9_mmseqs_input/all_samples_f1.ffn)
    n_after=$(grep -c '^>' 9_mmseqs/all_samples_f1_nr.ffn)
    n_removed=$((n_before - n_after))
    pct_removed=$(awk -v b="$n_before" -v a="$n_after" 'BEGIN{printf "%.2f", (b-a)/b*100}')
    echo "Sequences before clustering: $n_before" | tee -a 0_workflow_progress.txt
    echo "Sequences after clustering: $n_after" | tee -a 0_workflow_progress.txt
    echo "Removed redundant sequences: $n_removed ($pct_removed%)" | tee -a 0_workflow_progress.txt

    # Compress the output directory
    compressed_file="9_mmseqs.tar.gz"
    itens_to_compress=(9_mmseqs)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
##  9.3) MMseqs2 -> AMRFinderPlus

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) MMseqs2 -> AMRFinderPlus"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "9_mmseqs_amrfinder.tar.gz" ] && [ -f "9_mmseqs_amrfinder.tar.gz.md5" ] && md5sum -c "9_mmseqs_amrfinder.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (9_mmseqs_amrfinder.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo,
    # since this is a single call with no reliable partial resume.
    if [ -d "9_mmseqs_amrfinder" ]; then
        echo "${workflow_step}: Found incomplete 9_mmseqs_amrfinder directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "9_mmseqs_amrfinder"
    fi
    rm -f "9_mmseqs_amrfinder.tar.gz" "9_mmseqs_amrfinder.tar.gz.md5"

    # Create an output directory
    mkdir -p 9_mmseqs_amrfinder

    # Activate Conda environment
    conda activate amrfinder
    # Run main software
    amrfinder \
        --threads $(nproc --ignore=1) \
        --database /db/amrfinder/latest \
        --plus \
        --protein 9_mmseqs/all_samples_f1_nr.faa \
        --name "all_samples_f1_nr" \
        --protein_output "9_mmseqs_amrfinder/all_samples_f1_nr_amrfinder.faa" \
        --output "9_mmseqs_amrfinder/all_samples_f1_nr_amrfinder.tsv"
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="9_mmseqs_amrfinder.tar.gz"
    itens_to_compress=(9_mmseqs_amrfinder)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 9.4) MMseqs2 -> dbCAN

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) MMseqs2 -> dbCAN"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "9_mmseqs_dbcan.tar.gz" ] && [ -f "9_mmseqs_dbcan.tar.gz.md5" ] && md5sum -c "9_mmseqs_dbcan.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (9_mmseqs_dbcan.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create an output directory
    mkdir -p 9_mmseqs_dbcan/tmp

    # Split the input file (only if not already split from a previous run)
    if [ -z "$(ls -A 9_mmseqs_dbcan/tmp/*.faa 2>/dev/null)" ]; then
        echo "${workflow_step}: Spliting input files in parts" | tee -a 0_workflow_progress.txt
        conda activate seqkit 
        seqkit split2 -s 100000 9_mmseqs/all_samples_f1_nr.faa -O 9_mmseqs_dbcan/tmp
        conda deactivate
    else
        echo "${workflow_step}: Split input files already present from a previous run. Reusing them." | tee -a 0_workflow_progress.txt
    fi

    # Calculate the sample count to display loop progress
    i=1
    faa_files=(9_mmseqs_dbcan/tmp/*.faa)
    sample_count=${#faa_files[@]}

    # Activate Conda environment
    conda activate dbcan
    # Loop through a list of sample files
    for file in "${faa_files[@]}"; do
        # Extract file name
        filename=${file##*/}
        # Extract sample name
        sample=${filename%%.faa}

        # Skip part if output file already exists and looks complete (non-empty)
        output_file="9_mmseqs_dbcan/tmp/${sample}_dbcan/overview.tsv"
        if [ -s "$output_file" ]; then
            echo "${workflow_step} output file already exists and is valid for part: $sample ($output_file). Skipping part."
            i=$((i + 1))
            continue
        fi

        # Inform current part
        echo "${workflow_step} is processing part: ${sample} (${i}/${sample_count})"

        # Delete a previous incomplete run
        rm -rf "9_mmseqs_dbcan/tmp/${sample}_dbcan"
        # Create output directory
        mkdir -p "9_mmseqs_dbcan/tmp/${sample}_dbcan"

        # Run main software
        run_dbcan CAZyme_annotation \
            --threads $(nproc --ignore=1) \
            --db_dir /db/dbcan \
            --mode protein \
            --input_raw_data "$file" \
            --output_dir "9_mmseqs_dbcan/tmp/${sample}_dbcan"

        i=$((i + 1))

    done
    # Deactivate Conda environment
    conda deactivate

    # Concatenate output files
    echo "${workflow_step}: Concatenating result files" | tee -a 0_workflow_progress.txt
    # List of tabular files generated by dbCAN
    output_files=(
        "dbCAN_hmm_results.tsv"
        "dbCANsub_hmm_raw.tsv"
        "dbCANsub_hmm_results.tsv"
        "diamond.out"
        "overview.tsv"
    )
    for output_file in "${output_files[@]}"; do
        # Capture all matching file paths into a Bash array
        files=(9_mmseqs_dbcan/tmp/*_dbcan/"$output_file")
        # Check if at least one valid file exists
        if (( ${#files[@]} )); then
            target_out="9_mmseqs_dbcan/$output_file"
            # Write the header from the first available file
            if [[ -s "${files[0]}" ]]; then
                head -n 1 "${files[0]}" > "$target_out"
            fi
            # Append data file by file to avoid 'Argument list too long' with awk
            for file_path in "${files[@]}"; do
                tail -n +2 "$file_path" >> "$target_out"
            done
        fi
    done

    # Delete temporary directory
    rm -rf 9_mmseqs_dbcan/tmp

    # Compress the output directory
    compressed_file="9_mmseqs_dbcan.tar.gz"
    itens_to_compress=(9_mmseqs_dbcan)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
##  9.5) MMseqs2 -> eggNOG-mapper

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) MMseqs2 -> eggNOG-mapper"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "9_mmseqs_emapper.tar.gz" ] && [ -f "9_mmseqs_emapper.tar.gz.md5" ] && md5sum -c "9_mmseqs_emapper.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (9_mmseqs_emapper.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial directory or archive; wipe and redo,
    # since this is a single call with no reliable partial resume.
    if [ -d "9_mmseqs_emapper" ]; then
        echo "${workflow_step}: Found incomplete 9_mmseqs_emapper directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "9_mmseqs_emapper"
    fi
    rm -f "9_mmseqs_emapper.tar.gz" "9_mmseqs_emapper.tar.gz.md5"

    # Create an output directory
    mkdir -p 9_mmseqs_emapper

    # Activate Conda environment
    conda activate eggnog-mapper
    # Run main software
    emapper.py \
        --cpu $(nproc --ignore=1) \
        -i 9_mmseqs/all_samples_f1_nr.faa \
        --output all_samples_f1_nr \
        --output_dir 9_mmseqs_emapper
    # Deactivate Conda environment
    conda deactivate

    # Create result table without comments
    grep -v '^##' 9_mmseqs_emapper/all_samples_f1_nr.emapper.annotations \
        > 9_mmseqs_emapper/all_samples_f1_nr_emapper.tsv

    # Create a table of results comments
    grep '^##' 9_mmseqs_emapper/all_samples_f1_nr.emapper.annotations \
        > 9_mmseqs_emapper/all_samples_f1_nr_emapper_comments.tsv

    # Delete original table
    rm 9_mmseqs_emapper/all_samples_f1_nr.emapper.annotations

    # Compress files
    find 9_mmseqs_emapper/ -maxdepth 1 -type f ! -name '*.gz' -print0 | \
        xargs -0 pigz -p $(nproc --ignore=1)

    # Compress the output directory
    compressed_file="9_mmseqs_emapper.tar.gz"
    itens_to_compress=(9_mmseqs_emapper)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

# ############################################################
# ##  9.6) MMseqs2 -> VFDB (BLASTP)

# # Avoid literal glob pattern
# shopt -s nullglob

# # Software name for tracking progress in 0_workflow_progress.txt
# workflow_step="9) MMseqs2 -> VFDB (BLASTP)"
# # Update the file 0_workflow_progress.txt
# echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# # Start counting the running time
# start_time=$SECONDS

# # Skip if this step already completed successfully (final archive present and valid)
# if [ -f "9_mmseqs_vfdb.tar.gz" ] && [ -f "9_mmseqs_vfdb.tar.gz.md5" ] && md5sum -c "9_mmseqs_vfdb.tar.gz.md5" >/dev/null 2>&1; then
#     echo "${workflow_step} already completed successfully (9_mmseqs_vfdb.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
# else
#     # A previous interrupted run may have left a partial directory or archive; wipe and redo,
#     # since this is a single call with no reliable partial resume.
#     if [ -d "9_mmseqs_vfdb" ]; then
#         echo "${workflow_step}: Found incomplete 9_mmseqs_vfdb directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
#         rm -rf "9_mmseqs_vfdb"
#     fi
#     rm -f "9_mmseqs_vfdb.tar.gz" "9_mmseqs_vfdb.tar.gz.md5"

#     # Create an output directory
#     mkdir -p 9_mmseqs_vfdb

#     # Activate Conda environment
#     conda activate blastp
#     # Run main software
#     vfdb.py \
#         --cpu $(nproc --ignore=1) \
#         -i 9_mmseqs/all_samples_f1_nr.faa \
#         --output all_samples_f1_nr \
#         --output_dir 9_mmseqs_vfdb
#     # Deactivate Conda environment
#     conda deactivate

#     # Compress the output directory
#     compressed_file="9_mmseqs_vfdb.tar.gz"
#     itens_to_compress=(9_mmseqs_vfdb)
#     echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
#     tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
#     # Generate checksum file of compressed directory file
#     md5sum "${compressed_file}" > "${compressed_file}".md5
#     # Check file integrity
#     echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
#     md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
# fi

# # Stop counting the running time
# elapsed_time=$((SECONDS - $start_time))
# running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
# # Show the running time
# echo "${workflow_step}: running time ${running_time}" | tee -a 0_workflow_progress.txt
# # Update the file 0_workflow_progress.txt
# echo "${workflow_step} step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 9.7) MMseqs2 -> Salmon index
# conda create -n salmon -c bioconda -c conda-forge salmon -y

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) Salmon index"
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "9_mmseqs_salmon_index.tar.gz" ] && [ -f "9_mmseqs_salmon_index.tar.gz.md5" ] && md5sum -c "9_mmseqs_salmon_index.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (9_mmseqs_salmon_index.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # A previous interrupted run may have left a partial index directory; wipe and redo,
    # since this is a single call with no reliable partial resume.
    if [ -d "9_mmseqs_salmon_index" ]; then
        echo "${workflow_step}: Found incomplete 9_mmseqs_salmon_index directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "9_mmseqs_salmon_index"
    fi
    rm -f "9_mmseqs_salmon_index.tar.gz" "9_mmseqs_salmon_index.tar.gz.md5"

    # Create output directory
    mkdir -p 9_mmseqs_salmon_index

    # Inform current step
    echo "${workflow_step} is creating the index of non-redundant gene catalog"

    # Activate Conda environment
    conda activate salmon
    # Run main software
    salmon index \
        -p $(nproc --ignore=1) \
        -t 9_mmseqs/all_samples_f1_nr.ffn \
        -i 9_mmseqs_salmon_index/nr_index
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="9_mmseqs_salmon_index.tar.gz"
    itens_to_compress=(9_mmseqs_salmon_index)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 9.8) MMseqs2 -> Salmon quant

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) Salmon quantification"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip the whole step if it already completed successfully (both final archives present and valid).
# Necessary because the script deletes per-sample directories under 9_mmseqs_salmon_quant/ at the end.
if [ -f "9_mmseqs_salmon_quant.tar.gz" ] && [ -f "9_mmseqs_salmon_quant.tar.gz.md5" ] && md5sum -c "9_mmseqs_salmon_quant.tar.gz.md5" >/dev/null 2>&1 \
    && [ -f "9_mmseqs_salmon_quant_merge.tar.gz" ] && [ -f "9_mmseqs_salmon_quant_merge.tar.gz.md5" ] && md5sum -c "9_mmseqs_salmon_quant_merge.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (both archives verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create an output directory
    mkdir -p 9_mmseqs_salmon_quant

    # Calculate sample size
    i=1
    r1_files=(5_bwa_reads/*_1.fq.gz)
    sample_count=${#r1_files[@]}

    # Activate Conda environment
    conda activate salmon
    # Loop through a list of sample files
    for r1 in "${r1_files[@]}"; do

        # Obtain r2 path
        r2=${r1/_1.fq.gz/_2.fq.gz}
        # Extract r1 file name
        filename=${r1##*/}
        # Extract sample name
        sample=${filename%%_*}

        # Skip sample if output already exists and looks complete (quant.sf non-empty via symlink,
        # plus lib_format_counts.json which Salmon writes only after a successful run)
        output_file="9_mmseqs_salmon_quant/${sample}/${sample}_salmon.tsv"
        completion_marker="9_mmseqs_salmon_quant/${sample}/lib_format_counts.json"
        if [ -s "$output_file" ] && [ -f "$completion_marker" ]; then
            echo "9) Salmon output file already exists and is valid for sample: $sample ($output_file). Skipping sample."
            i=$((i + 1))
            continue
        elif [ -d "9_mmseqs_salmon_quant/${sample}" ]; then
            echo "9) Salmon found incomplete output for sample: $sample. Removing partial directory and reprocessing."
            rm -rf "9_mmseqs_salmon_quant/${sample}"
        fi

        # Inform current sample
        echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Run main software
        salmon quant \
            -p $(nproc --ignore=1) \
            -i 9_mmseqs_salmon_index/nr_index \
            -l A \
            --meta \
            -1 "${r1}" \
            -2 "${r2}" \
            -o "9_mmseqs_salmon_quant/${sample}"

        # Rename count file
        ln -s "quant.sf" "$output_file"

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

        i=$((i + 1))

    done

    # Merge results | tee -a 0_workflow_progress.txt
    echo "Merging outputs" | tee -a 0_workflow_progress.txt
    out_dir="9_mmseqs_salmon_quant_merge"
    mkdir -p ${out_dir}
    sample_dirs=(9_mmseqs_salmon_quant/*/)
    salmon quantmerge --quants "${sample_dirs[@]}" --column tpm      -o "${out_dir}/salmon_tpm.tsv"
    salmon quantmerge --quants "${sample_dirs[@]}" --column numreads -o "${out_dir}/salmon_numreads.tsv"
    salmon quantmerge --quants "${sample_dirs[@]}" --column len      -o "${out_dir}/salmon_length.tsv"

    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="9_mmseqs_salmon_quant.tar.gz"
    itens_to_compress=(9_mmseqs_salmon_quant)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

    # Compress the output directory
    compressed_file="9_mmseqs_salmon_quant_merge.tar.gz"
    itens_to_compress=(9_mmseqs_salmon_quant_merge)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

    # Clean temporary files
    rm -rf 9_mmseqs_salmon_quant/*/
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt


############################################################
# 10) Binning - Single/Multi-sample - Input files
############################################################

# This session uses information from file 5_metagenomes.tsv
# The information is used to decide to perform single binning (1 sample per isolation source / host) or multi-binning (>= 2 samples per isolation source / host) 
# The default method of binning in this session uses self-supervised mode
# Using a GPU is highly recommended. It can reduce the running time from hours to minutes

############################################################
## 10.1) SeqKit

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) SeqKit"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 10_seqkit

# Calculate the sample count to display loop progress
i=1
fasta_files=(7_megahit/*.fasta)
sample_count=${#fasta_files[@]}

# Defining minimum contig size
msize=1000

# Activate Conda environment
conda activate seqkit
# Loop through a list of sample files
for file in "${fasta_files[@]}"; do

    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Skip sample if output already exists and is a valid, complete gzip file
    output_fasta="10_seqkit/${sample}_megahit_filtered.fasta"
    output_gz="${output_fasta}.gz"
    if [ -f "$output_gz" ] && gzip -t "$output_gz" 2>/dev/null; then
        echo "${workflow_step} output file already exists and is valid for sample: $sample. Skipping sample."
        i=$((i + 1))
        continue
    elif [ -f "$output_fasta" ] || [ -f "$output_gz" ]; then
        echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
        rm -f "$output_fasta" "$output_gz" "${output_gz}.md5"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Run SeqKit
    seqkit seq \
        --min-len $msize \
        "$file" \
        > "$output_fasta"

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress filtered assembly files (skip ones already compressed in a previous run)
echo "Compressing files"
for file in 10_seqkit/*.fasta; do
    [ -f "${file}.gz" ] && gzip -t "${file}.gz" 2>/dev/null && continue
    pigz -k -p $(nproc --ignore=1) "${file}"
done

# Generate checksum files (skip files already verified against an existing .md5)
(cd 10_seqkit && for file in *.fasta.gz; do
    [ -f "${file}.md5" ] && md5sum -c "${file}.md5" >/dev/null 2>&1 && continue
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done) | tee -a 0_workflow_progress.txt
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
(cd 10_seqkit && md5sum -c *.md5) | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 10.2) Seqkit -> SemiBin concatenate_fasta

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) Semibin concatenate_fasta"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 10_seqkit_concat

# Activate Conda environment
conda activate semibin
# Inform input table
input_file="5_metagenomes.tsv"
# awk to group samples by source
tr -d '\r' < "$input_file" | awk 'BEGIN {FS="\t"} 
{
    # Concatenate Sample ID ($1) for each Isolaton Source ($4)
    if (samples_by_source[$4] == "") {
        samples_by_source[$4] = $1;
    } else {
        samples_by_source[$4] = samples_by_source[$4] "," $1;
    }
}
END {
    # Print the grouped results
    for (source in samples_by_source) {
        print source "\t" samples_by_source[source];
    }
}' |
# Loop through each group
while IFS=$'\t' read -r source sample_list; do
    # Split the comma-separated sample list into a shell array
    IFS=',' read -r -a samples_array <<< "$sample_list"
    # Number of samples
    num_samples=${#samples_array[@]}
    # Check if number of samples is 2 or more
    if [ "$num_samples" -ge 2 ]; then

        # Skip group if output already exists and is a valid, complete gzip file
        output_gz="10_seqkit_concat/${source}_concat/concatenated.fa.gz"
        if [ -f "$output_gz" ] && gzip -t "$output_gz" 2>/dev/null; then
            echo "${workflow_step}: output already exists and is valid for source ${source}. Skipping."
            continue
        elif [ -d "10_seqkit_concat/${source}_concat" ]; then
            echo "${workflow_step}: found incomplete/corrupted output for source ${source}. Removing partial directory and reprocessing."
            rm -rf "10_seqkit_concat/${source}_concat"
        fi

        # Declare empty array
        input_files_array=()
        # Iterate over all samples and apply the full prefix/suffix
        for sample in "${samples_array[@]}"; do
            # Input file path: 10_seqkit/sample_megahit_filtered.fasta
            input_files_array+=("10_seqkit/${sample}_megahit_filtered.fasta.gz")
        done
        # Start counting the running time
        loop_start_time=$SECONDS
        # Inform source and execute the command line
        echo "▶ ${workflow_step} — ${source} (${num_samples} samples) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Run the script
        SemiBin2 concatenate_fasta \
            --input-fasta "${input_files_array[@]}" \
            --output "10_seqkit_concat/${source}_concat"
        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — $source (${num_samples} samples) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    else
        echo "${workflow_step}: ⚠️ Skipping concatenation for ${source}: only ${num_samples} sample(s) found (minimum 2 required)."
    fi

done
# Deactivate Conda environment
conda deactivate

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 10.3) Seqkit -> SemiBin concatenate_fasta -> Minimap2 index

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) Minimap2 index"
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Activate conda environment
conda activate minimap2
# Inform input table
input_file="5_metagenomes.tsv"
# awk to group samples by source
tr -d '\r' < "$input_file" | awk 'BEGIN {FS="\t"} 
{
    # Concatenate Sample ID ($1) for each Isolaton Source ($4)
    if (samples_by_source[$4] == "") {
        samples_by_source[$4] = $1;
    } else {
        samples_by_source[$4] = samples_by_source[$4] "," $1;
    }
}
END {
    # Print the grouped results
    for (source in samples_by_source) {
        print source "\t" samples_by_source[source];
    }
}' |
# Loop through each group
while IFS=$'\t' read -r source sample_list; do

    # Split the comma-separated sample list into a shell array
    IFS=',' read -r -a samples_array <<< "$sample_list"
    # Number of samples
    num_samples=${#samples_array[@]}

    # Skip source if output file already exists and looks complete (non-empty)
    output_file="10_seqkit_concat_minimap_index/${source}/${source}.mmi"
    if [ -s "$output_file" ]; then
        echo "${workflow_step}: file already exists and is valid for source: $source ($output_file). Skipping index."
        continue
    elif [ -f "$output_file" ] || [ -d "10_seqkit_concat_minimap_index/${source}" ]; then
        echo "${workflow_step}: found incomplete/corrupted output for source: $source. Removing partial files and reprocessing."
        rm -rf "10_seqkit_concat_minimap_index/${source}"
    fi

    # Check if number of samples is 2 or more
    if [ "$num_samples" -ge 2 ]; then

        # Create output directory
        mkdir -p 10_seqkit_concat_minimap_index/${source}

        # Start counting the running time
        loop_start_time=$SECONDS

        # Inform source and execute the command line
        echo "▶ ${workflow_step} — ${source} (${num_samples} samples) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        
        minimap2 \
            -k21 -w11 \
            -d "10_seqkit_concat_minimap_index/${source}/${source}.mmi" \
            "10_seqkit_concat/${source}_concat/concatenated.fa.gz"

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — $source (${num_samples} samples) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    elif [ "$num_samples" -eq 1 ]; then

        # Create output directory
        mkdir -p 10_seqkit_concat_minimap_index/${source}

        # Start counting the running time
        loop_start_time=$SECONDS

        # Inform source and execute the command line
        echo "▶ ${workflow_step} — ${source} (${num_samples} samples) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

        # Run main software
        minimap2 \
        -k21 -w11 \
        -d "10_seqkit_concat_minimap_index/${source}/${source}.mmi" \
        "10_seqkit/${sample_list}_megahit_filtered.fasta.gz"

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — $source (${num_samples} samples) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    else
        echo "${workflow_step}: ⚠️ Skipping index for ${source}: no sample found (minimum 1 required)."
    fi

done
# Deactivate Conda environment
conda deactivate

# Generate checksum files for the index files (skip files already verified against an existing .md5)
for file in 10_seqkit_concat_minimap_index/*/*.mmi; do
    [ -f "${file}.md5" ] && md5sum -c "${file}.md5" >/dev/null 2>&1 && continue
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done

# Compress the output directory
compressed_file="10_seqkit_concat_minimap_index.tar.gz"
itens_to_compress=(10_seqkit_concat_minimap_index)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 10.4) Seqkit -> SemiBin concatenate_fasta -> Minimap2

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) Minimap2"
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create output directory
mkdir -p 10_seqkit_concat_minimap_mapping

# Calculate sample size
i=1
sample_count=$(awk 'END {print NR}' 5_metagenomes.tsv)

# Activate conda environment
conda activate minimap2
# Loop through a list of file lines
tr -d '\r' < 5_metagenomes.tsv | awk '1'| \
while IFS=$'\t' read -r sample ref_accession ref_name isolation_source others; do

    # Skip sample if output file already exists and passes a BAM structural integrity check
    output_file="10_seqkit_concat_minimap_mapping/${sample}.mapped.sorted.bam"
    if [ -f "$output_file" ] && samtools quickcheck "$output_file" 2>/dev/null; then
        echo "${workflow_step} mapping file already exists and is valid for sample: $sample ($output_file). Skipping sample."
        i=$((i + 1))
        continue
    elif [ -f "$output_file" ]; then
        echo "${workflow_step} mapping found incomplete/corrupted BAM for sample: $sample. Removing partial files and reprocessing."
        rm -f "$output_file" \
            "10_seqkit_concat_minimap_mapping/${sample}_alignment.log" \
            "10_seqkit_concat_minimap_mapping/${sample}_allreads_flagstat.txt" \
            "10_seqkit_concat_minimap_mapping/${sample}_mappedreads_flagstat.txt" \
            "${output_file}.md5"
    fi

    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Map the reads to contigs and generate bam file (no intermediate files)
    minimap2 \
        -t $(nproc --ignore=1) \
        -ax sr \
        "10_seqkit_concat_minimap_index/${isolation_source}/${isolation_source}.mmi" \
        5_bwa_reads/"${sample}"[^0-9]*_1.fq.gz \
        5_bwa_reads/"${sample}"[^0-9]*_2.fq.gz \
        2> "10_seqkit_concat_minimap_mapping/${sample}_alignment.log" \
        | tee >(samtools flagstat - > "10_seqkit_concat_minimap_mapping/${sample}_allreads_flagstat.txt") \
        | samtools view -b -h -F 4 -@ 8 - \
        | tee >(samtools flagstat - > "10_seqkit_concat_minimap_mapping/${sample}_mappedreads_flagstat.txt") \
        | samtools sort -@ 8 \
            -o "$output_file"

    # Validate the resulting BAM file; if corrupted, remove so the next run reprocesses this sample
    if ! samtools quickcheck "$output_file" 2>>10_seqkit_concat_minimap_mapping_corrupted.txt; then
        echo "10) Minimap2 output failed integrity check for sample: $sample" | tee -a 0_workflow_progress.txt
        rm -f "$output_file"
    fi

    # # Map the reads to contigs and generate bam file (with intermediate files)
    # minimap2 \
    #     -t $(nproc --ignore=1) \
    #     -ax sr \
    #     "10_seqkit_concat_minimap_index/${isolation_source}/${isolation_source}.mmi" \
    #     5_bwa_reads/"${sample}"[^0-9]*_1.fq.gz \
    #     5_bwa_reads/"${sample}"[^0-9]*_2.fq.gz \
    #     > "10_seqkit_concat_minimap_mapping/${sample}.sam" \
    #     2> "10_seqkit_concat_minimap_mapping/${sample}_alignment.log"
    # # Generate mapping report
    # samtools flagstat \
    #     "10_seqkit_concat_minimap_mapping/${sample}.sam" \
    #     > "10_seqkit_concat_minimap_mapping/${sample}_allreads_flagstat.txt"
    # # Convert to bam and keep the header and only mapped reads
    # samtools view -b -h -F 4 -@ "$(nproc --ignore=1)" \
    #     "10_seqkit_concat_minimap_mapping/${sample}.sam" \
    #     > "10_seqkit_concat_minimap_mapping/${sample}.mapped.bam"
    # # Generate mapping report
    # samtools flagstat \
    #     "10_seqkit_concat_minimap_mapping/${sample}.mapped.bam" \
    #     > "10_seqkit_concat_minimap_mapping/${sample}_mappedreads_flagstat.txt"
    # # Delete intermediary file
    # rm "10_seqkit_concat_minimap_mapping/${sample}.sam"
    # # Sort reads by coordinate
    # samtools sort -@ "$(nproc --ignore=1)" \
    #     -o "10_seqkit_concat_minimap_mapping/${sample}.mapped.sorted.bam" \
    #     "10_seqkit_concat_minimap_mapping/${sample}.mapped.bam"
    # # Delete intermediary file
    # rm "10_seqkit_concat_minimap_mapping/${sample}.mapped.bam"
    # # Create the index bai file for the sorted bam file
    # samtools index \
    #     "10_seqkit_concat_minimap_mapping/${sample}.mapped.sorted.bam"

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Generate checksum files for the reads (skip files already verified against an existing .md5)
(cd 10_seqkit_concat_minimap_mapping && for file in *.bam; do
    [ -f "${file}.md5" ] && md5sum -c "${file}.md5" >/dev/null 2>&1 && continue
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done) | tee -a 0_workflow_progress.txt
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
(cd 10_seqkit_concat_minimap_mapping && md5sum -c *.md5) | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
# 11) Binning - Single/Multi-sample (Self-supervised mode)
############################################################

############################################################
## 11.1) Semibin

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="11) Semibin"
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create output directory
mkdir -p 11_semibin

# Activate conda environment
conda activate semibin
input_file="5_metagenomes.tsv" 
# awk to group samples by source
tr -d '\r' < "$input_file" | awk 'BEGIN {FS="\t"} 
{
    # Concatenate Sample ID ($1) for each Isolaton Source ($4)
    if (samples_by_source[$4] == "") {
        samples_by_source[$4] = $1;
    } else {
        samples_by_source[$4] = samples_by_source[$4] "," $1;
    }
}
END {
    # Print the grouped results
    for (source in samples_by_source) {
        print source "\t" samples_by_source[source];
    }
}' |
# Loop through each group
while IFS=$'\t' read -r source sample_list; do

    # Split the comma-separated sample list into a shell array
    IFS=',' read -r -a samples_array <<< "$sample_list"
    # Number of samples
    num_samples=${#samples_array[@]}

    # Check if number of samples is 2 or more
    if [ "$num_samples" -ge 2 ]; then

        # Skip source if already completed (final per-source log present and no leftover intermediate dir)
        if [ -f "11_semibin/${source}_SemiBinRun.log" ] && [ ! -d "11_semibin/${source}_semibin" ]; then
            echo "${workflow_step}: source ${source} already binned and organized. Skipping."
            continue
        elif [ -d "11_semibin/${source}_semibin" ] || [ -f "11_semibin/${source}_SemiBinRun.log" ] || [ -f "11_semibin/${source}.mapped.sorted.bam_0_data_cov.csv" ]; then
            echo "${workflow_step}: found incomplete output for source ${source}. Removing partial files and reprocessing."
            rm -rf "11_semibin/${source}_semibin"
            rm -f "11_semibin/${source}_SemiBinRun.log" "11_semibin/${source}.mapped.sorted.bam_0_data_cov.csv"
            # Also remove any per-sample dirs this source may have partially populated
            for sample in "${samples_array[@]}"; do
                rm -rf "11_semibin/${sample}_semibin"
            done
        fi

        # Declare empty array
        input_files_array=()
        # Iterate over all samples and apply the full prefix/suffix
        for sample in "${samples_array[@]}"; do

           input_files_array+=("10_seqkit_concat_minimap_mapping/${sample}.mapped.sorted.bam")

        done

        # Inform source and execute the command line
        echo "▶ ${workflow_step} — ${source} (${num_samples} samples) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


        # Start counting the running time
        loop_start_time=$SECONDS

        # Run the script
        SemiBin2 multi_easy_bin \
            --threads $(nproc --ignore=1) \
            --input-fasta "10_seqkit_concat/${source}_concat/concatenated.fa.gz" \
            --input-bam "${input_files_array[@]}" \
            --output "11_semibin/${source}_semibin"

        # Organize output files per sample
        for sample in "${samples_array[@]}"; do

            # Create sample directory
            mkdir -p 11_semibin/"${sample}"_semibin

            # Rename, extract and organize bin files
            # Rename and move bin files
            for file in 11_semibin/"${source}"_semibin/samples/"${sample}"_megahit_filtered/output_bins/*.fa.gz; do

                [ -e "$file" ] || continue
                filename=$(basename "$file" .fa.gz)
                newname="${filename/#SemiBin_/${sample}Bin}"
                mv "$file" 11_semibin/"${sample}"_semibin/"$newname".fasta.gz

            done

            # Move other bin files to sample directory
            mv 11_semibin/"${source}"_semibin/samples/"${sample}"_megahit_filtered/*.csv \
            11_semibin/"${source}"_semibin/samples/"${sample}"_megahit_filtered/*.pt \
            11_semibin/"${source}"_semibin/samples/"${sample}"_megahit_filtered/*.tsv \
            11_semibin/"${sample}"_semibin

        done

        # Move contacenated coverage file to 11_semibin/
        for concatcovfile in 11_semibin/"${source}"_semibin/samples/*.mapped.sorted.bam_0_data_cov.csv; do

            [ -e "$concatcovfile" ] || continue
            mv "$concatcovfile" 11_semibin/"${source}".mapped.sorted.bam_0_data_cov.csv

        done

        # Move log file to 11_semibin/
        mv 11_semibin/"${source}"_semibin/SemiBinRun.log 11_semibin/"${source}"_SemiBinRun.log

        # Delete intermediate directory
        rm -r "11_semibin/${source}_semibin"

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — $source (${num_samples} samples) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    elif [ "$num_samples" -eq 1 ]; then

        # Skip source if already completed (bin files exist and no leftover output_bins/ subdir)
        if [ -d "11_semibin/${sample_list}_semibin" ] && [ ! -d "11_semibin/${sample_list}_semibin/output_bins" ] \
            && ls 11_semibin/"${sample_list}"_semibin/*Bin*.fasta.gz >/dev/null 2>&1; then
            echo "${workflow_step}: source ${source} (single sample ${sample_list}) already binned and organized. Skipping."
            continue
        elif [ -d "11_semibin/${sample_list}_semibin" ]; then
            echo "${workflow_step}: found incomplete output for source ${source} (single sample ${sample_list}). Removing partial directory and reprocessing."
            rm -rf "11_semibin/${sample_list}_semibin"
        fi

        # Inform sample and source and execute the command line
        echo "${workflow_step} is binning the sample from ${source} (${num_samples} sample):"
        # Start counting the running time
        loop_start_time=$SECONDS

        # Use GPU to reduce the required time to train the models
        SemiBin2 single_easy_bin \
            --threads $(nproc --ignore=1) \
            --self-supervised \
            --input-fasta "10_seqkit/${sample_list}_megahit_filtered.fasta" \
            --input-bam "10_seqkit_concat_minimap_mapping/${sample_list}.mapped.sorted.bam" \
            --output "11_semibin/${sample_list}_semibin"

        # Rename, extract and organize bin files
        # Rename and move bin files
        for file in 11_semibin/"${sample_list}"_semibin/output_bins/SemiBin_*.fa.gz; do

            [ -e "$file" ] || continue
            filename=$(basename "$file" .fa.gz)
            newname="${filename/#SemiBin_/${sample_list}Bin}"
            mv "$file" 11_semibin/"${sample_list}"_semibin/"$newname".fasta.gz

        done

        # Delete intermediate directory
        rm -r 11_semibin/"${sample_list}"_semibin/output_bins

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time and udate the file 0_workflow_progress.txt
        echo "✔ ${workflow_step} — $source (${num_samples} samples) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    else

        echo "⚠️ Skipping concatenation for ${source}: no sample found (minimum 1 required)."

    fi

done

# Deactivate Conda environment
conda deactivate

# Create a file listing all bins
echo -e "Sample\tBin" > 11_semibin/semibin_all.tsv
# Files loop 
for filepath in 11_semibin/*_semibin/*.fasta.gz; do
    dir_name=$(basename "$(dirname "$filepath")")
    sample="${dir_name%_semibin}"
    filename=$(basename "$filepath" .fasta.gz)
    echo -e "${sample}\t${filename}" >> 11_semibin/semibin_all.tsv
done

# Archive the output directory
compressed_file="11_semibin.tar"
itens_to_compress=(11_semibin)
echo "Archiving output directory" | tee -a 0_workflow_progress.txt
tar -cvf "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt


############################################################
# 12) Bin quality control
############################################################

############################################################
## 12.1) CheckM2

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="12) CheckM2"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip the whole step if it already completed successfully (final archive present and valid).
# Necessary because the script deletes 12_checkm2/ at the end.
if [ -f "12_checkm2.tar.gz" ] && [ -f "12_checkm2.tar.gz.md5" ] && md5sum -c "12_checkm2.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (12_checkm2.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create an output directory
    mkdir -p 12_checkm2

    # Calculate the sample count to display loop progress
    i=1
    sample_dirs=(11_semibin/*_semibin/)
    sample_count=${#sample_dirs[@]}

    # Activate Conda environment
    conda activate checkm2
    # Loop through a list of sample directories
    for dir in "${sample_dirs[@]}"; do

        # Extract directory name
        dirname=${dir#*/}
        # Extract sample name
        sample=${dirname%%_semibin/}

        # Skip sample if output file already exists and looks complete (header + at least one data line)
        output_file="12_checkm2/${sample}_checkm2.tsv"
        if [ -f "$output_file" ] && [ "$(wc -l < "$output_file")" -ge 2 ]; then
            echo "${workflow_step} output file already exists and is valid for sample: $sample ($output_file). Skipping sample."
            i=$((i + 1))
            continue
        elif [ -f "$output_file" ] || [ -d "12_checkm2/${sample}_checkm2" ]; then
            echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
            rm -f "$output_file"
            rm -rf "12_checkm2/${sample}_checkm2"
        fi

        # Inform current sample
        echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Run main software
        checkm2 predict \
            --threads $(nproc --ignore=1) \
            -x fasta.gz \
            --input "11_semibin/${sample}_semibin" \
            --output-directory "12_checkm2/${sample}_checkm2"

        # Copy and rename the output file
        cp "12_checkm2/${sample}_checkm2/quality_report.tsv" "12_checkm2/${sample}_checkm2.tsv"

        # Delete the samples output directories
        rm -r "12_checkm2/${sample}_checkm2"

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

        i=$((i + 1))

    done
    # Deactivate Conda environment
    conda deactivate

    # Concatenate results
    shopt -s nullglob
    files=(12_checkm2/*_checkm2.tsv)
    if [ ${#files[@]} -gt 0 ]; then
        first_file="${files[0]}"
        head -n 1 "$first_file" > 12_checkm2/checkm2_all.tsv
        for f in "${files[@]}"; do
            tail -n +2 "$f" >> 12_checkm2/checkm2_all.tsv
        done
    fi

    # Copy the concatenated results file to the main directory
    cp 12_checkm2/checkm2_all.tsv 12_checkm2.tsv

    # Compress the output directory
    compressed_file="12_checkm2.tar.gz"
    itens_to_compress=(12_checkm2 12_checkm2.tsv)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
    # Delete the output directory
    rm -r 12_checkm2
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 12.2) GUNC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="12) GUNC"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip the whole step if it already completed successfully (final archive present and valid).
if [ -f "12_gunc.tar.gz" ] && [ -f "12_gunc.tar.gz.md5" ] && md5sum -c "12_gunc.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (12_gunc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create an output directory
    mkdir -p 12_gunc

    # Calculate the sample count to display loop progress
    i=1
    sample_dirs=(11_semibin/*_semibin/)
    sample_count=${#sample_dirs[@]}

    # Activate Conda environment
    conda activate gunc
    # Loop through a list of sample directories
    for dir in "${sample_dirs[@]}"; do

        # Extract directory name
        dirname=${dir#*/}
        # Extract sample name
        sample=${dirname%%_semibin/}

        # Skip sample if output file already exists and looks complete (header + at least one data line)
        output_file="12_gunc/${sample}_gunc.tsv"
        if [ -f "$output_file" ] && [ "$(wc -l < "$output_file")" -ge 2 ]; then
            echo "${workflow_step} output file already exists and is valid for sample: $sample ($output_file). Skipping sample."
            i=$((i + 1))
            continue
        elif [ -f "$output_file" ] || [ -d "12_gunc/${sample}_gunc" ] || [ -d "12_gunc/${sample}_gunc_temp" ]; then
            echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
            rm -f "$output_file"
            rm -rf "12_gunc/${sample}_gunc" "12_gunc/${sample}_gunc_temp"
        fi

        # Inform current sample
        echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Create an output directory
        mkdir -p "12_gunc/${sample}_gunc" "12_gunc/${sample}_gunc_temp"

        # Run main software
        gunc run \
            --threads $(nproc --ignore=1) \
            --contig_taxonomy_output \
            --file_suffix .fasta.gz \
            --input_dir "11_semibin/${sample}_semibin" \
            --temp_dir "12_gunc/${sample}_gunc_temp" \
            --out_dir "12_gunc/${sample}_gunc"

        # Copy and rename the output file
        cp 12_gunc/${sample}_gunc/*maxCSS_level.tsv "12_gunc/${sample}_gunc.tsv"

        # # Plotting the data
        # # Create an output directory
        # mkdir -p 12_gunc/${sample}_gunc_plot
        # # Loop through a list of files
        # j=1
        # bin_count=$(ls -1 12_gunc/${sample}_gunc/diamond_output/*.out | wc -l)
        # for plotfile in 12_gunc/${sample}_gunc/diamond_output/*.out; do\
        #    plotfilename=${plotfile##*/}
        #    plotfilename=${plotfilename%%.diamond*}
        #    echo "GUNC is ploting bin ${plotfilename} from sample: ${sample} (${j}/${bin_count})"
        #    gunc plot \
        #    --contig_display_num 0 \
        #    --diamond_file $plotfile \
        #    --out_dir 12_gunc/${sample}_gunc_plot/;\
        #    j=$((j + 1))
        # done

        # Delete the intermediary directories
        rm -r 12_gunc/${sample}_gunc 12_gunc/${sample}_gunc_temp

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

        i=$((i + 1))

    done
    # Deactivate Conda environment
    conda deactivate

    # Concatenate results
    shopt -s nullglob
    files=(12_gunc/*_gunc.tsv)
    if [ ${#files[@]} -gt 0 ]; then
        first_file="${files[0]}"
        head -n 1 "$first_file" > 12_gunc/gunc_all.tsv
        for f in "${files[@]}"; do
            tail -n +2 "$f" >> 12_gunc/gunc_all.tsv
        done
    fi

    # Copy the concatenated results file to the main directory
    cp 12_gunc/gunc_all.tsv 12_gunc.tsv

    # Compress the output directory
    compressed_file="12_gunc.tar.gz"
    itens_to_compress=(12_gunc 12_gunc.tsv)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

    # Delete the output directory
    rm -r 12_gunc
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 12.3) GTDB-Tk

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="12) GTDB-Tk"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Check if this step already completed successfully (final archive present and valid)
if [ -f "12_gtdbtk.tar.gz" ] && [ -f "12_gtdbtk.tar.gz.md5" ] && md5sum -c "12_gtdbtk.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (12_gtdbtk.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # If 12_gtdbtk exists without a valid final archive, a previous run was interrupted.
    # GTDB-Tk has no per-genome resume, so the only safe option is to wipe and restart.
    if [ -d "12_gtdbtk" ]; then
        echo "${workflow_step}: Found incomplete 12_gtdbtk directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "12_gtdbtk"
    fi
    # Also remove any leftover partial archive/checksum from an interrupted compression step
    rm -f "12_gtdbtk.tar.gz" "12_gtdbtk.tar.gz.md5"

    # Create an output directory
    mkdir -p 12_gtdbtk

    # Batchfile listing every bin passing QC (path<TAB>genome_id), plus a
    # mapping file (genome_id<TAB>sample) used later to split results back per sample.
    batchfile="12_gtdbtk/all_bins_batchfile.tsv"
    mapfile_tsv="12_gtdbtk/genome_sample_map.tsv"
    passed_bins_list="12_gtdbtk/passed_qc_bins.txt"

    # Filter genomes based on CheckM2 quality thresholds:
    # Completeness >= 70%, Contamination <= 5%, Contig_N50 >= 5000 bp
    echo "${workflow_step}: Filtering genomes from 12_checkm2.tsv (Completeness >= 70%, Contamination <= 5%, N50 >= 5000 bp)" | tee -a 0_workflow_progress.txt
    awk -F'\t' 'NR>1 { if ($2 >= 70 && $3 <= 5 && $7 >= 5000) print $1 }' 12_checkm2.tsv > "$passed_bins_list"

    n_passed=$(wc -l < "$passed_bins_list")
    echo "${workflow_step}: ${n_passed} bins passed quality control filtering." | tee -a 0_workflow_progress.txt

    echo "${workflow_step}: Building batchfile from 11_semibin/*_semibin/"
    > "$batchfile"
    > "$mapfile_tsv"

    sample_dirs=(11_semibin/*_semibin/)
    for dir in "${sample_dirs[@]}"; do
        # Extract directory name
        dirname=${dir#*/}
        # Extract sample name
        sample=${dirname%%_semibin/}

        for bin in "${dir}"*.fasta.gz; do
            [ -f "$bin" ] || continue
            bin_filename=$(basename "$bin")
            bin_basename_raw="${bin_filename%.gz}"

            # Check if the current bin passed the CheckM2 quality filter (matching with or without .gz extension)
            if grep -qxF "$bin_filename" "$passed_bins_list" || grep -qxF "$bin_basename_raw" "$passed_bins_list"; then
                binname=$(basename "$bin" .fasta.gz)
                genome_id="${sample}_${binname}"
                printf '%s\t%s\n' "$(readlink -f "$bin")" "$genome_id" >> "$batchfile"
                printf '%s\t%s\n' "$genome_id" "$sample" >> "$mapfile_tsv"
            fi
        done
    done

    n_genomes=$(wc -l < "$batchfile")
    echo "${workflow_step} batchfile has ${n_genomes} quality-filtered genomes across $(ls -d 11_semibin/*_semibin/ | wc -l) samples"

    # Scratch directory
    scratch_dir="/tmp/gtdbtk_scratch_${USER}_$$"
    # Create scratch directory
    mkdir -p "$scratch_dir"

    # Activate Conda environment
    conda activate gtdbtk
    # Run main software
    gtdbtk classify_wf \
        --cpus $(nproc --ignore=1) \
        --pplacer_cpus 12 \
        --scratch_dir "$scratch_dir" \
        --batchfile "$batchfile" \
        --out_dir "12_gtdbtk/all_samples_gtdbtk"
    # Deactivate Conda environment
    conda deactivate

    # Delete scratch directory
    rm -rf "$scratch_dir"

    # Copy the combined summary files
    for map in "classify/gtdbtk.bac120.summary.tsv:gtdbtk_bacteria_all.tsv" \
               "classify/gtdbtk.ar53.summary.tsv:gtdbtk_archaea_all.tsv" \
               "gtdbtk.log:all_samples_gtdbtk.log" ; do
        src="12_gtdbtk/all_samples_gtdbtk/${map%%:*}"
        dst="12_gtdbtk/${map#*:}"
        [ -f "$src" ] && cp "$src" "$dst"
    done

    # Delete the raw gtdbtk output directory now that the summaries are copied
    rm -r "12_gtdbtk/all_samples_gtdbtk"

    # Combined bacteria+archaea file
    if [ -f "12_gtdbtk/gtdbtk_bacteria_all.tsv" ] || [ -f "12_gtdbtk/gtdbtk_archaea_all.tsv" ]; then
        first_file=""
        for f in "12_gtdbtk/gtdbtk_bacteria_all.tsv" "12_gtdbtk/gtdbtk_archaea_all.tsv"; do
            [ -f "$f" ] && first_file="$f" && break
        done
        head -n 1 "$first_file" > "12_gtdbtk/gtdbtk_archaea_bacteria_all.tsv"
        for f in "12_gtdbtk/gtdbtk_bacteria_all.tsv" "12_gtdbtk/gtdbtk_archaea_all.tsv"; do
            [ -f "$f" ] || continue
            tail -n +2 "$f" >> "12_gtdbtk/gtdbtk_archaea_bacteria_all.tsv"
        done
    fi

    # Split the combined summaries back into one file per sample
    for domain in bacteria archaea; do
        src="12_gtdbtk/gtdbtk_${domain}_all.tsv"
        [ -f "$src" ] || continue
        awk -F'\t' -v domain="$domain" '
            NR==FNR { sample[$1] = $2; next }
            FNR==1  { header = $0; next }
            {
                s = sample[$1]
                if (s != "") {
                    out = "12_gtdbtk/" s "_gtdbtk_" domain ".tsv"
                    if (!(out in seen)) { print header > out; seen[out] = 1 }
                    print $0 > out
                }
            }
        ' "$mapfile_tsv" "$src"
    done

    # Compress the output directory
    compressed_file="12_gtdbtk.tar.gz"
    itens_to_compress=(12_gtdbtk)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
echo "${workflow_step} running time: ${running_time}" | tee -a 0_workflow_progress.txt
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 12.4) QUAST

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="12) QUAST"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip the whole step if it already completed successfully
if [ -f "12_quast.tar.gz" ] && [ -f "12_quast.tar.gz.md5" ] && md5sum -c "12_quast.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (12_quast.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create an output directory
    mkdir -p 12_quast

    # Calculate the sample count to display loop progress
    i=1
    sample_dirs=(11_semibin/*_semibin/)
    sample_count=${#sample_dirs[@]}

    # Activate Conda environment
    conda activate quast
    # Loop through a list of sample directories
    for dir in "${sample_dirs[@]}"; do
        # Extract directory name
        dirname=${dir#*/}
        # Extract sample name
        sample=${dirname%%_semibin/}

        # Skip sample if output file already exists and looks complete (header + at least one data line)
        output_file="12_quast/${sample}_quast.tsv"
        if [ -f "$output_file" ] && [ "$(wc -l < "$output_file")" -ge 2 ]; then
            echo "${workflow_step} output file already exists and is valid for sample: $sample ($output_file). Skipping sample."
            i=$((i + 1))
            continue
        elif [ -f "$output_file" ] || [ -d "12_quast/${sample}_quast" ]; then
            echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
            rm -f "$output_file"
            rm -rf "12_quast/${sample}_quast"
        fi

        # Inform current sample
        echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS
        # Create an output directory
        mkdir -p "12_quast/${sample}_quast"
        # Run main software
        quast.py -t $(nproc --ignore=1) -m 0 -o \
            "12_quast/${sample}_quast" \
            11_semibin/${sample}_semibin/*.fasta.gz
        # Copy transposed report
        cp "12_quast/${sample}_quast/transposed_report.tsv" "12_quast/${sample}_quast.tsv"
        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        i=$((i + 1))
    done
    # Deactivate Conda environment
    conda deactivate

    # Concatenate results
    shopt -s nullglob
    files=(12_quast/*_quast.tsv)
    if [ ${#files[@]} -gt 0 ]; then
        first_file="${files[0]}"
        head -n 1 "$first_file" > 12_quast/quast_all.tsv
        for f in "${files[@]}"; do
            tail -n +2 "$f" >> 12_quast/quast_all.tsv
        done
    fi
    # Compress the output directory
    compressed_file="12_quast.tar.gz"
    itens_to_compress=(12_quast)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
    # Delete the output directory
    rm -r 12_quast
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt


############################################################
# 13) Bin functional abundance profile
############################################################

############################################################
## 13.1) Aragorn

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="13) Aragorn"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "13_aragorn.tar.gz" ] && [ -f "13_aragorn.tar.gz.md5" ] && md5sum -c "13_aragorn.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (13_aragorn.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create an output directory
    mkdir -p 13_aragorn

    # Calculate the sample count to display loop progress
    i=1
    sample_dirs=(11_semibin/*_semibin/)
    sample_count=${#sample_dirs[@]}

    # Activate Conda environment
    conda activate aragorn
    # Loop through a list of sample directories
    for dir in "${sample_dirs[@]}"; do

        # Extract directory name
        dirname=${dir#*/}
        # Extract sample name
        sample=${dirname%%_semibin/}

        # Inform current sample
        echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Create an output directory
        mkdir -p "13_aragorn/${sample}_aragorn"

        # Loop through a list of files
        for file in 11_semibin/${sample}_semibin/*.fasta.gz; do

            # Extract sample name
            prefix=$(basename ${file} .fasta.gz)

            # Skip bin if output already exists and looks complete (non-empty)
            bin_output="13_aragorn/${sample}_aragorn/${prefix}_aragorn.txt"
            if [ -s "$bin_output" ]; then
                echo "${workflow_step} output already exists and is valid for bin: $prefix. Skipping bin."
                continue
            fi

            # Run software
            aragorn \
                "$file" \
                > "$bin_output"

        done
        
        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

        i=$((i + 1))

    done
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="13_aragorn.tar.gz"
    itens_to_compress=(13_aragorn)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

    # Delete the output directory
    rm -r 13_aragorn
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 13.2) Pybarrnap

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="13) Pybarrnap"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "13_pybarrnap.tar.gz" ] && [ -f "13_pybarrnap.tar.gz.md5" ] && md5sum -c "13_pybarrnap.tar.gz.md5" >/dev/null 2>&1; then
    echo "${workflow_step} already completed successfully (13_pybarrnap.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create an output directory
    mkdir -p 13_pybarrnap

    # Calculate the sample count to display loop progress
    i=1
    sample_dirs=(11_semibin/*_semibin/)
    sample_count=${#sample_dirs[@]}

    # Activate Conda environment
    conda activate pybarrnap
    # Loop through a list of sample directories
    for dir in "${sample_dirs[@]}"; do

        # Extract directory name
        dirname=${dir#*/}
        # Extract sample name
        sample=${dirname%%_semibin/}

        # Inform current sample
        echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Create an output directory
        mkdir -p "13_pybarrnap/${sample}_pybarrnap"

        # Loop through a list of files
        for file in 11_semibin/${sample}_semibin/*.fasta.gz; do

            # Extract sample name
            prefix=$(basename ${file} .fasta.gz)

            # Skip bin if both bac and arc outputs already exist and look complete
            bac_output="13_pybarrnap/${sample}_pybarrnap/${prefix}_bac_pybarrnap.fasta"
            arc_output="13_pybarrnap/${sample}_pybarrnap/${prefix}_arc_pybarrnap.fasta"
            if [ -s "$bac_output" ] && [ -s "$arc_output" ]; then
                echo "${workflow_step} outputs already exist and are valid for bin: $prefix. Skipping bin."
                continue
            fi

            # Run barrnap for bacteria
            zcat $file | pybarrnap \
                --threads $(nproc --ignore=1) \
                --quiet \
                --kingdom bac \
                 -o "$bac_output" \
                > "13_pybarrnap/${sample}_pybarrnap/${prefix}_bac_pybarrnap.gff" \
                2> "13_pybarrnap/${sample}_pybarrnap/${prefix}_bac_pybarrnap.log"

            # Run the program for archaea
            zcat $file | pybarrnap \
                --threads $(nproc --ignore=1) \
                --quiet \
                --kingdom arc \
                -o "$arc_output" \
                > "13_pybarrnap/${sample}_pybarrnap/${prefix}_arc_pybarrnap.gff" \
                2> "13_pybarrnap/${sample}_pybarrnap/${prefix}_arc_pybarrnap.log"

        done
        
        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

        i=$((i + 1))

    done
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    compressed_file="13_pybarrnap.tar.gz"
    itens_to_compress=(13_pybarrnap)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

    # Delete the output directory
    rm -r 13_pybarrnap
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 13.3) Pyrodigal

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="13) Pyrodigal"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 13_pyrodigal

# Calculate the sample count to display loop progress
i=1
sample_dirs=(11_semibin/*_semibin/)
sample_count=${#sample_dirs[@]}

# Activate Conda environment
conda activate pyrodigal
# Loop through a list of sample directories
for dir in "${sample_dirs[@]}"; do

    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_semibin/}
    
    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create an output directory
    mkdir -p "13_pyrodigal/${sample}_pyrodigal"

    # Loop through a list of files
    for file in 11_semibin/${sample}_semibin/*.fasta.gz; do

        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.fasta.gz}

        # Skip bin if all output files already exist and look complete (non-empty)
        bin_dir="13_pyrodigal/${sample}_pyrodigal/${binname}"
        if [ -s "${bin_dir}/${binname}.fsa" ] && \
           [ -s "${bin_dir}/${binname}.ffn" ] && \
           [ -s "${bin_dir}/${binname}.faa" ] && \
           [ -s "${bin_dir}/${binname}.gff" ]; then
           echo "${workflow_step} output files (.fsa, .ffn, .faa, .gff) already exist and are valid for sample $sample -> bin $binname. Skipping bin." | tee -a 0_workflow_progress.txt
           continue
        elif [ -d "$bin_dir" ]; then
           echo "${workflow_step} found incomplete output for sample $sample -> bin $binname. Removing partial directory and reprocessing." | tee -a 0_workflow_progress.txt
           rm -rf "$bin_dir"
        fi

        # Create output file
        mkdir -p "13_pyrodigal/${sample}_pyrodigal/${binname}"

        # Extract input file
        zcat "${file}" > "13_pyrodigal/${sample}_pyrodigal/${binname}/${binname}.fsa"

        # Run Pyrodigal
        pyrodigal \
            -j $(nproc --ignore=1) \
            -m \
            -p single \
            --no-stop-codon \
            -f gff \
            -i "13_pyrodigal/${sample}_pyrodigal/${binname}/${binname}.fsa" \
            -d "13_pyrodigal/${sample}_pyrodigal/${binname}/${binname}.ffn" \
            -a "13_pyrodigal/${sample}_pyrodigal/${binname}/${binname}.faa" \
            -o "13_pyrodigal/${sample}_pyrodigal/${binname}/${binname}.gff"
    
    done

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Gene count per bin
echo -e "Sample\tBin\tGene_count" > 13_pyrodigal/pyrodigal_all.tsv
# Files loop
shopt -s nullglob
for filepath in 13_pyrodigal/*_pyrodigal/*/*.faa; do
    dir_sample=$(basename "$(dirname "$(dirname "$filepath")")")
    sample="${dir_sample%_pyrodigal}"
    filename=$(basename "$filepath" .faa)
    gene_count=$(grep -c ">" "$filepath")
    echo -e "${sample}\t${filename}\t${gene_count}" >> 13_pyrodigal/pyrodigal_all.tsv
done

# Compress the output directory
compressed_file="13_pyrodigal.tar.gz"
itens_to_compress=(13_pyrodigal)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 13.4) Pyrodigal -> AMRFinderPlus

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="13) AMRFinderPlus"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 13_pyrodigal_amrfinder

# Calculate the sample count to display loop progress
i=1
sample_dirs=(13_pyrodigal/*_pyrodigal/)
sample_count=${#sample_dirs[@]}

# Activate Conda environment
conda activate amrfinder
# Loop through a list of sample directories
for dir in "${sample_dirs[@]}"; do

    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_pyrodigal/}
    
    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create an output directory
    mkdir -p "13_pyrodigal_amrfinder/${sample}_amrfinder"

    # Loop through a list of files
    for file in ${dir}/*/*.faa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.faa}

        # Skip bin if output already exists and looks complete (non-empty, no leftover format file)
        output_tsv="13_pyrodigal_amrfinder/${sample}_amrfinder/${binname}_amrfinder.tsv"
        format_faa="13_pyrodigal_amrfinder/${sample}_amrfinder/${binname}_amrfinder_format.faa"
        nuc_out="13_pyrodigal_amrfinder/${sample}_amrfinder/${binname}_amrfinder.fasta"
        prot_out="13_pyrodigal_amrfinder/${sample}_amrfinder/${binname}_amrfinder.faa"
        if [ -s "$output_tsv" ] && [ ! -f "$format_faa" ]; then
            echo "${workflow_step} output already exists and is valid for bin: $binname from sample $sample. Skipping bin."
            continue
        elif [ -f "$output_tsv" ] || [ -f "$format_faa" ] || [ -f "$nuc_out" ] || [ -f "$prot_out" ]; then
            echo "${workflow_step} found incomplete/corrupted output for bin: $binname from sample $sample. Removing partial files and reprocessing."
            rm -f "$output_tsv" "$format_faa" "$nuc_out" "$prot_out"
        fi

        # Format faa files for AMRFinderPlus
        awk '{
          if ($0 ~ /^>/) {
            match($0, /^>([^ ]+)/, a)
            id=a[1]
            sub(/ID=[^;]+/, "ID="id)
            print
          } else print
        }' "${file}" > "$format_faa"

        #Inform bin
        echo "${workflow_step} is processing bin: ${binname} from sample ${sample}"
        # Run main software
        amrfinder \
            --threads $(nproc --ignore=1) \
            --database /db/amrfinder/latest \
            --plus \
            --annotation_format prodigal \
            --nucleotide "13_pyrodigal/${sample}_pyrodigal/${binname}/${binname}.fsa" \
            --protein "$format_faa" \
            --gff "13_pyrodigal/${sample}_pyrodigal/${binname}/${binname}.gff" \
            --name "${binname}" \
            --nucleotide_output "$nuc_out" \
            --protein_output "$prot_out" \
            --output "$output_tsv"
        
        # Remove intermediate file
        rm "$format_faa"
    done

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="13_pyrodigal_amrfinder.tar.gz"
itens_to_compress=(13_pyrodigal_amrfinder)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 13.5) Pyrodigal -> dbCAN

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="13) dbCAN"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 13_pyrodigal_dbcan

# Calculate the sample count to display loop progress
i=1
sample_dirs=(13_pyrodigal/*_pyrodigal/)
sample_count=${#sample_dirs[@]}

# Activate Conda environment
conda activate dbcan
# Loop through a list of sample directories
for dir in "${sample_dirs[@]}"; do

    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_pyrodigal/}
    
    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create an output directory
    mkdir -p "13_pyrodigal_dbcan/${sample}_dbcan"

    # Loop through a list of files
    for file in ${dir}/*/*.faa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.faa}

        # Skip bin if output already exists and looks complete (non-empty)
        output_file="13_pyrodigal_dbcan/${sample}_dbcan/${binname}/overview.tsv"
        if [ -s "$output_file" ]; then
            echo "${workflow_step} output already exists and is valid for bin: $binname from sample $sample. Skipping bin."
            continue
        elif [ -d "13_pyrodigal_dbcan/${sample}_dbcan/${binname}" ]; then
            echo "${workflow_step} found incomplete output for bin: $binname from sample $sample. Removing partial directory and reprocessing."
            rm -rf "13_pyrodigal_dbcan/${sample}_dbcan/${binname}"
        fi

        # Create an output directory
        mkdir -p "13_pyrodigal_dbcan/${sample}_dbcan/${binname}"

        #Inform bin
        echo "13) dbCAN is processing bin: ${binname} from sample ${sample}"
        # Run main software
        run_dbcan CAZyme_annotation \
            --threads $(nproc --ignore=1) \
            --db_dir /db/dbcan \
            --mode protein \
            --input_raw_data $file \
            --output_dir "13_pyrodigal_dbcan/${sample}_dbcan/${binname}"

        # Delete uniInput.faa
        rm -f "13_pyrodigal_dbcan/${sample}_dbcan/${binname}/uniInput.faa"
        # Delete temporary directory
        rm -rf "13_pyrodigal_dbcan/${sample}_dbcan/${binname}/tmp"
    done

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="13_pyrodigal_dbcan.tar.gz"
itens_to_compress=(13_pyrodigal_dbcan)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 13.6) Pyrodigal -> eggNOG-mapper

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="13) eggNOG-mapper"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 13_pyrodigal_emapper

# Calculate the sample count to display loop progress
i=1
sample_dirs=(13_pyrodigal/*_pyrodigal/)
sample_count=${#sample_dirs[@]}

# Activate Conda environment
conda activate eggnog-mapper
# Loop through a list of sample directories
for dir in "${sample_dirs[@]}"; do

    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_pyrodigal/}
    
    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    
    # Start counting the running time
    loop_start_time=$SECONDS
    # Loop through a list of files
    for file in ${dir}/*/*.faa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.faa}

        # Skip bin if output already exists and looks complete (tsv present, everything compressed, nothing left uncompressed)
        bin_dir="13_pyrodigal_emapper/${sample}_emapper/${binname}_emapper"
        output_tsv="${bin_dir}/${binname}_emapper.tsv"
        annotations_gz="${bin_dir}/${binname}.emapper.annotations.gz"
        hits_gz="${bin_dir}/${binname}.emapper.hits.gz"
        seed_gz="${bin_dir}/${binname}.emapper.seed_orthologs.gz"
        if [ -s "$output_tsv" ] && [ -f "$annotations_gz" ] && [ -f "$hits_gz" ] && [ -f "$seed_gz" ] \
            && [ ! -f "${bin_dir}/${binname}.emapper.annotations" ]; then
            echo "${workflow_step} output already exists and is valid for bin: $binname from sample $sample. Skipping bin."
            continue
        elif [ -d "$bin_dir" ]; then
            echo "${workflow_step} found incomplete output for bin: $binname from sample $sample. Removing partial directory and reprocessing."
            rm -rf "$bin_dir"
        fi

        # Create an output directory
        mkdir -p "$bin_dir"

        #Inform bin
        echo "${workflow_step} is processing bin: ${binname} from sample ${sample}"
        # Run main software
        emapper.py \
            --cpu $(nproc --ignore=1) \
            -i $file \
            --output "${binname}" \
            --output_dir "$bin_dir"

        # Adjust output table
        grep -v '^##' "${bin_dir}/${binname}.emapper.annotations" \
        > "$output_tsv"

        # Compress original annotation files
        pigz -f "${bin_dir}/${binname}.emapper."{annotations,hits,seed_orthologs}
    
    done
    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt   
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="13_pyrodigal_emapper.tar.gz"
itens_to_compress=(13_pyrodigal_emapper)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 13.7) Pyrodigal -> VFDB (BLASTP)

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="13) VFDB (BLASTP)"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 13_pyrodigal_vfdb

# Calculate the sample count to display loop progress
i=1
sample_dirs=(13_pyrodigal/*_pyrodigal/)
sample_count=${#sample_dirs[@]}

# Activate Conda environment
conda activate blast
# Loop through a list of sample directories
for dir in "${sample_dirs[@]}"; do

    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_pyrodigal/}
    
    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    
    # Start counting the running time
    loop_start_time=$SECONDS
    # Loop through a list of files
    for file in ${dir}/*/*.faa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.faa}

        # Create an output directory
        mkdir -p "$bin_dir"

        #Inform bin
        echo "${workflow_step} is processing bin: ${binname} from sample ${sample}"
        # Run main software
        vfdb.py \
            --cpu $(nproc --ignore=1) \
            -i $file \
            --output "${binname}" \
            --output_dir "$bin_dir"
   
    done
    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt   
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="13_pyrodigal_vfdb.tar.gz"
itens_to_compress=(13_pyrodigal_vfdb)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt


############################################################
# 14) Bin mobile genetic elements
############################################################

############################################################
## 14.1) MOB-suite

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="14) MOB-suite"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 14_mobsuite

# Calculate the sample count to display loop progress
i=1
sample_dirs=(13_pyrodigal/*_pyrodigal/)
sample_count=${#sample_dirs[@]}

# Activate Conda environment
conda activate mob_suite
# Loop through a list of sample directories
for dir in "${sample_dirs[@]}"; do

    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_pyrodigal/}
    
    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create an output directory
    mkdir -p "14_mobsuite/${sample}_mobsuite"

    # Loop through a list of files
    for file in "${dir}"/*/*.fsa; do

        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.fsa}

        # Skip bin if output already exists and looks complete (contig_report.txt present)
        outdir="14_mobsuite/${sample}_mobsuite/${binname}_mobsuite"
        if [ -s "${outdir}/contig_report.txt" ]; then
            echo "${workflow_step} output already exists and is valid for bin: $binname from sample $sample. Skipping bin."
            continue
        elif [ -d "$outdir" ]; then
            echo "${workflow_step} found incomplete output for bin: $binname from sample $sample. Removing partial directory and reprocessing."
            rm -rf "$outdir"
        fi

        #Inform bin
        echo "${workflow_step} is processing bin: ${binname} from sample ${sample}"
        # Run main software
        mob_recon -n $(nproc --ignore=1) \
            --infile "${file}" \
            --outdir "$outdir"

        # Rename MOB-suite fasta files
        (
            cd "$outdir" || exit
            [ -f chromosome.fasta ] && mv chromosome.fasta "${binname}"_chromosome.fasta
            for p in plasmid*; do
                [ -e "$p" ] || continue
                mv "$p" "${binname}_$p"
            done
        )
    done

    # Merge contig_report.txt files
    # Delete merged file if it exists
    if [ -f 14_mobsuite/${sample}_mobsuite/contig_report_all.tsv ]; then
        rm 14_mobsuite/${sample}_mobsuite/contig_report_all.tsv
    fi
    # Initialize control variable to check in the header was printed
    header_printed=0
    # Check if any contig_report.txt files exist
    if find 14_mobsuite/${sample}_mobsuite/ -maxdepth 2 -type f -name "contig_report.txt" | grep -q .; then
        # Iterate over all found contig_report.txt files
        for file in 14_mobsuite/${sample}_mobsuite/*/contig_report.txt; do
            # Test if the header was not printed yet
            if [ $header_printed -eq 0 ]; then
                # If not, contatenate the entire file
                cat "$file" >> 14_mobsuite/${sample}_mobsuite/contig_report_all.tsv
                # Mark that the header was printed
                header_printed=1
            else
            # The header was printed, so concatenate file and ignore its header
            tail -n +2 "$file" >> 14_mobsuite/${sample}_mobsuite/contig_report_all.tsv  
            fi
            # Add a new line to separate the results of each sample
            # echo >> 14_mobsuite/${sample}_mobsuite/contig_report_all.tsv
        done
    fi
        # Merge mobtyper_results.txt files
    # Delete merged file if it exists
    if [ -f 14_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv ]; then
        rm 14_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv
    fi
    # Initialize control variable to check in the header was printed
    header_printed=0
    # Check if any mobtyper_results.txt files exist
    if find 14_mobsuite/${sample}_mobsuite/ -maxdepth 2 -type f -name "mobtyper_results.txt" | grep -q .; then
        # Iterate over all found mobtyper_results.txt files
        for file in 14_mobsuite/${sample}_mobsuite/*/mobtyper_results.txt; do
            # Test if the header was not printed yet
            if [ $header_printed -eq 0 ]; then
                # If not, contatenate the entire file
                cat "$file" >> 14_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv
                # Mark that the header was printed
                header_printed=1
            else
            # The header was printed, so concatenate file and ignore its header
            tail -n +2 "$file" >> 14_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv  
            fi
            # Add a new line to separate the results of each sample
            # echo >> 14_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv
        done
    fi

    # Merge mge.report.txt files
    # Delete merged file if it exists
    if [ -f 14_mobsuite/${sample}_mobsuite/mge.report_all.tsv ]; then
        rm 14_mobsuite/${sample}_mobsuite/mge.report_all.tsv
    fi
    # Initialize control variable to check in the header was printed
    header_printed=0
    # Check if any mge.report.txt files exist
    if find 14_mobsuite/${sample}_mobsuite/ -maxdepth 2 -type f -name "mge.report.txt" | grep -q .; then
        # Iterate over all found mge.report.txt files
        for file in 14_mobsuite/${sample}_mobsuite/*/mge.report.txt; do
            # Test if the header was not printed yet
            if [ $header_printed -eq 0 ]; then
                # If not, contatenate the entire file
                cat "$file" >> 14_mobsuite/${sample}_mobsuite/mge.report_all.tsv
                # Mark that the header was printed
                header_printed=1
            else
            # The header was printed, so concatenate file and ignore its header
            tail -n +2 "$file" >> 14_mobsuite/${sample}_mobsuite/mge.report_all.tsv  
            fi
            # Add a new line to separate the results of each sample
            # echo >> 14_mobsuite/${sample}_mobsuite/mge.report_all.tsv
        done
    fi

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt   
    
    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="14_mobsuite.tar.gz"
itens_to_compress=(14_mobsuite)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt

############################################################
## 14.2) VIBRANT

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="14) VIBRANT"
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ■■■" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 14_vibrant

# Calculate the sample count to display loop progress
i=1
sample_dirs=(13_pyrodigal/*/)
sample_count=${#sample_dirs[@]}

# Activate Conda environment
conda activate vibrant
# Loop through a list of sample directories
for dir in "${sample_dirs[@]}"; do

    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_pyrodigal/}
    
    # Inform current sample
    echo "▶ ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    # Start counting the running time
    loop_start_time=$SECONDS

    # Create an output directory
    mkdir -p "14_vibrant/${sample}_vibrant"

    # Loop through a list of files
    for file in ${dir}/*/*.fsa; do

        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.fsa}

        # Skip bin if output already exists and looks complete (final quality summary present)
        bin_outdir="14_vibrant/${sample}_vibrant/${binname}_vibrant"
        quality_file="${bin_outdir}/VIBRANT_${binname}/VIBRANT_results_${binname}/VIBRANT_genome_quality_${binname}.tsv"
        if [ -s "$quality_file" ]; then
            echo "${workflow_step} output already exists and is valid for bin: $binname from sample $sample. Skipping bin."
            continue
        elif [ -d "$bin_outdir" ]; then
            echo "${workflow_step} found incomplete output for bin: $binname from sample $sample. Removing partial directory and reprocessing."
            rm -rf "$bin_outdir"
        fi

        #Inform bin
        echo "14) VIBRANT is processing bin: ${binname} from sample ${sample}"

        # Create an output directory
        mkdir -p "$bin_outdir"

        # Run main software
        VIBRANT_run.py \
            -t $(nproc --ignore=1) \
            -f nucl \
            -d /db/vibrant/databases/ \
            -i ${file} \
            -folder "$bin_outdir"
    done

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔ ${workflow_step} — ${sample} (${i}/${sample_count}) — ${loop_running_time} @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt    
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
compressed_file="14_vibrant.tar.gz"
itens_to_compress=(14_vibrant)
echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"
# Generate checksum file of compressed directory file
md5sum "${compressed_file}" > "${compressed_file}".md5
# Check file integrity
echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo "■■■ ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — total: ${running_time} ■■■" | tee -a 0_workflow_progress.txt