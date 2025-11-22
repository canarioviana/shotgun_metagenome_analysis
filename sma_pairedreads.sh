# Bash script for shotgun metagenome analysis from short-read sequencing data
#
# Author: Marcus Vinicius Canário Viana
# Date: 21/11/2025
# More info: see README.md in the repository
#
# Instructions:
# ⚠️ DO NOT execute this script entirely at once!
# Copy and paste individual command lines into the Linux terminal as needed.
# This file uses the .sh extension only to enable Bash syntax highlighting in text editors.


############################################################
## SUMMARY OF SHOTGUN METAGENOME ANAYLSIS WORKFLOW (FROM SHORT-READS)
############################################################

## 1) Reads files and renaming
    # 1.1) Reads stored as local files
    # 1.2) Reads from NCBI SRA
## 2) Raw reads quality assessment
    # 2.1) FastQC
    # 2.2) MultiQC
## 3) Raw reads trimming
    # 3.1) Fastp
## 4) Trimmed reads quality assessment
    # 4.1) FastQC
    # 4.2) MultiQC
## 5) Host decontamination (optional)
    # 5.1) NCBI Datasets
    # 5.2) Bwa-mem2 index
    # 5.3) Bwa-mem2 mapping
    # 5.4) FastQC
    # 5.5) MultiQC
## 6) Taxonomic abundance profile and prophages
    # 6.1) Kraken
    # 6.2) Kraken -> Bracken
    # 6.3) Kraken -> Bracken -> Krona
    # 6.4) Kraken -> Bracken -> Comparison
    # 6.5) MetaPhlAn
    # 6.6) MetaPhlAn -> Comparison
## 7) Metagenome assembly
    # 7.1) MEGAHIT
    # 7.2) QUAST
## 8) Functional abundance profile
    # 8.1) Barrnap
    # 8.2) Aragorn
    # 8.3) Pyrodigal
    # 8.4) eggNOG-mapper
    # 8.5) dbCAN
    # 8.6) AMRFinderPlus
    # 8.7) VIBRANT
## 9A) Binning (Single-sample only)
    # 9.1) Bwa-mem2 index
    # 9.2) Bwa-mem2 mapping
    # 9.3) SemiBin
## 9B) Binning (Single or Multi-sample)
    # 9.1) SemiBin concatenate_fasta
    # 9.2) Bwa-mem2 index
    # 9.3) Bwa-mem2 mapping
    # 9.4) Semibin binning
## 10) Bin quality control and taxonomy
    # 10.1) QUAST
    # 10.2) CheckM2
    # 10.3) GUNC
    # 10.4) Barrnap
    # 10.5) GTDB-Tk
## 11) Bin functional abundance profile
    # 11.1) Prokka
    # 11.2) eggNOG-mapper
    # 11.3) dbCAN
    # 11.4) DeepGOPlus
    # 11.5) AMRFinderPlus
## 12) Bin mobile genetic elements
    # 12.1) MOB-suite
    # 12.2) VIBRANT


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
## 1.2) Reads from NCBI SRA

# 1. Create a **tab-separated file** named **"1_reads_accessions.tsv"**.
# 2. This file **must contain** the NCBI SRA **accession number** in the first column and the **sample name** in the second column. Other columns will be ignored.
# 3. **Do not use** special characters in the sample names.
# 4. Place the **"1_reads_accessions.tsv"** file in the working directory.

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
        if [ -f "1_reads/${sample}_1.fq.gz" ] && [ -f "1_reads/${sample}_2.fq.gz" ]; then
            echo "1) Sample $sample paired files found. Skipping download."
        elif [ -f "1_reads/${sample}.fq.gz" ]; then
            echo "1) Sample $sample single-end file found. Skipping download."
        else
            echo "1) Downloading sample: $sample (accession: $accession)"

            # Remove any incomplete files
            rm -f "1_reads/${sample}_1.fq.gz" "1_reads/${sample}_2.fq.gz" "1_reads/${sample}.fq.gz"    
            
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
            rm -r "${accession}"
            # Compress files
            echo "Compressing fastq files."
            pigz -p $(nproc --ignore=1) ${accession}*.fastq

            # Check if the reads are paired-end
            if [ -f "${accession}_1.fastq.gz" ] && [ -f "${accession}_2.fastq.gz" ];  then
                echo "1) Sample ${sample} has paired-end reads."
                # Rename files
                echo "1) Renaming files."
                mv "${accession}_1.fastq.gz" "${sample}_1.fq.gz"
                mv "${accession}_2.fastq.gz" "${sample}_2.fq.gz"
                # Go back to main directory
                cd ..
            # Check if the reads are not paired-end
            elif [ -f "${accession}.fastq.gz" ]; then
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
            fi
        fi
    done
    echo "1) Download process complete. Deactivating the environment."
    # Deactivate Conda environment
    conda deactivate
else
    echo "1) The file 1_reads_accessions.tsv was not found. Proceeding using local files."
fi


############################################################
# 2) Raw reads quality assessment
############################################################

############################################################
## 2.1) FastQC

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="2) FastQC"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 2_fastqc
# Activate Conda environment
conda activate fastqc
# Run main software
fastqc -t $(nproc --ignore=1) 1_reads/*.gz -o 2_fastqc
# Deactivate Conda environment
conda deactivate
# Compress the output directory
zip -q -r 2_fastqc.zip 2_fastqc
# Generate checksum file of compressed directory file
md5sum 2_fastqc.zip > 2_fastqc.zip.md5

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
# Show the running time
echo "$workflow_step running time ${running_time}" | tee -a 0_workflow_progress.txt

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 2.2) MultiQC

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="2) MultiQC"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Activate Conda environment
conda activate multiqc
# Run main software
multiqc 2_fastqc/*_fastqc.zip -o 2_fastqc_multiqc
# Deactivate Conda environment
conda deactivate
# Compress the output directory
zip -q -r 2_fastqc_multiqc.zip 2_fastqc_multiqc
# Generate checksum file of compressed directory file
md5sum 2_fastqc_multiqc.zip > 2_fastqc_multiqc.zip.md5
# Delete the output directory
rm -r 2_fastqc 2_fastqc_multiqc

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 3) Raw reads trimming
############################################################

############################################################
## 3.1) Fastp

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="3) Fastp"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 3_fastp

# Calculate sample size
i=1
sample_count=$(ls -1 1_reads/*_1.fq.gz | wc -l)

# Activate Conda environment
conda activate fastp
# Loop through a list of sample files
for r1 in 1_reads/*_1.fq.gz; do
    # Obtain r2 path
    r2="${r1/_1.fq.gz/_2.fq.gz}"
    # Extract r1 file name
    r1filename=${r1##*/}
    # Extract sample name
    sample=${r1filename%%_*}
   
    # Inform current sample
    echo "3) Fastp is processing sample: ${sample} (${i}/${sample_count})"
    echo "3) Read 1 file: ${r1}"
    echo "3) Read 2 file: ${r2}"
    # Start counting the running time
    start_time=$SECONDS

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
        --out1 "3_fastp/${sample}_trimmed_1.fq.gz" \
        --out2 "3_fastp/${sample}_trimmed_2.fq.gz" \
        --html "3_fastp/${sample}_trimmed_fastp.html" \
        --json "3_fastp/${sample}_trimmed_fastp.json"
    
    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress output files
echo "Compressing output directory"
zip -q -r 3_fastp.zip 3_fastp/*.json 3_fastp/*.html
# Generate checksum file of compressed directory file
md5sum 3_fastp.zip > 3_fastp.zip.md5

# Generate checksum files for the reads
cd 3_fastp
for file in *.gz; do
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done    
cd ..

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 4) Trimmed reads quality assessment
############################################################

############################################################
## 4.1) FastQC

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="4) FastQC"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 4_fastqc
# Activate Conda environment
conda activate fastqc
# Run main software
fastqc -t $(nproc --ignore=1) 3_fastp/*.gz -o 4_fastqc
# Deactivate Conda environment
conda deactivate
# Compress the output directory
zip -q -r 4_fastqc.zip 4_fastqc
# Generate checksum file of compressed directory file
md5sum 4_fastqc.zip > 4_fastqc.zip.md5

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
# Show the running time
echo "$workflow_step running time ${running_time}" | tee -a 0_workflow_progress.txt

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 4.2) MultiQC

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="4) MultiQC"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Activate Conda environment
conda activate multiqc
# Run main software
multiqc 4_fastqc/*_fastqc.zip -o 4_fastqc_multiqc
# Deactivate Conda environment
conda deactivate
# Compress the output directory
echo "Compressing output directory"
zip -q -r 4_fastqc_multiqc.zip 4_fastqc_multiqc
# Generate checksum file of compressed directory file
md5sum 4_fastqc_multiqc.zip > 4_fastqc_multiqc.zip.md5
# Delete the output directory
rm -r 4_fastqc 4_fastqc_multiqc

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 5) Host decontamination (optional)
############################################################

# If you want to skip host decontamination, create a soft link that simulates that this step was done
ln -s 3_fastp 5_bwa_reads
ls -ld 5_bwa_reads

############################################################
## 5.1) NCBI Datasets (Download reference genomes)

# Create the tab-separated text file named "5_ref_ids.tsv"
# Column 1: The GenBank genome assembly ID of the reference genome. It will be used to download the genome.
# Column 2: The species name of the reference genome. It will name the Bwa-mem2 index. Do not use spaces or special characters.
# Place the file in the working directory.

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="5) NCBI Datasets"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Start counting the running time
start_time=$SECONDS

# Create output directory
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
# Unzip genome_data.zip
unzip genome_data.zip
# Rehydrate directory
datasets rehydrate --directory .
# Move assembly files to 5_ref_genomes
mv ncbi_dataset/data/*/*.fna 5_ref_genomes
# Remove sufix from file names and change file extension
cd 5_ref_genomes
rename 's/(.*?_.*?)_.*/$1.fasta/' *.fna
cd ..
# Delete temporary files and directory
rm -rf genome_data.zip ncbi_dataset README.md md5sum.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
# Show the running time
echo "$workflow_step running time ${running_time}" | tee -a 0_workflow_progress.txt

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 5.2) Bwa-mem2 index (Create Bwa-mem2 indexes)

# Software name for tracking progress in progress.txt
workflow_step="5) Bwa-mem2 index"
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create output directory
mkdir -p 5_bwa_index

# Activate Conda environment
conda activate bwa-mem2
# Loop through a list of file lines (Species as index name)
tr -d '\r' < 5_ref_ids.tsv | awk '1'| \
while IFS=$'\t' read -r ref_accession ref_name others; do
    # Create output directory
    mkdir -p 5_bwa_index/${ref_name}

    # Inform current sample
    echo "5) Bwa-mem2 index is processing reference: ${ref_name}"
    # Start counting the running time
    start_time=$SECONDS

    # Run main software
    bwa-mem2 index \
    -p "5_bwa_index/${ref_name}/${ref_name}" \
    "5_ref_genomes/${ref_accession}.fasta"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for ${ref_name}: running time ${running_time}" | tee -a 0_workflow_progress.txt
done
# Deactivate Conda environment
conda deactivate

# Generate a checksum file for the indexes
find 5_bwa_index -type f -exec md5sum {} + > 5_bwa_index.md5

# Compress original reference genome files
echo "Compressing reference genomes files"
for file in 5_ref_genomes/*.fasta; do
    pigz -p $(nproc --ignore=1) "${file}"
done

# Generate a checksum file for the original reference genome files
cd 5_ref_genomes
for file in *.fasta.gz; do
    md5sum "${file}" > ${file}.md5
done
cd ..

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 5.3) Bwa-mem2 mapping (Filtering out host reads)

# Software name for tracking progress in progress.txt
workflow_step="5) Bwa-mem2 mapping"
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create the tab-separated text file named "5_metagenomes.tsv":
# Column 1: The sample name.
# Column 2: The GenBank genome assembly ID of the reference genome. It will be used to register the specific version of the reference genome assembly.
# Column 3: The species name of the reference genome. The same name you used to create the Bwa-mem2 index.
# Column 4: The species name of the sample host (or isolation sorce). It will be used in the binning step.
# Place the file in the working directory.

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
    echo "5) Bwa-mem2 is processing sample ${sample} from ${ref_name}:
    3 (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Run main software (without intermediate files)
    bwa-mem2 mem \
    -t $(nproc --ignore=1) \
    "5_bwa_index/${ref_name}/${ref_name}" \
    "3_fastp/${sample}_trimmed_1.fq.gz" \
    "3_fastp/${sample}_trimmed_2.fq.gz" \
    2> "5_bwa_reads/${sample}_alignment.log" \
    | tee >(samtools flagstat - > "5_bwa_reads/${sample}_nofilter_flagstat.txt") \
    | samtools view -h -b -f 12 -@ $(nproc --ignore=1) - \
    | tee >(samtools flagstat - > "5_bwa_reads/${sample}_hostfilter_flagstat.txt") \
    | samtools sort -n -@ $(nproc --ignore=1) - \
    | samtools fastq \
    -1 "5_bwa_reads/${sample}_trimmed_nohost_1.fq.gz" \
    -2 "5_bwa_reads/${sample}_trimmed_nohost_2.fq.gz" \
    -0 /dev/null \
    -s /dev/null \
    -@ $(nproc --ignore=1) \
    -

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Generate checksum files for the reads
cd 5_bwa_reads
for file in *.gz; do
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done    
cd ..

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 5.4) FastQC

# Software name for tracking progress in progress.txt
workflow_step="5) FastQC"
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 5_bwa_reads_fastqc
# Activate Conda environment
conda activate fastqc
# Run main software
fastqc -t $(nproc --ignore=1) 5_bwa_reads/*.gz -o 5_bwa_reads_fastqc
# Deactivate Conda environment
conda deactivate
# Compress the output directory
zip -q -r 5_bwa_reads_fastqc.zip 5_bwa_reads_fastqc
# Generate checksum file of compressed directory file
md5sum 5_bwa_reads_fastqc.zip > 5_bwa_reads_fastqc.zip.md5

# Show the running time
echo "$workflow_step running time ${running_time}" | tee -a 0_workflow_progress.txt
# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 5.5) MultiQC

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="5) MultiQC"
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Activate Conda environment
conda activate multiqc
# Run main software
multiqc 5_bwa_reads_fastqc/*_fastqc.zip -o 5_bwa_reads_fastqc_multiqc
# Deactivate Conda environment
conda deactivate
# Compress the output directory
zip -q -r 5_bwa_reads_fastqc_multiqc.zip 5_bwa_reads_fastqc_multiqc
# Generate checksum file of compressed directory file
md5sum 5_bwa_reads_fastqc_multiqc.zip > 5_bwa_reads_fastqc_multiqc.zip.md5
# Delete the output directory
rm -r 5_bwa_reads_fastqc 5_bwa_reads_fastqc_multiqc 

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 6) Taxonomic abundance profile
############################################################

############################################################
## 6.1) Kraken

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) Kraken"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 6_kraken_report
mkdir -p 6_kraken_output

# Calculate sample size
i=1
sample_count=$(ls -1 5_bwa_reads/*_1.fq.gz | wc -l)

# Activate Conda environment
conda activate kraken2
# Loop through a list of sample files
for r1 in 5_bwa_reads/*_1.fq.gz; do
    # Obtain r2 path
    r2=${r1/_1.fq.gz/_2.fq.gz}
    # Extract r1 file name
    filename=${r1##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Inform current sample
    echo "6) Kraken is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS
    
    # Run main software
    kraken2 \
    --threads $(nproc --ignore=1) \
    --db /db/kraken2/k2_pluspf_20251015/ \
    --paired "${r1}" "${r2}" \
    --report "6_kraken_report/${sample}_kreport.tsv" \
    --report-minimizer-data \
    --minimum-hit-groups 3 \
    --report-zero-counts \
    | pigz -p $(nproc --ignore=1) -c > "6_kraken_output/${sample}.kraken.gz"
    
    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compressing the report directory
tar -c --use-compress-program=pigz -f 6_kraken_output.tar.gz 6_kraken_output
# Generate checksum of file
echo "Processing checksum of compressed report file."
md5sum 6_kraken_output.tar.gz > 6_kraken_output.tar.gz.md5

# Generate checksum of Kraken output (.kraken.gz) files 
# Go to target directory
cd 6_kraken_output
# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 *.kraken.gz | wc -l)
# Loop through a list of sample files
for file in *.kraken.gz; do
    echo "Processing checksum of file: ${file} (${i}/${sample_count})"
    md5sum ${file} > "${file}.md5"
    i=$((i + 1))
done
# Leave target directory
cd ..

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 6.2) Bracken

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) Kraken -> Bracken"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 6_bracken_report

# Calculate sample size
i=1
sample_count=$(ls -1 6_kraken_report/*_kreport.tsv | wc -l)

# Activate Conda environment
conda activate kraken2
# Loop through a list of sample files
for file in 6_kraken_report/*_kreport.tsv; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Inform current sample
    echo "6) Bracken is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Run main software
    bracken \
    -d /db/kraken2/k2_pluspf_20251015/ \
    -i "6_kraken_report/${sample}_kreport.tsv" \
    -w "6_bracken_report/${sample}_breport.tsv" \
    -o "6_bracken_report/${sample}.bracken" \
    -r 100 \
    -l S \
    -t 0

    # Convert report to MetaPhlAn format
    kreport2mpa.py \
    -r "6_bracken_report/${sample}_breport.tsv" \
    -o "6_bracken_report/${sample}_breport_mpa.tsv"

    # Calculate alfa diversity
    > "6_bracken_report/${sample}_breport_alfadiversity.txt"
    alpha_diversity.py -f "6_bracken_report/${sample}.bracken" -a BP >> "6_bracken_report/${sample}_breport_alfadiversity.txt"
    alpha_diversity.py -f "6_bracken_report/${sample}.bracken" -a Sh >> "6_bracken_report/${sample}_breport_alfadiversity.txt"
    alpha_diversity.py -f "6_bracken_report/${sample}.bracken" -a F >> "6_bracken_report/${sample}_breport_alfadiversity.txt"
    alpha_diversity.py -f "6_bracken_report/${sample}.bracken" -a Si >> "6_bracken_report/${sample}_breport_alfadiversity.txt"
    alpha_diversity.py -f "6_bracken_report/${sample}.bracken" -a ISi >> "6_bracken_report/${sample}_breport_alfadiversity.txt"
    
    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))"
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compressing the report directory
echo "Compressing output directory"
zip -q -r 6_bracken_report.zip 6_bracken_report
# Generate checksum of file
echo "Processing checksum of compressed report file."
md5sum 6_bracken_report.zip > 6_bracken_report.zip.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 6.3) Kraken -> Bracken -> Krona

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) Kraken -> Bracken -> Krona"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 6_bracken_krona

# Calculate sample size
i=1
sample_count=$(ls -1 6_bracken_report/*_breport.tsv | wc -l)

# Activate Conda environment
conda activate kraken2
# Loop through a list of sample files
for file in 6_bracken_report/*_breport.tsv; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Inform current sample
    echo "6) Generating Krona graphs for sample ${sample} (${i}/${sample_count})"

    # Generate .txt file
    kreport2krona.py \
    -r "${file}" \
    -o "6_bracken_krona/${sample}_bkrona.txt" \
    --no-intermediate-ranks

    # Generate .html file
    ktImportText "6_bracken_krona/${sample}_bkrona.txt" \
    -o "6_bracken_krona/${sample}_bkrona.html"
done
# Deactivate Conda environment
conda deactivate

# Compressing the report directory
echo "Compressing output directory"
zip -q -r 6_bracken_krona.zip 6_bracken_krona
# Generate checksum of file
echo "Processing checksum of compressed report file."
md5sum 6_bracken_krona.zip > 6_bracken_krona.zip.md5
# Delete output directory
rm -r 6_bracken_krona

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 6.4) Kraken -> Bracken -> Comparison

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) Kraken -> Bracken -> Comparison"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 6_bracken_comparison

# Start counting the running time
start_time=$SECONDS

# Activate Conda environment
conda activate kraken2

# Combine Bracken outputs of all samples
combine_bracken_outputs.py \
--files 6_bracken_report/*.bracken \
-o 6_bracken_comparison/combined_allsamples_breport.tsv

# Calculate beta diversity for all samples
beta_diversity.py \
--input-files 6_bracken_report/*.bracken \
--type bracken \
> 6_bracken_comparison/combined_allsamples_breport_betadiversity.txt

# Calculate beta diversity per isolaton source using information from file 5_metagenomes.tsv
input_file="5_metagenomes.tsv"
# awk to group samples by isolaton source
awk 'BEGIN {FS="\t"} 
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
}' "$input_file" | tr -d '\r' |
# Loop through each group
while IFS=$'\t' read -r source sample_list; do
    # Split the comma-separated sample list into a shell array
    IFS=',' read -r -a samples_array <<< "$sample_list"
    # Number of samples
    num_samples=${#samples_array[@]}
    # Check if number of samples is 2 or more
    if [ "$num_samples" -ge 2 ]; then
        # Declare empty array
        input_files_array=()
        # Iterate over all samples and apply the full prefix/suffix
        for sample in "${samples_array[@]}"; do
            input_files_array+=("6_bracken_report/${sample}.bracken")
        done
        # Run the script
        echo "6) Calculating beta diversity for samples from ${source} (${num_samples} samples):"
        beta_diversity.py \
        --input-files "${input_files_array[@]}" \
        --type bracken \
        > "6_bracken_comparison/combined_${source}_breport_betadiversity.txt"
    else
        # Skip host if there is only one sample
        echo "6) ⚠️ Skipping ${source}: only ${num_samples} sample(s) found (minimum 2 required)."
    fi
done

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
# Show the running time
echo "$workflow_step running time ${running_time}" | tee -a 0_workflow_progress.txt

# Deactivate Conda environment
conda deactivate

# Compressing the report directory
echo "Compressing output directory"
zip -q -r 6_bracken_comparison.zip 6_bracken_comparison
# Generate checksum of file
echo "Processing checksum of compressed report file"
md5sum 6_bracken_comparison.zip > 6_bracken_comparison.zip.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 6.5) MetaPhlAn

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) MetaPhlAn"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 6_metaphlan

# Calculate sample size
i=1
sample_count=$(ls -1 5_bwa_reads/*_1.fq.gz | wc -l)

# Activate Conda environment
conda activate metaphlan
# Loop through a list of sample files
for r1 in 5_bwa_reads/*_1.fq.gz; do
    # Obtain r2 path
    r2=${r1/_1.fq.gz/_2.fq.gz}
    # Extract r1 file name
    filename=${r1##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Inform current sample
    echo "6) MetaPhlAn is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Run main software
    metaphlan \
     "${r1}","${r2}" \
    --input_type fastq \
    --nproc $(nproc --ignore=1) \
    --verbose \
    --db_dir /db/metaphlan/ \
    --mapout "6_metaphlan/${sample}_metaphlan_bwa2.bz2" \
    -o "6_metaphlan/${sample}_mprofile.txt"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 6_metaphlan.tar.gz 6_metaphlan
# Create checksum file
md5sum 6_metaphlan.tar.gz > 6_metaphlan.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 6.6) MetaPhlAn -> Comparison

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) MetaPhlAn -> Comparison"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 6_metaphlan_comparison

# Activate Conda environment
conda activate metaphlan

# Merge abundance tables of all samples
merge_metaphlan_tables.py 6_metaphlan/*_mprofile.txt \
> 6_metaphlan_comparison/merged_allsamples_mprofile.txt

# Merge abundance tables per source using information from file 5_metagenomes.tsv
input_file="5_metagenomes.tsv" 
# awk to group samples by source
awk 'BEGIN {FS="\t"} 
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
}' "$input_file" | tr -d '\r' |
# Loop through each group
while IFS=$'\t' read -r source sample_list; do
    # Split the comma-separated sample list into a shell array
    IFS=',' read -r -a samples_array <<< "$sample_list"
    # Number of samples
    num_samples=${#samples_array[@]}
    # Check if number of samples is 2 or more
    if [ "$num_samples" -ge 2 ]; then
        # Declare empty array
        input_files_array=()
        # Iterate over all samples and apply the full prefix/suffix
        for sample in "${samples_array[@]}"; do
            # Metaphlan input file path: 6_metaphlan/sample_mprofile.txt
            input_files_array+=("6_metaphlan/${sample}_mprofile.txt")
        done
        # Run the script
        echo "6) Merging Metaphlan abundance tables for ${source} (${num_samples} samples):"
        merge_metaphlan_tables.py "${input_files_array[@]}" \
        > "6_metaphlan_comparison/merged_${source}_mprofile.txt"
    else
        echo "6) ⚠️ Skipping Metaphlan merge for ${source}: only ${num_samples} sample(s) found (minimum 2 required)."
        echo "${source}" > 6_metaphlan_comparison_skiped_soruces.tsv
    fi
done

# Calculate beta diversity
# Find path of the script "calculate_diversity.R" in the conda enviroment directory 
script=$(find $CONDA_PREFIX -name "calculate_diversity.R" 2>/dev/null)
# Run the script
for file in 6_metaphlan_comparison/merged_*_mprofile.txt; do
    output=${file%.txt}
    Rscript $script \
    -f ${file} \
    -d beta \
    -m bray-curtis \
    > "${output}_betadiversity.txt"
done

# Heatmap visualization
for file in 6_metaphlan_comparison/merged_*_mprofile.txt; do
    output=${file%.txt}
    hclust2.py \
      -i ${file} \
      -o "${output}_sqrt_scale.png" \
      --skip_rows 1 \
      --ftop 50 \
      --f_dist_f correlation \
      --s_dist_f braycurtis \
      --cell_aspect_ratio 9 \
      -s --fperc 99 \
      --flabel_size 4 \
      --metadata_rows 2,3,4 \
      --legend_file "${output}_sqrt_scale_legend.png" \
      --max_flabel_len 100 \
      --metadata_height 0.075 \
      --minv 0.01 \
      --no_slabels \
      --dpi 300 \
      --slinkage complete
done

# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 6_metaphlan_comparison.tar.gz 6_metaphlan_comparison
# Create checksum file
md5sum 6_metaphlan_comparison.tar.gz > 6_metaphlan_comparison.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 7) Metagenome assembly
############################################################

############################################################
## 7.1) MEGAHIT

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="7) MEGAHIT"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 7_megahit

# Calculate sample size
i=1
sample_count=$(ls -1 5_bwa_reads/*_1.fq.gz | wc -l)

# Activate Conda environment
conda activate megahit
# Loop through a list of sample files
for r1 in 5_bwa_reads/*_1.fq.gz; do
    # Obtain r2 path
    r2=${r1/_1.fq.gz/_2.fq.gz}
    # Extract r1 file name
    filename=${r1##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Inform current sample
    echo "7) MEGAHIT is assembling sample: $sample (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Run main software
    megahit \
    -t $(nproc --ignore=1) \
    -1 "${r1}" \
    -2 "${r2}" \
    -o "7_megahit/${sample}_megahit" \
    --min-contig-len 200

    # Move and rename assembly file
    mv "7_megahit/${sample}_megahit/final.contigs.fa" "7_megahit/${sample}_megahit.fasta"

    # Move and rename log file
    mv "7_megahit/${sample}_megahit/log" "7_megahit/${sample}_megahit.log"

    # Delete the sample directory
    rm -r "7_megahit/${sample}_megahit/"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 7_megahit.tar.gz 7_megahit
# Create checksum file
md5sum 7_megahit.tar.gz > 7_megahit.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 7.2) QUAST

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="7) QUAST"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 7_megahit_quast

# Activate Conda environment
conda activate quast
#Loop (usar --out-prefix "${sample}" na proxima vez)
for file in 7_megahit/*.fasta; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

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

# Compress directory
echo "Compressing output directory"
zip -q -r 7_megahit_quast.zip 7_megahit_quast
# Generate checksum file of compressed directory file
md5sum 7_megahit_quast.zip > 7_megahit_quast.zip.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 8) Functional abundance profile
############################################################

############################################################
## 8.1) Barrnap

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) Barrnap"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 8_barrnap

# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 7_megahit/*.fasta | wc | awk '{print $1}')

# Activate Conda environment
conda activate barrnap

# Loop through a list of sample files
for file in 7_megahit/*.fasta; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Inform current sample
    echo "8) Barrnap is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create output directory
    mkdir "8_barrnap/${sample}_barrnap"

    # Run barrnap for bacteria
    barrnap \
    --threads $(nproc --ignore=1) \
    --quiet \
    --kingdom bac \
    -o "8_barrnap/${sample}_barrnap/${sample}_bac_barrnap.fasta" \
    < $file \
    > "8_barrnap/${sample}_barrnap/${sample}_bac_barrnap.gff"
    # # Run the program for archaea
    # barrnap \
    # --threads $(nproc --ignore=1) \
    # --quiet \
    # --kingdom arc \
    # -o "8_barrnap/${sample}_barrnap/${sample}_arc_barrnap.fasta" \
    # < $file \
    # > "8_barrnap/${sample}_barrnap/${sample}_arc_barrnap.gff"
    
    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 8_barrnap.tar.gz 8_barrnap
# Create checksum file
md5sum 8_barrnap.tar.gz > 8_barrnap.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 8.2) Aragorn

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) Aragorn"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 8_aragorn

# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 7_megahit/*.fasta | wc | awk '{print $1}')

# Activate Conda environment
conda activate aragorn

# Loop through a list of sample files
for file in 7_megahit/*.fasta; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Inform current sample
    echo "8) Aragorn is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create output directory
    mkdir -p "8_aragorn/${sample}_aragorn"

    # Run main software
    aragorn \
    $file \
    > "8_aragorn/${sample}_aragorn/${sample}_aragorn.txt"
    
    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 8_aragorn.tar.gz 8_aragorn
# Create checksum file
md5sum 8_aragorn.tar.gz > 8_aragorn.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 8.3) Pyrodigal

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) Pyrodigal"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 8_pyrodigal

# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 7_megahit/*.fasta | wc | awk '{print $1}')

# Activate Conda environment
conda activate pyrodigal
# Loop through a list of sample files
for file in 7_megahit/*.fasta; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Inform current sample
    echo "8) Pyrodigal is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create output directory
    mkdir -p "8_pyrodigal/${sample}_pyrodigal"
    # Create link for nucleotide file
    # ln -s "$PWD/$file" "$PWD/8_pyrodigal/${sample}_pyrodigal/${sample}.fasta"

    # Run main software
    pyrodigal \
    -j $(nproc --ignore=1) \
    -m \
    -p meta \
    --no-stop-codon \
    -f gff \
    -i ${file} \
    -d "8_pyrodigal/${sample}_pyrodigal/${sample}.ffn" \
    -a "8_pyrodigal/${sample}_pyrodigal/${sample}.faa" \
    -o "8_pyrodigal/${sample}_pyrodigal/${sample}.gff"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))"
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 8_pyrodigal.tar.gz 8_pyrodigal
# Generate checksum file of compressed dire
md5sum 8_pyrodigal.tar.gz > 8_pyrodigal.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 8.4) eggNOG-mapper

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) eggNOG-mapper"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 8_emapper

# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 8_pyrodigal/*/*.faa | wc | awk '{print $1}')

# Activate Conda environment
conda activate eggnog-mapper
# Loop through a list of sample files
for file in 8_pyrodigal/*/*.faa; do
    # Extract sample name
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.faa}
    
    # Inform current sample
    echo "8) eggNOG-mapper is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create output directory
    mkdir -p "8_emapper/${sample}_emapper"

    # Run main software
    emapper.py \
    --cpu $(nproc --ignore=1) \
    -i $file \
    --output "${sample}" \
    --output_dir "8_emapper/${sample}_emapper"

    # Adjust output table
    tail +5 "8_emapper/${sample}_emapper/${sample}.emapper.annotations" \
    > "8_emapper/${sample}_emapper/${sample}_emapper.tsv"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 8_emapper.tar.gz 8_emapper
# Generate checksum file of compressed directory file
md5sum 8_emapper.tar.gz > 8_emapper.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 8.5) dbCAN

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) dbCAN"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 8_dbcan

# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 8_pyrodigal/*/*.faa | wc | awk '{print $1}')

# Activate Conda environment
conda activate dbcan
# Loop through a list of sample files
for file in 8_pyrodigal/*/*.faa; do
    # Extract sample name
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.faa}
    
    # Inform current sample
    echo "8) dbCAN is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create output directory
    mkdir -p "8_dbcan/${sample}_dbcan"

    # Run main software
    run_dbcan CAZyme_annotation \
    --threads $(nproc --ignore=1) \
    --db_dir /db/dbcan \
    --mode protein \
    --input_raw_data $file \
    --output_dir "8_dbcan/${sample}_dbcan"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 8_dbcan.tar.gz 8_dbcan
# Create checksum file
md5sum 8_dbcan.tar.gz > 8_dbcan.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
## 8.6) AMRFinderPlus

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) AMRFinderPlus"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 8_amrfinder

# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 8_pyrodigal/*/*.faa | wc | awk '{print $1}')

# Activate Conda environment
conda activate amrfinder
# Loop through a list of sample files
for file in 8_pyrodigal/*/*.faa; do
    # Extract sample name
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.faa}
    
    # Inform current sample
    echo "8) AMRFinderPlus is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Method 1 - Run AMRFinderPlus using nucleotide and proteins sequences
    # Format faa files for AMRFinderPlus
    awk '{
      if ($0 ~ /^>/) {
        match($0, /^>([^ ]+)/, a)
        id=a[1]
        sub(/ID=[^;]+/, "ID="id)
        print
      } else print
    }' "${file}" > "8_amrfinder/${sample}_amrfinder_format.faa" 
    # Run main software
    amrfinder \
    --threads $(nproc --ignore=1) \
    --database /db/amrfinder/latest \
    --plus \
    --annotation_format prodigal \
    --nucleotide "7_megahit/${sample}_megahit.fasta" \
    --protein "8_amrfinder/${sample}_amrfinder_format.faa"  \
    --gff "8_pyrodigal/${sample}_pyrodigal/${sample}.gff" \
    --name "${sample}" \
    --nucleotide_output "8_amrfinder/${sample}_amrfinder.fasta"\
    --protein_output "8_amrfinder/${sample}_amrfinder.faa"\
    --output "8_amrfinder/${sample}_amrfinder.tsv"
    # Remove intermediate file
    rm "8_amrfinder/${sample}_amrfinder_format.faa"

    # # Method 2 - Run AMRFinderPlus using only proteins sequences
    # # Run main software
    # amrfinder \
    # --threads $(nproc --ignore=1) \
    # --database /db/amrfinder/latest \
    # --plus \
    # --protein "${file}"  \
    # --name "${sample}" \
    # --protein_output "8_amrfinder/${sample}_amrfinder.faa"\
    # --output "8_amrfinder/${sample}_amrfinder.tsv"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 8_amrfinder.tar.gz 8_amrfinder
# Create checksum file
md5sum 8_amrfinder.tar.gz > 8_amrfinder.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 8.7) VIBRANT

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) VIBRANT"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 8_vibrant

# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 7_megahit/*.fasta | wc | awk '{print $1}')

# Activate Conda environment
conda activate vibrant
# Loop through a list of sample files
for file in 7_megahit/*.fasta; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}
    
    # Inform current sample
    echo "8) VIBRANT is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create output directory
    mkdir -p "8_vibrant/${sample}_vibrant"

    # Run main software
    VIBRANT_run.py \
    -t $(nproc --ignore=1) \
    -f nucl \
    -d /db/vibrant/databases/ \
    -i "${file}" \
    -folder "8_vibrant/${sample}_vibrant"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 8_vibrant.tar.gz 8_vibrant
# Create checksum file
md5sum 8_vibrant.tar.gz > 8_vibrant.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 9A) Binning (Single-sample only)
############################################################

# The default method of binning in this session uses human gut model
# You can change the model in step 9.3
# You can change the SemiBin parameter "--environment" to one of the following
# human_gut
# dog_gut
# ocean
# soil
# cat_gut
# human_oral
# mouse_gut
# pig_gut
# built_environment
# wastewater
# chicken_caecum
# global

############################################################
## 9.1) Bwa-mem2 index

# Software name for tracking progress in progress.txt
workflow_step="9) Bwa-mem2 index"
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 7_megahit/*.fasta | wc | awk '{print $1}')

# Create output directory
mkdir -p 9_bwa_index

# Activate Conda environment
conda activate bwa-mem2
# Loop through a list of sample files
for file in 7_megahit/*.fasta; do
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}

    # Create output directory
    mkdir -p "9_bwa_index/${sample}"
    # Inform current sample
    echo "9) Bwa-mem2 index is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Run main software
    bwa-mem2 index \
    -p "9_bwa_index/${sample}/${sample}" \
    ${file}

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing directory: 9_bwa_index"
tar -c --use-compress-program=pigz -f 9_bwa_index.tar.gz 9_bwa_index

# Create checksum file
md5sum 9_bwa_index.tar.gz > 9_bwa_index.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 9.2) Bwa-mem2 mapping

# Software name for tracking progress in progress.txt
workflow_step="9) Bwa-mem2 mapping"
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create output directory
mkdir -p 9_bwa_mapping

# Calculate the sample count to display loop progress
i=1
dir=(9_bwa_index/*/)
sample_count=${#dir[@]}

# Activate conda environment
conda activate bwa-mem2
# Loop through a list of sample files
for dir in 9_bwa_index/*/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%/}

    # Inform current sample
    echo "9) Bwa-mem2 is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Map the reads to contigs and generate bam file (Faster, no intermediate files)
    bwa-mem2 mem \
        -t $(nproc --ignore=1) \
        "9_bwa_index/${sample}/${sample}" \
        5_bwa_reads/"${sample}"*_1.fq.gz \
        5_bwa_reads/"${sample}"*_2.fq.gz \
        2> "9_bwa_mapping/${sample}_alignment.log" \
    | tee >(samtools flagstat - > "9_bwa_mapping/${sample}_allreads_flagstat.txt") \
    | samtools view -b -h -F 4 -@ $(nproc --ignore=1) - \
    | tee >(samtools flagstat - > "9_bwa_mapping/${sample}_mappedreads_flagstat.txt") \
    | samtools sort -@ $(nproc --ignore=1) \
        -o "9_bwa_mapping/${sample}.mapped.sorted.bam"

    # Create the index bai file for the sorted bam file
    samtools index \
        "9_bwa_mapping/${sample}.mapped.sorted.bam"

    # # Map the reads to contigs and generate bam file (with intermediate files)
    # bwa-mem2 mem \
    #     -t "$(nproc --ignore=1)" \
    #     "9_bwa_index/${sample}/${sample}" \
    #     "5_bwa_reads/${sample}_trimmed_nohost_1.fq.gz" \
    #     "5_bwa_reads/${sample}_trimmed_nohost_2.fq.gz" \
    #     > "9_bwa_mapping/${sample}.sam" \
    #     2> "9_bwa_mapping/${sample}_alignment.log"
    # # Generate mapping report
    # samtools flagstat \
    #     "9_bwa_mapping/${sample}.sam" \
    #     > "9_bwa_mapping/${sample}_allreads_flagstat.txt"
    # # Convert to bam and keep the header and only mapped reads
    # samtools view -b -h -F 4 -@ "$(nproc --ignore=1)" \
    #     "9_bwa_mapping/${sample}.sam" \
    #     > "9_bwa_mapping/${sample}.mapped.bam"
    # # Generate mapping report
    # samtools flagstat \
    #     "9_bwa_mapping/${sample}.mapped.bam" \
    #     > "9_bwa_mapping/${sample}_mappedreads_flagstat.txt"
    # # Delete intermediary file
    # rm "9_bwa_mapping/${sample}.sam"
    # # Sort reads by coordinate
    # samtools sort -@ "$(nproc --ignore=1)" \
    #     -o "9_bwa_mapping/${sample}.mapped.sorted.bam" \
    #     "9_bwa_mapping/${sample}.mapped.bam"
    # # Delete intermediary file
    # rm "9_bwa_mapping/${sample}.mapped.bam"
    # # Create the index bai file for the sorted bam file
    # samtools index \
    #     "9_bwa_mapping/${sample}.mapped.sorted.bam"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Generate checksum files for the reads
cd 9_bwa_mapping
for file in *.mapped.sorted.bam; do
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done    
cd ..

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 9.3) SemiBin binning (single-sample binning)

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) SemiBin binning (single-sample binning)"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 9_semibin

# Calculate the sample count to display loop progress
i=1
sample_count=$(ls -1 7_megahit/*.fasta | wc | awk '{print $1}')

# Activate Conda environment
conda activate semibin
# Loop through a list of sample files
for file in 7_megahit/*.fasta; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%_*}
    
    # Inform current sample
    echo "9) SemiBin is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create output directory
    mkdir -p "9_semibin/${sample}_semibin"

    # SemiBin (Easy mode) using short reads mapping, single binning AND HUMAN GUT MODEL (--environment human_gut)
    SemiBin2 single_easy_bin \
    --threads $(nproc --ignore=1) \
    --environment human_gut \
    --input-fasta ${file} \
    --input-bam "9_bwa_mapping/${sample}.mapped.sorted.bam" \
    --output "9_semibin/${sample}_semibin"

    # # SemiBin (Easy mode) using short reads mapping, single binning AND SELF-SUPERVISED MODEL (--self-supervised)
    # # Use GPU to reduce the required time to train the models
    # SemiBin2 single_easy_bin \
    # --threads $(nproc --ignore=1) \
    # --self-supervised \
    # --input-fasta ${file} \
    # --input-bam "9_bwa_mapping/${sample}.mapped.sorted.bam" \
    # --output "9_semibin/${sample}_semibin"

    # # SemiBin (Step by step) using short reads mapping and single binning
    # # Generate features (mandatory)
    # SemiBin2 generate_sequence_features_single \
    # --threads $(nproc --ignore=1) \
    # --input-fasta ${file} \
    # --input-bam "9_bwa_mapping/${sample}.mapped.sorted.bam" \
    # --output "9_semibin/${sample}_semibin"
    # # Self-train model (optional)
    # SemiBin2 train_self \
    # --threads $(nproc --ignore=1) \
    # --data "9_semibin/${sample}_semibin/data.csv" \
    # --data-split "9_semibin/${sample}_semibin/data_split.csv" \
    # -o "9_semibin/${sample}_semibin"
    # # Binning from short reads and self-trained model (optional)
    # # You can use the self trained model or other model
    # SemiBin2 bin_short \
    # --threads $(nproc --ignore=1) \
    # --input-fasta ${file} \
    # --model "9_semibin/${sample}_semibin/model.pt" \
    # --data "9_semibin/${sample}_semibin/data.csv" \
    # --output "9_semibin/${sample}_semibin"


    # Rename, extract and organize bin files
    # Rename and move bin files
    for file in 9_semibin/"${sample}"_semibin/output_bins/SemiBin_*.fa.gz; do
        [ -e "$file" ] || continue
        filename=$(basename "$file" .fa.gz)
        newname="${filename/#SemiBin_/${sample}Bin}"
        mv "$file" 9_semibin/"${sample}"_semibin/"$newname".fasta.gz
    done
    # Delete intermediate directory
    rm -r 9_semibin/"${sample}"_semibin/output_bins

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "9) $workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing directory: 9_semibin"
tar -c --use-compress-program=pigz -f 9_semibin.tar.gz 9_semibin
# Create checksum file
md5sum 9_semibin.tar.gz > 9_semibin.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 9B) Binning - Single/Multi-sample (Self-supervised mode)
############################################################

# This session uses information from file 5_metagenomes.tsv
# The information is used to decide to perform single binning (1 sample per isolation source / host) or multi-binning (>= 2 samples per isolation source / host) 
# The default method of binning in this session uses self-supervised mode
# Using a GPU is highly recommended. It can reduce the running time from hours to minutes

############################################################
## 9.1) SemiBin concatenate_fasta

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) Semibin concatenate_fasta"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 9_semibin_concatenate

# Activate Conda environment
conda activate semibin
# Inform input table
input_file="5_metagenomes.tsv"
# awk to group samples by source
awk 'BEGIN {FS="\t"} 
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
}' "$input_file" | tr -d '\r' |
# Loop through each group
while IFS=$'\t' read -r source sample_list; do
    # Split the comma-separated sample list into a shell array
    IFS=',' read -r -a samples_array <<< "$sample_list"
    # Number of samples
    num_samples=${#samples_array[@]}
    # Check if number of samples is 2 or more
    if [ "$num_samples" -ge 2 ]; then
        # Declare empty array
        input_files_array=()
        # Iterate over all samples and apply the full prefix/suffix
        for sample in "${samples_array[@]}"; do
            # Metaphlan input file path: 6_metaphlan/sample_mprofile.txt
            input_files_array+=("7_megahit/${sample}_megahit.fasta")
        done

        # Start counting the running time
        start_time=$SECONDS

        # Inform source and execute the command line
        echo "9) Concatenating metagenomes from source ${source} (${num_samples} samples):"
        # Run the script
        SemiBin2 concatenate_fasta \
        --input-fasta "${input_files_array[@]}" \
        --output "9_semibin_concatenate/${source}_concat"

        # Stop counting the running time
        elapsed_time=$((SECONDS - $start_time))
        running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
        # Show the running time
        echo "9) Concatenation of metagenomes from $source: running time ${running_time}" | tee -a 0_workflow_progress.txt

    else
        echo "9) ⚠️ Skipping concatenation for ${source}: only ${num_samples} sample(s) found (minimum 2 required)."
    fi
done

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 9.2) Bwa-mem2 index

# Software name for tracking progress in progress.txt
workflow_step="9) Bwa-mem2 index"
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Activate conda environment
conda activate bwa-mem2
# Inform input table
input_file="5_metagenomes.tsv"
# awk to group samples by source
awk 'BEGIN {FS="\t"} 
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
}' "$input_file" | tr -d '\r' |
# Loop through each group
while IFS=$'\t' read -r source sample_list; do
    # Split the comma-separated sample list into a shell array
    IFS=',' read -r -a samples_array <<< "$sample_list"
    # Number of samples
    num_samples=${#samples_array[@]}
    # Check if number of samples is 2 or more
    if [ "$num_samples" -ge 2 ]; then
        # Create output directory
        mkdir -p 9_bwa_index/${source}

        # Start counting the running time
        start_time=$SECONDS

        # Inform source and execute the command line
        echo "9) Creating index for concatenated metagenomes from source ${source} (${num_samples} samples):"
        
        # Uncompress input file
        gunzip -k 9_semibin_concatenate/${source}_concat/concatenated.fa.gz
        # Run main software
        bwa-mem2 index \
        -p "9_bwa_index/${source}/${source}" \
        9_semibin_concatenate/${source}_concat/concatenated.fa
        # Delete uncompressed input file
        rm 9_semibin_concatenate/${source}_concat/concatenated.fa

        # Stop counting the running time
        elapsed_time=$((SECONDS - $start_time))
        running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
        # Show the running time
        echo "9) Create index of concatenated metagenomes from $source: running time ${running_time}" | tee -a 0_workflow_progress.txt
    elif [ "$num_samples" -eq 1 ]; then
        # Create output directory
        mkdir -p 9_bwa_index/${source}

        # Start counting the running time
        start_time=$SECONDS

        # Inform source and execute the command line
        echo "9) Creating index for the metagenome from source ${source} (${num_samples} sample):"

        # Run main software
        bwa-mem2 index \
        -p "9_bwa_index/${source}/${source}" \
        "7_megahit/${sample_list}_megahit.fasta"

        # Stop counting the running time
        elapsed_time=$((SECONDS - $start_time))
        running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
        # Show the running time
        echo "9) Create index of concatenated metagenomes from $source: running time ${running_time}" | tee -a 0_workflow_progress.txt

    else
        echo "9) ⚠️ Skipping index for ${source}: no sample found (minimum 1 required)."
    fi
done

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 9.3) Bwa-mem2 mapping

# Software name for tracking progress in progress.txt
workflow_step="9) Bwa-mem2 mapping"
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create output directory
mkdir -p 9_bwa_mapping

# Calculate sample size
i=1
sample_count=$(awk 'END {print NR}' 5_metagenomes.tsv)

# Activate conda environment
conda activate bwa-mem2
# Loop through a list of file lines
tr -d '\r' < 5_metagenomes.tsv | awk '1'| \
while IFS=$'\t' read -r sample ref_accession ref_name isolation_source others; do
    
    # Inform current sample
    echo "9) Bwa-mem2 is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Map the reads to contigs and generate bam file (no intermediate files)
    bwa-mem2 mem \
        -t $(nproc --ignore=1) \
        "9_bwa_index/${isolation_source}/${isolation_source}" \
        5_bwa_reads/"${sample}"*_1.fq.gz \
        5_bwa_reads/"${sample}"*_2.fq.gz \
        2> "9_bwa_mapping/${sample}_alignment.log" \
    | tee >(samtools flagstat - > "9_bwa_mapping/${sample}_allreads_flagstat.txt") \
    | samtools view -b -h -F 4 -@ $(nproc --ignore=1) - \
    | tee >(samtools flagstat - > "9_bwa_mapping/${sample}_mappedreads_flagstat.txt") \
    | samtools sort -@ $(nproc --ignore=1) \
        -o "9_bwa_mapping/${sample}.mapped.sorted.bam"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Generate checksum files for the reads
cd 9_bwa_mapping
for file in *.bam; do
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done    
cd ..

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 9.4) Semibin binning

# Software name for tracking progress in progress.txt
workflow_step="9) Semibin binning"
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create output directory
mkdir -p 9_semibin

# Activate conda environment
conda activate semibin
input_file="5_metagenomes.tsv" 
# awk to group samples by source
awk 'BEGIN {FS="\t"} 
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
}' "$input_file" | tr -d '\r' |
# Loop through each group
while IFS=$'\t' read -r source sample_list; do
    # Split the comma-separated sample list into a shell array
    IFS=',' read -r -a samples_array <<< "$sample_list"
    # Number of samples
    num_samples=${#samples_array[@]}
    # Check if number of samples is 2 or more
    if [ "$num_samples" -ge 2 ]; then
        # Declare empty array
        input_files_array=()
        # Iterate over all samples and apply the full prefix/suffix
        for sample in "${samples_array[@]}"; do
           input_files_array+=("9_bwa_mapping/${sample}.mapped.sorted.bam")
        done

        # Inform source and execute the command line
        echo "9) SemiBin binning samples from ${source} (${num_samples} samples):"
        # Start counting the running time
        start_time=$SECONDS

        # Run the script
        SemiBin2 multi_easy_bin \
        --threads $(nproc --ignore=1) \
        --input-fasta "9_semibin_concatenate/${source}_concat/concatenated.fa.gz" \
        --input-bam ${input_files_array} \
        --output "9_semibin/${source}_semibin"

        # Organize output files per sampe
        for sample in "${samples_array[@]}"; do
            # Create sample directory
            mkdir 9_semibin/"${sample}"_semibin

            # Rename, extract and organize bin files
            # Rename and move bin files
            for file in 9_semibin/"${source}"_semibin/samples/"${sample}"_megahit/output_bins/*.fa.gz; do
                [ -e "$file" ] || continue
                filename=$(basename "$file" .fa.gz)
                newname="${filename/#SemiBin_/${sample}Bin}"
                mv "$file" 9_semibin/"${sample}"_semibin/"$newname".fasta.gz
            done
            # Move other bin files to sample directory
            mv 9_semibin/"${source}"_semibin/samples/"${sample}"_megahit/*.csv \
            9_semibin/"${source}"_semibin/samples/"${sample}"_megahit/*.pt \
            9_semibin/"${source}"_semibin/samples/"${sample}"_megahit/*.tsv \
            9_semibin/"${sample}"_semibin
        done

        # Move contacenated coverage file to 9_semibin/
        for concatcovfile in 9_semibin/"${source}"_semibin/samples/*.mapped.sorted.bam_0_data_cov.csv; do
            [ -e "$concatcovfile" ] || continue
            mv "$concatcovfile" 9_semibin/"${source}".mapped.sorted.bam_0_data_cov.csv
        done

        # Delete intermediate directory
        rm -r "9_semibin/${source}_semibin"

        # Stop counting the running time
        elapsed_time=$((SECONDS - $start_time))
        running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
        # Show the running time
        echo "9) Binning of metagenomes from ${source}: running time ${running_time}" | tee -a 0_workflow_progress.txt
    elif [ "$num_samples" -eq 1 ]; then
        # Inform sample and source and execute the command line
        echo "9) SemiBin binning sample from ${source} (${num_samples} sample):"
        # Start counting the running time
        start_time=$SECONDS

        # Use GPU to reduce the required time to train the models
        SemiBin2 single_easy_bin \
        --threads $(nproc --ignore=1) \
        --self-supervised \
        --input-fasta "7_megahit/${sample_list}_megahit.fasta" \
        --input-bam "9_bwa_mapping/${sample_list}.mapped.sorted.bam" \
        --output "9_semibin/${sample_list}_semibin"

        # Rename, extract and organize bin files
        # Rename and move bin files
        for file in 9_semibin/"${sample_list}"_semibin/output_bins/SemiBin_*.fa.gz; do
            [ -e "$file" ] || continue
            filename=$(basename "$file" .fa.gz)
            newname="${filename/#SemiBin_/${sample_list}Bin}"
            mv "$file" 9_semibin/"${sample_list}"_semibin/"$newname".fasta.gz
        done
        # Delete intermediate directory
        rm -r 9_semibin/"${sample_list}"_semibin/output_bins

        # Stop counting the running time
        elapsed_time=$((SECONDS - $start_time))
        running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
        # Show the running time
        echo "9) Binning of metagenomes from ${source}: running time ${running_time}" | tee -a 0_workflow_progress.txt
    else
        echo "⚠️ Skipping concatenation for ${source}: no sample found (minimum 1 required)."
    fi
done

# Deactivate Conda environment
conda activate base

# Compress directory
echo "Compressing directory: 9_semibin"
tar -c --use-compress-program=pigz -f 9_semibin.tar.gz 9_semibin
# Create checksum file
md5sum 9_semibin.tar.gz > 9_semibin.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 10) Bin quality control
############################################################

############################################################
## 10.1) QUAST

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) QUAST"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 10_quast

# Calculate the sample count to display loop progress
i=1
dir=(9_semibin/*_semibin/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate quast

# Loop through a list of sample directories
for dir in 9_semibin/*_semibin/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_semibin/}

    # Inform current sample
    echo "10) QUAST is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create an output directory
    mkdir -p "10_quast/${sample}_quast"

    # Run main software
    quast.py -t $(nproc --ignore=1) -m 0 -o "10_quast/${sample}_quast" 9_semibin/${sample}_semibin/*.fasta.gz

    # Copy transposed report
    cp "10_quast/${sample}_quast/transposed_report.tsv" "10_quast/${sample}_quast.tsv"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done

# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 10_quast.tar.gz 10_quast
# Create checksum file
md5sum 10_quast.tar.gz > 10_quast.tar.gz.md5
# Delete the output directory
rm -r 10_quast

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 10.2) CheckM2

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) CheckM2"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 10_checkm2

# Calculate the sample count to display loop progress
i=1
dir=(9_semibin/*_semibin/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate checkm2
# Loop through a list of sample directories
for dir in 9_semibin/*_semibin/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_semibin/}
    
    # Inform current sample
    echo "10) CheckM2 is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Run main software
    checkm2 predict \
    --threads $(nproc --ignore=1) \
    -x fasta.gz \
    --input "9_semibin/${sample}_semibin" \
    --output-directory "10_checkm2/${sample}_checkm2"

    # Copy and rename the output file
    cp "10_checkm2/${sample}_checkm2/quality_report.tsv" "10_checkm2/${sample}_checkm2.tsv"

    # Delete the samples output directories
    rm -r "10_checkm2/${sample}_checkm2"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 10_checkm2.tar.gz 10_checkm2
# Create checksum file
md5sum 10_checkm2.tar.gz > 10_checkm2.tar.gz.md5
# Delete the output directory
rm -r 10_checkm2

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 10.3) GUNC

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) GUNC"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 10_gunc

# Calculate the sample count to display loop progress
i=1
dir=(9_semibin/*_semibin/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate gunc

# Loop through a list of sample directories
for dir in 9_semibin/*_semibin/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_semibin/}

    # Inform current sample
    echo "10) GUNC is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create an output directory
    mkdir -p "10_gunc/${sample}_gunc" "10_gunc/${sample}_gunc_temp"

    # Run main software
    gunc run \
    --threads $(nproc --ignore=1) \
    --contig_taxonomy_output \
    --file_suffix .fasta.gz \
    --input_dir "9_semibin/${sample}_semibin" \
    --temp_dir "10_gunc/${sample}_gunc_temp" \
    --out_dir "10_gunc/${sample}_gunc"
    # Copy and rename the output file
    cp 10_gunc/${sample}_gunc/*maxCSS_level.tsv "10_gunc/${sample}_gunc.tsv"

    # # Plotting the data
    # # Create an output directory
    # mkdir -p 10_gunc/${sample}_gunc_plot
    # # Loop through a list of files
    # j=1
    # bin_count=$(ls -1 10_gunc/${sample}_gunc/diamond_output/*.out | wc | awk '{print $1}')
    # for plotfile in 10_gunc/${sample}_gunc/diamond_output/*.out; do\
    #    plotfilename=${plotfile##*/}
    #    plotfilename=${plotfilename%%.diamond*}
    #    echo "GUNC is ploting bin ${plotfilename} from sample: ${sample} (${j}/${bin_count})"
    #    gunc plot \
    #    --contig_display_num 0 \
    #    --diamond_file $plotfile \
    #    --out_dir 10_gunc/${sample}_gunc_plot/;\
    #    j=$((j + 1))
    # done

    # Delete the intermediary directories
    rm -r 10_gunc/${sample}_gunc 10_gunc/${sample}_gunc_temp

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 10_gunc.tar.gz 10_gunc
# Create checksum file
md5sum 10_gunc.tar.gz > 10_gunc.tar.gz.md5
# Delete the output directory
rm -r 10_gunc

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 10.4) Barrnap

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) Barrnap"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 10_barrnap

# Calculate the sample count to display loop progress
i=1
dir=(9_semibin/*_semibin/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate barrnap

# Loop through a list of sample directories
for dir in 9_semibin/*_semibin/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_semibin/}

    # Inform current sample
    echo "10) Barrnap is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS


    # Create an output directory
    mkdir -p "10_barrnap/${sample}_barrnap"

    # Loop through a list of files
    for file in 9_semibin/${sample}_semibin/*.fasta.gz; do
        # Extract sample name
        prefix=$(basename ${file} .fasta.gz)

        # Run barrnap for bacteria
        zcat $file | barrnap \
        --threads $(nproc --ignore=1) \
        --quiet \
        --kingdom bac \
        -o "10_barrnap/${sample}_barrnap/${prefix}_bac_barrnap.fasta" \
        > "10_barrnap/${sample}_barrnap/${prefix}_bac_barrnap.gff"
        # Run the program for archaea
        zcat $file | barrnap \
        --threads $(nproc --ignore=1) \
        --quiet \
        --kingdom arc \
        -o "10_barrnap/${sample}_barrnap/${prefix}_arc_barrnap.fasta" \
        > "10_barrnap/${sample}_barrnap/${prefix}_arc_barrnap.gff"
    done
    
    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 10_barrnap.tar.gz 10_barrnap
# Create checksum file
md5sum 10_barrnap.tar.gz > 10_barrnap.tar.gz.md5
# Delete the output directory
rm -r 10_barrnap

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 10.5) GTDB-Tk

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) GTDB-Tk"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Requires 64GB of RAM if the species is not identified by the ANI screening step

# Create an output directory
mkdir -p 10_gtdbtk

# Calculate the sample count to display loop progress
i=1
dir=(9_semibin/*_semibin/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate gtdbtk

# Loop through a list of sample directories
for dir in 9_semibin/*_semibin/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_semibin/}

    # Inform current sample
    echo "10) GTDB-Tk is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create an output directory
    mkdir -p "10_gtdbtk/${sample}_gtdbtk"

    # Run main software
    gtdbtk classify_wf \
    --cpus $(nproc --ignore=1) \
    --extension .fasta.gz \
    --mash_db /db/gtdbtk \
    --genome_dir "9_semibin/${sample}_semibin" \
    --out_dir "10_gtdbtk/${sample}_gtdbtk"

    # Copy and rename output files
    if [ -f "10_gtdbtk/${sample}_gtdbtk/gtdbtk.bac120.summary.tsv" ]; then
        cp "10_gtdbtk/${sample}_gtdbtk/gtdbtk.bac120.summary.tsv" "10_gtdbtk/${sample}_gtdbtk_bacteria.tsv"
    fi
    if [ -f "10_gtdbtk/${sample}_gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
        cp "10_gtdbtk/${sample}_gtdbtk/gtdbtk.ar53.summary.tsv" "10_gtdbtk/${sample}_gtdbtk_archaea.tsv"
    fi

    # Delete the temporary directorys
    rm -r "10_gtdbtk/${sample}_gtdbtk"

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    i=$((i + 1))
done

# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 10_gtdbtk.tar.gz 10_gtdbtk
# Create checksum file
md5sum 10_gtdbtk.tar.gz > 10_gtdbtk.tar.gz.md5
# Delete the output directory
rm -r 10_gtdbtk

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 11) Bin functional abundance profile
############################################################

############################################################
## 11.1) Prokka 

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="11) Prokka"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 11_prokka

# Calculate the sample count to display loop progress
i=1
dir=(9_semibin/*_semibin/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate prokka
# Loop through a list of sample directories
for dir in 9_semibin/*_semibin/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_semibin/}
    
    # Inform current sample
    echo "11) Prokka is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create an output directory
    mkdir -p "11_prokka/${sample}_prokka"

    # Loop through a list of files
    for file in 9_semibin/${sample}_semibin/*.fasta.gz; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.fasta.gz}

        # Delete a uncompresed input file in case it already exists
        [ -f "9_semibin/${sample}_semibin/${binname}.fasta" ] && rm -f "9_semibin/${sample}_semibin/${binname}.fasta"

        # Extract input file
        gunzip -k "${file}" # 2>/dev/null

        # Run Prokka
        prokka \
        --cpus $(nproc --ignore=1) \
        --fast \
        --metagenome \
        --addgenes \
        --prefix ${binname} \
        --force \
        --outdir "11_prokka/${sample}_prokka/${binname}" \
        "9_semibin/${sample}_semibin/${binname}.fasta"

        # Delete uncompressed input file
        rm -f "9_semibin/${sample}_semibin/${binname}.fasta"
    done

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 11_prokka.tar.gz 11_prokka
# Create checksum file
md5sum 11_prokka.tar.gz > 11_prokka.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 11.2) eggNOG-mapper

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="11) eggNOG-mapper"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 11_emapper

# Calculate the sample count to display loop progress
i=1
dir=(11_prokka/*/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate eggnog-mapper
# Loop through a list of sample directories
for dir in 11_prokka/*_prokka/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_prokka/}
    
    # Inform current sample
    echo "11) eggNOG-mapper is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Loop through a list of files
    for file in ${dir}/*/*.faa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.faa}

        # Create an output directory
        mkdir -p "11_emapper/${sample}_emapper/${binname}_emapper"

        # Run main software
        emapper.py \
        --cpu $(nproc --ignore=1) \
        -i $file \
        --output "${binname}" \
        --output_dir "11_emapper/${sample}_emapper/${binname}_emapper"

        # Adjust output table
        tail +5 "11_emapper/${sample}_emapper/${binname}_emapper/${binname}.emapper.annotations" \
        > "11_emapper/${sample}_emapper/${binname}_emapper/${binname}_emapper.tsv"
    done

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
zip -q -r 11_emapper.zip 11_emapper
# Generate checksum file of compressed directory file
md5sum 11_emapper.zip > 11_emapper.zip.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
## 11.3) dbCAN

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="11) dbCAN"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 11_dbcan

# Calculate the sample count to display loop progress
i=1
dir=(11_prokka/*/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate dbcan
# Loop through a list of sample directories
for dir in 11_prokka/*_prokka/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_prokka/}
    
    # Inform current sample
    echo "11) dbCAN is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create an output directory
    mkdir -p "11_dbcan/${sample}_dbcan"

    # Loop through a list of files
    for file in ${dir}/*/*.faa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.faa}

        # Create an output directory
        mkdir -p "11_dbcan/${sample}_dbcan/${binname}"

        # Run main software
        run_dbcan CAZyme_annotation \
        --threads $(nproc --ignore=1) \
        --db_dir /db/dbcan \
        --mode protein \
        --input_raw_data $file \
        --output_dir "11_dbcan/${sample}_dbcan/${binname}"
    done

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
zip -q -r 11_dbcan.zip 11_dbcan
# Generate checksum file of compressed directory file
md5sum 11_dbcan.zip > 11_dbcan.zip.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 11.4) DeepGOPlus

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="11) DeepGOPlus"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 11_deepgoplus

# Calculate the sample count to display loop progress
i=1
dir=(11_prokka/*/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate deepgoplus
# Loop through a list of sample directories
for dir in 11_prokka/*_prokka/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_prokka/}
    
    # Inform current sample
    echo "11) DeepGOPlus is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create an output directory
    mkdir -p "11_deepgoplus/${sample}_deepgoplus"

    # Loop through a list of files
    for file in ${dir}/*/*.faa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.faa}

        # Run main software
        deepgoplus \
        --data-root /db/deepgoplus/data/ \
        --in-file $file \
        --out-file "11_deepgoplus/${sample}_deepgoplus/${binname}_deepgoplus.tsv"
    done

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
zip -q -r 11_deepgoplus.zip 11_deepgoplus
# Generate checksum file of compressed directory file
md5sum 11_deepgoplus.zip > 11_deepgoplus.zip.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 11.5) AMRFinderPlus

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="11) AMRFinderPlus"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 11_amrfinder

# Calculate the sample count to display loop progress
i=1
dir=(11_prokka/*/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate amrfinder
# Loop through a list of sample directories
for dir in 11_prokka/*_prokka/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_prokka/}
    
    # Inform current sample
    echo "11) AMRFinderPlus is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create an output directory
    mkdir -p "11_amrfinder/${sample}_amrfinder"

    # Loop through a list of files
    for file in ${dir}/*/*.faa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.faa}

        # Run main software
        amrfinder \
        --threads $(nproc --ignore=1) \
        --database /db/amrfinder/latest \
        --plus \
        --annotation_format prokka \
        --nucleotide "11_prokka/${sample}_prokka/${binname}/${binname}.fsa" \
        --protein $file \
        --gff "11_prokka/${sample}_prokka/${binname}/${binname}.gff" \
        --name "${binname}" \
        --nucleotide_output "11_amrfinder/${sample}_amrfinder/${binname}_amrfinder.fasta"\
        --protein_output "11_amrfinder/${sample}_amrfinder/${binname}_amrfinder.faa"\
        --output "11_amrfinder/${sample}_amrfinder/${binname}_amrfinder.tsv"
    done

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
zip -q -r 11_amrfinder.zip 11_amrfinder
# Generate checksum file of compressed directory file
md5sum 11_amrfinder.zip > 11_amrfinder.zip.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt


############################################################
# 12) Bin mobile elements
############################################################

############################################################
## 12.1) MOB-suite

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="12) MOB-suite"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 12_mobsuite

# Calculate the sample count to display loop progress
i=1
dir=(11_prokka/*_prokka/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate mob_suite
# Loop through a list of sample directories
for dir in 11_prokka/*_prokka/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_prokka/}
    
    # Inform current sample
    echo "12) MOB-suite is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create an output directory
    mkdir -p "12_mobsuite/${sample}_mobsuite"

    # Loop through a list of files
    for file in "${dir}"/*/*.fsa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.fsa}

        # Run main software
        mob_recon -n $(nproc --ignore=1) \
        --infile "${file}" \
        --outdir "12_mobsuite/${sample}_mobsuite/${binname}_mobsuite"

        # Rename MOB-suite fasta files
        outdir="12_mobsuite/${sample}_mobsuite/${binname}_mobsuite"
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
    if [ -f 12_mobsuite/${sample}_mobsuite/contig_report_all.tsv ]; then
        rm 12_mobsuite/${sample}_mobsuite/contig_report_all.tsv
    fi
    # Initialize control variable to check in the header was printed
    header_printed=0
    # Check if any contig_report.txt files exist
    if find 12_mobsuite/${sample}_mobsuite/ -maxdepth 2 -type f -name "contig_report.txt" | grep -q .; then
        # Iterate over all found contig_report.txt files
        for file in 12_mobsuite/${sample}_mobsuite/*/contig_report.txt; do
            # Test if the header was not printed yet
            if [ $header_printed -eq 0 ]; then
                # If not, contatenate the entire file
                cat "$file" >> 12_mobsuite/${sample}_mobsuite/contig_report_all.tsv
                # Mark that the header was printed
                header_printed=1
            else
            # The header was printed, so concatenate file and ignore its header
            tail -n +2 "$file" >> 12_mobsuite/${sample}_mobsuite/contig_report_all.tsv  
            fi
            # Add a new line to separate the results of each sample
            # echo >> 12_mobsuite/${sample}_mobsuite/contig_report_all.tsv
        done
    fi
        # Merge mobtyper_results.txt files
    # Delete merged file if it exists
    if [ -f 12_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv ]; then
        rm 12_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv
    fi
    # Initialize control variable to check in the header was printed
    header_printed=0
    # Check if any mobtyper_results.txt files exist
    if find 12_mobsuite/${sample}_mobsuite/ -maxdepth 2 -type f -name "mobtyper_results.txt" | grep -q .; then
        # Iterate over all found mobtyper_results.txt files
        for file in 12_mobsuite/${sample}_mobsuite/*/mobtyper_results.txt; do
            # Test if the header was not printed yet
            if [ $header_printed -eq 0 ]; then
                # If not, contatenate the entire file
                cat "$file" >> 12_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv
                # Mark that the header was printed
                header_printed=1
            else
            # The header was printed, so concatenate file and ignore its header
            tail -n +2 "$file" >> 12_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv  
            fi
            # Add a new line to separate the results of each sample
            # echo >> 12_mobsuite/${sample}_mobsuite/mobtyper_results_all.tsv
        done
    fi

    # Merge mge.report.txt files
    # Delete merged file if it exists
    if [ -f 12_mobsuite/${sample}_mobsuite/mge.report_all.tsv ]; then
        rm 12_mobsuite/${sample}_mobsuite/mge.report_all.tsv
    fi
    # Initialize control variable to check in the header was printed
    header_printed=0
    # Check if any mge.report.txt files exist
    if find 12_mobsuite/${sample}_mobsuite/ -maxdepth 2 -type f -name "mge.report.txt" | grep -q .; then
        # Iterate over all found mge.report.txt files
        for file in 12_mobsuite/${sample}_mobsuite/*/mge.report.txt; do
            # Test if the header was not printed yet
            if [ $header_printed -eq 0 ]; then
                # If not, contatenate the entire file
                cat "$file" >> 12_mobsuite/${sample}_mobsuite/mge.report_all.tsv
                # Mark that the header was printed
                header_printed=1
            else
            # The header was printed, so concatenate file and ignore its header
            tail -n +2 "$file" >> 12_mobsuite/${sample}_mobsuite/mge.report_all.tsv  
            fi
            # Add a new line to separate the results of each sample
            # echo >> 12_mobsuite/${sample}_mobsuite/mge.report_all.tsv
        done
    fi

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 12_mobsuite.tar.gz 12_mobsuite
# Create checksum file
md5sum 12_mobsuite.tar.gz > 12_mobsuite.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

############################################################
## 12.2) VIBRANT

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="12) VIBRANT"
# Update the file 0_workflow_progress.txt
echo "$workflow_step step started at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

# Create an output directory
mkdir -p 12_vibrant

# Calculate the sample count to display loop progress
i=1
dir=(11_prokka/*/)
sample_count=${#dir[@]}

# Activate Conda environment
conda activate vibrant
# Loop through a list of sample directories
for dir in 11_prokka/*_prokka/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_prokka/}
    
    # Inform current sample
    echo "12) VIBRANT is processing sample: ${sample} (${i}/${sample_count})"
    # Start counting the running time
    start_time=$SECONDS

    # Create an output directory
    mkdir -p "12_vibrant/${sample}_vibrant"

    # Loop through a list of files
    for file in ${dir}/*/*.fsa; do
        # Extract bin file name
        filename=${file##*/}
        # Extract bin name
        binname=${filename%%.fsa}

        # Create an output directory
        mkdir -p "12_vibrant/${sample}_vibrant/${binname}_vibrant"

        # Run main software
        VIBRANT_run.py \
        -t $(nproc --ignore=1) \
        -f nucl \
        -d /db/vibrant/databases/ \
        -i ${file} \
        -folder "12_vibrant/${sample}_vibrant/${binname}_vibrant"
    done

    # Stop counting the running time
    elapsed_time=$((SECONDS - $start_time))
    running_time=$(date -u -d "@$elapsed_time" +"%H:%M:%S")
    # Show the running time
    echo "$workflow_step for sample $sample: running time ${running_time} (Finished at $(date +'%Y-%m-%d %H:%M:%S'))" | tee -a 0_workflow_progress.txt
    
    i=$((i + 1))
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
tar -c --use-compress-program=pigz -f 12_vibrant.tar.gz 12_vibrant
# Create checksum file
md5sum 12_vibrant.tar.gz > 12_vibrant.tar.gz.md5

# Update the file 0_workflow_progress.txt
echo "$workflow_step step finished at $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
