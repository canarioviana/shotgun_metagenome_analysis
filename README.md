# Bash Scripts for Shotgun Metagenome Analysis

This repository contains a collection of command-lines and **Bash scripts** essential for metagenome **analysis from short paired-end reads (Illumina)**.

It also includes detailed instructions for the installation of all necessary software.

---

## Usage Guidelines

* The scripts provide comprehensive **software installation instructions** and **a worfkflow for metagenomic analysis**
* **_end2end.sh Files (Workflow):** These scripts encapsulate the complete workflow and **can be executed at once** (end-to-end). To execute them:
    
    * **Grant execution permission:**
        ```bash
        chmod +x XXX_end2end.sh
        ```
    * **Run the script:**
        ```bash
        ./XXX_end2end.sh
        ```

* **_script.sh Files (Modular Commands):** These files are collections of commands grouped by function (e.g., QC only, assembly only).

    * **⚠️ IMPORTANT: These scripts SHOULD NOT be executed in their entirety.**

    * Instead, you should copy (or modify) and paste the relevant command lines directly into your Linux terminal as needed for modular use.

---

* Metadata for sequencing reads from NCBI SRA (reads_accessions.tsv)

1. Create a **tab-separated file** named **"reads_accessions.tsv"**.
2. This file **must contain** the NCBI SRA **accession number** in the first column and the **sample name** in the second column. Other columns will be ignored.
3. **Do not use** special characters in the sample names.
4. Place the **"0_reads_accessions.tsv"** file in the working directory.

---

* Sequencing reads as local files (*_1.fq.gz and *_2.fq.gz)

1. The sequencing reads must be in FASTQ format and compressed, with the suffixes `_1.fq.gz` and `_2.fq.gz`, or `_1.fastq.gz` and `_2.fastq.gz` or `_R1_001.fastq.gz` and `_R2_001.fastq.gz`
2. In the working directory, create the directory `1_reads/` and place the read files inside it.

---

* Metadata for samples (metagenomes.tsv)

1. Create a tab-separated text file named `metagenomes.tsv` in the **working directory**, containing four columns in the following order. Any subsequent columns will be ignored:

| Column | Description |
| :--- | :--- |
| **`sample`** | The sample name. |
| **`ref_accession`** | The GenBank genome assembly ID of the reference genome. Used to register the specific version of the reference genome assembly. If the sample is not host-associated, use `NA`, `none`, or leave empty. |
| **`ref_name`** | The species name of the reference genome. Use the same name as in the `ref_name` column of `ref_genomes_ids.tsv`. If the sample is not host-associated, use `NA`, `none`, or leave empty. |
| **`isolation_source`** | The species name of the sample host (or isolation source). This will be used in the binning step. |

---

* Metadata for reference genomes from GenBank (ref_genomes_ids.tsv)

1. Create a tab-separated text file named `ref_genomes_ids.tsv` containing the following columns in this order. Any subsequent columns will be ignored:

| Column | Description |
| :--- | :--- |
| **`ref_accession`** | The GenBank genome assembly ID of the reference genome. Used to download the specific version of the reference genome assembly. |
| **`ref_name`** | The species name of the reference genome. Must match the name used to create the BWA-MEM2 index. |

---

## The Shotgun Metagenome Analysis Workflow

1) Reads files and renaming
    * Reads from NCBI SRA (sra-tools)
    * Reads stored as local files
2) Raw reads quality assessment
    * FastQC
    * FastQC -> MultiQC
3) Raw reads trimming
    * Fastp
4) Trimmed reads quality assessment
    * FastQC
    * FastQC -> MultiQC
5) Host decontamination (optional)
    * NCBI Datasets
    * Bwa-mem2 index
    * Bwa-mem2 mem
    * Bwa-mem2 reads
    * Bwa-mem2 -> FastQC
    * Bwa-mem2 -> FastQC -> MultiQC
6) Taxonomic abundance profile
    * Kraken
    * Kraken -> Bracken
    * Kraken -> Bracken -> Comparison
    * Kraken -> Bracken -> Krona
    * MetaPhlAn
    * MetaPhlAn -> Comparison
7) Metagenome assembly
    * MEGAHIT
    * MEGAHIT -> QUAST
8) Functional abundance profile and prophages
    * Pybarrnap
    * Aragorn
    * Pyrodigal
    * Pyrodigal -> AMRFinderPlus
    * Pyrodigal -> dbCAN
    * Pyrodigal -> eggNOG-mapper
    * Pyrodigal -> VFDB (BLASTP)
9) Functional abundance profile of gene catalog
    * MMseqs2 input (SeqKit)
    * MMseqs2 easy-linclust
    * MMseqs2 -> AMRFinderPlus
    * MMseqs2 -> dbCAN
    * MMseqs2 -> eggNOG-mapper
    * MMseqs2 -> Pyrodigal -> VFDB (BLASTP)
    * MMseqs2 -> Salmon index
    * MMseqs2 -> Salmon quant
10) Binning - Single/Multi-sample - Input files
    * Seqkit
    * Seqkit -> SemiBin concatenate_fasta
    * Seqkit -> SemiBin concatenate_fasta -> Minimap2 index
    * Seqkit -> SemiBin concatenate_fasta -> Minimap2
11) Binning - Single/Multi-sample (Self-supervised mode)
    * SemiBin
12) Bin quality control and taxonomy
    * CheckM2
    * GUNC
    * GTDB-Tk
    * QUAST
13) Bin functional abundance profile
    * Aragorn
    * Pybarrnap
    * Pyrodigal
    * Pyrodigal -> AMRFinderPlus
    * Pyrodigal -> dbCAN
    * Pyrodigal -> eggNOG-mapper
    * Pyrodigal -> VFDB (BLASTP)
14) Bin mobile genetic elements
    * MOB-suite
    * VIBRANT
