# Bash Scripts for Shotgun Metagenome Analysis

This repository contains a collection of command-lines and **Bash scripts** essential for metagenome **analysis from short paired-end reads (Illumina)**.

It also includes detailed instructions for the installation of all necessary software.

---

## 📋 Usage Guidelines

* The scripts provide comprehensive **software installation instructions** and **a worfkflow for metagenome analysis**.

* **_script.sh Files (Modular Commands):** These files are collections of commands grouped by function (e.g., QC only, assembly only).

    * **⚠️ IMPORTANT: These scripts SHOULD NOT be executed in their entirety.**

    * Instead, you should copy (or modify) and paste the relevant command lines directly into your Linux terminal as needed for modular use.

---

## 🛠️ The Shotgun Metagenome Analysis Workflow

1) Reads files and renaming
    * Reads stored as local files
    * Reads from NCBI SRA
2) Raw reads quality assessment
    * FastQC
    * MultiQC
3) Raw reads trimming
    * Fastp
4) Trimmed reads quality assessment
    * FastQC
    * MultiQC
5) Host decontamination (optional)
    * NCBI Datasets
    * Bwa-mem2 index
    * Bwa-mem2 mapping
    * Bwa-mem2 reads
    * FastQC
    * MultiQC
6) Taxonomic abundance profile
    * Kraken
    * Kraken -> Bracken
    * Bracken -> Krona
    * Bracken -> Comparison
    * MetaPhlAn
    * MetaPhlAn -> Comparison
7) Metagenome assembly
    * MEGAHIT
    * QUAST
8) Functional abundance profile
    * Barrnap
    * Aragorn
    * Pyrodigal
    * eggNOG-mapper
    * dbCAN
    * AMRFinderPlus
    * VIBRANT
9) Binning (Single-sample only)
    * Bwa-mem2 index
    * Bwa-mem2 mapping
    * SemiBin
10) Bin quality control
    * QUAST
    * CheckM2
    * GUNC
    * GTDB-Tk
    * Barrnap
11) Bin functional abundance profile
    * Prokka
    * eggNOG-mapper
    * dbCAN
    * DeepGOPlus
    * AMRFinderPlus
12) Bin mobile elements
    * MOB-suite 
    * VIBRANT






