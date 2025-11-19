# Bash script for installation of software used for shotgun metagenome analysis from short-read sequencing
#
# ⚠️ DO NOT execute this script entirely at once!
# Copy and paste individual command lines into the Linux terminal as needed.
# This file uses the .sh extension only to enable Bash syntax highlighting in text editors.
#
# Author: Marcus Vinicius Canário Viana
# Date: 19/11/2025
# More info: see README.md in the repository


# Summary
# A) System requirements
# B) Software installation
# B.1) Software installation - Miniconda
# B.2) Software installation - Metagenomic workflow
# C) Connecting to a server and using Screen


##########################################################################
## A) System requirements
##########################################################################
# RAM:
# 16GB is the minimum required
# 64GB if you run taxonomic analysis with GTDB-Tk
# 128GB if you want to use the PlusPF database from Kraken 2
#
# Static IP address (For use on a computer server)

##########################################################################
## B) Software installation
##########################################################################

##########################################################################
## B.1) Software installation - Miniconda
##########################################################################

##########################################################################
# Miniconda (For a single user. Local computer.)
# Download the installer
cd
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x87_64.sh
sh Miniconda3-latest-Linux-x87_64.sh
# Follow the instructions to install Miniconda in /home/$USER/miniconda3

# Do you accept the license terms? [yes|no]
# >>> yes

# Miniconda3 will now be installed into this location:
# /root/miniconda3
#
#   - Press ENTER to confirm the location
#   - Press CTRL-C to abort the installation
#   - Or specify a different location below
#
# [/home/user_name/miniconda3] >>>

# Do you wish to update your shell profile to automatically initialize conda?
# This will activate conda on startup and change the command prompt when activated.
# If you'd prefer that conda's base environment not be activated on startup,
#    run the following command when conda is activated:

# conda config --set auto_activate_base false

# You can undo this by running `conda init --reverse $SHELL`? [yes|no]
# [no] >>> yes

# Apply changes in current session
source ~/.bashrc
# Now, "(base)" should appear in the terminal before the user name
# Add Conda channels
conda config --show channels
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
conda config --add channels anaconda
conda config --add channels r
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
# Set libmamba as enviroment solver. It is faster than the classic solver
conda config --set solver libmamba
# Update Miniconda
conda update conda -y # If the message "Terms of Service have not been accepted" appears, execute the command lines shown in the message and try this command again.
# Show Conda configuration
conda config --show-sources

# Directory for databases
# Create the directory "db" for the databases in the disk root
sudo mkdir /db
# OR create a softlink "/db" for the databases directory located in anoter place
cd /
sudo ln -s /path/to/real/database/directory/db db
# Change the directory owner and group to the current user. This is required to download databases later.
sudo chown $USER /db
sudo chgrp $USER /db

##########################################################################
# Miniconda (For all users. Server computer.)
# Log in as root
sudo su
cd
# Download the installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x87_64.sh
# Follow the instructions and install Miniconda in /usr/local/miniconda3

# Do you accept the license terms? [yes|no]
# >>> yes

# Miniconda3 will now be installed into this location:
# /root/miniconda3
#
#   - Press ENTER to confirm the location
#   - Press CTRL-C to abort the installation
#   - Or specify a different location below
#
# [/root/miniconda3] >>> /usr/local/miniconda3

# Do you wish to update your shell profile to automatically initialize conda?
# This will activate conda on startup and change the command prompt when activated.
# If you'd prefer that conda's base environment not be activated on startup,
#    run the following command when conda is activated:

# conda config --set auto_activate_base false

# You can undo this by running `conda init --reverse $SHELL`? [yes|no]
# [no] >>> yes

# Apply changes in the current session
source ~/.bashrc
# Now, "(base)" should appear in the terminal before "root¨
# Add Conda channels
conda config --show channels
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels defaults
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
# Set libmamba as enviroment solver. It is faster than the classic solver
conda config --set solver libmamba
# Update Miniconda
conda update conda -y
# Show Conda configuration
conda config --show-sources
# Create a simbolic link for all users
ln -s /usr/local/miniconda3/bin/conda /usr/bin/conda
# Log out as root
exit

# For each user, execute this command only once to add Conda.
conda init # To write Conda information to ~/.bashrc
# Apply changes in the current session
source ~/.bashrc
# Now, "(base)" should appear in the terminal before the username.

# ⚠️ IMPORTANT! Log in as root before creating a Conda environment to make it available to all users.

# Directory for databases
# Create the directory "db" for the databases in the disk root
sudo mkdir /db
# OR create a softlink "/db" for the databases directory located in anoter place
cd /
sudo ln -s /path/to/real/database/directory/db db

# After download a database in /db/, give all users the permission to access the directories
# Give to all users the permission to read databases directories (+x is necessary for directories) 
sudo find /db/ -type d -exec chmod a+rx {} \;
# Give to all users the permission to read databases files
sudo find /db/ -type f -exec chmod a+r {} \;
# OR give to all users the permission to read databases directories and files
sudo chmod -R o+rx /db

##########################################################################
## B.2) Software installation - Metagenomics
##########################################################################

# Check the software manuals for database installation

# ⚠️ IMPORTANT! If you installed Miniconda as root, log in as root before creating a Conda environment to make it available to all users.
# To log in as root you should execute the command line below and type the password:
# sudo su
# To leave, execute the command line below
# exit
# After download a database in /db/, give all users the permission to access the directories
# Give to all users the permission to read databases directories (+x is necessary for directories) 
sudo find /db/ -type d -exec chmod a+rx {} \;
# Give to all users the permission to read databases files
sudo find /db/ -type f -exec chmod a+r {} \;
# OR give to all users the permission to read databases directories and files
sudo chmod -R o+rx /db

############################################################
# Essential software
sudo apt-get install bioperl -y
sudo apt-get install build-essential -y
sudo apt-get install cmake -y
sudo apt-get install git -y
sudo apt-get install glibc-source -y
sudo apt-get install gzip -y
sudo apt-get install libgd-dev -y
sudo apt-get install libgsl-dev -y
sudo apt-get install htop -y
sudo apt-get install moreutils -y
sudo apt-get install pigz -y
sudo apt-get install rename -y
sudo apt-get install screen -y
sudo apt-get install zip -y
sudo apt-get install nvidia-cuda-toolkit -y # For NVDIA GPUs
sudo ubuntu-drivers autoinstall -y
# OpenSSH to access computer remotely:
sudo apt install openssh-server openssh-client # To access computer remotely
sudo ufw allow ssh
sudo ufw enable
sudo ufw status
sudo echo "ClientAliveInterval 60" >> /etc/ssh/sshd_config
sudo systemctl restart ssh.service

############################################################
# AMRFinderPlus (tested version: 4.0.23)
# https://github.com/ncbi/amr
conda create -c conda-forge -c bioconda -n amrfinder --strict-channel-priority ncbi-amrfinderplus -y
conda activate amrfinder
mkdir /db/amrfinder
amrfinder_update -d /db/amrfinder
conda deactivate

############################################################
# Aragorn
# https://github.com/morloclib/aragorn
conda create -n aragorn -c bioconda aragorn -y

############################################################
# Bwa-mem2 (tested version: Bwa-mem2-v2.3)
# https://github.com/bwa-mem2/bwa-mem2
conda create -n bwa-mem2 -c bioconda -c conda-forge bwa-mem2 samtools -y

############################################################
# Barrnap
# https://github.com/tseemann/barrnap
conda create -n barrnap -c bioconda -c conda-forge barrnap gsl=2.5 -y

############################################################
# CheckM2 (tested version: 1.0.2, database uniref100.KO.1.dmnd)
# https://github.com/chklovski/CheckM2
cd
git clone --recursive https://github.com/chklovski/checkm2.git
cd checkm2
conda env create -n checkm2 -f checkm2.yml -y
conda activate checkm2
python setup.py install
checkm2 -h
cd ..
rm -r checkm2
# Database in standad path
# checkm2 database --download # install in /home/user/databases
# Database in custom path
mkdir -p /db/checkm2
checkm2 database --download --path /db/checkm2
conda env config vars set CHECKM2DB="/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
conda deactivate
conda activate checkm2
echo $CHECKM2DB
conda deactivate

############################################################
# DeepGOPlus (tested version: 1.0.2)
# https://github.com/bio-ontology-research-group/deepgoplus
git clone https://github.com/bio-ontology-research-group/deepgoplus.git
cd deepgoplus
conda create -n deepgoplus -c bioconda python=3.7 diamond -y
conda activate deepgoplus
pip install -r requirements.txt
pip install deepgoplus
cd ..
mkdir -p /db/deepgoplus
wget http://deepgoplus.bio2vec.net/data/data.tar.gz -P /db/deepgoplus
tar -zxvf /db/deepgoplus/data.tar.gz -C /db/deepgoplus/
mv /db/deepgoplus/metadata /db/deepgoplus/data/
rm /db/deepgoplus/data.tar.gz

############################################################
# dbCAN (tested version: 5.2.1)
# https://github.com/bcb-unl/run_dbcan
conda create -n dbcan -c bioconda -c conda-forge dbcan -y
mkdir /db/dbcan
# wget https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/CAZy.dmnd -P /db/dbcan
# wget https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN_sub.hmm -P /db/dbcan
conda activate dbcan
run_dbcan database --db_dir /db/dbcan
conda activate dbcan
cd /db/dbcan
for file in *.hmm; do
	hmmpress $file
done

############################################################
# eggNOG-mapper (tested version: 2.1.9, database v5.0.2)
# https://github.com/eggnogdb/eggnog-mapper
conda create -n eggnog-mapper -c bioconda eggnog-mapper -y
conda activate eggnog-mapper
conda env config vars set EGGNOG_DATA_DIR="/db/eggnog/"
conda deactivate
conda activate eggnog-mapper
echo $EGGNOG_DATA_DIR
conda deactivate
# Database directory
mkdir /db/eggnog
cd /db/eggnog
# Download main annotation database (v5.0.2, 39 GB)
wget -c --tries=0 http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
gunzip eggnog.db.gz
# Download taxa database (v5.0.2, 272 MB)
wget -c --tries=0 http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
tar -xzvf eggnog.taxa.tar.gz
rm eggnog.taxa.tar.gz
# Download pfam database (v5.0.2, 2.8 GB)
wget -c --tries=0 http://eggnog5.embl.de/download/emapperdb-5.0.2/pfam.tar.gz
tar -xzvf pfam.tar.gz
rm pfam.tar.gz
# Download diamond database (v5.0.2, 8.7 GB)
wget -c --tries=0 http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
gunzip eggnog_proteins.dmnd.gz
# # Download MMseqs2 database (v5.0.2, 11 GB)
# wget -c --tries=0 http://eggnog5.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz
# tar -xzvf mmseqs.tar.gz
# rm mmseqs.tar.gz
# # Index the MMseqs2 database (v5.0.2, 76 GB)
# mmseqs createindex mmseqs/mmseqs.db tmp
# rm -r tmp
# # Pad the MMseqs2 database for GPU support
# mmseqs makepaddedseqdb /db/eggnog/mmseqs/mmseqs.db /db/eggnog/mmseqs/mmseqs_padded.db
# #mmseqs easy-search query.fasta /db/eggnog/mmseqs/mmseqs_padded.db result.txt tmp --gpu 1

############################################################
# FastQC (tested version: 0.12.1)
# https://github.com/s-andrews/FastQC
conda create -n fastqc -c bioconda fastqc -y

############################################################
# Fastp (tested version: 1.0.1)
# https://github.com/OpenGene/fastp
conda create -n fastp -c bioconda fastp -y

############################################################
# GTDB-Tk (tested version: 2.4.1, database r226)
# https://github.com/Ecogenomics/GTDBTk
mkdir -p /db/gtdbtk/
cd /db
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
#In case it fails you can continue the download with the command below
#wget -c --progress=bar https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
echo '24b476ea5a4ef30519d461e56cc4a27f gtdbtk_data.tar.gz' | md5sum -c #md5sum for version r226
tar -xvzf gtdbtk_data.tar.gz -C "/db/gtdbtk" --strip 1 > /dev/null
conda create -n gtdbtk -c bioconda python=3.10 gtdbtk=2.4.1 -y
conda activate gtdbtk
conda env config vars set GTDBTK_DATA_PATH="/db/gtdbtk"
conda deactivate
conda activate gtdbtk
echo $GTDBTK_DATA_PATH
gtdbtk check_install
conda deactivate

############################################################
# GUNC (tested version: 1.0.6, database v2.1)
# https://github.com/grp-bork/gunc
mkdir /db/gunc
conda create -n gunc -c bioconda gunc -y
conda activate gunc
conda install setuptools -y
gunc download_db /db/gunc/
conda env config vars set GUNC_DB=$(find /db/gunc/ -type f -name "*.dmnd" -print -quit)
conda deactivate
conda activate gunc
echo $GUNC_DB
conda deactivate

############################################################
# Kraken 2 (tested version: 2.1.6, database k2_pluspf_20251015)
# https://github.com/DerrickWood/kraken2
# https://benlangmead.github.io/aws-indexes/k2
conda create -n kraken2 -c conda-forge -c bioconda kraken2 krakentools bracken krona r bowtie2 samtools -y
# Database PlusPF-16 (14.9GB)
cd /db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_20251015.tar.gz
mkdir -p /db/kraken2/k2_pluspf_16_GB_20251015
tar -xzf k2_pluspf_16_GB_20251015.tar.gz -C kraken2/k2_pluspf_16_GB_20251015
# rm k2_pluspf_16_GB_20251015.tar.gz
cd
# Database PlusPF (100. 6GB)
cd /db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20251015.tar.gz
mkdir -p kraken2/k2_pluspf_20251015
tar -xzf /k2_pluspf_20251015.tar.gz -C kraken2/k2_pluspf_20251015
# rm /db/k2_pluspf_20251015.tar.gz
cd

############################################################
# MEGAHIT (tested version: 1.2.9)
# https://github.com/voutcn/megahit
conda create -n megahit -c bioconda megahit -y

############################################################
# MetaPhlAn (tested version: 4.2.4, database mpa_vJan25_CHOCOPhlAnSGB_202503)
# https://github.com/biobakery/MetaPhlAn
conda create -n metaphlan -c conda-forge -c bioconda metaphlan matplotlib=3.2.2 -y
conda activate metaphlan
mkdir /db/metaphlan/
metaphlan --install --db_dir /db/metaphlan/
# conda env config vars set METAPHLAN_DB_DIR=/db/metaphlan
conda deactivate
conda create -n graphlan -c biobakery graphlan -y

##########################################################################
# MOB-suite (Plasmid identification)
conda create -n mob_suite -c bioconda mob_suite -y

############################################################
# MultiQC (tested version: 1.23)
# https://github.com/MultiQC/MultiQC
conda create -n multiqc -c bioconda multiqc -y

############################################################
# NCBI Datasets (tested version: 18.0.2)
# https://github.com/ncbi/datasets
conda create -n datasets ncbi-datasets-cli -y

##########################################################################
# Prokka (tested version: 1.14.6)
conda create -n prokka -c conda-forge prokka -y

############################################################
# Pyrodigal (tested version: 3.6.3)
# https://github.com/althonos/pyrodigal
conda create -n pyrodigal -c bioconda pyrodigal -y

############################################################
# QUAST (tested version: 5.2.0)
# https://github.com/ablab/quast
conda create -n quast -c bioconda quast -y

############################################################
# SemiBin (tested version: 2.2.0)
# https://github.com/BigDataBiology/SemiBin
# No support to GPU
conda create -n semibin -c conda-forge -c bioconda semibin -y
# OR with support to GPU (tested GPU: RTX 3060 12GB on Ubuntu 24.04.3)
sudo apt-get purge -y 'nvidia-*'
sudo apt-get autoremove -y
ubuntu-drivers devices
sudo apt-get install nvidia-driver-580 -y
sudo reboot now
conda create -n semibin python=3.13 pip -y # SemiBin 2.2.0 supports Python 3.7-3.13
conda activate semibin
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available(), torch.cuda.device_count()); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'no GPU')"
conda install -c conda-forge -c bioconda semibin=2.2.0 -y
conda deactivate

##########################################################################
# SRA Tools (tested version: 3.2.1)
conda create -n sra-tools -c bioconda sra-tools -y

############################################################
# VIBRANT (tested version: 1.2.1)
# https://github.com/AnantharamanLab/VIBRANT
conda create -n vibrant -c bioconda vibrant -y
conda activate vibrant
mkdir /db/vibrant
download-db.sh /db/vibrant
conda deactivate

##########################################################################
## C) Connecting to a server and using Screen
##########################################################################

# Connection to server
ssh user_name@server_ip

# Create screen
screen -S screen_name

# Leave screen without closing it (Detach)
CTRL+a+d

# List screens
screen -ls

# Retrive a screen (Attach)
screen -rd screen_name

# Close a screen (When attached to it)
exit

# Kill a screen (When a process inside the screen freezes)
screen -ls #To list screens
screen -XS screen_name kill