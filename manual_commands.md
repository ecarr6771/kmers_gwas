## Pipeline Setup & Installation
- Login to the HPCC & enter your NetID password

    `ssh carrele1@hpcc.msu.edu`

- `ssh` to a development node

    `ssh dev-amd20`

- Load conda for the HPCC login  

    `module load conda`

- Create a conda environment for this project  

    `conda create css844`

    - or activate the environment you've already created  
    
    `conda activate css844`
Note: I installed Python 2.7 in the `conda` environment before running any code, but did not install R or any other packages. 

- Make a directory for the project

    `/mnt/scratch/carrele1/kmersGWAS_dir`

- Download the source code from GitHub into a folder for the project

    `wget https://github.com/voichek/kmersGWAS/releases/download/v0.2-beta/v0_2_beta.zip`

- Unzip the file
    - You should get a bunch of folders and some files, including a makefile

        `unzip v0_2_beta.zip`

- Run `make` in the folder where the Makefile is

## Directory Setup
The following code is written into `gwas_run.sh`:    
- Create directories for each individual you will count the k-mers for
    
    `mkdir "${wd}"/"${INDIVIDUAL}"`

- Copy all (2) sequence files (the forward and reverse) into the directory for the accession

    `cp *"${INDIVIDUAL}"*.fastq.gz "${wd}"/"$INDIVIDUAL"/`  

- make a file `input_files.txt` inside the directory with the absolute paths to each sequence file in the folder

    `echo "${wd}"/"${str2}" >> "${wd}"/"${INDIVIDUAL}/input_files.txt"`

## _k_-mer Table Formation
- Get the list of k-mers using `kmc` 
    - Do the canonized version first (with filtering)

    `/mnt/scratch/carrele1/kmersGWAS_dir/external_programs/kmc_v3 -t2 -k31 -ci2 @input_files.txt "$INDIVIDUAL"_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2`

    - Then get the list of all the k-mers

    `/mnt/scratch/carrele1/kmersGWAS_dir/external_programs/kmc_v3 -t2 -k31 -ci0 -b @input_files.txt "$INDIVIDUAL"_kmc_all ./ 1> kmc_all.1 2> kmc_all.2`

- Combine the 2 lists 

    `/mnt/scratch/carrele1/kmersGWAS_dir/bin/kmers_add_strand_information -c output_kmc_canon -n "$INDIVIDUAL"_kmc_all -k 31 -o kmers_with_strand`

- Put the kmer lists out in textual format

    `./kmc_tools transform "$INDIVIDUAL" dump "$INDIVIDUAL".txt` 
    
- Remove .kmc files because they're really big
    
    `rm *.kmc*`

- Make a file with the absolute paths of the kmers_with_strand files (one for each individual)

    `echo -e "$(pwd)/$INDIVIDUAL/kmers_with_strand\t$INDIVIDUAL" >> "kmers_list_paths.txt""`

- Combine and filter the k-mer lists

    `/mnt/scratch/carrele1/kmersGWAS_dir/bin/list_kmers_found_in_multiple_samples -l kmers_list_paths.txt -k 31 --mac 5 -p 0.2 -o kmers_to_use`

    This is where the flag output happens, so if you want to save it somewhere, put in `1> flag_stats.txt` on the end of the above command

- You may have to install additional packages at this point  

    `conda install libgcc=5.2.0`

- Make the k-mer table

    `/mnt/scratch/carrele1/kmersGWAS_dir/bin/build_kmers_table -l kmers_list_paths.txt -k 31 -a kmers_to_use -o kmers_table`

- Put the k-mer table out into a readable format
    - Note: still working on where kmers_list.txt is generated, but I think it's from the `bowtie2` output.

    `/mnt/scratch/carrele1/kmersGWAS_dir/bin/filter_kmers -t kmers_table -k kmers_list.txt -o output.txt`

- Generate a kinship matrix:  

    `/mnt/scratch/carrele1/kmersGWAS_dir/bin/emma_kinship_kmers -t kmers_table -k 31 --maf 0.05 > kmers_table.xqkinship`

- Convert the kmer table to PLINK binary format

    `/mnt/scratch/carrele1/kmersGWAS_dir/bin/kmers_table_to_bed -t /mnt/scratch/carrele1/kmersGWAS_dir/kmers_count_stats/kmers_table -k 31 -p phenotype.pheno --maf 0.05 --mac 11 -b 10000000 -o plink_table/`

- Tabular output  

    `plink --bfile plink_table --recode vcf --allow-no-sex --out kmertable`

- Setup for the SLURM job:

   ```bash
    #SBATCH --time=24:00:00
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=2
    #SBATCH --mem-per-cpu=1G
    #SBATCH --job-name make_kmer_table
    #SBATCH --nodelist=amr-048
   ```
- Run the GWAS

    `/mnt/scratch/carrele1/kmersGWAS_dir/kmers_gwas.py --pheno phenotype.pheno --kmers_table kmers_table -l 31 -p 8 --outdir output_dir`

- Convert PLINK binary files to .`vcf` format

    `plink --bfile plink_table."${NUMBER}" --recode vcf --out kmertable."${NUMBER}"`

- Compress the `.vcf` files

    `bgzip vcftable/kmertable."${NUMBER}".vcf`

- Index the `.vcf` files

    `bcftools index vcftable/kmertable."${NUMBER}".vcf.gz`

- Merge the resulting .vcf files

    `bcftools concat -f vcf-list.txt`

## Downstream Analyses
- Once we run the GWAS, we'll need to align the kmers to the reference genome to generate a Manhattan plot

    `bowtie2 -p 2 -x ref.txt -f /mnt/home/carrele1/css844/Period_kmers_list.fa -S potato_output.sam 2> /mnt/home/carrele1/css844/bowtie2_log.txt`


for (i in seq_along(phenos$accession_id)) {
  if (!grepl('_', phenos$accession_id[i], fixed = TRUE)) {
    phenos$accession_id[i] <- paste(phenos$accession_id[i], "_", phenos$accession_id[i], sep = "")
  }
}

phenos$accession_id <- as.factor(phenos$accession_id)

write.table (phenos, "phenotype.pheno", sep = "\t", row.names = FALSE)