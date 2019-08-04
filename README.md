# RNA_Seq-Analysis-workflow

RNA-Seq analyis using STAR two pass mode for alliging the raw reads and RSEM for transcript quantification using snakemake workflow. 

Required Tools

1 - FastQC (Quality control)

2 - Trim-galore ( Adapter trimming)

3 - STAR (Allign the reads to reference genome) 

4 - Picard (Remove Duplicates)

5 - RNA-SeQC(Qualiy Control of RNA-Seq data)

6 - RSEM(Transcript Quantification)

# Setting upconda environment for tools and their dependencies

Install anaconda or load it if its already on your server

conda create --name rnaseq-env

source activate rnaseq-env

conda install -c bioconda star

conda install -c bioconda fastqc

conda install -c bioconda rsem


# generate the .txt report of all fastqc "run fastqc-summary script in fastqc output directory"

python3 fastqc-summary -s $INDIR > "QC_Report.txt"

# Run the pipeline on cluster using this command 'modify parameters according to your cluster configuration

snakemake -j 999 --configfile config.yaml --use-conda --nolock --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition}  -N {cluster.N} -n {cluster.n}  -t {cluster.time} --mem {cluster.mem}"
