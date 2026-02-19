GBS-Genomics-Pipeline: User Manual and Complete Guide

GBS-Genomics-Pipeline is an open-source, modular, and containerized workflow for whole genome analysis of Streptococcus agalactiae (Group B Streptococcus).
It is built with Nextflow and designed to run on laptops, workstations, servers, and HPC clusters using Docker, Singularity, or Conda.

The pipeline automates the full genomic workflow from raw sequencing reads or curated genome assemblies to phylogenetic tree construction and integrated visualization.

1. What This Pipeline Does
This workflow performs complete bacterial genome analysis in an automated and reproducible way.
It performs:

• Quality control of raw reads
• Removal of human contamination
• Genome assembly
• Assembly quality assessment
• Taxonomic confirmation
• Background genome selection
• Functional annotation
• Virulence factor detection
• MLST typing
• Core genome analysis
• Phylogenetic tree reconstruction
• Integrated AMR and virulence visualization

The workflow has two input pathways:
Curated Pathway
For assembled genome FASTA files (.fa, .fna, .fasta)

Raw Pathway
For paired-end FASTQ files (_1.fastq and _2.fastq)
Note: The raw read pathway is currently under scientific validation and optimization.

2. Why This Pipeline Uses Containers
Bioinformatics tools are difficult to install and often conflict with each other.
This pipeline is fully containerized using the Docker image:
kizitodevbio/strepto-pipeline:latest

This means:

• No manual installation of SPAdes, Prokka, BLAST, etc.
• Identical results on any system
• Clean and reproducible execution

3. What You Need Before Running

You must install:

Nextflow (workflow manager)

One of the following:
Docker (recommended)
Singularity (for HPC systems)
Conda (alternative method)

You only need ONE container system.

4. Install Nextflow

Run this on Linux or macOS:
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/


Check installation:
nextflow -version

5. Install Docker (Recommended)
Ubuntu:
sudo apt-get update
sudo apt-get install docker.io
sudo systemctl start docker
sudo systemctl enable docker


Verify Docker:
docker --version

Add your user to Docker group (so you do not need sudo):
sudo usermod -aG docker $USER

Then log out and log back in.

6. Pull the Docker Image Manually (Optional but Recommended)
You can manually download the container image:
docker pull kizitodevbio/strepto-pipeline:latest

If you use Singularity:
singularity pull docker://kizitodevbio/strepto-pipeline:latest

Nextflow will automatically pull the image if it is not already present.

7. Download the Pipeline
Move to the directory where you want to keep the project:
cd ~/Desktop

Clone the repository:
git clone https://github.com/kizito-devbio/GBS-Genomics-Pipeline.git


Enter the folder:
cd GBS-Genomics-Pipeline

You are now inside the pipeline directory.
8. Preparing Your Data
You must organize your data before running.

For curated genomes:
Create a folder:
mkdir curated_data

Place your FASTA files inside:
example:
sample1.fasta
sample2.fna

For raw reads:
Create a folder:
mkdir raw_data

Files must follow this naming format:
sample1_1.fastq
sample1_2.fastq

Paired-end naming is required.

9. Running the Pipeline

You must always run the command from inside the GBS-Genomics-Pipeline folder.

The parameters:
--curated_dir
or
--raw_dir

are simply names that tell the pipeline which pathway to use. They point to your folder.

Option A: Run with Docker (Recommended)
Curated genome pathway:
nextflow run pipeline.nf -profile docker --curated_dir curated_data --outdir results


Raw reads pathway:
nextflow run pipeline.nf -profile docker --raw_dir raw_data --outdir results

Option B: Run with Singularity (HPC environments)
Curated genome pathway:
nextflow run pipeline.nf -profile singularity --curated_dir curated_data --outdir results

Raw reads pathway:
nextflow run pipeline.nf -profile singularity --raw_dir raw_data --outdir results

Option C: Run with Conda(Recommended)
Conda builds the environment locally.

Curated genome pathway:
nextflow run pipeline.nf -profile conda --curated_dir curated_data --outdir results


Raw reads pathway:
nextflow run pipeline.nf -profile conda --raw_dir raw_data --outdir results

10. Resume an Interrupted Run
If your system shuts down or the job stops: by adding -resume as seen in the example below
nextflow run pipeline.nf -profile docker --curated_dir curated_data --outdir results -resume

Nextflow will continue from the last successful step.


12. Output Structure

After completion, the results folder contains:

Annotation: move into the fuctional annotation
fuctional annotation and AMR files

Virulence: move into virulence
Virulence factor reports

MLST: move into mlst
Sequence typing results

Core Genome: move into core gemone
Alignment files

Phylogeny: move into phylogeny
Newick tree files

Logs
Execution timeline
Performance report
Workflow DAG graph

Reports are automatically generated in:

results/logs/

12. System Resource Configuration
The pipeline automatically detects available CPUs and assigns resources safely.

Default behavior:

• Uses available CPUs minus one
• Uses up to 6 GB RAM unless modified

For cluster systems (SLURM), use:

nextflow run pipeline.nf -profile cluster --curated_dir curated_data --outdir results


You can adjust memory and CPU in the configuration file if needed.

13. Important Notes

• Raw read pathway is still under validation
• FASTQ files must be paired-end
• FASTA files must be properly formatted
• Always run from inside the project directory
• Docker abd conda is strongly recommended for reproducibility

14. Contact
For issues, suggestions, or collaboration:

Email: kizitosylvester@gmail.com
