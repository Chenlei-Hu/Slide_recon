# Slide_recon

Reconstruction for imaging-free spatial genomics.


## Getting Started

### Environments
```
conda env create -f environment.yml
```
activate the conda environment with:
```
conda activate recon
```

## Demo 
The concept of diffusion-based recontruction is implemented in '/demo' with example image and code.

## Reconstruction
For experimental data of reconstruction with Slide-seq or Slide-tags, the following two files are used to parse the sequencing results and reconstruct bead spatial locations.

### Generating diffusion matrix with fastq 
Need to change 'fq_dir' and 'out_dir' in 'fiducial_seq_blind_whitelist.py'.\

* Run in command-line
```shell
python fiducial_seq_blind_whitelist.py  -d DATE -s SAMPLE
```
* The usage of this command
If using V15 (including both V15A and V15T) for fiducial beads, need to define
```shell
usage: fiducial_seq_blind_whitelist.py -d DATE -s SAMPLE [-r2 READ2TYPE]

optional arguments:
  -d DATE, --date DATE  string, name of folder within which contains a '/fastq' folder
  -s SAMPLE, --sample SAMPLE
                        string, output folder name that will be created under the date folder
  -r2  READ2TYPE, --read2type READ2TYPE
                        string, optional, defaul='V9', type of fiducial beads used. Choose 'V15'
                        if V15A or V15T beads are used
```
The expected return includes:
* 'SAMPLE_QC.pdf' file showing the bead barcode rank plots and read distribution plots.
* 'SAMPLE_blind_raw_reads_filtered.csv.gz' file containing the diffusion matrix in sparse format.
* 'SAMPLE_blind_statistics_filtered.csv' file with read parsing statistics summary.


### Reconstruction with diffusion matrix
Need to change 'sample_folder' in 'reconstruction_blind.py'.\
Using multicore can largely reduce the running time.
* Run in command-line
```shell
python3 reconstruction_blind.py -d DATE -s SAMPLE -a CAPTURE_BEAD -t FIDUCIAL_BEAD -e EXPERIMENT
```
* The usage of this command
Slide-seq or Slide-tags use different flag.
```shell
usage: reconstruction_blind.py -d DATE -s SAMPLE -a CAPTURE_BEAD -t FIDUCIAL_BEAD -e EXPTYPE
                               [-c CORE]   

optional arguments:
  -d DATE, --date DATE  string, name of folder within which contains a '/fastq' folder
  -s SAMPLE, --sample SAMPLE
                        string, output folder name that will be created under the date folder
  -a  CAPTURE_BEAD, --anchor CAPTURE_BEAD
                        string, type of capture beads used. Only for naming files
  -t  FIDUCIAL_BEAD, --target FIDUCIAL_BEAD
                        string, type of fiducial beads used. Only for naming files
  -e EXPTYPE, --exptype EXPTYPE
                        string, 'seq' for Slide-seq and 'tags' for Slide-tags
  -c CORE, --core CORE
                        string, optional, default='CPU'. Choose 'GPU' if GPU is available
```
The expected return includes:
* Four png files showing the bead barcode UMI and covered bead distribution.
* 'CAPTURE_BEAD_recon_loc.csv' file containing reconstructed coordinates for capture beads.
* 'SAMPLE_UMAP.png' showing the locations. 'SAMPLE_UMAP_density.png' showing the density of reconstruction locations. 'SAMPLE_UMAP_convex.png' showing the shape.

## Notes
* I usually save the file with the structure: /.../date/. I will copy fastq file under /.../date/fastq. The output will be in /.../date/sample.
* 'date' is the folder name of one experiment. 'sample_name' is the name of the sample in the front of the fastq file. 'capture_bead_name' and 'fiducial_bead_name' is name of capture and fiducial beads, but they are just names for saving files and won't affect results.
* The code has been tested on Python 3.11.7
