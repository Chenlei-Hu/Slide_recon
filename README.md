# Slide_recon

Reconstruction for imaging-free spatial transcriptomics


## Getting Started

### Environments
```
conda env create -f environment.yml
```
activate the conda environment with
```
conda activate recon
```

### Generating diffusion matrix with fastq 
Need to change 'fq_dir' and 'out_dir' in 'fiducial_seq_blind_whitelist.py'
Also change information in 'qsub_seq_blind.sh'

* I usually do it on UGER. On broad cluster
```
use UGER
conda activate recon
qsub qsub_seq_blind.sh
```

* Can also excute the python file
```
python fiducial_seq_blind_whitelist.py  -d <date> -s <sample_name>
```
If using V15 (including both V15A and V15TSO) for fiducial beads, need to define
```
python fiducial_seq_blind_whitelist.py  -d <date> -s <sample_name> -r2 V15
```

### Reconstruction with diffusion matrix
Need to change 'sample_folder' in 'reconstruction_blind.py'
I suggest using multicore to excute
* If doing reconstruction with Slide-seq
```
python3 reconstruction_blind.py -d <date> -s <sample_name> -a <capture_bead_name> -t <fiducial_bead_name> -e seq
```

* If doing reconstruction with Slide-tags
```
python3 reconstruction_blind.py -d <date> -s <sample_name> -a <capture_bead_name> -t <fiducial_bead_name> -e tags
```

* Default is runing with CPU, but if have GPU
```
python3 reconstruction_blind.py -d <date> -s <sample_name> -a <capture_bead_name> -t <fiducial_bead_name> -e <seq_or_tags> -c GPU
```

## Notes
* I usually save the file with the structure: /.../date/. I will copy fastq file under /.../date/fastq. The output will be in /.../date/sample.
* 'date' is the folder name of one experiment. 'sample_name' is the name of the sample in the front of the fastq file. 'capture_bead_name' and 'fiducial_bead_name' is name of capture and fiducial beads, but they are just names for saving files and won't affect results. 
