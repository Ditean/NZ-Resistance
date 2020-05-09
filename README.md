# Drug Resistance Profiler
New Zealand Tuberculosis Drug Resistance Pipeline

This pipeline is intended to be used within New Zealand for clinical purposes.

This pipeline will perform variant calling on raw sequences, and identify drug-conferring SNPs. Processed samples will be automatically stored within a repository, as well as produce a clinical report on drug-susceptibility.

## Getting started
The resistance pipeline is designed to be flexible for a range of situations. Primarily, the pipeline has been designed so that clinicians will have a user-friendly experience.

### ***Running sequences***
There are two methods for the pipeline to detect resistance mutations

**Placing sequences within the input directory**\
If no sequences are provided, the pipeline will search the input directory for sequences.
```
resistance.sh
```
**Providing sequencing files**\
The user can specify forward and reverse reads for processing
```
resistance.sh --forward /sequences/sample_1_R1.fastq.gz --reverse /sequences/sample_1_R1.fastq.gz
```

### Prerequisites

Software needed to run
* **Python 3.7.1**
* **Fastq-mcf** *v1.04.807*
* **BWA** *v0.7.17-r1188*
* **Samtools** *v1.9* (htslib 1.9)
* **Bcftools** *v1.9* (htslib 1.9)
* **Freebayes** *v1.3.1-dirty*

### Installing

It is recommended that you perform a dry-run when setting up for the first time:

**Dry Run**
```
resistance.sh --setup
```
After completing the dry-run, the pipeline will exit and be ready for use. The set-up log files can also be consulted for issues with the setup.

Alternatively, the pipeline can detect if the pipeline has not been set-up, and will configure the environment before running.
## Development

## Built with

## Authors
* **Jordan Taylor** - *Bioinformatician* - [University of Otago](https://micro.otago.ac.nz/our-people/other-research-staff/tom-devine-2/)

## Acknowledgements

* Person 1
* Person 2
* Reason 1
