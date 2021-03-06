# ZWA

The Zero Waste Algorithm (ZWA) is a pipeline comprised of various bioinformatic tools for context-based read trimming. The goal of this pipeline is to assist in viral genome assembly by retaining the viral part of reads and discarding any ribosomal hybrid bases the reads might contain. 



_DEPENDECIES:_

- [BIOAWK](https://github.com/lh3/bioawk)

- [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

- [BBMAP reformat.sh](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)

- [faSOMERECORDS](https://github.com/santiagosnchez/faSomeRecords)

- [bwa](http://bio-bwa.sourceforge.net/)

- [samtools](http://www.htslib.org/)

- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)



### Before running
Make sure you edit the script to specify paths and system related variables


### Typical run

` ./zwa.sh RAW_READS.fastq reference.fasta `




_!!!**LIMITATIONS**!!!_

> 1)This is an unoptimized proof of concept pipeline. It is quite slow

> 2)If your raw reads contain the letter 'S' in the contig name you should replace it before running the pipeline (eg. use sed)

> 3)This script has been tested on single-end raw reads only. Nothing suggests that it will not work for paired-end reads after pairing them.


For ease of use and to maximize performance, a Docker image with all the required software preinstalled has been created at https://hub.docker.com/repository/docker/ndovroli/zwa/. It is ready to deploy on any system/cloud with Docker.


