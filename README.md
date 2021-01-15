# ZWA

The Zero Waste Algorithm (ZWA) is a pipeline comprised of various bioinformatic tools for context-based read trimming. The goal of this pipeline is to assist in viral genome assembly by retaining hte viral part of reads and discarding any ribosomal hybrid bases the reads might contain.


!!!LIMITATIONS!!!
1)This is an unoptimized proof of concept pipeline. It is quite slow

2)The algorithm is unable to "clean" the hybrid ribosomal part if its in the middle of the read

3)If your raw reads contain the letter 'S' in the contig name you should replace it before running the pipeline (eg. use sed) 
