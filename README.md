# Ribozyme RNA-seq

This repository contains the scripts used to analyze sequencing files to generate results from RNA-seq experiments for measuring the activity of ribozymes and ribozyme switches in mammalian cell transient transfection assays. The scripts need to be modified with user input for sequencing barcode information and sequencing run names, and only represent a semi-automated workflow. 

1. Process FASTQ files
After downloading gzipped fastq data files, which may have been demultiplexed, paired-end sequencing reads are joined using default parameters in PEAR. The version of pear that is used is v0.9.10 [May 30, 2016]. Pear can be downloaded from https://sco.h-its.org/exelixis/web/software/pear/
The workflow from paired-end joining to custom barcode demultiplexing are implemented in the jupyter notebook processFASTQ.
