# Ribozyme RNA-seq

This repository contains the scripts used to analyze demultiplexed sequencing file in the form of a .csv file containing sequences in the first column, and read counts in subsequent columns from each indexed run, to generate results from RNA-seq experiments for measuring the activity of ribozymes and ribozyme switches in mammalian cell transient transfection assays. 

Example input file provided is RSIX_combinedcounts.csv. 
Run analyzeRSIX.m as a demo. This MatLab script calls RNAseq.m, which is a class containing methods to create a RNAseq object with computed parameters, such as normalized RNA levels as RDratio, and RNA/DNA ratio of spiked-in control ribozymes as ctrls_RDratio. 
To use, analyzeRSIX.m needs to be modified to indicate which index columns correspond to RNA and DNA read counts, which are replicates and which contain different ligand concentrations. Currently, only with and without ligand concentrations are allowed, and two replicate experiments are required. At least one DNA read count is needed for normalization of RNA read counts.
