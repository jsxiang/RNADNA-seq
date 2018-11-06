# Ribozyme RNA-seq

This repository contains the scripts used to analyze demultiplexed sequencing file in the form of a .csv file containing sequences in the first column, and read counts in subsequent columns from each indexed run, to generate results from RNA-seq experiments for measuring the activity of ribozymes and ribozyme switches in mammalian cell transient transfection assays. 

Example input file provided is RSIX_combinedcounts.csv. 
Run analyzeRSIX.m as a demo. This MatLab script calls RNAseq.m, which is a class containing methods to create a RNAseq object with computed parameters:

             seqs: 1xm cell array of unique sequences
          RDratio: mx2 vector of RNA/DNA ratios, 1st column in no ligand condition, and 2nd column in with ligand condition; values are mean of replicates
    ctrls_RDratio: 5x2 vector of RDratio for control ribozymes
      val_RDratio: RDratio for validation candidates
            minus: above fields, for two replicates, in no ligand condition
             plus: above fields, for two replicates, in with ligand condition
           isdiff: binary vector results from unpaired t-test, for significance, for specified alpha, and fold change >1.5
       foldchange: mx1 vector for fold change of mean RNA/DNA in the different ligand conditions
                H: binary vector results from unpaired t-test, for significance, for specified alpha
                P: p-value results from unpaired t-test, for significance, for specified alpha
      ctrlsisdiff: binary vector results from unpaired t-test, for significance, for specified alpha, and fold change >1.5, for controls
      totalcounts: mx1 vector of total sequence read count in all the experiment conditions
       ctrlcounts: vector of total sequence read count in all the experiment conditions, for control ribozymes
         ctrlseqs: cell array of control ribozyme sequences
        valcounts: vector of total sequence read count in all the experiment conditions, for validation candidates
          valseqs: cell array of validation candidate sequences

To use, analyzeRSIX.m needs to be modified to indicate which index columns correspond to RNA and DNA read counts, which are replicates and which contain different ligand concentrations. Currently, only with and without ligand concentrations are allowed, and two replicate experiments are required. At least one DNA read count is needed for normalization of RNA read counts.
