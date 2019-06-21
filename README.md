# ConsHMM
ConsHMM provides tools for parsing a multiple species alignment and training a Hidden Markov Model (HMM) to learn a conservation state annotation of the reference genome in the alignment, at single nucleotide resolution. The HMM is learned using an updated version of the [ChromHMM](http://compbio.mit.edu/ChromHMM/) software, which is included in this repository. Tools for visualizing and interpreting ConsHMM output are also provided. 

## Requirements
1. `Python 3`
2. `Java 1.5 or later`
2. [Biopython](http://biopython.org/wiki/Download)
3. [Numpy](http://www.numpy.org/)

If you are in a conda environment, the following lines will install the necessary python libraries
```
conda install -c conda-forge biopython
conda install -c anaconda numpy
```

## Getting started
The segmentation and browser files mentioned in the paper are available at https://ernst.cass.idre.ucla.edu/public/ConsHMM/. Model files and segmentations for hg38 are also available. **Note that the hg38 states are based on a different multiple sequence alignment and have not been annotated yet**. The link provides the intermediate files produced by the pipeline using the hg19 Multiz 100-way alignment.

The [Wiki](https://github.com/ernstlab/ConsHMM/wiki) contains useful tutorials, including how to reproduce the model and segmentation from the original ConsHMM paper or create your own based on a different reference species and/or multiple-sequence alignment.

## Citation

For any use of the ConsHMM software or ConsHMM state annotations, please cite:

Sperlea A, Ernst J. 2018. Systematic Discovery of Conservation States for Single-Nucleotide Annotation of the Human Genome. bioRxiv doi: https://doi.org/10.1101/262097

## Authors

Adriana Sperlea (University of California, Los Angeles)

Jason Ernst (University of California, Los Angeles)

## Collaborators

Bruins In Genomics students Brooke Felsheim (Washington University in St. Louis) and Jennifer Chien (Wellesley College) helped test the pipeline during the summer of 2018 and implemented several additional features.
