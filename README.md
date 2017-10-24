# ConsHMM
ConsHMM provides Python code for parsing a multiple species alignment into the necessary format for training a Hidden Markov Model(HMM) to learn conservation states. The HMM is learned using an updated version of the [ChromHMM](http://compbio.mit.edu/ChromHMM/) software, which is included in this repository. 

## Requirements
1. `Python (2.7 <= version < 3)`
2. `Java 1.5 or later`
2. [Biopython](http://biopython.org/wiki/Download)
3. [Numpy](http://www.numpy.org/)

## Getting started
The segmentation and browser files mentioned in the paper are available at https://ernst.cass.idre.ucla.edu/public/ConsHMM/.
The [Wiki](https://github.com/ernstlab/ConsHMM/wiki) contains more information about the code behind the model and segmentation, as well as tutorials for how to learn your own model and create a segmentation from a multiple sequence alignment.

