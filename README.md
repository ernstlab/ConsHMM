# ConsHMM

<<<<<<< HEAD
## H2
https://ernst.cass.idre.ucla.edu/public/ConsHMM/.

'''
createMAFFeatures.py takes in an alignment file in MAF format and extracts just the sequence information
creates a feature file in which
each line corresponds to a nucleotide in the human genome, and contains columns explaining whether
the other species have the same nucleotide at that position or not.
----The stuff below is not accurate----
The output format is:
pos,species1,species2,...,species99
pos = position in the human genome
speciesX = 0 if there is an aligned nucleotide in that particular species
         = 1 if there is no aligned nucleotide in that particular species
The order of the species is the one specified in the first line of the file.
=======
The segmentation and browser files mentioned in the paper are available at https://ernst.cass.idre.ucla.edu/public/ConsHMM/.
The Wiki contains more information about the code behind the model and segmentation, as well as tutorials for how to learn your own model and create a segmentation from a multiple sequence alignment.
>>>>>>> fb86868e1b2b3e2ebecba2951b1cf3f5ec2750a0
