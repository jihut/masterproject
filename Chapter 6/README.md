# This folder contains different files for the proposed CUSUM chart method. 

The `test em` folder contains the function for calculating the proposed chart based on the EM-based models (Perme et al.,[2012](https://pubmed.ncbi.nlm.nih.gov/21689081/)).
A couple of files to produce some of the results in the text are also included. 

The `test piecewise 1` folder contains the function for calculating the proposed chart based on the piecewise constant baseline models. 
Some examples to reproduce the results from the text are also given here. 

The `test piecewise 2` folder contains a specific example of a piecewise baseline used right before Chapter 6.4 in the text. If comments are needed here, we refer to the files in the `test piecewise 1` as the only differences are the choice of baseline parameters and follow-up partition.

NOTE! In the text, we also tested our methods against the threshold values obtained from the Markov Chain method (Gandy et al.,[2010](https://academic.oup.com/biomet/article/97/2/375/219017)). 
The file `Nevents.R` contains the functions to calculate the threshold $c$ using this method and is obtained from the authors. We do not own the copyrights of this file!  
