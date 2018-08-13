# Amino Acid Typing by Sequence Probability for Protein NMR Backbone Assignments
aka 'seqprob' 
by Matt Stetz circa 2009

## Introduction
This is a modified version of the original seqprob approach described by Grzesiek and Bax in 1993 (JBNMR 1993). It is based on calculating a Guassian probability density function using mean and standard deviation values from the BMRB database

Amino Acid Typing is a necessary procedure in the assignment of NMR backbone resonances for applications in biophysics.

## Requirements
This program is run as a GUI. You must have the following libraries:
* wxPython

## Accuracy
The accuracy of the program was assessed using reference data for a high molecular weight protein (56 kDa) which was assigned manually. The results are shown in the image final accuracy.png. The accuracy is better than 90% if 6 sequential amino acids can be linked by inter-residue CA and CB connections
