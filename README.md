# FEMGravSA
Fortran-based code to invert gravity anomaly data using the simulated annealing method. 
The forward calculation of the synthetic gravity is done using a technique based on 
the finite-element method. 

Compilation:
To compile the code using a command-line prompt - "make FemGrav".

Data:
The data folder contains input files needed to prepare the input data file for 
FemGrav. The input files include raw gravity data available from (http://icgem.gfz-potsdam.de/home). 
The open-source Python library of Fatiando e terra is also needed, and is avaliable at:
https://github.com/fatiando/fatiando
