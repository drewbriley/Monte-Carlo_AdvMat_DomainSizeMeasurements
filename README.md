# Monte-Carlo_DomainSizeMeasurements
This is the code used to create the simulations for the domainsize measurement technique.

The file consists of the src and included files in the relevant folders alongside the parameters folder which contains a number of exmaple files that this program will run on. Note that these are the ones used on the Sunbird computer, to increase the speed of the computation you can lower the "lowerExcitonLimit" value as well as make the "domainSizeStart" and "domainSizeStop" values the same.  

The make file can be called with make -f Makefile rho=21-22, where the rho=21-22 is appended to the .exe file. 

The program can then be run with ./3D-MV-DomainSizeExp_21-22 parameters/DomainSizeMeasurements_21-22.txt
