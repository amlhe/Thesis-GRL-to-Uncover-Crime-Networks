# Using-Graph-Representation-Learning-to-Uncover-Structure-in-Crime-Networks
GitHub repository containing all code used in the work process of completing the thesis "Using Graph Representation Learning to Uncover Structure in Crime Networks" at the Technical University of Denmark. 
The folder "Original datasets" contains the files of the network data as they have been downloaded directly from their respective sources. Within the SBM folder, three different R scripts can be found for the pre-processing of the datasets for further analyses. 
The structure of files and folders are kept in such a way, that each script should call the files used in their proper directories. 

### ESBM
Within the SBM folder, the implementation and scripts used for the different SBM models with a variety of priors are included in the Source sub-folder. Much of the esbm.R file stems from the code supplied by "ESBM: Extended stochastic block models" from danieldurante on GitHub. A licence for the free use of and allowance for editing can be found within the SBM folder. Extensions have been made to the code for the SBM on multilayer networks and for performing link prediction and subsequent performance assessment. 

### HM-LDM
The code builds upon the work made in "Hybrid-Membership Latent Distance Model (HM-LDM)" available on Nicknakis Github. It uses the functions and methods implemented by Nicknakis, but the evaluation, assessment and outputs have been modified and changed to fit the purpose of the thesis with the considered network datasets. 
