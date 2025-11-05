# Onco Microbiome GEM Machine Learning (OMG-ML)
Welcome to OMG-ML, a systems biology approach to predicting colon cancer treatments. This guide walks you through setting up and running OMG-ML in MATLAB. 
This tool allows for predicting the effect of multi-drug combinations on colorectal cancer (CRC) and the effect of microorganisms on treatment responses.


<img width="862" height="875" alt="image" src="https://github.com/user-attachments/assets/7a16019c-750f-400c-aa9e-bc555e6ccc27" />


This model can be used along with the included drug combination files (sourced from [SynergyxDB]([url](https://www.synergxdb.ca/))). Additional drug or microbial treatments can also be added to make novel predictions by adding log fold changes of gene expression in treatment compared to control to the included gene expression file. 


The basic approach of this model was adapted from CARAMeL (Condition-specific Antibiotic Regimen Assessment using Mechanistic Learning), which was used for antibioitic combination predictions (source:  Carolina H Chung, Sriram Chandrasekaran, A flux-based machine learning model to simulate the impact of pathogen metabolic heterogeneity on drug interactions, PNAS Nexus, Volume 1, Issue 3, July 2022, pgac132, https://doi.org/10.1093/pnasnexus/pgac132). 

# License 
Released via GPL GNU License
© 2025 The Regents of the University of Michigan
Chandrasekaran Research Group - https://systemsbiologylab.org/
Contact: csriram@umich.edu

# Dependencies 
Using CRC-CARAMeL requires the installation of the following packages/programs. 
1. [MATLAB]([url](https://www.mathworks.com/products/matlab.html)) 2019b or higher
2. [Gurobi]([url](https://www.gurobi.com/)) optimization solver
3. [COBRA Toolbox]([url](https://github.com/opencobra/cobratoolbox))

# Installation 
Install the folders and packages included by cloning this repository to your local computer

git clone https://github.com/anniejs02/CRC_CARAMeL.git

# Instructions 
See CRC_CARAMEL_main.m or CRC_CARAMEL_main.mlx for step-by-step instructions in applying CRC_CARAMeL. To run those files in MATLAB, you must download the contents of 'data' and 'empty_prediction_files' and have an appropriate human metabolic model. 
