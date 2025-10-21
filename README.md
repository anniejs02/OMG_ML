# CRC_CARAMeL
Welcome to CRC-CARAMeL: Colorectal Cancer Condition-specific Anticancer Regimen Assessment using Mechanistic Learning. This guide walks you through setting up and running CRC-CARAMeL in MATLAB. 
This tool allows for predicting the effect of multi-drug combinations on colorectal cancer (CRC) and the effect of microorganisms on treatment responses.

<img width="3000" height="2100" alt="Untitled (91)" src="https://github.com/user-attachments/assets/f41fd422-2d44-44ba-87d0-c2910f3ef078" />

This model can be used along with the included drug combination files (sourced from [SynergyxDB]([url](https://www.synergxdb.ca/))). Additional drug or microbial treatments can also be added to make novel predictions. 


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
