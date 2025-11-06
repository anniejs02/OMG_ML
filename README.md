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

# Files 
CARAMeL
└───code: contains all relevant code files
│   └───CARAMeL_suite: contains all code files created to implement the CARAMeL approach
│   │   │   caramel.m:                      main script the constructs a CARAMeL model
│   │   │   caramel_anova.m:                conducts a one-way ANOVA to determine GEM reactions with differential flux activity 
│   │   │   caramel_assessPerf.m:           assesses the predictive performance of a CARAMeL model
│   │   │   caramel_auroc:                  calculates the area under the receiver operating curve (AUROC) for a CARAMeL model 
│   │   │   caramel_classify:               classifies drug interactions as synergistic, additive, or antagonistic
│   │   │   caramel_crossValidation.m:      conducts a k-fold cross-validation for a CARAMeL model
│   │   │   caramel_features2gem.m:         extracts GEM reactions associated with top CARAMeL model features
│   │   │   caramel_featurize.m:            implements featurization of raw input data into ML compatible format
│   │   │   caramel_leaveOut.m:             conducts a leave-out analysis for a CARAMeL model
│   │   │   caramel_plot.m:                 generates plots visualizing CARAMeL model performance
│   │   │   caramel_processInteractions.m:  processes drug interaction data into a standardized format
│   │   │   caramel_rankSubsystems.m:       determines metabolic pathways enriched by GEM reactions tied to top CARAMeL features
│   │   │   caramel_screen.m:               screens all possible drug combination outcomes given a list of drugs
│   │   │   caramel_topFeatures.m:          extracts the top CARAMeL model features 
│   └───GEM_functions: contains all code files relevant to simulating metabolism using GEMs
│   │   │   change_media.m:                 changes the simulated media condition for a given GEM
│   │   │   constrain_flux_regulation.m:    simulates reaction fluxes based on omics data constraints
│   │   │   derive_flux.m:                  derives flux simulations for a specified list of conditions
│   │   │   process_flux.m:                 processes simulated flux data into a format compatible with CARAMeL model construction
└───data: contains all data files used for OMG-ML model development
│   │   all_training_combinations.csv:      data sourced from SynergyxDB on drug combination synergies for HT29 and HCT11 cells 
│   │   hct116_fn_single_drug.xlsx:         single drug + f. nucleatum empty file to use for prediction
│   │   ht29_fuso_1drug.xlsx:               single drug + f. nucleatum empty file to use for prediction
│   │   logfc_geneexp.csv:                  Log fold change gene expression values for all treatments, microbes and drugs 
│   OMGML_main.mlx:             MATLAB livescript file that generates predictions used in OMGML paper
│   OMGML_main.m:               MATLAB file that generates predictions used in OMGML paper
