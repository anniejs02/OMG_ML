# Onco Microbiome GEM Machine Learning (OMG-ML)
Welcome to OMG-ML, a systems biology approach to predicting colon cancer treatments. This guide walks you through setting up and running OMG-ML in MATLAB. 
This tool allows for predicting the effect of multi-drug combinations on colorectal cancer (CRC) and the effect of microorganisms on treatment responses.



<p align="center">
  <img width="601" height="598" alt="image" src="https://github.com/user-attachments/assets/169f4028-60c8-49ad-a38c-097626af8598" />
</p>


This model can be used along with the included drug combination files (sourced from [SynergyxDB]([url](https://www.synergxdb.ca/))). Additional drug or microbial treatments can also be added to make novel predictions by adding log fold changes of gene expression in treatment compared to control to the included gene expression file. 



# License 
Released via GPL GNU License
© 2025 The Regents of the University of Michigan
Chandrasekaran Research Group - https://systemsbiologylab.org/
Contact: csriram@umich.edu

# Dependencies 
Using OMG-ML requires the installation of the following packages/programs. 
1. [MATLAB]([url](https://www.mathworks.com/products/matlab.html)) 2019b or higher
2. [Gurobi]([url](https://www.gurobi.com/)) optimization solver
3. [COBRA Toolbox]([url](https://github.com/opencobra/cobratoolbox))

# Installation 
Install the folders and packages included by cloning this repository to your local computer

git clone https://github.com/anniejs02/OMG_ML.git

# Instructions 
See OMGML_main.m or OMGML_main.mlx for step-by-step instructions in applying OMG-ML. The code outlines the following steps: 
1. Load in the RECON1 human metabolic model (Recon1_Duarte.mat)
2. Constrain RECON1 using transcriptomics from drug treatment and co-cultures
3. Generate features from constrained models
4. Train OMG-ML model using known drug combination synergies
5. Employ model

# Files


<pre><code>
OMGML
├── functions/                        # Contains all relevant code files
│   ├── CARAMeL_suite/           # Scripts implementing the CARAMeL approach
│   │   ├── caramel.m                        # Main script that constructs a CARAMeL model
│   │   ├── caramel_anova.m                  # One-way ANOVA to identify differential flux reactions
│   │   ├── caramel_assessPerf.m             # Assesses model predictive performance
│   │   ├── caramel_auroc.m                  # Calculates AUROC
│   │   ├── caramel_classify.m               # Classifies drug interactions
│   │   ├── caramel_crossValidation.m        # k-fold cross-validation
│   │   ├── caramel_features2gem.m           # Extracts GEM reactions linked to top features
│   │   ├── caramel_featurize.m              # Converts raw input data into ML format
│   │   ├── caramel_leaveOut.m               # Leave-one-out analysis
│   │   ├── caramel_plot.m                   # Plots model performance
│   │   ├── caramel_processInteractions.m    # Standardizes drug interaction data
│   │   ├── caramel_rankSubsystems.m         # Finds enriched pathways among top features
│   │   ├── caramel_screen.m                 # Screens all possible drug combinations
│   │   └── caramel_topFeatures.m            # Extracts top model features
│   └── GEM_functions/          # Code related to GEM-based metabolic simulations
│       ├── change_media.m                   # Changes simulated media for GEM
│       ├── constrain_flux_regulation.m      # Simulates fluxes with omics constraints
│       ├── derive_flux.m                    # Derives flux simulations for given conditions
│       └── process_flux.m                   # Processes flux data for CARAMeL model use
├── predictions/                  # All prediction output files and folders
│   ├── prediction_1microbe_1drug/           # Original folder moved here
│   ├── prediction_2drugs_1microbe/          # Original folder moved here
│   ├── hct116_2d_pred.xlsx                    # Renamed prediction file
│   └── MOREMICROBES_averaged_10predictions_HCT116.xlsx  # Newly added predictions file
├── OMGML_main.mlx               # MATLAB Live Script for OMGML paper predictions
├── OMGML_main.m                 # MATLAB script version for OMGML paper predictions
└── figure_clean.rmd             # R code for all figures included in manuscript

</code></pre>




