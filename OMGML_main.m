%% OMG-ML
% OMG-ML is a Flux Balance Analysis + Machine learning framework to predict 
% effective drug combination and microbiome based therapies for colorectal cancer. 
% This begins by importing the training data (known drug combination synergies) 
% and generating flux data via flux balance analysis (derive_flux). Term 'bug' 
% throughout this code is used as shorthand to represent microorganism. 
% 
% 
% 
% 
% 
% Figure 1. Graphical abstract of the CRC-CARAMeL model
% 
% Original CARAMeL paper citation: Carolina H Chung, Sriram Chandrasekaran, 
% A flux-based machine learning model to simulate the impact of pathogen metabolic 
% heterogeneity on drug interactions, PNAS Nexus, Volume 1, Issue 3, July 2022, 
% pgac132, <https://doi.org/10.1093/pnasnexus/pgac132 https://doi.org/10.1093/pnasnexus/pgac132>
% 
% 
% 
% Initialize enviornment

%Initialize CobraToolbox & Gurobi 

% Clear environment (if prompted)
sound(sin(1:1000))
clear_env = input('Clear environment? (true/false) ');
if clear_env
    clearvars           % clear workspace
    close all hidden    % close all figures
    clc                 % clear command window
end

% Initialize COBRA Toolbox
initCobraToolbox
changeCobraSolver('gurobi')

% Define current working directory
folder = erase(which('CRC_CARAMeL_main.mlx'), 'CRC_CARAMeL_main.mlx'); 
%% 
% 
% 
% Load in known drug synergy data (sourced from <https://www.synergxdb.ca/ SynergyxDB>)
% 
% Description: two drug names per combinations with a 3rd column for the synergy 
% score (Loewe)

train_all = readtable("all_training_combinations.csv");

%% 
% 
%% Generate Flux
% Necessary files: 
%% 
% # Human metabolic model (RECON1, RECON3D, iHuman, etc)
% # Logfold change gene expression data from drug and/or microbial treatments 
%% 
% Genes will be mapped to human metabolic model to derive flux using the 'derive_flux' 
% function. Thresholds for differentially expressed genes may need to be adapted 
% for your specific data. 
% 
% 
% 
% Load in gene expression file

ht29_hct116.added_hct_fuso = readtable('logfc_geneexp.csv');
head(ht29_hct116.added_hct_fuso)
%% 
% Normalize


% Extract gene names 
genes = ht29_hct116.added_hct_fuso{:, 1};

% Extract drug treatments (excluding the first column with gene names)
conditions = ht29_hct116.added_hct_fuso.Properties.VariableNames(2:end);

% Grab all columns except first (gene names)
idx = 2:width(ht29_hct116.added_hct_fuso); 
fold_data = table2array(ht29_hct116.added_hct_fuso(:, idx));

% Normalize 
med = median(fold_data, 2, 'omitnan');
mdev = mad(fold_data, 1, 2);            % MATLAB's MAD, normalized by 1/n
scale = 1.4826 * mdev;                  % ~consistent with SD under normal dist
scale(scale == 0 | isnan(scale)) = 1;   % prevent divide-by-zero

fold_data_robust = (fold_data - med) ./ scale;
fold_data_robust = round(fold_data_robust, 2);

% --- Explicitly define bug columns ---
bugNames = {'sgg', 'entero', 'fuso_hct116', 'ht29_fuso', 'ecoli', 'averonii'};
bugCols = ismember(conditions, bugNames);

% Everything else is a drug column ---
drugCols = ~bugCols;

% First chunk (drugs)
gene_struct_drugs = struct( ...
    'genes', {genes}, ...
    'conditions', {conditions(drugCols)}, ...
    'fold_data', {fold_data_robust(:, drugCols)} );

% Second chunk (bugs)
gene_struct_bugs = struct( ...
    'genes', {genes}, ...
    'conditions', {conditions(bugCols)}, ...
    'fold_data', {fold_data_robust(:, bugCols)} );

%% 
% Derive Flux Values for Bugs and Drugs seperately 

% Bugs 
specs = [1e-2 1e-2 1e-2]; 
thresh = [-0.2, 0.2];
flux_bugs = derive_flux(ht29_GEMmodel, gene_struct_bugs, 'ModelSpecs', specs, 'Threshold', thresh);

% Drugs 
specs = [1e-2 1e-2 1e-2]; 
thresh = [-1, 1];
flux_drugs = derive_flux(ht29_GEMmodel, gene_struct_drugs, 'ModelSpecs', specs, 'Threshold', thresh);

metadata = flux_bugs.flux(:, 1:3);
%% 
% Combine flux into one table 

% Extract numeric flux columns
flux_bugs_data  = flux_bugs.flux(:, 4:end);
flux_drugs_data = flux_drugs.flux(:, 4:end);


% Combine into one table
flux_combined = [metadata, flux_bugs_data, flux_drugs_data];

% Process for input to ML model 
combined_ps = process_flux(flux_combined, ht29_GEMmodel);
%% Assess Perfomance 
% 10 fold cross validation 

[cv, best, idx] = ...
    caramel_crossValidation(combined_ps, train_all, 'Kfold', 10, 'Key', key, 'MLtype', 'fitrensemble', 'Method', 'Bag');
%% Train Model 

[training_data, model] = caramel(combined_ps, train_all_filtered, ...
    'train', 'MLtype', 'fitrensemble', 'Method', 'Bag', 'Key', key);
%% Make Predictions - No Microbes/Bugs
% Description: Predicting novel combinations (not included in training) of 2 
% drug combinations. 

%% ------------------------------------------------------------
% 2-Drug Combination Predictions — HT29
%% ------------------------------------------------------------

% Step 1: Get all drug names (excluding 'gene')
allVars = ht29_hct116.added_hct_fuso.Properties.VariableNames;
drugNames = setdiff(allVars, {'gene'});

% Normalize naming (lowercase, underscores)
normalizeNames = @(x) lower(strrep(strtrim(x), '-', '_'));
drugNames = normalizeNames(drugNames);

% Normalize names in training data for consistency
ht29.Drug1 = normalizeNames(string(ht29.Drug1));
ht29.Drug2 = normalizeNames(string(ht29.Drug2));

% Step 2: Generate all unique 2-drug combinations
allCombos = nchoosek(drugNames, 2);  % Nx2 cell array

% Step 3: Build a new table of test combinations
Drug1 = string(allCombos(:, 1));
Drug2 = string(allCombos(:, 2));
LoeweScores = NaN(size(Drug1));  % Placeholder column

testComboTable = table(Drug1, Drug2, LoeweScores, ...
    'VariableNames', {'Drug1', 'Drug2', 'Loewe'});

% Step 4: Format test data for model input
interaction_ht29_struct = struct( ...
    'names', {table2array(testComboTable(:, 1:2))}, ...
    'scores', {table2array(testComboTable(:, 3))} ...
);

% Step 5: Run predictions over multiple iterations
numIterations = 10;
numTestCombinations = size(interaction_ht29_struct.names, 1);

allInteractionNames = [];
cumulativePredScores = [];

fprintf('Running CARAMeL predictions for HT29 (%d iterations)...\n', numIterations);

for iteration = 1:numIterations
    fprintf('Iteration %d of %d...\n', iteration, numIterations);

    % Prepare test struct for current iteration
    testStruct = struct( ...
        'names', {interaction_ht29_struct.names(:, 1:2)}, ...
        'scores', {interaction_ht29_struct.scores(:, 1)} ...
    );

    % Run prediction
    iterationPredictions = caramel(combined_ps, testStruct, 'predict', ...
        'MLModel', all_model2, ...
        'MLtype', 'fitrensemble', ...
        'Method', 'Bag', ...
        'Key', key ...
    );

    % Aggregate predictions
    for j = 1:numel(iterationPredictions)
        if iteration == 1
            allInteractionNames = [allInteractionNames; iterationPredictions(j).interactionNames];
            cumulativePredScores = [cumulativePredScores; iterationPredictions(j).predScores];
        else
            cumulativePredScores = cumulativePredScores + iterationPredictions(j).predScores;
        end
    end
end

% Step 6: Average across all iterations
averagedPredScores = cumulativePredScores / numIterations;

% Step 7: Combine into final table
interactionTable = array2table(allInteractionNames, 'VariableNames', {'Drug1', 'Drug2'});
scoreTable = array2table(averagedPredScores, 'VariableNames', {'PredScores'});
combinedResults = [interactionTable, scoreTable];

% Step 8: Save results
outputFile = 'HT29_2drug_predictions_avg10.xlsx';
writetable(combinedResults, outputFile, 'Sheet', 'Sheet1', 'WriteMode', 'overwrite');
fprintf('Averaged predictions saved to %s\n\n', outputFile);



%% ------------------------------------------------------------
% 2-Drug Combination Predictions — HCT116
%% ------------------------------------------------------------

% Step 1: Get all drug names (excluding 'gene')
allVars = ht29_hct116.added_hct_fuso.Properties.VariableNames;
drugNames = setdiff(allVars, {'gene'});

normalizeNames = @(x) lower(strrep(strtrim(x), '-', '_'));
drugNames = normalizeNames(drugNames);

% Normalize names in HCT116 training data
hct116.Drug1 = normalizeNames(string(hct116.Drug1));
hct116.Drug2 = normalizeNames(string(hct116.Drug2));

% Step 2: Generate all possible 2-drug pairs
allCombos = nchoosek(drugNames, 2);

% Step 3: Build test table
Drug1 = string(allCombos(:, 1));
Drug2 = string(allCombos(:, 2));
CellLine = repmat("hct116", size(Drug1));
LoeweScores = NaN(size(Drug1));

testComboTable = table(Drug1, Drug2, CellLine, LoeweScores, ...
    'VariableNames', {'Drug1', 'Drug2', 'CellLine', 'Loewe'});

% Step 4: Format as struct for prediction
interactionStruct = struct( ...
    'names', {table2array(testComboTable(:, 1:3))}, ...
    'scores', {table2array(testComboTable(:, 4))} ...
);

% Step 5: Run predictions
numIterations = 10;
numTestCombos = size(interactionStruct.names, 1);

allInteractionNames = [];
cumulativePredScores = [];

fprintf('Running CARAMeL predictions for HCT116 (%d iterations)...\n', numIterations);

for iter = 1:numIterations
    fprintf('Iteration %d of %d...\n', iter, numIterations);

    testStruct = struct( ...
        'names', {interactionStruct.names}, ...
        'scores', {interactionStruct.scores} ...
    );

    predictionBatch = caramel(combined_ps, testStruct, 'predict', ...
        'MLModel', all_model2, ...
        'MLtype', 'fitrensemble', ...
        'Method', 'Bag', ...
        'Key', key ...
    );

    for j = 1:numel(predictionBatch)
        if iter == 1
            allInteractionNames = [allInteractionNames; predictionBatch(j).interactionNames];
            cumulativePredScores = [cumulativePredScores; predictionBatch(j).predScores];
        else
            cumulativePredScores = cumulativePredScores + predictionBatch(j).predScores;
        end
    end
end

% Step 6: Average predictions and save
averagedPredScores = cumulativePredScores / numIterations;

interactionNameTable = array2table(allInteractionNames, ...
    'VariableNames', {'Drug1', 'Drug2', 'CellLine'});
averagedScoreTable = array2table(averagedPredScores, ...
    'VariableNames', {'PredScores'});

predictionResults = [interactionNameTable, averagedScoreTable];

outputFile = 'HCT116_2drug_predictions_avg10.xlsx';
writetable(predictionResults, outputFile, ...
    'Sheet', 'Sheet1', 'WriteMode', 'overwrite');
fprintf('Averaged predictions saved to "%s"\n', outputFile);
%% Make Predictions for Drug + Bug
% The below code can be repeated for the other bugs 

% Step 1: Load and prepare HT29 single-drug data
inputFile = "home/anniejs/ht29_fuso_1drug.xlsx";
singleDrugData = readtable(inputFile);

% Add empty 'hct' column for format consistency
singleDrugData.hct = repmat({''}, height(singleDrugData), 1);

% Rename column if necessary
if any(strcmp(singleDrugData.Properties.VariableNames, 'log2foldchange'))
    singleDrugData.Properties.VariableNames{'log2foldchange'} = 'ht29_fuso';
end

% Convert to struct for CARAMeL
singleDrugStruct = struct( ...
    'names', {table2array(singleDrugData(:, 1:3))}, ...
    'scores', {table2array(singleDrugData(:, 4))} ...
);

% Replace any remaining 'log2foldchange' entries in names
for i = 1:numel(singleDrugStruct.names)
    if ischar(singleDrugStruct.names{i}) || isstring(singleDrugStruct.names{i})
        if strcmp(singleDrugStruct.names{i}, 'log2foldchange')
            singleDrugStruct.names{i} = 'ht29_fuso';
        end
    end
end

% Step 2: Predict across multiple iterations
numIterations = 10;
aggregatedResults = [];

for iter = 1:numIterations
    fprintf('Running prediction iteration %d/%d...\n', iter, numIterations);

    pred = caramel(combined_ps, singleDrugStruct, ...
        'predict', ...
        'MLmodel', all_model2, ...
        'MLType', 'fitrensemble', ...
        'Key', key ...
    );

    % Format predictions
    predNamesTable = cell2table(pred.interactionNames, ...
        'VariableNames', {'Drug', 'Cell', 'Fuso'});
    predScoresTable = array2table(pred.predScores, 'VariableNames', {'PredScores'});
    iterResults = [predNamesTable, predScoresTable];

    % Aggregate
    if iter == 1
        aggregatedResults = iterResults;
    else
        aggregatedResults.PredScores = aggregatedResults.PredScores + iterResults.PredScores;
    end
end

% Step 3: Average and export
aggregatedResults.PredScores = aggregatedResults.PredScores / numIterations;
outputFile = 'HT29_fuso_1d.xlsx';
writetable(aggregatedResults, outputFile, 'Sheet', 'Sheet1');
fprintf('Averaged HT29+Fuso single-drug predictions saved to "%s"\n', outputFile);


%% Drug + F. nucleatum Prediction - HCT116
% Step 1: Load and prepare HCT116 single-drug data
inputFile = "home/anniejs/hct116_fn_single_drug.xlsx";
singleDrugData = readtable(inputFile);

% Add 'fuso' column
singleDrugData.fuso = repmat({'fuso_hct116'}, height(singleDrugData), 1);

% Remove placeholder/duplicate column
singleDrugData(:,3) = [];

% Reorder columns
desiredColumns = {'gene', 'hct', 'fuso', 'score'};
singleDrugData = singleDrugData(:, desiredColumns);

% Convert to struct for CARAMeL
singleDrugStruct = struct( ...
    'names', {table2array(singleDrugData(:,1:3))}, ...
    'scores', {table2array(singleDrugData(:,4))} ...
);

% Step 2: Predict across multiple iterations
numIterations = 10;
aggregatedResults = [];

for iter = 1:numIterations
    fprintf('Running prediction iteration %d/%d...\n', iter, numIterations);

    pred = caramel(combined_ps, singleDrugStruct, ...
        'predict', ...
        'MLmodel', all_model2, ...
        'MLType', 'fitrensemble', ...
        'Method', 'Bag', ...
        'Key', key ...
    );

    predNamesTable = cell2table(pred.interactionNames, ...
        'VariableNames', {'Drug', 'Cell', 'Fuso'});
    predScoresTable = array2table(pred.predScores, 'VariableNames', {'PredScores'});
    iterResults = [predNamesTable, predScoresTable];

    if isempty(aggregatedResults)
        aggregatedResults = iterResults;
    else
        aggregatedResults.PredScores = aggregatedResults.PredScores + iterResults.PredScores;
    end
end

% Step 3: Average and export
aggregatedResults.PredScores = aggregatedResults.PredScores / numIterations;
outputFile = 'HT29_fuso_1d.xlsx';
writetable(aggregatedResults, outputFile, 'Sheet', 'Sheet1');
fprintf('Averaged HCT116+Fuso single-drug predictions saved to "%s"\n', outputFile);
%% Extract Top Features 


numTopFeatures = 50;
topFeaturesTable = caramel_topFeatures(model, 'N', numTopFeatures);

%% Step 2: Featurize Training Set Using Flux Profiles
[phenotypeData, jointProfiles, interactionData] = caramel_featurize( ...
    combined_ps, ...
    train_all_filtered, ...
    'Key', key ...
);

% Package interaction scores and joint profiles into a struct
data_struct = struct( ...
    'interactionScores', interactionData.scores, ...
    'jointProfiles', jointProfiles(2:end, :) ...
);

%% Step 3: Display Feature Information
disp('Sample of joint feature profiles (first 5 rows, 5 columns):');
disp(data_struct.jointProfiles(1:5, 1:5));

disp('Feature variable names:');
disp(data_struct.jointProfiles.Properties.VariableNames);

disp('Top 10 features:');
disp(jointProfiles.Feature(1:10));

fprintf('Number of unique features: %d\n', ...
    numel(unique(data_struct.jointProfiles.Feature)));

%% Step 4: Map Top Features to GEM Model
gemFeatureTable = caramel_features2gem( ...
    topFeaturesTable, ...
    ht29_GEMmodel, ...
    data_struct ...
);

%% Step 5: Export to CSV
outputFile = '=topfeatures_fusomodel.csv';
writetable(gemFeatureTable, outputFile);
%% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%
