%% rawToDICOM
% This is the main script for the pipeline in converting raw Bruker files to DICOM files.
% This also includes processing CS data.

% Run whole script for the complete pipeline or each section as necessary

%% Preparation module - user set parameters
% Project name - e.g. AGORA
project = 'AGORA';
% Path to cohort inside project - e.g. AG_9\cohort1\week43
cohortPath = 'AG_9\cohort1\week43';

%% 1) sortRawData
% Copies data from the project's cohort path in R:\DataTransfer to Paravision into R:\Preprocessed data from Paravision
% Sorts only data into folders based on keywords = {'FLASH','TPM', 't1', 'MRE', 'LGE', 'tagged', 'CINE'}
sortRawData(project, cohortPath)

%% 2) Reconstruct CS


%% 3) creating DICOM files