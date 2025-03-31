%% rawToDICOM
% This is the main script for the pipeline in converting raw Bruker files to DICOM files.
% This also includes processing CS data.

% Run whole script for the complete pipeline or each section as necessary

%% Preparation module - user set parameters

% Set up pathStruct for easy navigating
% Project name - e.g. AGORA
pathStruct.project         = 'AGORA';
% Path to cohort inside project - e.g. AG_9\cohort1\week43
pathStruct.cohort          = 'AG_9\cohort1\week43';
% Root paths
pathStruct.oldRoot         = 'R:\Henrik Elias'; % 'R:\Preprocessed from Paravision'
pathStruct.newRoot         = 'R:\Henrik Elias\DICOM_test'; % 'R:\Projects'

% adding Bruker functions for reading raw files
addpath('R:\Felles_PCRTP\functions\BrukerFiles');

% adding common_utils which should be in parallell folder to current directory
addpath(fullfile(fileparts(pwd), 'common_utils'));

%% 1 - sortRawData
% Copies data from the project's cohort path in R:\DataTransfer to Paravision into R:\Preprocessed data from Paravision
% Sorts only data into folders based on keywords = {'FLASH','TPM', 't1', 'MRE', 'LGE', 'tagged', 'CINE'}

sortRawData(pathStruct.project, pathStruct.cohort);

%% 2 - Create DICOM files of CINE images
% Finds all scans in the CINE folder
subjectStruct           = dir(fullfile(pathStruct.oldRoot, pathStruct.project, 'CINE', pathStruct.cohort));
subjectStruct           = subjectStruct(~ismember({subjectStruct.name},{'..', '.'}));

%% 2.1 - For each scan
% Sort kspace into [x, y, slice, frame, MEG, coil]
% Reconstructs CS data if undersampled
% Converts into DICOM and saves in corresponding project folder under R:\Projects


%for scan = 1:length(sortedStruct)
    pathStruct.subjName     = subjectStruct(1).name;
    createDICOMCine(pathStruct)
%end