%% rawToDICOM
% This is the main script for the pipeline in converting raw Bruker files to DICOM files.
% This also includes processing CS data.
 
% Run whole script for the complete pipeline or each section as necessary

%% User set parameters - project and cohort names

% Set up pathStruct for easy navigating
% Project name - e.g. AGORA
pathStruct.project         = 'AGORA';
% Path to cohort inside project - e.g. AG_9\cohort1\week43
pathStruct.cohort          = 'MI_9\cohort 1\Week 4';

if isempty(pathStruct.project) || isempty(pathStruct.cohort)
    error('Please make sure to fill out project field and cohort field correctly before proceeding.')
end

%% Preparation module - path settings
% Root paths
% This is where the raw data is collected
pathStruct.rawRoot            = 'R:\DataTransfer from ParaVision';
% This is where the raw data is sorted into
pathStruct.sortedRoot         = 'R:\Preprocessed data from Paravision';
% This is where the DICOM files are saved
pathStruct.DICOMRoot          = 'R:\Projects';

% adding Bruker functions for reading raw files
addpath('R:\Felles_PCRTP\functions\BrukerFiles');

% adding functions folder and common_utils which should be in parallell folder to current directory
addpath('functions\')
addpath('helpers\');

%% 1 - sortRawData
% Copies data from the project's cohort path in R:\DataTransfer to Paravision into R:\Preprocessed data from Paravision
% Sorts only data into folders based on keywords = {'FLASH','TPM', 't1', 'MRE', 'LGE', 'tagged', 'CINE'}

sortRawData(pathStruct);

%% 2 - Create DICOM files of CINE images
% Finds all scans in the CINE folder
subjectStruct              = dir(fullfile(pathStruct.sortedRoot, pathStruct.project, 'CINE', pathStruct.cohort));
subjectStruct              = subjectStruct(~ismember({subjectStruct.name},{'..', '.'}));

if isempty(subjectStruct)
    warning('No CINE scans found for %s. Make sure project and cohort name is correct.', pathStruct.cohort)
    return
end

%% 2.1 - Perfrom operation for each scan
% Sort kspace into [x, y, slice, frame, MEG, coil]
% Reconstructs CS data if undersampled
% Converts into DICOM and saves in corresponding project folder under R:\Projects

for scan = 1:length(subjectStruct)
    pathStruct.subjName     = subjectStruct(scan).name;
    disp('-------------------------------')
    disp(['Creating DICOM files for ', pathStruct.subjName])
    createDICOMCine(pathStruct)
    disp('Completed.')
end

disp('-------------------------------')
disp(['DICOM files stored in ', pathStruct.DICOMRoot,'\', pathStruct.project,'\', pathStruct.cohort])