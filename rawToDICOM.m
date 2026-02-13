%% rawToDICOM
% This is the main script for the pipeline in converting raw Bruker files to DICOM files.
% This also includes processing
    % Compressed sensing
    % Self gated images
    % Partial echo
 
% Run whole script for the complete pipeline or each section as necessary

%% User set parameters - project/cohort names and settings

% Set up parameter struct (pm) for easy navigating
% Project name - e.g. AGORA
pm.project         = 'AGORA';
% Path to cohort inside project - e.g. AG_9\cohort1\week43
pm.cohort          = 'AB_24\cohort1\week 6';

% Some flags for tailoring
pm.skipSort        = false; % skips sortRawData
pm.forceRecon      = false; % forces to do reconstruction even if imageData.mat exist
pm.forceDICOM      = false; % overwrites existing DICOM files

%% Preparation module - path settings
% Root paths
% This is where the raw data is collected
pm.rawRoot            = 'R:\DataTransfer from ParaVision';
% This is where the raw data is sorted into
pm.sortedRoot         = 'R:\Preprocessed data from Paravision';
% This is where the DICOM files are saved
pm.DICOMRoot          = 'R:\Projects';

% adding Bruker functions for reading raw files
addpath('R:\Felles_PCRTP\functions\BrukerFiles_2019\pvtools');

% adding functions/helpers folders which is in the same directory as this file
addpath('functions\');
addpath('helpers\');
addpath(genpath('self-gating'));

%% 1 - sortRawData
% Copies data from the project's cohort path in R:\DataTransfer to Paravision into R:\Preprocessed data from Paravision
% Sorts only data into folders based on keywords = {'FLASH','TPM', 't1', 'MRE', 'LGE', 'tagged', 'CINE'}

sortRawData(pm);

%% 2 - Locate CINE folder
% Finds all scans in the CINE folder
subjectStruct              = dir(fullfile(pm.sortedRoot, pm.project, 'CINE', pm.cohort));
subjectStruct              = subjectStruct(~ismember({subjectStruct.name},{'..', '.'}));

if isempty(subjectStruct)
    warning('No CINE scans found for %s. Make sure project and cohort name is correct.', pm.cohort)
    return
end

%% 3 - Perfom reconstruction and DICOM conversion for each subject
% Sort kspace into [x, y, slice, frame, MEG, coil]
% Reconstructs CS data if undersampled
% Converts into DICOM and saves in corresponding project folder under R:\Projects

for subj = 1:length(subjectStruct)
    pm.subjName     = subjectStruct(subj).name;
    disp('--------')
    disp(['Creating DICOM files for ', pm.subjName])
    createDICOMCine(pm)
    disp('Completed.')
end

disp('--------')
disp(['DICOM files stored in ', pm.DICOMRoot,'\', pm.project,'\', pm.cohort])