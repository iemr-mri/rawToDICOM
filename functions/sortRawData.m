function sortRawData(project, cohortPath)
    %% sortRawData.mat
    % project [string] = name of project folder inside R:\DataTransfer from ParaVision
    % cohortPath [string] = path of folders to cohort
    % Should be on the form [group]\[cohort] e.g. cohortPath = '\AG_9\cohort1\week43'

    % 1) Preparation
    % adds common_utils path and prepares before the loop
    % 2) Loop through each subject in the cohort

    %% 1) Preparation
    
    % Data location
    originalPath = fullfile('R:\DataTransfer from ParaVision', project, cohortPath);
    disp(['Locating data in ', originalPath])
    disp('-------------------------------')

    newRoot = 'R:\Henrik Elias';
    %newRoot = 'R:\Preproccessed data from Paravision';
    
    % Struct of the cohort folder (originalPath) that contains each subject
    cohortStruct            = dir(originalPath);
    cohortStruct            = cohortStruct(~ismember({cohortStruct.name},{'.', '..'}));

    %What are the protocol name keywords to search for?
    keywords                 = {'FLASH','TPM', 't1', 'MRE', 'LGE', 'tagged', 'CINE'};
    
    %% 2) Loop through each subject
    for subject=1:length(cohortStruct)
        disp(['Sorting data for: ', cohortStruct(subject).name, ' into R:\Preprocessed data from Paravision.'])
        subjectStruct           = dir(fullfile(originalPath, cohortStruct(subject).name));
        subjectStruct           = subjectStruct(~ismember({subjectStruct.name},{'.', '..'}));
        
        % Loop through each scan in the subjectStruct
        for scan = 1:length(subjectStruct)
            % Construct the full file path for 'acqp' in the current folder
            scanFolder = fullfile(subjectStruct(scan).folder, subjectStruct(scan).name);
            acqpPath = fullfile(scanFolder, 'acqp');

            % Check if the 'acqp' file exists in the current folder
            if isfile(acqpPath)
                % Read the scan name
                acqpStruct = readBrukerParamFile(acqpPath);
                scan_name = acqpStruct.ACQ_scan_name;

                % If the scan name matches one of the keywords, we perform copy the data into the new folder, sorted by keyword 
                for key = 1:length(keywords)
                    if contains(lower(scan_name), lower(keywords(key)))
                        if contains(keywords(key),'FLASH')
                            folderKey = 'CINE';
                        else
                            folderKey = upper(keywords{key});
                        end

                        newPath = fullfile(newRoot, project, folderKey, cohortPath, cohortStruct(subject).name, scan_name);
                        if ~isfile([newPath,'\rawdata.job0']) || ~isfile([newPath,'\acqp'])
                            copyfile(scanFolder, newPath);
                            disp([scan_name, ' copied.'])
                        else
                            disp([scan_name, ' already exist.'])
                        end
                    end
                end
            end
        end
        disp('-------------------------------')
    end
close all
end