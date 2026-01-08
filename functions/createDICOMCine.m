function createDICOMCine(pathStruct)
    % Input:
        % pathStruct: struct containing various path strings for folder structure

    %% 1 - Locating all CINE files
    % struct with all cine scans for the subject
    scansCINE       = dir(fullfile(pathStruct.sortedRoot, pathStruct.project,'CINE',pathStruct.cohort, pathStruct.subjName));
    scansCINE       = scansCINE(~ismember({scansCINE.name},{'..', '.'}));

    % identify any self-gated scans and segregate into two lists
    scansSG         = scansCINE;
    scansCINE       = scansCINE(~contains({scansCINE.name},'SG'));
    scansSG         = scansSG(contains({scansSG.name},'SG'));
    if ~isempty(scansSG)
        scansSG = slicePositionSort(scansSG);
    end

    %% 2 - Perform reconstruction and DICOM conversion of each scan on scansCINE
    for scan = 1:length(scansCINE)
        %% 2.1 Check existence of dir and DICOM file
        destination = fullfile(pathStruct.DICOMRoot, pathStruct.project, pathStruct.cohort, 'CINE_DICOM', pathStruct.subjName, scansCINE(scan).name);
        [dirPath]                           = fileparts(destination);
        % "7" specifically checks if dirPath is a folder
        if exist(dirPath, 'dir') ~= 7
            mkdir(dirPath)
        end
        
        if exist([destination,'.dcm'], 'file')
            disp(['DICOM file ', destination, '.dcm already exist.'])
            continue
        end

        imagePath       = fullfile(scansCINE(scan).folder,scansCINE(scan).name);
        try
            rawObj          = RawDataObject(imagePath, 'dataPrecision', 'double');
        catch
            warning('rawObj not found for %s', imagePath)
            continue
        end
        
        %% 2.2 Check existence of .mat file
        if ~exist(fullfile(imagePath, 'imageData.mat'), 'file')

            %% 2.3 - Rearrange kspace data to [x, y, slices, movieFrames, flowEncDir, coils]
            try
                kspaceSorted    = kspaceSort(rawObj);
            catch ME
                warning('Problem with kspace for %s', imagePath)
                fprintf('Error: %s\n', ME.message);
                continue
            end
        
             %% 2.4 - Performing CS reconstruction if CS file
            if contains(imagePath, 'CS_191021')
                disp('-------------------------------')
                disp(['Reconstructing CS data for ', scansCINE(scan).name])
                final_kspace = reconstructCS(kspaceSorted);
            else
                final_kspace = kspaceSorted;
            end
        
            %% 2.5 - Combine coils
            combined_im = combineCoils(final_kspace);

            %% 2.6 - Image corrections
            imageData    = imageCorrections(combined_im, rawObj);
    
            %% 2.7 - Save image data
            save(fullfile(imagePath, 'imageData.mat'), "imageData")
        else
            disp(['Using previously processed data for ', scansCINE(scan).name])
        end
        %% 2.8 - Convert into DICOM and save in new root
        convertToDICOM(imagePath, rawObj, destination)
    end
    
    %% 3 - Perform reconstruction and DICOM conversion of each scan on scansSG with self-gating module
    if ~isempty(scansSG)
        %% 3.1 Perform self-gating module on the whole stack of slices to process them together
        SGcine(scansSG);

        %% 3.2 Set DICOM destination and check if DICOM already exist
        for scan = 1:length(scansSG)
            destination = fullfile(pathStruct.DICOMRoot, pathStruct.project, pathStruct.cohort, 'CINE_DICOM', pathStruct.subjName, scansSG(scan).name);
            [dirPath]   = fileparts(destination);
            % "7" specifically checks if dirPath is a folder
            if exist(dirPath, 'dir') ~= 7
                mkdir(dirPath)
            end
            %check if DICOM file already exist
            if exist([destination,'.dcm'], 'file')
                disp(['DICOM file ', destination, '.dcm already exist.'])
                continue
            end
    
            imagePath       = fullfile(scansSG(scan).folder,scansSG(scan).name);
            try
                rawObj          = RawDataObject(imagePath, 'dataPrecision', 'double');
            catch
                warning('rawObj not found for %s', imagePath)
                continue
            end

            %% 3.3 Make DICOM files for each individual scan
            convertToDICOM(imagePath, rawObj, destination)
        end
    end

    fclose('all');
end