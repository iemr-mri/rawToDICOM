function createDICOMCine(pathStruct)
    % Input:
        % pathStruct: struct containing various path strings for folder structure

    %% 1 - Locating all CINE and SG files
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
        destination = fullfile(pathStruct.DICOMRoot, pathStruct.project, pathStruct.cohort, 'CINE_DICOM', pathStruct.subjName, scansCINE(scan).name);
        scanName    = scansCINE(scan).name;

        %% 2.1 Check existence of DICOM file
        existD      = existDICOM(destination, scanName);
        if existD
            continue
        else
            disp('--------')
            disp(['Reconstructing ', scanName])
        end
        
        %% 2.2 Locate rawObj
        imagePath       = fullfile(scansCINE(scan).folder,scanName);
        try
            rawObj          = RawDataObject(imagePath, 'dataPrecision', 'double');
        catch
            warning('rawObj not found for %s', imagePath)
            continue
        end
        
        %% 2.3 Check existence of .mat file
        if ~exist(fullfile(imagePath, 'imageData.mat'), 'file')

            %% 2.3.1 - Rearrange kspace data to [x, y, slices, movieFrames, flowEncDir, coils]
            try
                kspaceSorted    = kspaceSort(rawObj);
            catch ME
                warning('Problem with kspace for %s', imagePath)
                fprintf('Error: %s\n', ME.message);
                continue
            end
        
             %% 2.3.2 - Performing CS reconstruction if CS file
            if contains(imagePath, 'CS_191021')
                disp('--------')
                disp(['Reconstructing CS data for ', scanName])
                final_kspace = reconstructCS(kspaceSorted);
            else
                final_kspace = kspaceSorted;
            end
        
            %% 2.3.3 - Combine coils
            imageData = combineCoils(final_kspace);
    
            %% 2.3.4 - Save image data
            save(fullfile(imagePath, 'imageData.mat'), "imageData")

        else % use previously processed imageData file if exists
            disp(['Using previously processed data for ', scanName])
        end

        %% 2.4 - Convert into DICOM and save in new root
        disp('--------')
        disp(['Converting ', scanName, ' to DICOM'])
        convertToDICOM(imagePath, rawObj, destination)
    end
    
    %% 3 - Perform reconstruction and DICOM conversion of each scan on scansSG with self-gating module
    if ~isempty(scansSG)
        existSG = 0;
        for scan = 1:length(scansSG)
            destination = fullfile(pathStruct.DICOMRoot, pathStruct.project, pathStruct.cohort, 'CINE_DICOM', pathStruct.subjName, scansSG(scan).name);
            scanName    = scansCINE(scan).name;

            %% 3.1 Check existence of DICOM file
            existD      = existDICOM(destination, scanName);
            if existD
                disp(['DICOM file already exist for ', scanName])
                existSG = 1;
            end
        end
        
        if existSG
            recon = reconCheck();
        else
            recon = 1;
        end
        
        %% 3.1 Perform self-gating module on the whole stack of slices to process them together
        if recon
            SGcine(scansSG);
    
            for scan = 1:length(scansSG)
                imagePath       = fullfile(scansSG(scan).folder,scansSG(scan).name);
                try
                    rawObj          = RawDataObject(imagePath, 'dataPrecision', 'double');
                catch
                    warning('rawObj not found for %s', imagePath)
                    continue
                end
    
                %% 3.3 Make DICOM files for each individual scan
                disp('--------')
                disp(['Converting ', scanName, ' to DICOM'])
                convertToDICOM(imagePath, rawObj, destination)
            end
        end
    end

    fclose('all');
end