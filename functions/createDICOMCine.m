function createDICOMCine(pathStruct)
    % Input:
        % pathStruct: struct containing various path strings for folder structure

    %% 1 - Locating all CINE files
    addpath('R:\Felles_PCRTP\functions\BrukerFiles');
    
    % struct with all cine scans for the subject
    scansCINE       = dir(fullfile(pathStruct.sortedRoot, pathStruct.project,'CINE',pathStruct.cohort, pathStruct.subjName));
    scansCINE       = scansCINE(~ismember({scansCINE.name},{'..', '.'}));

    %% 2 - Perform reconstruction and DICOM conversion of each scan
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
            visuParam       = readBrukerParamFile(fullfile(rawObj.Filespath.auto,'\pdata\1\visu_pars'));
        catch ME
            warning('rawObj or visuParam not found for %s', imagePath)
            fprintf('Error: %s\n', ME.message);
            continue
        end
        
        %% 2.2 Check existence of .mat file
        if ~exist(fullfile(imagePath, 'imageData.mat'), 'file')

            %% 2.3 - Rearrange kspace data to [x, y, slices, movieFrames, flowEncDir, coils]
            kspaceSorted    = kspaceSort(rawObj);
        
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
            final_im    = imageCorrections(combined_im, rawObj);
    
            %% 2.7 - Save image data
            save(fullfile(imagePath, 'imageData.mat'), "final_im")
        else
            disp(['Using previously processed data for ', scansCINE(scan).name])
        end
        %% 2.8 - Convert into DICOM and save in new root
        convertToDICOM(imagePath, rawObj, destination)
    end

    fclose('all');
end