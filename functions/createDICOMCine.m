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
        %% Check existence of dir and DICOM file
        destination = fullfile(pathStruct.DICOMRoot, pathStruct.project, pathStruct.cohort, pathStruct.subjName, 'CINE_DICOM', scansCINE(scan).name);
        [dirPath]                           = fileparts(destination);
        % "7" specifically checks if dirPath is a folder
        if exist(dirPath) ~= 7
            mkdir(dirPath)
        end
        
        if exist([destination,'.dcm'])
            disp(['DICOM file ', destination, '.dcm already exist.'])
            continue
        end

        %% 2.1 - Rearrange kspace data to [x, y, slices, movieFrames, flowEncDir, coils]
        imagePath       = fullfile(scansCINE(scan).folder,scansCINE(scan).name);
        try
            rawObj          = RawDataObject(imagePath, 'dataPrecision', 'double');
            visuParam       = readBrukerParamFile(fullfile(rawObj.Filespath.auto,'\pdata\1\visu_pars'));
        catch
            warning('rawObj or visuParam not found for %s', imagePath)
            return
        end
        kspaceSorted    = kspaceSort(rawObj);
    
         %% 2.2 - Performing CS reconstruction if CS file
        if contains(imagePath, 'CS_191021')
            disp('-------------------------------')
            disp(['Reconstructing CS data for ', scansCINE(scan).name])
            final_kspace = reconstructCS(kspaceSorted);
        else
            final_kspace = kspaceSorted;
        end
    
        %% 2.3 - Combine coils
        combined_im = combineCoils(final_kspace);

        %% 2.4 - Image corrections
        final_im    = imageCorrections(combined_im, rawObj, visuParam);
    
        %% 2.5 - Convert into DICOM and save in new root
        convertToDICOM(final_im, rawObj, destination)
    end

    fclose('all');
end