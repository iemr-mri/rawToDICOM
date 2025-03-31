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
        %% 2.1 - Rearrange kspace data to [x, y, slices, movieFrames, flowEncDir, coils]
        imagePath       = fullfile(scansCINE(scan).folder,scansCINE(scan).name);
        rawObj          = RawDataObject(imagePath, 'dataPrecision', 'double');
        
        kspaceSorted    = kspaceSort(rawObj);
    
         %% 2.2 - Performing CS reconstruction if CS file
        % CS files contains the 'CS_191021' string
        if ismember(imagePath, 'CS_191021')
            final_kspace = reconstructCS(kspaceSorted);
        else
            final_kspace = kspaceSorted;
        end
    
        %% 2.3 - Combine coils
        final_im = combineCoils(final_kspace);
    
        %% 2.4 - Convert all scans into DICOM and saving in new root
        destination = fullfile(pathStruct.DICOMRoot, pathStruct.project, pathStruct.cohort, pathStruct.subjName, 'CINE_DICOM', scansCINE(scan).name);
        convertToDICOM(final_im, rawObj, destination)
    end

    close all
end