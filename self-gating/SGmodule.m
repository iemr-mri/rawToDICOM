function SGmodule(scansSG, pm)
    
    %% Check existence of DICOM file
    n_SG = 0;
    for scan = 1:length(scansSG)
        destination = fullfile(pm.DICOMRoot, pm.project, pm.cohort, 'CINE_DICOM', pm.subjName, scansSG(scan).name);
        scanName    = scansSG(scan).name;
    
        existD      = existDICOM(destination, scanName);
        if existD
            disp(['DICOM file already exist for ', scanName])
            n_SG = n_SG+1;
        end
    end
    
    %% Perform self-gating module on the whole stack of slices to process them together
    % self-gating module will run if number of DICOM files are less than number of SG scans or forceSG is set to true
    if (n_SG<length(scansSG)) || pm.forceSG
        if pm.forceSG
            disp('--------')
            disp('Note: forceSG set to true.')
        end
        
        % Call self-gating module
        SGcine(scansSG);
    
    
        for scan = 1:length(scansSG)
            destination = fullfile(pm.DICOMRoot, pm.project, pm.cohort, 'CINE_DICOM', pm.subjName, scansSG(scan).name);
            scanName    = scansSG(scan).name;
            imagePath       = fullfile(scansSG(scan).folder,scansSG(scan).name);
            
            try
                rawObj          = RawDataObject(imagePath, 'dataPrecision', 'double');
            catch
                warning('rawObj not found for %s', imagePath)
                continue
            end
        
            %% Make DICOM files for each individual scan
            disp('--------')
            disp(['Converting ', scanName, ' to DICOM'])
            convertToDICOM(imagePath, rawObj, destination)
        end

    end
end