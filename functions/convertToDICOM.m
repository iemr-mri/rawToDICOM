function convertToDICOM(imagePath,rawObj,destination)
    % Converts image data (kspace) into DICOM with metainfo from rawObj and saves at a specified destination.
    % Input:
        % imageData - [xData, yData, slices, frames]
        % rawObj - meta data object with structs
        % destination - name of destination path 
    
    data = load(fullfile(imagePath, 'imageData.mat'));
    imageData = data.final_im;
    
    %% Initializing metadata structs
    visuParam                           = readBrukerParamFile(fullfile(rawObj.Filespath.auto,'\pdata\1\visu_pars'));
    acqp                                = rawObj.Acqp;
    method                              = rawObj.Method;

    %% Sorting slices in correct order if slice
    sliceNum = size(imageData,3);
    
    if sliceNum > 1
        imageData = sliceShuffler(imageData);
    end

    %% Making DICOM files per slice in imageData
    for slice=1:sliceNum
        

        %% Orientation fix
        sliceData = orientRotation(imageData(:,:,slice,:), rawObj, visuParam);

        %% Filename
        if sliceNum > 1
            sliceName                           = append(destination, '_', num2str(slice));
        else 
            sliceName                          = destination;
        end

        
        %% Initializing file and metadata
        try 
            dicomwrite(sliceData,[sliceName,'.dcm'])
        catch
            disp('-----------------')
            disp(['Problem initializing DICOM file for ', sliceName, '.'])
            disp('Possibly corrupted file.')
            disp('-----------------')
            continue
        end
        info                                = dicominfo([sliceName,'.dcm']);
    
        %% Geometrical information
        info.SliceThickness                 = method.PVM_SliceThick;
        position                            = method.PVM_EffSliceOffset(slice);
        info.SliceLocation                  = position;
       
        matrixFOV                           = [method.PVM_DefMatrix(1), method.PVM_DefMatrix(2)]; % seems to be only parameter that is consistent for both CS undersampled images and partial echo images
        sizeFOV                             = visuParam.VisuCoreExtent;
        spatialResolution                   = sizeFOV ./ matrixFOV;
        info.PixelSpacing                   = spatialResolution;   
        
        info.ImagePositionPatient           = visuParam.VisuCorePosition(slice,1:3);
        info.ImageOrientationPatient        = visuParam.VisuCoreOrientation(slice,1:6);
        
        %% Multi-frame metadata
        info.NumberOfFrames                 = size(sliceData,4);
    
        %% General metadata
        info.PatientID                      = visuParam.VisuSubjectId;
        info.PatientName.FamilyName         = visuParam.VisuSubjectId;
        info.HeartRate                      = 60/((acqp.ACQ_repetition_time)/1000); % estimation
        info.ImageType                      = 'ORIGINAL\PRIMARY\OTHER';
        info.Modality                       = 'MR';
        info.ScanningSequence               = 'RM\GR';
        info.SequenceVariant                = 'SP';
        info.MRAcquisitionType              = '2D';
        info.InPlanePhaseEncodingDirection  = 'ROW';
        info.ProtocolName                   = visuParam.VisuAcquisitionProtocol;
        info.AcquisitionMatrix              = [0; matrixFOV(1); matrixFOV(2); 0];
        info.AnatomicalOrientation          = 'QUADRUPED';
    
        %% Saving DICOM file with info
        dicomwrite(sliceData,[sliceName,'.dcm'], info,'CreateMode','Copy');
        
    end
end