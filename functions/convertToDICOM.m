function convertToDICOM(imagePath,rawObj,destination)
    % Converts image data (kspace) into DICOM with metainfo from rawObj and saves at a specified destination.
    % Input:
        % imageData - [xData, yData, slices, frames]
        % rawObj - meta data object with structs
        % destination - name of destination path 
    
    data = load(fullfile(imagePath, 'imageData.mat'));
    imageData = data.imageData;
    %% Image corrections
    imageData    = imageCorrections(imageData, rawObj);
    
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
        sliceData = orientRotation(imageData(:,:,slice,:), rawObj);

        %% Filename
        if sliceNum > 1
            sliceName                           = append(destination, '_', num2str(slice));
        else 
            sliceName                          = destination;
        end

        
        %% Initializing file and metadata
        try
            dicomwrite(sliceData, [sliceName, '.dcm']);
            info  = dicominfo([sliceName,'.dcm']);
        catch ME
            fprintf('-----------------\n');
            fprintf('Problem initializing DICOM file for %s.\n', sliceName);
            fprintf('Error: %s\n', ME.message);
            fprintf('-----------------\n');
            continue
        end
        
    
        %% Geometrical information
        info.SliceThickness                 = method.PVM_SliceThick;
        position                            = method.PVM_EffSliceOffset(slice);
        info.SliceLocation                  = position;
       
        % matrixFOV                           = method.PVM_EncMatrix;
        % matrixFOV(1)                        = matrixFOV(1)*method.PVM_EncPft(1); %partial echo adjustment
        % if isfield(rawObj.Method, "CSPhaseEncList")
        %     matrixFOV(2)                        = matrixFOV(2)*method.CSacceleration; % compressed sensing adjustment
        % end
        % imageSize                           = [size(imageData,1), size(imageData,2)];
        % if isequal(imageSize, 2*matrixFOV)
        %     matrixFOV = imageSize;
        % end

        imageSize                           = [size(imageData,1), size(imageData,2)];
        matrixFOV                           = imageSize;
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

        try
            info.HeartRate = data.heartRate; % estimation of heart rate from self-gating
        catch
            info.HeartRate                  = 60/((acqp.ACQ_repetition_time)/1000); % estimation based on TR
        end

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