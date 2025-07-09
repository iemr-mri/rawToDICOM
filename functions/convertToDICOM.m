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

    %% Initializing DICOM file and info struct
    try 
        dicomwrite(imageData,[destination,'.dcm'])
    catch
        disp('-----------------')
        disp(['Problem initializing DICOM file for ', destination, '.'])
        disp('Possibly corrupted file.')
        disp('-----------------')
        return
    end
    info                                = dicominfo([destination,'.dcm']);

    %% Geometrical information
    info.SliceThickness                 = method.PVM_SliceThick;
    position                            = method.PVM_SPackArrSliceOffset;
    info.SliceLocation                  = position;
   
    matrixFOV                           = [size(imageData,1), size(imageData,2)];
    sizeFOV                             = visuParam.VisuCoreExtent;
    spatialResolution                   = sizeFOV ./ matrixFOV;
    info.PixelSpacing                   = spatialResolution;   
    
    % Affine matrix transformation
    affineMatrix                        = build_affine(visuParam, method, spatialResolution);
    [imageMat, imagePos]                = to_matvec(affineMatrix);
    imageMat                            = imageMat * diag((1/spatialResolution(1))*[1;1;1]);
    imageOrientation                    = reshape(imageMat,1,[]);

    info.ImagePositionPatient           = imagePos;
    info.ImageOrientationPatient        = imageOrientation(1:6);
    
    %% Multi-frame metadata
    info.NumberOfFrames                 = size(imageData,4);

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
    info.AcquisitionMatrix              = [0; 128; 128; 0];
    info.AnatomicalOrientation          = 'QUADRUPED';

    %% Saving DICOM file with info
    dicomwrite(imageData,[destination,'.dcm'], info,'CreateMode','Copy');
    
end