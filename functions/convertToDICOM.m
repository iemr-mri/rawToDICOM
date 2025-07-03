function convertToDICOM(imageData,rawObj,destination)
    % Converts image data (kspace) into DICOM with metainfo from rawObj and saves at a specified destination.
    % Input:
        % imageData - [xData, yData, slices, frames]
        % rawObj - meta data object with structs
        % destination - name of destination path 
    
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
    
    % Only use second dimension for matrixFOV since CS has half phase steps
    matrixFOV                           = [visuParam.VisuCoreSize(2), visuParam.VisuCoreSize(2)];
    sizeFOV                             = visuParam.VisuCoreExtent;
    spatialResolution                   = sizeFOV ./ matrixFOV;
    info.PixelSpacing                   = spatialResolution;   
    %voxelResolution                    = [spatialResolution, visuParam.VisuCoreFrameThickness];
        
    affineMatrix                        = build_affine(visuParam, method, spatialResolution);
    [imageMat, imagePos]                = to_matvec(affineMatrix);
    imageMat                            = imageMat * diag((1/spatialResolution(1))*[1;1;1]);
    imageOrientation                    = reshape(imageMat,1,[]);

    % sliceGeo                            = method.PVM_SliceGeo;
    % [imagePos, imageOrientation]        = sliceGeometryParser(sliceGeo);


    info.ImagePositionPatient           = imagePos;

    % info.ImagePositionPatient           = visuParam.VisuCorePosition(1:3);

    if contains(visuParam.VisuAcquisitionProtocol, 'LAX')
        %disp('LAX orientation correction.')
        info.ImageOrientationPatient        = [imageOrientation(1:3),imageOrientation(4:6)];
    else
        info.ImageOrientationPatient        = imageOrientation(1:6);
    end
    

    % zdir = cross(visuParam.VisuCoreOrientation(1:3),visuParam.VisuCoreOrientation(4:6))';
    % pos = visuParam.VisuCorePosition(1:3)' + ...
    %     spatialResolution(1)*visuParam.VisuCoreOrientation(4:6)'+...
    %     spatialResolution(1)*visuParam.VisuCoreOrientation(1:3)'-...
    %     acqp.ACQ_slice_thick*zdir;

    info.InPlanePhaseEncodingDirection  = 'ROW';

    
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
    info.ProtocolName                   = 'SegFLASH';
    info.AcquisitionMatrix              = [0; 128; 128; 0];
    info.AnatomicalOrientation          = 'QUADRUPED';

    %% Saving DICOM file with info
    dicomwrite(imageData,[destination,'.dcm'], info,'CreateMode','Copy');
    
end