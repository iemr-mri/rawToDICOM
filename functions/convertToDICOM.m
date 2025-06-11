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

    %% Normalize image
    % Calculate normalization factor to scale maximum intensity to 30,000
    normFactor                          = 30000/max(imageData,[],'all');
    % Normalize image data
    Inorm                               = normFactor .* imageData;
    Inorm                               = int16(Inorm);
    
    % if contains(visuParam.VisuAcqSequenceName, 'LAX')
    %     Inorm                           = flip(Inorm,2);
    % end

    %% Initializing DICOM file and info struct
    try 
        dicomwrite(Inorm,[destination,'.dcm'])
    catch
        disp('-----------------')
        disp(['Problem initializing DICOM file for ', destination, '.'])
        disp('Possibly corrupted file.')
        disp('-----------------')
        return
    end
    info                                = dicominfo([destination,'.dcm']);

    %% Geometrical information
    info.SliceThickness                 = acqp.ACQ_slice_thick;
    info.SliceLocation                  = acqp.ACQ_slice_offset;
    
    % Only use second dimension for matrixFOV since CS has half phase steps
    matrixFOV                           = [visuParam.VisuCoreSize(2), visuParam.VisuCoreSize(2)];
    sizeFOV                             = visuParam.VisuCoreExtent;
    spatialResolution                   = sizeFOV ./ matrixFOV;
    info.PixelSpacing                   = spatialResolution;   
    %voxelResolution                    = [spatialResolution, visuParam.VisuCoreFrameThickness];
        
    affineMatrix                        = build_affine(visuParam, spatialResolution);
    [imageMat, imagePos]                = to_matvec(affineMatrix);
    imageMat                            = imageMat * diag((1/spatialResolution(1))*[1;1;1]);
    imageOrientation                    = imageMat(:);
    info.ImagePositionPatient           = imagePos;
    info.ImageOrientationPatient        = imageOrientation(1:6);

    info.InPlanePhaseEncodingDirection  = 'ROW';

    
    %% Multi-frame metadata
    info.NumberOfFrames                 = size(imageData,4);

    %% General metadata
    info.PatientID                      = visuParam.VisuSubjectId;
    info.PatientName.FamilyName         = visuParam.VisuSubjectId;
    info.HeartRate                      = 60/((acqp.ACQ_repetition_time)/1000);
    info.ImageType                      = 'ORIGINAL\PRIMARY\OTHER';
    info.Modality                       = 'MR';
    info.ScanningSequence               = 'RM\GR';
    info.SequenceVariant                = 'SP';
    info.MRAcquisitionType              = '2D';
    info.InPlanePhaseEncodingDirection  = 'ROW';
    info.AcquisitionMatrix              = [0; spatialResolution(1); spatialResolution(2); 0];
    info.AnatomicalOrientation          = 'QUADRUPED';

    %% Saving DICOM file with info
    dicomwrite(Inorm,[destination,'.dcm'], info,'CreateMode','Copy');
    
end