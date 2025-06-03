function convertToDICOM(imageData,rawObj,destination)
    % Converts image data (kspace) into DICOM with metainfo from rawObj and saves at a specified destination.
    % Input:
        % imageData - [xData, yData, slices, frames]
        % rawObj - meta data object with structs
        % destination - name of destination path 

    %% Normalize image
    % Calculate normalization factor to scale maximum intensity to 30,000
    normFactor                          = 30000/max(imageData,[],'all');
    % Normalize image data
    Inorm                               = normFactor .* imageData;
    Inorm                               = int16(Inorm);
    
    % if contains(destination, 'LAX')
    %     Inorm                           = flip(Inorm,2);
    % else
    %     Inorm                           = flip(Inorm,1);
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
    acqp                                = rawObj.Acqp;

    %% Initializing metadata structs
    visuParam                           = readBrukerParamFile(fullfile(rawObj.Filespath.auto,'\pdata\1\visu_pars'));

    %% LAX vs SAX orientation
    % if contains(destination, 'LAX')
    %    orientationVector = [visuParam.VisuCoreOrientation(1), visuParam.VisuCoreOrientation(2), visuParam.VisuCoreOrientation(3),...
    %         visuParam.VisuCoreOrientation(4), visuParam.VisuCoreOrientation(5), visuParam.VisuCoreOrientation(6)];
    % else
    %     orientationVector = [visuParam.VisuCoreOrientation(4), visuParam.VisuCoreOrientation(5), visuParam.VisuCoreOrientation(6),...
    %         visuParam.VisuCoreOrientation(1), visuParam.VisuCoreOrientation(2), visuParam.VisuCoreOrientation(3)];
    % end

    %% Geometrical information
    info.SliceThickness                 = acqp.ACQ_slice_thick;
    info.SliceLocation                  = acqp.ACQ_slice_offset;

    pixels                              = rawObj.Method.PVM_EncMatrix;
    roFOV                               = acqp.ACQ_fov(1);
    info.PixelSpacing(1, 1)             = roFOV*10/pixels(1);
    info.PixelSpacing(2, 1)             = roFOV*10/pixels(1);
    
    ImagePos                            = [visuParam.VisuCorePosition(1), visuParam.VisuCorePosition(2), visuParam.VisuCorePosition(3)];
    orientationVector                   = [visuParam.VisuCoreOrientation(1), visuParam.VisuCoreOrientation(2), visuParam.VisuCoreOrientation(3),...
                                           visuParam.VisuCoreOrientation(4), visuParam.VisuCoreOrientation(5), visuParam.VisuCoreOrientation(6)];
    
    info.ImagePositionPatient           = ImagePos(:);
    info.ImageOrientationPatient        = orientationVector(:);

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
    info.AcquisitionMatrix              = [0; pixels(1); pixels(1); 0];
    info.ProtocolName                   = 'SegFLASH';
    info.AnatomicalOrientation          = 'QUADRUPED';
    info.PatientOrientation             = 'F\A';

    %% Saving DICOM file with info
    dicomwrite(Inorm,[destination,'.dcm'], info,'CreateMode','Copy');
    
end