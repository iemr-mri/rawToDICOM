function cineToDicom(magnitude, destination, expName, scanNumbers)
% Takes in our 4D magnitude cine images [x y cardiacPhase slice] and the
% desired location for dicoms and creates dicoms.

% K: This function is mostly GPT, as the dicomwrite is a bit complicated.

% K: Importantly, GPT does not know out how to avoid the "photometric
% interpretation" error I initially got when sending [x y cardiacPhase]
% slices into the dicomwriter. I found that reshaping to 
% [x y 1 cardiacPhase] avoids this, as the dicomwriter uses the third 
% dimension to check something, and therefore crashes if size(..., 3) = 40 
% and not 1.

% K: Also, the dicomwrite gives these warnings that clog the commmand
% window

tic;
disp('Section 12: Cine to dicom')
% Input: cine4D [rows cols phases slices]
[rows, cols, numSlices, cardiacPhases] = size(magnitude);

% scale/convert to integer dtype suited for DICOM (example: uint16)
% adjust scaling to your data range; here we map min..max to 0..4095 (12-bit)
minv = min(magnitude(:));
maxv = max(magnitude(:));
if maxv > minv
    scaled = uint16( round( (magnitude - minv) / (maxv - minv) * 4095 ) );
else
    scaled = uint16( magnitude ); % all constant
end

% Create base metadata (use a real template DICOM if available)
info                        = struct();
info.PatientName.FamilyName = 'Rat';
info.PatientID              = 'P001';
info.Modality               = 'MR';
info.StudyInstanceUID       = dicomuid;
info.SeriesInstanceUID      = dicomuid; % same for this series
info.SOPClassUID            = '1.2.840.10008.5.1.4.1.1.4'; % MR Image Storage (non-enhanced)

% Pixel/format tags
info.Rows                   = rows;
info.Columns                = cols;
info.SamplesPerPixel        = 1;
info.PhotometricInterpretation = 'MONOCHROME2';
info.BitsAllocated          = 16;
info.BitsStored             = 12;   % adjust to actual bit depth used
info.HighBit                = 11;
info.PixelRepresentation    = 0; % unsigned
info.NumberOfFrames         = cardiacPhases;
% Spatial tags
info.PixelSpacing = [1.25; 1.25];   % [row mm; col mm] — set to your acquisition values
info.SliceThickness = 8.0;          % mm
info.ImageOrientationPatient = [1 0 0 0 1 0]; % row/col direction cosines

% Choose starting slice position and z-step (example)
startZ = 0.0;   % mm
zstep  = info.SliceThickness;  % if contiguous

% Check if the folder is there
DICOMfolder = [destination expName];
if ~exist(DICOMfolder,'dir') == true
    mkdir(DICOMfolder);
end

% Prevent annoying warning
stuipdWarning   = 'images:dicomwrite:inconsistentIODAndCreateModeOptions';
oldWarnState    = warning('query', stuipdWarning);
warning('off', stuipdWarning);

% Loop over slices and write one multiframe file per slice
for slice = 1:numSlices
    frames = squeeze(scaled(:,:,slice,:)); % rows x cols x nphases
    frames = reshape(frames, [rows cols 1 cardiacPhases]);
    info.InstanceNumber = slice;
    % set unique SOPInstanceUID per file
    info.SOPInstanceUID = dicomuid;
    % ImagePositionPatient: X,Y,Z location of top-left pixel (3 values)
    info.ImagePositionPatient = [0; 0; startZ + (slice-1)*zstep];
    % Optionally set a timing tag for frames: FrameTime (ms) or TriggerTime
    % Many viewers expect timing info elsewhere; you can set:
    info.FrameTime = 25; % milliseconds between frames (if known)
    % Write multiframe DICOM: dicomwrite accepts 3-D array as multiframe
    fname = fullfile(destination, expName, sprintf('%d.dcm', scanNumbers(slice)));
    dicomwrite(frames, fname, info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);
    % dicomwrite(frames,fname,'MultiframeSingleFile',true)
end

% Restore usual warnings
warning(oldWarnState.state, stuipdWarning);

disp(['Section 12: Cine to dicom completed in ' num2str(toc) ' seconds.'])