function [pm, midlines, rawWithMid, CSvector] = CSSGcineReader(folderName, expName, scanName, scanPartioner, selectedPortion, visualSwitch)
% Reads compress sensing, self-gated rawdata, creates mag videos and returns midlines for later PCA
% processing.

tic;

% Load scan, setup parameters, matrices and animation folder
pathName            = [folderName expName '\' char(scanName)];
rawScan             = RawDataObject(pathName, 'dataPrecision', 'double');

pm.folderName       = folderName;
pm.expName          = expName;
pm.scanNumber       = num2str(scanName);
pm.saveFolder       = ['CSSG reconstructed data\' pm.expName '\' pm.scanNumber];
pm.CSacceleration   = rawScan.Method.CSacceleration;
pm.midlineRate      = rawScan.Method.MidlineRate;
pm.sliceNum         = rawScan.Method.PVM_SPackArrNSlices;
pm.coilNum          = rawScan.Method.PVM_EncNReceivers;
pm.movieFrames      = rawScan.Method.PVM_NMovieFrames; % number of cardiac phases
pm.scanPartioner    = scanPartioner;
pm.selectedPortion  = selectedPortion;
if scanPartioner == 1
    pm.repetitions  = rawScan.Method.PVM_NRepetitions;
else
    pm.repetitions  = rawScan.Method.PVM_NRepetitions/pm.scanPartioner * (selectedPortion(2)-selectedPortion(1)); % Use scanTimeReducer to artificially shorten dataset.
end
pm.acquiredMatrix   = rawScan.Method.PVM_EncMatrix;
pm.partialEcho      = rawScan.Method.PVM_EncPft(1);
pm.TR               = rawScan.Method.FrameRepTime;
pm.kxPixels         = pm.acquiredMatrix(1);
pm.kyLines          = pm.acquiredMatrix(2);
pm.midlineNum       = floor(pm.movieFrames/pm.midlineRate) * pm.kyLines * pm.repetitions; % Number of midlines recorded in total
pm.actualMatrix     = [pm.kxPixels * pm.partialEcho, pm.CSacceleration * pm.kyLines];
pm.cardiacFreq      = 'Unset'; % Will be set during curveCleaner.m
pm.breathingFreq    = 'Unset'; % Will be set during curveCleaner.m
pm.cutoffSTD        = 1.5; % Used to set cutoff peak height in the peak finder. Low number throws out more peaks.
pm.phaseOffset      = rawScan.Acqp.ACQ_phase1_offset;
pm.FoV              = rawScan.Method.PVM_Fov(1); % Assumed square. In millimeters
pm.resolution       = pm.FoV/pm.actualMatrix(1);
pm.MechanicalPhases = 1; % Artificially set MechanicalPhases to 1 to facilitate use of same functions as MRE self-gating.

CSvector        = round(rawScan.Method.CSPhaseEncList * pm.actualMatrix(2)/(2*pm.CSacceleration)); % Go from [-accel, accel] to [-47, 48]. Round to force int
CSvector        = repmat(CSvector,1,pm.repetitions);

% Set ratOrMouse
if pm.coilNum == 4
    pm.ratOrMouse  = "rat";
elseif pm.coilNum == 2
    pm.ratOrMouse = "mouse";
end

% Setup animation folders
pm.animationFolder = 'Animations/';
if ~(exist(pm.animationFolder, 'dir') == 7)
    mkdir(pm.animationFolder);
end

if ~(exist([pm.animationFolder pm.expName], 'dir') == 7)
    mkdir([pm.animationFolder pm.expName]);
end

% Partition scan
rawWithMid      = squeeze(rawScan.data{1});
if scanPartioner > 1
    portionLength   = size(rawWithMid,3)/pm.scanPartioner;
    portionStart    = 1 + pm.selectedPortion(1) * portionLength;
    portionEnd      = pm.selectedPortion(2) * portionLength;
    rawWithMid      = squeeze(rawWithMid(:,:,portionStart:portionEnd));
end

coilRawWithMid  = zeros(pm.coilNum, pm.kxPixels, pm.kyLines, pm.movieFrames, pm.repetitions);
coilRawNoMid    = zeros(pm.coilNum, pm.kxPixels, pm.kyLines, pm.movieFrames, pm.repetitions);

% Show scan time and midlinerate
pm.scanTime         = pm.movieFrames * pm.kyLines * pm.repetitions * pm.TR / 1000;
disp(['Total scan time: ' num2str(pm.scanTime) 's.']);
pm.midlinePeriod    = pm.midlineRate*pm.TR;
disp(['Midlinerate: ' num2str(pm.midlineRate) ' (' num2str(pm.midlinePeriod) ' ms period)']);

% Calculate in what positions we will find midlines for each set of frames
midlinePositions    = [];
countUp             = pm.midlineRate;
while countUp <= pm.movieFrames 
    midlinePositions    = [midlinePositions; countUp];
    countUp             = countUp + pm.midlineRate;
end

% Visualize rawWithMid
if visualSwitch == true
    figure(1);
    imagesc(abs(squeeze(rawWithMid(1,:,1:pm.movieFrames*pm.kyLines)))'); % Arbitarily chose coil 1
    title("Rawdata unshuffled");
end

% Shuffle data into frame-by-frame and rep-by-rep images
tempRaw = reshape(rawWithMid, [pm.coilNum, pm.kxPixels, pm.movieFrames, pm.kyLines, pm.repetitions]); % [pm.kxPixels, pm.movieFrames, pm.kyLines, pm.repetitions]);
tempRaw = permute(tempRaw, [1,2,4,3,5]); % [coils, Pixels, lines, movieFrames, repetitions]
coilRawWithMid = tempRaw;

% Sort midlines and other lines into separate matrices
midlines        = zeros(pm.coilNum, pm.kxPixels, pm.midlineNum);
midlinesCounter = 1;

for rep = 1:pm.repetitions
    for line = 1:pm.kyLines
        for frame = 1:pm.movieFrames
            if ismember(frame,midlinePositions)
                midlines(:,:,midlinesCounter)   = squeeze(coilRawWithMid(:,:,line,frame,rep));
                midlinesCounter                 = midlinesCounter + 1;
            else
                coilRawNoMid(:,:,:,frame,rep) = coilRawWithMid(:,:,:,frame,rep);
            end
        end
    end
end

if visualSwitch == true
    for frame = 1:pm.movieFrames
        for rep = 1:pm.repetitions            
            figure(2);
            imagesc(abs(squeeze(coilRawWithMid(1,:,:,frame,rep)))');
            title(['CoilRawWithMid - Frame: ' num2str(frame) ' Rep: ' num2str(rep)]);
        
            figure(3);
            imagesc(abs(squeeze(midlines(1,:,:)))');
            title("midlines");

            % pause(0.5);
        end
    end
end

coilMagImages   = zeros(pm.coilNum, pm.kxPixels, pm.kyLines, pm.movieFrames, pm.repetitions);
totalMagImage   = zeros(pm.kxPixels, pm.kyLines, pm.movieFrames);

% Create images
for rep = 1:pm.repetitions
    for frame = 1:pm.movieFrames
        for coil = 1:pm.coilNum
            coilMagImages(coil,:,:,frame,rep)   = abs(ifftshift(ifft2(fftshift(squeeze(coilRawNoMid(coil,:,:,frame,rep))))));
            totalMagImage(:,:,frame)            = squeeze(totalMagImage(:,:,frame)) + squeeze(coilMagImages(coil,:,:,frame,rep));
    
            % if visualSwitch == true
            %     figure(4);
            %     imagesc(squeeze(coilMagImages(coil,:,:,frame,rep)));
            %     titleString = ["Coil: " num2str(coil) " Frame: " num2str(frame) " Rep: " num2str(rep)];
            %     title(titleString);
            %     pause(0.05);
            % end
        end
    end
end

if visualSwitch == true
    figure(5);
    imagesc(squeeze(totalMagImage(:,:,1)));
    title("Combined");
end

disp(['Section 1: SGreader - Finished in ' num2str(toc) ' seconds.'])
end % Function end