% Self-gating of MRI images using PCA on k-space centerlines. No ECG.

% 2024 November,    Kasper: PPU prototype
% 2025 February,    Kasper: Actual self-gating
% 2025 March,       Kasper: CompressedSensing built on Gary and Emil's old script
% 2025 July,        Kasper: Merged Cine and MRE SG pipelines
% 2026 January,     Kasper: Added slice synchronizer and DICOM export, LAX
%                           and SAX flexibility

% Section 1: Read file
% -Loads specified file
% -Retrieves parameters of interest from metadata
% -Sorts midlines and actual acquisitions

% Section 2: PCA analysis
% -For each acquired midline it combines imag and real for all coils into
%   1 long vector
% -Runs PCA on that vector using each acquisition as a timepoint.
% -Retrieves first and second prinicpal component, which should represent 
% cardiac and breathing motion.

% Section 3: Time correction
% -Expands the timeline and interpolates between points of the PCA curve to
%  avoid post-processing issues related to having few timepoints per
%  heartbeat

% Section 4: PCA curve clean up
% -Takes in the curves of the (e.g.) 10 first principal components.
% -Selects cardiac and breathing PCs.
% -Uses butterworth bandpass filter centered around the cardiac and 
%   breathing frequencies to clean up signal.

% Section 5: Peak detection
% -Takes in cleaned PC curves 
% -Calculates peak detection parameters
% -Finds "rough" cardiac peaks by looking for rises followed by drops and 
%   finding the max value within that area
% -Finds "fine" cardiac peaks by using the rough peaks as initial guesses 
%   for a curve fit. Curve fit uses skewed gaussian shape.
% -Finds breathing peaks and determines what is the start of the breath 
%   cycle.

% Section 6: Rolling window visualization
% -Takes in cardiac and breathing curves.
% -Rolls through the data rendering them on top of each other

% Section 7: Data shuffler
% -Takes in rawdata iterates through every timepoint in the scan.
% -Sorts out timepoints that are midlines, during breathing or not within 
%   good cardiac peaks.
% -Bins remaining acquisitions into appropriate frames
% -Averages the incomplete k-spaces and shows the resulting images.

% Section 8: CS reconstruction
% -Takes in correctly shuffled data
% -Loops through all videos [x y movieFrames]
% -Temporal FFT. Take abs(). Threshold out lower values, which will likely 
%   be noise.
% -iFFT. Set non-missing lines back to their acquired values. Compare 
%   difference before and after reintroduction.
% -Repeat until the difference is low.

% Section 9: Zero-filling interpolation
% -Takes in reconstructed data from CS.
% -Increases "resolution" by expanding k-space with zeros
% -Effectively interpolates between points

% Section 10: Data saver
% -Takes in the ZIPed data
% -Saves it in a local folder to skip CS recon in future analysis

% Section 11: Slice synchronizer
% -Takes in ZIPed SAX data
% -Has the user select the heart area and some myocardium
% -For the middle slice, it determines diastole by looking for the
%  timepoint with the least "air" around the heart
% -Other slices are synchronized from the middle slice and out by checking
%  which circshift to the cardiacPhase gives the least difference between
%  the slice in question and the neighbour slice

function magnitude = SGcine(folderStructs)
%% Section 0: Basic parameters
disp('________ Self-Gating - New Run ________');

% State scan(s) to analyze
[folderName, expName]       = fileparts(folderStructs(1).folder); % All folderStructs have same folder and expName
folderName(end+1)           = '\'; % Script expects folderName to end on '\'
[SAXnames, LAXnames]        = scanNameSplitter(folderStructs);
scanNames                   = [SAXnames, LAXnames];

% Set pipeline parameters
reconFinished   = false; % If we've already reconstructed the files, just load them and skip to slice synchronizer
visualSwitch    = false;
showRolling     = false;
redoSync        = false;

%% Main section
if reconFinished == false
    for scan = 1:length(scanNames)
        % Check if this slice is already reconstructed
        % Commented out: currently just making this decision for all slices beforehand 
        % saveName = [folderName, expName, '\', char(scanNames(scan)), '\imageData.mat'];
        % if exist(saveName, 'file')
        %     continue; % Skip to next scan
        % end
        
        %% Section 1: Read data
        longTic         = tic;      % Count the whole reconstruction process
        scanPartioner   = 1;        % How many sections should the scan be divided into. Useful for comparison of fully vs not-fully sampled
        selectedPortion = [0 1];    % Starts from 0
        [pm, midlines, rawVector, CSvector] = CSSGcineReader(folderName, expName, scanNames(scan), scanPartioner, selectedPortion, false);
        
        %% Section 2: PCA analysis
        [PCvecs, PCscores, PCexplained] = PCArunner(midlines, pm, true);
        
        %% Section 3: Time correction
        timeStamps      = (0:(length(rawVector)-1))*pm.TR*1e-3; % We create some fake timestamps just to get the interpolated timeline
        pm.PCselection  = 1:10;
        selectPCs       = PCscores(pm.PCselection,:);
        [interpolated, timeline, pm] = timeCorrecter(selectPCs, timeStamps, pm, false);
        
        %% Section 4: PCA curve clean up
        pm.showPCs = false;
        [PCcardiac, cardiacClean, PCbreath, breathClean, pm] = curveCleaner(interpolated, PCexplained, pm, false);
        
        %% Section 5: Peak detection
        pm.diffOrWidth      = "width"; % Set whether breath peaks are found by sharpest differential or lowest width.
        pm.breathFlipBool   = false;   % Set whether breath peaks are flipped or not 
        pm.cardiacFlipBool  = false;   % Set whether cardiac peaks are flipped or not
        pm.performCurveFit  = false;
        [roughPeakList, finePeakList, breathData, breathStarts, pm] = PCApeakFinder(cardiacClean, breathClean, timeStamps, pm, false);
        
        %% Section 6: Rolling window visualization
        pm.breathRange = [0.25,0.75];
        if showRolling == true
            showRaw = true;
            rollingWindow(PCcardiac, cardiacClean, PCbreath, breathClean, roughPeakList, finePeakList, breathStarts, timeline, pm, showRaw);
        end
        
        %% Section 7: Data shuffling
        pm.breathRange  = [0.25,0.75];
        pm.newFrameNum  = 40;
        pm.beatTolerance = 0.05; % Discards heart beats with durations deviating from the median
        pm.binTolerance = 1; % Set percentage acceptance around each frame bin
        pm.fillBlanks   = false;
        pm.showFilling  = false;
        pm.showCSdist   = false;
        pm.absOrRel     = "abs";
        [rawShuffled, averagedShuffled] = dataShuffler(rawVector, finePeakList, breathStarts(1,:), CSvector, timeStamps, pm, false);
        
        %% Section 8: CS reconstruction
        pm.maxIter          = 50;
        pm.CSprctThresh     = 50;
        pm.CSdiffThresh     = 0.01;
        pm.selectFrame      = 3;
        pm.showCSprocess    = false;
        [reconkSpace, combReconkSpace, combRealSpace] = CSreconstructor(averagedShuffled, pm, false);
        
        %% Section 9: Zero-filling interpolation
        if ~exist('magnitude','var')
            magnitude = zeros(2*pm.actualMatrix(1), 2*pm.actualMatrix(2), length([SAXnames, LAXnames]), pm.newFrameNum);
        end
        
        [~, magnitude(:,:,scan,:)] = zipper(combRealSpace, pm, false);
        
        %% Section 10: Datasaver
        imageData = squeeze(magnitude(:,:,scan,:));
        imageData = reshape(imageData, [size(imageData,1) size(imageData,2) 1 size(imageData,3)]);
        heartRate = pm.cardiacFreq*60;
        save(saveName,"imageData", "heartRate");
        disp(['Reconstruction of ' pm.expName ' E' pm.scanNumber ' completed in ' num2str(toc(longTic)) ' seconds.'])
    end
end

%% Section 11: Slice synchronizer
% if reconFinished == true || redoSync % We load stored data instead
%     unsyncSAX = false;
% else
%     unsyncSAX = squeeze(magnitude(:,:,1:length(SAXnames),:));
% end

%always load imageData at sliceSynchronizer
unsyncSAX = false;

if (~isempty(SAXnames) & exist('magnitude', 'var')) || redoSync % Synchronization only performed on SAX
    [syncSAX, rotDiffs, minDiff]            = sliceSynchronizer(unsyncSAX, expName, SAXnames, folderName, true);
    % overwrite imageData after slice synchronizer
    for SAXscan = 1:length(SAXnames)
        saveName  = [folderName, expName, '\', char(SAXnames(SAXscan)), '\imageData.mat'];
        imageData = squeeze(syncSAX(:,:,SAXscan,:));
        imageData = reshape(imageData, [size(imageData,1) size(imageData,2) 1 size(imageData,3)]);
        save(saveName,"imageData", "-append");
    end
end

