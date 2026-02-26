function [rawShuffledData, averagedShuffledData] = dataShuffler(rawData, finePeakList, breathIndices, CSvector, timeStamps, pm, visualSwitch)
% Purpose
% -Uses detected cardiac and breathing peaks in the PCA analysis of 
%   ungated MRI acquisition to temporally sort the rawdata.
%
% Input
% 
% -rawData:         [coils xPixels yLines frames repetitions]
%   The shuffled (but not temporally sorted) raw MRI acquisition.
% 
% -finePeakList:    [index&height numMidlines]
%   List of all peaks detected in the PC-curve. Curvefitted for
%   sub-pixel temporal resolution.
%
% -breathTimes:    [numMidlines]
%   List of all breath starts (peak directly before breath).
%
% -pm:
%   Struct holding the parameters of the scan in question.
%
% -visualSwitch:
%   Switching visualizing parts of the function on and off.
%
% Procedure
% -Loop through each acquired line in the rawdata.
% -Skip timepoints with bad breathing.
% -Skip timepoints between missed cardiac peaks.
% -Find absolute or relative time of acquisition in relation to the
% previous peak (abs) or both neighbouring peaks (rel). Set new frame.
% -Sort midlines by frame and store for averaging.
% -Sort all other lines by frame into the shuffledData.
% -Average the image in each frame over repetitions.
% -Add averaged midlines to each frame.
% -Combine resulting coil kspaces into full, combined images.
% 
% -Also: Tons of visualizing and diagnostic statistics through the function.

disp('Section 7: dataShuffler');
tic;

% Calculate relative binSize
binSize         = 1/pm.newFrameNum;
cardiacPeriod   = 1/pm.cardiacFreq;
framePeriod     = cardiacPeriod*1e3/pm.newFrameNum; % In milliseconds

% Convert breath indices to times
breathTimes     = breathIndices*pm.newTempRes; % In milliseconds

% Filter out peaks that are too tight or widely spaced (missed peaks or breath)
finePeakTimes   = finePeakList(:,1) * pm.newTempRes;
diffList        = diff(finePeakTimes);
upThreshold     = median(diffList) * (1 + pm.beatTolerance);
downThreshold   = median(diffList) * (1 - pm.beatTolerance);
notTooHigh      = diffList < upThreshold;
notTooLow       = diffList > downThreshold;
goodPeaks       = notTooHigh .* notTooLow;

% Setup for looping through acquisitions
totalAcqs       = pm.movieFrames * pm.kyLines * pm.MechanicalPhases * pm.repetitions;
scanTime        = timeStamps(end);
rawShuffledData = zeros(pm.coilNum, pm.actualMatrix(1), pm.actualMatrix(2), pm.newFrameNum, pm.MechanicalPhases, pm.repetitions);
shuffledMids    = zeros(pm.coilNum, pm.actualMatrix(1), pm.midlineNum, pm.newFrameNum, pm.MechanicalPhases);
kxEchoStart     = pm.actualMatrix(1) * (1 - 1/pm.partialEcho) + 1; % Where to start filling in data due to partial echo
midCounter      = zeros(pm.newFrameNum, pm.MechanicalPhases);
frameCounter    = zeros(pm.newFrameNum, pm.MechanicalPhases);

% Diagnostic counters
goodAcqs        = 0;
startEndBreath  = 0;
badBreath       = 0;
startEndCardiac = 0;
badBeat         = 0;
badTolerance    = 0;
highFrame       = 0;
midSum          = 0;
doubleData      = 0;

% Setup figures for visualization of frame filling
frameFillFigNum = 300;
if pm.showFilling == true
    visualizeFrameFilling(frameFillFigNum);
end
CSGIFname = [pm.expName '/scanTime ' num2str(round(pm.scanTime))];
CSGIFname = [CSGIFname ' ' char(pm.absOrRel) ' breath[' num2str(pm.breathRange(1)) ',' num2str(pm.breathRange(2)) ']'];
CSGIFname = [CSGIFname ' frames ' num2str(pm.newFrameNum) ' CSdist.gif'];

% Go through all acquisitions, line by line, discard according to rules.
% Sort approved acquisitions into correct new frames.
for acq = 1:totalAcqs
    % Find coordinates and data of acquisition in question
    % With CS we have to retrieve the true "line" value of the acq from the
    % CS vector.
    [frame, ~, mechPhase, rep] = MREcoordinateCalculator(acq, pm);
    acqData         = squeeze(rawData(:,:,acq));
    line            = CSvector(acq) + pm.actualMatrix(2)/2; % Go from [-47,48] to [1,96]
    acqTime         = timeStamps(acq)*1e3; % Convert to milliseconds

    % Skip timepoints with bad breath. Exclude timepoints before the
    % first or after the last peak. Then find relative position of
    % acquisition between peaks and filter for accepted range.
    prevBreath = max(breathTimes(breathTimes < acqTime));
    if isempty(prevBreath) || prevBreath == max(breathTimes)
        startEndBreath = startEndBreath + 1;
        continue;
    else
        prevBreathIndex = find(breathTimes == prevBreath,1);
        nextBreath      = breathTimes(prevBreathIndex + 1);
        breathPos       = (acqTime - prevBreath)/(nextBreath-prevBreath);
        if breathPos < pm.breathRange(1) || breathPos > pm.breathRange(2)
            badBreath = badBreath + 1;
            continue;
        else
            % Acquisition is in good breath range.
        end
    end

    % Find the previous cardiac peak to the point
    earlierPeaks    = finePeakTimes(finePeakTimes < acqTime);
    lastIndex       = find(earlierPeaks, 1, 'last');
    
    % Skip if before first peak or after last peak
    if isempty(lastIndex) || lastIndex > length(goodPeaks)
        startEndCardiac = startEndCardiac + 1;
        continue;
    end

    % Check if previous cardiac peak is a good peak. Skip if not.
    if goodPeaks(lastIndex) == false
        badBeat = badBeat + 1;
        continue;
    end
    lastPeak        = finePeakTimes(lastIndex);

    % Use either absolute or relative frame binning
    if pm.absOrRel == "abs"
        % Absolute binning only relates to previous peak
        absAcqTime      = acqTime - lastPeak;
        % Check if the timepoint is outside bin tolerance
        binDiff         = mod(absAcqTime,framePeriod);
        if binDiff > framePeriod*pm.binTolerance
            badTolerance = badTolerance + 1;
            continue;
        end
        newFrame        = ceil(absAcqTime/framePeriod);
    elseif pm.absOrRel == "rel"
        % Relative binning finds the relative position between two peaks
        peakDistance    = diffList(lastIndex);
        relAcqTime      = (acqTime - lastPeak)/peakDistance;
        % Check if the timepoint is outside bin tolerance
        binDiff         = mod(relAcqTime,binSize);
        if binDiff > binSize*pm.binTolerance
            badTolerance = badTolerance + 1;
            continue;
        end
        newFrame        = ceil(relAcqTime/binSize);
    else
        disp("absOrRel was not set to either 'abs' or 'rel'. Exiting.")
        return;
    end

    % If newFrame is outside desired number of frames, skip it
    if newFrame > pm.newFrameNum
        highFrame = highFrame + 1;
        continue;
    end

    % Midlines are stored and counted separately
    if mod(acq,pm.midlineRate) == 0
        curPos                                                      = midCounter(newFrame,mechPhase);
        shuffledMids(:,kxEchoStart:end,curPos+1,newFrame,mechPhase) = acqData;
        midCounter(newFrame,mechPhase)                              = midCounter(newFrame,mechPhase) + 1;
        midSum                                                      = midSum + 1;
        continue;
    end
    
    % If the acquisition is still not discarded, use it
    % If this dataline is already filled, average the existing data with
    % the new. K: Doesn't average correctly for n > 2 in one bin, but
    % shouldn't be relevant for real recordings.
    if ~isequal(squeeze(rawShuffledData(:,:,line,newFrame,mechPhase,rep)), zeros(pm.coilNum,pm.actualMatrix(1)))
        rawShuffledData(:,kxEchoStart:end,line,newFrame,mechPhase,rep) = (rawShuffledData(:,kxEchoStart:end,line,newFrame,mechPhase,rep) + acqData)/2;
        doubleData = doubleData + 1;
    else
        rawShuffledData(:,kxEchoStart:end,line,newFrame,mechPhase,rep) = acqData;
        frameCounter(newFrame,mechPhase) = frameCounter(newFrame,mechPhase) + 1;
    end
    goodAcqs = goodAcqs + 1;
    
    % If we want to show the filling of frames
    if pm.showFilling == true
        disp([line frame newFrame mechPhase rep]);
        figure(frameFillFigNum+newFrame-1);
        imagesc(abs(squeeze(rawShuffledData(1,:,:,newFrame,mechPhase,rep)))');
        title(['Coil 1: Frame ' num2str(newFrame) ' mechPhase ' num2str(mechPhase) ' Rep ' num2str(rep)]);
    end

end % End acq loop
goodAcqs = goodAcqs - doubleData;

% Calculate occupancy
estimOccupied  = MREoccupancyEstimator(pm.actualMatrix(2), mean(frameCounter, 'all')) / pm.actualMatrix(2) * 100; % Assumes random distribution. Is quite accurate for low filling rates, then the Laplacian distribution skews it.
actualOccupied = MREoccupancyCalculator(rawShuffledData, pm.showCSdist, frameFillFigNum, CSGIFname) / pm.actualMatrix(2) * 100; % Occupancy in percentage

% Display diagnostic counters
fprintf('Shuffling complete with the following counters of loop outcome:\n \n');
fprintf('Total scan time:                   %.2fs\n', scanTime);
fprintf('Total acquisitions:                %d\n', totalAcqs);
fprintf('Unique accepted acquisitions:      %d (%.2f%%)\n', goodAcqs, 100*goodAcqs/totalAcqs);
fprintf('Before first or after last breath: %d\n', startEndBreath);
fprintf('Bad breath:                        %d\n', badBreath);
fprintf('Before first or after last beat:   %d\n', startEndCardiac);
fprintf('Bad beats:                         %d\n', badBeat);
fprintf('Too high newFrame:                 %d\n', highFrame);
fprintf('Midlines used in images:           %d\n', midSum);
fprintf('DoubleData:                        %d\n', doubleData);
fprintf('Acqs per frame:                    Mean: %.2f - STD: %.2f - Max: %.2f - Min: %.2f\n', mean(frameCounter(:)), std(frameCounter(:)), max(frameCounter(:)), min(frameCounter(:)));
fprintf('Theoretical filling of each frame: %.2f%% \n', estimOccupied);
fprintf('Actual filling per frame:          Mean: %.2f%% - STD: %.2f%% - Max: %.2f%% - Min: %.2f%%\n', mean(actualOccupied(:)), std(actualOccupied(:)), max(actualOccupied(:)), min(actualOccupied(:)));

% Average over repetitions
averagedkSpaces = zeros(pm.coilNum, pm.actualMatrix(1), pm.actualMatrix(2), pm.newFrameNum, pm.MechanicalPhases);
counterMask     = zeros(pm.coilNum, pm.actualMatrix(1), pm.actualMatrix(2), pm.newFrameNum, pm.MechanicalPhases);
for rep = 1:pm.repetitions
    for mechPhase = 1:pm.MechanicalPhases
        for frame = 1:pm.newFrameNum
            temp                                    = squeeze(rawShuffledData(:,:,:,frame,mechPhase,rep));
            averagedkSpaces(:,:,:,frame,mechPhase)  = averagedkSpaces(:,:,:,frame,mechPhase) + temp;
            counterMask(:,:,:,frame,mechPhase)      = squeeze(counterMask(:,:,:,frame,mechPhase)) + (temp ~= 0);
        end
    end
end
nonZeroMask = counterMask ~= 0;
averagedkSpaces(nonZeroMask) = averagedkSpaces(nonZeroMask) ./ counterMask(nonZeroMask);

% Include midlines
summedMids = squeeze(sum(shuffledMids, 3)); % Sums over duplicate midlines
for mechPhase = 1:pm.MechanicalPhases
    for frame = 1:pm.newFrameNum
        % Skip midlines that weren't found to avoid /0 and NaNs
        if midCounter(frame,mechPhase) == 0
            continue
        end
        averagedkSpaces(:,:,pm.actualMatrix(2)/2,frame,mechPhase) = squeeze(summedMids(:,:,frame,mechPhase)/midCounter(frame,mechPhase));
    end
end

% Set averaged return value
averagedShuffledData = averagedkSpaces;

% Combine coils for mag image
combinedMag = zeros(pm.actualMatrix(1), pm.actualMatrix(2), pm.newFrameNum);
for coil = 1:pm.coilNum
    for frame = 1:pm.newFrameNum
        for mechPhase = 1:pm.MechanicalPhases
            combinedMag(:,:,frame) = combinedMag(:,:,frame) + abs(squeeze(fftshift(ifftn(averagedkSpaces(coil,:,:,frame,mechPhase))))).^2;
        end
    end
end
combinedMag = sqrt(combinedMag);

disp(['Section 7: Datashuffler finished in ' num2str(toc) ' seconds.' newline])

% Visualize averaged data
imageVisFigNum     = 200;
combinedGIF = [pm.expName '/scanTime ' num2str(round(pm.scanTime))];
combinedGIF = [combinedGIF ' ' char(pm.absOrRel)];
combinedGIF = [combinedGIF ' breath[' num2str(pm.breathRange(1)) ',' num2str(pm.breathRange(2)) ']'];
combinedGIF = [combinedGIF ' frames ' num2str(pm.newFrameNum)];
combinedGIF = [combinedGIF '.gif'];

if visualSwitch == true && exist(combinedGIF, 'file')
    delete(combinedGIF);
    disp([combinedGIF ' deleted for new gif.'])
end

if visualSwitch == true
    averagedVisualizer(pm, imageVisFigNum, averagedkSpaces, combinedMag, combinedGIF);
end

end % End of dataShuffler


%% Support functions

function visualizeFrameFilling(startFigNum)
% Set up figures for filling of frames
% Created for up to 10 frames, but 10+ still runs
handle1 = figure(startFigNum);
set(handle1, 'Position', [0,575,475,400]);

handle2 = figure(startFigNum + 1);
set(handle2, 'Position', [475,575,475,400]);

handle3 = figure(startFigNum + 2);
set(handle3, 'Position', [950,575,475,400]);

handle4 = figure(startFigNum + 3);
set(handle4, 'Position', [1425,575,475,400]);

handle5 = figure(startFigNum + 4);
set(handle5, 'Position', [0,100,475,400]);

handle6 = figure(startFigNum + 5);
set(handle6, 'Position', [475,100,475,400]);

handle7 = figure(startFigNum + 6);
set(handle7, 'Position', [950,100,475,400]);

handle8 = figure(startFigNum + 7);
set(handle8, 'Position', [1425,100,475,400]);

handle9 = figure(startFigNum + 8);
set(handle9, 'Position', [950,1200,475,400]);

handle10 = figure(startFigNum + 9);
set(handle10, 'Position', [1425,1200,475,400]);
end % End showFrameFilling


function [frame, line, mechPhase, rep] = MREcoordinateCalculator(acq, pm)
% Calculate frame, line, mechPhase, and rep based on acquisition number 
acqPerRep   = pm.movieFrames * pm.kyLines * pm.MechanicalPhases;
acqPerMech  = pm.movieFrames * pm.kyLines;

% Handle reps
rep         = ceil(acq/acqPerRep);
acqNoReps   = mod(acq,acqPerRep);

% Special case: acq == acqPerRep. Don't roll over to acqNoReps == 0.
if acqNoReps == 0
    acqNoReps = acqPerRep;
end

% Handle mechPhases
mechPhase   = ceil(acqNoReps/acqPerMech);
acqNoMech   = mod(acqNoReps,acqPerMech);

if acqNoMech == 0
    acqNoMech = acqPerMech;
end

% Handle lines (doesn't make sense for CS, but idc)
line   = ceil(acqNoMech/pm.movieFrames);

% Special case: frame == pm.movieFrames. Don't roll over to frame == 0
if mod(acqNoReps,pm.movieFrames) == 0
    frame  = pm.movieFrames;
else
    frame  = mod(acqNoReps,pm.movieFrames);
end
end % End of coordinateCalculator



function unaveragedVisualizer(pm, startFigNum, shuffledData)
% Sets up 4x2 figures showing the kspaces of the 4 coils and the mag images
% of the same coils. All reps, not averaged.
timestep = 7.5/size(shuffledData,4); % Let cycle last 7.5s
for rep = 1:pm.repetitions
    for frame = 1:pm.newFrameNum
        handle1 = figure(startFigNum);
        set(handle1, 'Position', [0,575,475,400]);
        imagesc(abs(fftshift(ifftn(squeeze(shuffledData(1,:,:,frame,rep))))));
        title(['Coil 1: Frame ' num2str(frame) ' Rep ' num2str(rep)]);
        
        handle2 = figure(startFigNum + 1);
        set(handle2, 'Position', [475,575,475,400]);
        imagesc(abs(fftshift(ifftn(squeeze(shuffledData(2,:,:,frame,rep))))));
        title(['Coil 2: Frame ' num2str(frame) ' Rep ' num2str(rep)]);
        
        if pm.ratOrMouse == "rat"
            handle3 = figure(startFigNum + 2);
            set(handle3, 'Position', [950,575,475,400]);
            imagesc(abs(fftshift(ifftn(squeeze(shuffledData(3,:,:,frame,rep))))));
            title(['Coil 3: Frame ' num2str(frame) ' Rep ' num2str(rep)]);
            
            handle4 = figure(startFigNum + 3);
            set(handle4, 'Position', [1425,575,475,400]);
            imagesc(abs(fftshift(ifftn(squeeze(shuffledData(4,:,:,frame,rep))))));
            title(['Coil 4: Frame ' num2str(frame) ' Rep ' num2str(rep)]);
        end
        
        handle5 = figure(startFigNum + 4);
        set(handle5, 'Position', [0,100,475,400]);
        imagesc(abs(squeeze(shuffledData(1,:,:,frame,rep)))');
        title(['Coil 1: Frame ' num2str(frame) ' Rep ' num2str(rep)]);
        
        handle6 = figure(startFigNum + 5);
        set(handle6, 'Position', [475,100,475,400]);
        imagesc(abs(squeeze(shuffledData(2,:,:,frame,rep)))');
        title(['Coil 2: Frame ' num2str(frame) ' Rep ' num2str(rep)]);
        
        if pm.ratOrMouse == "rat"
            handle7 = figure(startFigNum + 6);
            set(handle7, 'Position', [950,100,475,400]);
            imagesc(abs(squeeze(shuffledData(3,:,:,frame,rep)))');
            title(['Coil 3: Frame ' num2str(frame) ' Rep ' num2str(rep)]);
            
            handle8 = figure(startFigNum + 7);
            set(handle8, 'Position', [1425,100,475,400]);
            imagesc(abs(squeeze(shuffledData(4,:,:,frame,rep)))');
            title(['Coil 4: Frame ' num2str(frame) ' Rep ' num2str(rep)]);
        end
        
        pause(timestep);
    end
end
end % End of unaveragedVisualizer


function averagedVisualizer(pm, startFigNum, averagedkSpaces, combinedMag, combinedGIF)
% Sets up 4x2 figures showing the kspaces of the 4 coils and the mag images
% of the same coils. Also shows the total combined mag. Averaged data.
timestep = 7.5/size(averagedkSpaces,4);
for frame = 1:pm.newFrameNum
    handle1 = figure(startFigNum);
    set(handle1, 'Position', [0,575,475,400]);
    imagesc(abs(fftshift(ifftn(squeeze(averagedkSpaces(1,:,:,frame))))));
    title(['Coil 1: Frame ' num2str(frame)]);
    
    handle2 = figure(startFigNum + 1);
    set(handle2, 'Position', [475,575,475,400]);
    imagesc(abs(fftshift(ifftn(squeeze(averagedkSpaces(2,:,:,frame))))));
    title(['Coil 2: Frame ' num2str(frame)]);

    if pm.ratOrMouse == "rat"
        handle3 = figure(startFigNum + 2);
        set(handle3, 'Position', [950,575,475,400]);
        imagesc(abs(fftshift(ifftn(squeeze(averagedkSpaces(3,:,:,frame))))));
        title(['Coil 3: Frame ' num2str(frame)]);
        
        handle4 = figure(startFigNum + 3);
        set(handle4, 'Position', [1425,575,475,400]);
        imagesc(abs(fftshift(ifftn(squeeze(averagedkSpaces(4,:,:,frame))))));
        title(['Coil 4: Frame ' num2str(frame)]);
    end
    
    handle5 = figure(startFigNum + 4);
    set(handle5, 'Position', [0,100,475,400]);
    imagesc(abs(squeeze(averagedkSpaces(1,:,:,frame)))');
    title(['Coil 1: Frame ' num2str(frame)]);
    
    handle6 = figure(startFigNum + 5);
    set(handle6, 'Position', [475,100,475,400]);
    imagesc(abs(squeeze(averagedkSpaces(2,:,:,frame)))');
    title(['Coil 2: Frame ' num2str(frame)]);
    
    if pm.ratOrMouse == "rat"
        handle7 = figure(startFigNum + 6);
        set(handle7, 'Position', [950,100,475,400]);
        imagesc(abs(squeeze(averagedkSpaces(3,:,:,frame)))');
        title(['Coil 3: Frame ' num2str(frame)]);
        
        handle8 = figure(startFigNum + 7);
        set(handle8, 'Position', [1425,100,475,400]);
        imagesc(abs(squeeze(averagedkSpaces(4,:,:,frame)))');
        title(['Coil 4: Frame ' num2str(frame)]);
    end
    
    handle9 = figure(startFigNum + 8);
    set(handle9, 'Position', [712.5,1200,800,700]);
    imagesc(squeeze(combinedMag(:,:,frame)));
    title(['Combined coils: Frame ' num2str(frame)]);
    colormap gray;
    
    GIFmaker(handle9, combinedGIF, timestep);
    
    pause(timestep);
end
end % End averagedvisualizer

function occupied = MREoccupancyEstimator(numLines, acqPerFrame)
% Solve the occupancy problem. This source does it complicated:
% (https://probabilityandstats.wordpress.com/2010/03/27/the-occupancy-problem/)
% GPT does it simpler (but seems right based on my testing):
occupied = numLines - (numLines * (1-(1/numLines))^acqPerFrame);
end % End occupancyEstimator

function [occupied, lineCounter] = MREoccupancyCalculator(shuffledData, showCSdist, figNum, CSGIFname)
% Calculate actual occupancy of the frames and mechPhases after shuffling.
dims = size(shuffledData); % [coils x y frames mechphases reps]
dims = [dims 1 1 1];
lineCounter = zeros(dims(3:5));
occupied    = zeros(dims(4:5));
nullLine    = zeros(dims(1:2));
distCells   = cell(dims(4),dims(5));
for rep = 1:dims(6)
    for mechPhase = 1:dims(5)
        for frame = 1:dims(4)
            for line = 1:dims(3)
                % Skip empty lines
                if ~isequal(squeeze(shuffledData(:,:,line,frame,mechPhase,rep)), nullLine)
                    lineCounter(line,frame,mechPhase)   = squeeze(lineCounter(line,frame,mechPhase)) + 1;
                    distCells{frame,mechPhase}          = [distCells{frame,mechPhase}; line-dims(3)/2];
                end
            end
        end
    end
end

% Sum up over frames and mechPhases
for mechPhase = 1:dims(5)
    for frame = 1:dims(4)
        tempFrame       = squeeze(lineCounter(:,frame,mechPhase));
        occupied(frame,mechPhase) = sum(tempFrame ~= 0);
    end
end
if showCSdist == true
    % Visualize filling of first mechPhase
    figure(figNum + 20);
    imagesc(squeeze(lineCounter(:,:,1)));
    title('Line distribution per frame in MechPhase 1');
    colorbar;
    
    % Check line distribution in histograms
    edges       = (-dims(3)/2):(dims(3)/2) + 0.5;
    totalDist   = vertcat(distCells{:});
    
    figure(figNum + 21);
    histogram(totalDist,edges);
    title('Line distribution - All frames and MechPhases')
    
    figure(figNum + 22);
    frameDist = [];
    for frame = 1:dims(4)
        frameDist = [frameDist; distCells{frame,1}];
        histogram(frameDist,edges);
        title(['Line distribution - Frame: ' num2str(frame) ' - MechPhase: ' num2str(1)]);
        GIFmaker(figNum + 22, CSGIFname, 5/dims(4)); % 5 seconds divided by movieFrames
        pause(5/dims(4));
    end
end

end % End occupancyCalculator

