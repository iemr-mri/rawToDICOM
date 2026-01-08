function [roughPeakList, finePeakList, breathData, breathStarts, pm] = PCApeakFinder(cardiacData, breathData, timeStamps, pm, visualSwitch)
% Looks for areas with a certain rise and drop, finds a rough local peak by
% taking the max value. Uses that max value as a guess for curve fitting
% using a skewed gaussian fit. Curve fitted peak values are stored as fine
% peaks.

disp('Section 5: PCApeakFinder');
tic;

% Calculate rise and falltimes
cardiacPeriod   = 1/pm.cardiacFreq;
breathingPeriod = 1/pm.breathingFreq;
pointsPerHB     = cardiacPeriod*1000/pm.newTempRes;   % Compensate for milliseconds
pointsPerBreath = breathingPeriod*1000/pm.newTempRes; % Compensate for milliseconds
riseTime        = floor(0.45*pointsPerHB);
fallTime        = floor(0.45*pointsPerHB);


minHeight = 5;

% Find the peaks of the breath movement and determine the breath starts.
% A certain range of the area between peaks is discarded in a later function.

% Find the heights, locations, and widths of all breath peaks and valleys
[breathPeakHeights,breathPeaks,breathPeakWidths,~]          = findpeaks(breathData,'MinPeakProminence', std(breathData));
[breathValleyHeights,breathValleys,breathValleyWidths,~]    = findpeaks(-breathData,'MinPeakProminence', std(breathData));
meanBreathDur = mean(diff(breathPeaks));

% We investigate the differential to find where the curve changes fastest.
% We assume the most rapid change marks the early phase of breaths.
breathDiffs             = diff(breathData);
breathDiffs(end+1)      = 0; % Pad the end to ensure same length

% Find breathStarts using either the differential of the curve or the width
% of peaks.
if pm.diffOrWidth == "diff"
    % Find points of max change. Set semi-arbitrary MinPeakProminence to filter 
    % out small variations.
    [diffPeakHeights,diffPeaks,~,~]       = findpeaks(breathDiffs, 'MinPeakProminence', std(breathDiffs));
    [diffValleyHeights,diffValleys,~,~]   = findpeaks(-breathDiffs, 'MinPeakProminence', std(breathDiffs));
    breathStarts = diffBreathHandler(diffPeakHeights, diffValleyHeights, diffPeaks, diffValleys, breathPeaks, meanBreathDur, breathValleys, breathPeakHeights, breathValleyHeights, pm.breathFlipBool);
elseif pm.diffOrWidth == "width"
    breathStarts = widthBreathHandler(breathPeaks, breathPeakHeights, breathPeakWidths, breathValleys, breathValleyHeights, breathValleyWidths, pm.breathFlipBool);
else
    disp("Breath peak handling was not set to either 'diff' or 'width'. Exiting.")
    return;
end



finePoints  = 10000; % How many points the skewedPeakFinder should segment the truncChunk into to find its peak
figNumStart = 40; 

% Flip cardiac curve if desired
if pm.cardiacFlipBool == true
    cardiacData = -cardiacData;
end

% Loop through all timepoints and check for peaks
roughPeakList   = [];
time = 1;

while time <= length(cardiacData) - (riseTime + fallTime)

    % Check if the current point is a start of an increase
    if (cardiacData(time + riseTime) - cardiacData(time) > minHeight)
        % Check if the subsequent decrease meets the criteria
        if (cardiacData(time + riseTime) - cardiacData(time + riseTime + fallTime) > minHeight)
            % Make sure we don't exceed the end of the y vector
            if (time + riseTime + fallTime) > length(cardiacData)
                dataChunk = cardiacData(time:end);
            else
                dataChunk = cardiacData(time:time + riseTime + fallTime);
            end

            % Find the highest value within the detected peak range
            [peakValue, peakIndex]      = max(dataChunk);
            globalPeakIndex             = time - 1 + peakIndex;
            
            % Add the detected peak to the list
            roughPeakList   = [roughPeakList; globalPeakIndex, peakValue];
            
            % Visualize
            if visualSwitch == true
                figure(figNumStart);
                plot(dataChunk);
                hold on;
                plot(peakIndex, peakValue, 'ro', 'MarkerSize', 5, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
                hold off;
            end

            time = time + riseTime + fallTime;
        else
            time = time +1;
        end
    else
        time = time + 1;
    end
end

% Improve the roughPeakList to sub-sample rate temporal resolution using
% curve fitting to find a better peak
finePeakList = zeros(size(roughPeakList));
if pm.performCurveFit == true
    for peak = 1:length(roughPeakList)
        tempPeak                = roughPeakList(peak,:);
        truncStart              = tempPeak(1) - riseTime;
        truncEnd                = tempPeak(1) + fallTime;
        % Handle border cases
        if truncStart <= 0 || truncEnd >= length(cardiacData)
            continue;
        end
        peakChunk               = cardiacData(truncStart:truncEnd);
        maxPeak                 = [tempPeak(1) - truncStart, tempPeak(2)];
    
        % PCAfinePeakFinder takes the current chunk and fits it to a skewed
        % gaussian curve before finding the peak of that curve
        fitPeak                 = PCAfinePeakFinder(peakChunk, finePoints, maxPeak, visualSwitch);
        finePeakList(peak,1)    = truncStart + fitPeak(1);
        finePeakList(peak,2)    = cardiacData(round(finePeakList(peak,1))); % fitPeak(2); 
        
        % Display the differences between rough and fine peaks
        if visualSwitch == true
        disp(['Original peak position: [' num2str(roughPeakList(peak,1)) ',' num2str(roughPeakList(peak,2)) '] - Adjusted: [' ...
            num2str(finePeakList(peak,1)) ',' num2str(finePeakList(peak,2)) '] - Difference: [' ...
            num2str(roughPeakList(peak,1) - finePeakList(peak,1)) ',' num2str(roughPeakList(peak,2) - finePeakList(peak,2)) ']']);
        end
        % pause(1);
    end
elseif pm.performCurveFit == false
    finePeakList = roughPeakList;
end

% Cardiac data was only temporarily flipped, so the points need to be
% permanently flipped
if pm.cardiacFlipBool == true
    roughPeakList(:,2) = -roughPeakList(:,2);
    finePeakList(:,2) = -finePeakList(:,2);
end

% Display number of detected peaks
detectedCardiacPeaks    = length(finePeakList);
estimatedCardiacPeaks   = floor(pm.cardiacFreq*(timeStamps(end)-timeStamps(1))); % Heartbeats per second mulitplied by duration of recording.
detectedBreathPeaks     = length(breathStarts);
estimatedBreathPeaks    = floor(pm.breathingFreq*(timeStamps(end)-timeStamps(1)));
disp(['Accepted ' num2str(detectedCardiacPeaks) ' cardiac peaks out of an estimated total of ' num2str(estimatedCardiacPeaks)]);
disp(['Detected ' num2str(detectedBreathPeaks) ' breath peaks out of an estimated total of ' num2str(estimatedBreathPeaks)])

disp(['Section 5: PCApeakFinder - Finished in ' num2str(toc) ' seconds.' newline]);
end % End of PCApeakFinder.m

%% Support functions

function fitPeak = PCAfinePeakFinder(dataChunk, finePoints, maxPeak, visualSwitch)

fitPeak = [0,0];

% Compensate for counting issue during segment truncation
maxPeak(1) = maxPeak(1) + 1; 

% Curvefitting depends on data being barely above zero. Shift data
chunkShifter = min(dataChunk);
dataChunk = dataChunk - chunkShifter;

% Define x values for the truncated data
x = 1:length(dataChunk);

% Define a custom Skewed Gaussian model for curve fitting
skewedGaussianModel = fittype('a*exp(-((x-b)/c)^2).*(1+erf(d*(x-b)))', ...   % Skewed Gaussian function form
    'independent', 'x', 'coefficients', {'a', 'b', 'c', 'd'});

% Perform the fitting with an initial guess for the coefficients
guessPos = find(dataChunk == max(dataChunk));

if length(guessPos) > 1  % Avoid multiple guessPos
    guessPos = guessPos(1);
end

initialGuess        = [max(dataChunk), guessPos, std(x), 0.1]; % This STD might be problematic...?
fittedCurve         = fit(x(:), dataChunk(:), skewedGaussianModel, 'Start', initialGuess);

% Find peak of fittedCurve
fineRange           = linspace(1, length(dataChunk), finePoints); % More points for a finer resolution
fineFit             = fittedCurve(fineRange);
[~, finePeakIndex]  = max(fineFit);

fitPeak(1)          = fineRange(finePeakIndex);
fitPeak(2)          = fittedCurve(fitPeak(1)) + chunkShifter;

if visualSwitch == true
    % Visualize the raw and fitted data
    figure(100);
    plot(x, dataChunk, '-', 'Color', 'blue');
    hold on;
    plot(fittedCurve, '-'); % Plot the fitting result with a red line
    plot(maxPeak(1), dataChunk(maxPeak(1)), 'rx', 'MarkerSize', 15, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'none'); 
    plot(fitPeak(1), fittedCurve(fitPeak(1)), 'rx', 'MarkerSize', 15, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'none');
    grid on;
    hold off;
    % pause(0.5);
end

% Compensate for counting issue during truncation
fitPeak(1) = fitPeak(1) - 1;
end % End of PCAfinePeakFinder

function breathStarts = diffBreathHandler(diffPeakHeights, diffValleyHeights, diffPeaks, diffValleys, breathPeaks, meanBreathDur, breathValleys, breathPeakHeights, breathValleyHeights)
% Finds the breathStarts using the point of highest differential in 
% relation to the peaks. It is assumed that the point with the most change
% is the aspiration (breath in), which means that the peak before it is the
% breath start.

% If the diff peaks have stronger changes than valleys, the peaks are the
% points of max change. Else it is the valleys.
if mean(diffPeakHeights) > mean(diffValleyHeights)
    maxDiffs = diffPeaks;
else
    maxDiffs = diffValleys;
end

% Go through all peaks and find their distance to the previous maxDiff
peaksToMaxDiffs = [];
for peak = 1:length(breathPeaks)
    curPeak     = breathPeaks(peak);
    
    % Take the highest maxDiff that is before the current peak
    prevMaxDiff = max(maxDiffs(maxDiffs < curPeak));
    
    % Filter out empty peaks and peaks that are spaced too far (assumed
    % missed peaks)
    % Add distance from peak to previous maxDiff to the list.
    if ~isempty(prevMaxDiff) && (curPeak-prevMaxDiff) < meanBreathDur*1.5 
        peaksToMaxDiffs = [peaksToMaxDiffs; curPeak-prevMaxDiff];
    end
end

% Go through all valleys and find their distance to the previous maxDiff
valleysToMaxDiffs = [];
for peak = 1:length(breathValleys)
    curPeak     = breathValleys(peak);
    
    % Take the highest maxDiff that is before the current valley
    prevMaxDiff = max(maxDiffs(maxDiffs < curPeak));
    
    % Filter out empty valleys and valleys that are spaced too far (assumed
    % missed valleys)
    % Add distance from valley to previous maxDiff to the list.
    if ~isempty(prevMaxDiff) && (curPeak-prevMaxDiff) < meanBreathDur*1.5 
        valleysToMaxDiffs = [valleysToMaxDiffs; curPeak-prevMaxDiff];
    end
end

% Whichever is furthest away from the previous maxDiff is closest to the
% next. Whichever is closest to the next is also the start of the breath.
if mean(peaksToMaxDiffs) > mean(valleysToMaxDiffs)
    breathStarts(1,:) = breathPeaks;
    breathStarts(2,:) = breathPeakHeights;
else
    % We also flip the breath data to have the correct side up
    breathStarts(1,:) = breathValleys;
    breathStarts(2,:) = -breathValleyHeights;
end
end

function breathStarts = widthBreathHandler(breathPeaks, breathPeakHeights, breathPeakWidths, breathValleys, breathValleyHeights, breathValleyWidths, flipBool)
% Check the width of peaks. Assumed that the sharper peak represents
% breathing.
meanPeakWidth   = mean(breathPeakWidths);
meanValleyWidth = mean(breathValleyWidths);

% The thinner peaks win
if meanPeakWidth < meanValleyWidth
    if flipBool == false
        breathStarts(1,:) = breathPeaks; 
        breathStarts(2,:) = breathPeakHeights; 
    else
        % Flip the breath data to have the correct side up (peaks)
        breathStarts(1,:) = breathValleys; 
        breathStarts(2,:) = -breathValleyHeights; 
    end
else
    if flipBool == false
        breathStarts(1,:) = breathValleys; 
        breathStarts(2,:) = -breathValleyHeights; 
    else
        % Flip the breath data to have the correct side up (valleys)
        breathStarts(1,:) = breathPeaks; 
        breathStarts(2,:) = breathPeakHeights; 
    end
end
end