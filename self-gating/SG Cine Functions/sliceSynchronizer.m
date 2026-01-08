function [magnitude, rotDiffs, minDiff] = sliceSynchronizer(magnitude, expName, scanNumbers, folderName, showProcess)
% Self-gated slices from the same session are generally not synchronized as
% they are not R-peak triggered, but instead synchronized to the peaks of
% the navigator-PCA data.
% Takes in the slices, retrieves the center slice, and looks for diastole
% by having the user set a ROI around the heart and evaluating the
% development of air-pixels vs myocardium-pixels. The 7-timestep window with
% the least air-pixels is assumed to be end diastole. The neighbours are
% then synchronized by lowest RMS difference to the center slice, then 
% their neighbours and so on.

% Input:
% -magnitude        Either the pre-loaded slices [x, y, cardiacPhases, slice]
%                   or boolean false, indicating that we need to load them
% -expName          The name of the session so that we can load the data
% -scanNumbers      The scans/slices in question, so that we can load them

% Output:
% -magnitude        The synchronized version of the input data
% -origMagnitude    The unsynchronized version of the input data
% -rotDiffs         The RMS differences for each temporal rotation for all
%                   slices
% -minDiff          The minima of rotDiffs. One for each slice. Determines
%                   how much each slice should be rotated

tic;
disp('Section 11: Slice synchronizer')

% If the images weren't processed during this run, retrieve them from saved files
numSlices = length(scanNumbers);
if magnitude == false
    for slice = 1:numSlices
        temp = load([folderName expName '\' char(scanNumbers(slice)) '\imageData.mat']);
        tempMag = abs(squeeze(temp.imageData));
        if slice == 1
            magnitude = zeros([size(tempMag,1), size(tempMag,2), numSlices, size(tempMag,3)]);
        end
        magnitude(:,:,slice,:) = tempMag;
    end
end

% Set origMagnitude for comparison to sorted version
origMagnitude = magnitude;

% Retrieve dimensions
[xpixels, ypixels, numSlices, cardiacPhases] = size(magnitude);

% Brightness levels for visualizations
upperBright = prctile(magnitude, 99, 'all');
lowerBright = prctile(magnitude, 25, 'all');

% Crop area around heart for improved diastole detection
midIndex = ceil(numSlices/2);
midSlice = squeeze(magnitude(:,:,midIndex,1)); % Assume that the middle slice is best for diastole detection. Arbitrarily take first timepoint.
maskFigHandle = figure();
imagesc(midSlice);
colormap gray;
clim([lowerBright, upperBright]);
xlim([xpixels*1/4, xpixels*3/4]);
ylim([ypixels*1/4, ypixels*3/4]);
disp("Select a ROI that should include both the heart and the air around it in all cardiac phases and all slices for diastole detection and slice synchronization.")
heartPoly = drawpolygon();
heartMask = createMask(heartPoly);
heartMask = repmat(heartMask,[1,1,cardiacPhases]);
% heartMask = ones(size(heartMask));

% Select piece of myocardium to detect air
figure(maskFigHandle)
imagesc(midSlice);
colormap gray;
clim([lowerBright, upperBright]);
xlim([xpixels*1/4, xpixels*3/4]);
ylim([ypixels*1/4, ypixels*3/4]);
disp("Select piece of myocardium to help separate myocardium and air.")
myoPoly     = drawpolygon();
myoMask     = createMask(myoPoly);
myoVal      = mean(midSlice(myoMask),"all");
airThresh   = 0.6*myoVal; % Pixels with less than a certain intensity of the myocardium are considered air

% Detect diastole for central slice by looking for the least amount of air near the heart
% K: It might be easier to find peak systole in noisy acquisitions. Might
% need to change approach later
midMasked       = squeeze(magnitude(:,:,midIndex,:)).*heartMask;
midThreshed     = midMasked < airThresh;
airTrend        = squeeze(sum(sum(midThreshed,1),2));
paddedAirTrend  = repmat(airTrend,[1,3]); % Pad to enable averaging at borders too
smoothAirTrend  = zeros(1,cardiacPhases);
numNeighbours   = floor(cardiacPhases/10); % How many neighbouring timepoints to include in the average. /10 => 20% of the cardiac cycle is considered
for cardiacPhase = 1:cardiacPhases
    curRange = cardiacPhases + cardiacPhase - numNeighbours:cardiacPhases + cardiacPhase + numNeighbours;
    smoothAirTrend(cardiacPhase) = mean(paddedAirTrend(curRange));
end
[~,diastOffset] = min(smoothAirTrend); % Minimal amount of air = diastole
magnitude(:,:,midIndex,:) = circshift(squeeze(magnitude(:,:,midIndex,:)),1-diastOffset,3);

% Visualize the air trend. Should have 1 period of sinusoid behavior.
if showProcess == true
    figStartNum = 3000;
    figure(figStartNum)
    plot(airTrend, 'DisplayName','Raw air')
    hold on
    plot(smoothAirTrend, 'DisplayName','Smoothed air')
    plot(diastOffset,min(smoothAirTrend), 'o', 'MarkerSize', 13, 'MarkerEdgeColor', 'black', 'LineWidth', 2, 'DisplayName','Minimal air')
    xlabel("Cardiac phase")
    ylabel("Air pixels in the ROI")
    title("Trend of air pixels near heart area")
    legend('Location','best')
    hold off
end

% For single slice, we just exit
if numSlices == 1
    rotDiffs    = 0;
    minDiff     = 0;
    return;
end

% Synchronize each slice with the previous (assuming centerslice in diastole). 
% Rotate cardiac phases until they align (like a Rubik's cube).
% Start from the center and count up, then start from center and count down
% so that all slices are synchronized to the center slice.
sliceSequence   = [midIndex+1:numSlices, midIndex-1:-1:1];
rotDiffs        = zeros(length(sliceSequence),cardiacPhases);
minDiff         = zeros(1,length(sliceSequence));
upOrDown        = -1; % We compare each slice to either the slice above or below
for slice = 1:numSlices-1
    curIndex    = sliceSequence(slice);
    if curIndex == midIndex - 1 % When we start counting down, we start looking at the neighbour above
        upOrDown = +1;
    end
    neighbour   = curIndex + upOrDown;
    curSlice    = squeeze(magnitude(:,:,curIndex,:)).*heartMask;
    neighSlice  = squeeze(magnitude(:,:,neighbour,:)).*heartMask;

    % Visualize
    if showProcess == true
        curFigHand      = figure(figStartNum + 1);
        imagesc(squeeze(magnitude(:,:,curIndex,1)));
        xlim([xpixels*1/4, xpixels*3/4]);
        ylim([ypixels*1/4, ypixels*3/4]);
        clim([lowerBright, upperBright]);
        colormap gray
        title("Current slice")

        neighFigHand    = figure(figStartNum + 2);
        imagesc(squeeze(magnitude(:,:,neighbour,1)));
        xlim([xpixels*1/4, xpixels*3/4]);
        ylim([ypixels*1/4, ypixels*3/4]);
        clim([lowerBright, upperBright]);
        colormap gray
        title("Neighbour slice")

        stackDiffHand   = figure(figStartNum + 3);
        colormap gray

        diffPlotHand    = figure(figStartNum + 4);
    end

    % Check all cardiac offsets to find minimal difference
    for cardiacPhase = 1:cardiacPhases
        curOffset   = cardiacPhase - 1;
        curRotSlice = circshift(curSlice,curOffset,3); % Shift in cardiac dimension
        diffSlice   = neighSlice - curRotSlice;
        rotDiffs(slice,cardiacPhase) = rms(diffSlice(:));

        % Visualize
        if showProcess == true
            figure(curFigHand)
            imagesc(squeeze(curRotSlice(:,:,1)))
            xlim([xpixels*1/4, xpixels*3/4]);
            ylim([ypixels*1/4, ypixels*3/4]);
            title(['Current slice - Rotated ' num2str(curOffset) ' steps'])
            
            figure(stackDiffHand)
            avgDiffSlice = mean(sqrt(diffSlice.^2),3);
            % imagesc(squeeze(diffSlice(:,:,1)));
            colorbar
            imagesc(avgDiffSlice)
            title(['Visual difference - Rotated ' num2str(curOffset) ' steps']);

            figure(diffPlotHand)
            plot(squeeze(rotDiffs(slice,:)))
            title(['RMS difference - Rotated ' num2str(cardiacPhase) ' steps']);
            xlim([1,cardiacPhases])

            pause(0.1);
        end
    end
    % Find the point where the difference was minimal
    [~,minDiff(slice)] = min(squeeze(rotDiffs(slice,:)));

    % Rotate current slice until it matches with the neighbour
    magnitude(:,:,curIndex,:) = circshift(squeeze(magnitude(:,:,curIndex,:)), minDiff(slice),3);
end

disp(['Section 11: Slice synchronizer completed in ' num2str(toc) ' seconds.']);
close all
end
