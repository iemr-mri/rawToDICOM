function [interpolated, timeline, pm] = timeCorrecter(PCcurves, timestamps, pm, visualSwitch)
% MRE acquisitions have varying TR due to waiting for the shaker. At 400 Hz
% this wait can be up to 2 x shaker period = 5 ms.
% Takes in the PCA-curves, expands the temporal resolution, and uses the
% timestamps to place every acquired line in its correct position on the
% new timeline. 

disp('Section 3: timeCorrecter');
tic;
% Adjustable parameters
chunkSelection  = [1,3];
visualPC        = 3;

% Setup timepoints etc
timestamps      = timestamps * 1000;                                % Use milliseconds
midstamps       = timestamps(pm.midlineRate:pm.midlineRate:end);    % Select the midlines. Start from e.g. 3rd, and take every e.g. 3rd
pm.newTempRes   = pm.TR/100;                                        % Resolution set to 1/100 of TR. Use milliseconds as unit. [FIKS] Could maybe use physiology-defined resolution? Maybe 1/1000th of 8 Hz = 125 us
totalTime       = (timestamps(end) - timestamps(1));                % Time from end of first read-out to end of last read-out. Millisecond units
timeline        = timestamps(1):pm.newTempRes:(totalTime+1);        % Partition totalTime based on chosen resolution. Add 1 ms buffer at end to ensure we cover last stamp.

% Find new temporal positions of midlines given new temporal resolution
roundedMids         = round(midstamps/pm.newTempRes);
incorrectTimeline   = (1:length(PCcurves))*(pm.TR*(pm.midlineRate)) - pm.TR; % Setup incorrectTime to visualize temporal shift.

% Setup NaN-filled, empty timeline, compatible with fillmissing function.
padded          = nan(1,length(timeline));

% Setup interpolated matrix
PCnum           = size(PCcurves,1);
interpolated    = zeros(PCnum,length(timeline));

% Insert original values at each roundedMid timepoint, then interpolate, for each PC
for PC = 1:PCnum
    curPC       = squeeze(PCcurves(PC,:));
    curPadded   = padded;
    for midline = 1:pm.midlineNum
        curIndex            = roundedMids(midline);
        curPadded(curIndex) = curPC(midline);
    end
    interpolated(PC,:) = fillmissing(curPadded,'linear');
end

% Scale original signal to new timeline
origTempRes     = pm.TR*pm.midlineRate;

% Visualize interpolation vs original data
if visualSwitch == true
    chunkTime       = 5000; % 5000 milliseconds per chunk (~= 5 breaths)
    longTime        = floor(chunkTime/pm.newTempRes);
    shortTime       = floor(chunkTime/origTempRes);
    
    interpRange     = (longTime*chunkSelection(1) + 1):longTime*chunkSelection(2);
    interpTime      = timeline(interpRange);
    origRange       = (shortTime*chunkSelection(1) + 1):shortTime*(chunkSelection(2));
    origTime        = incorrectTimeline(origRange);
    
    axisLims        = [interpTime(1)-1000,interpTime(end)+1000];
    
    figure(501);
    plot(interpTime, squeeze(interpolated(visualPC,interpRange)), 'b-', 'DisplayName', 'Interpolated Signal');
    hold on;
    plot(origTime, squeeze(PCcurves(visualPC,origRange)), 'r-', 'DisplayName', 'Original Signal');
    xlabel('Time');
    ylabel('Signal Value');
    xlim(axisLims);
    title('Temporally corrected data vs uncorrected');
    legend;
    hold off;
end

% FFT
interpFFT   = fftshift(fft(squeeze(interpolated(visualPC,:))));
origFFT     = fftshift(fft(squeeze(PCcurves(visualPC,:))));
 
interpNyqFreq   = 1/(2*pm.newTempRes) * 1000; % Division by millisecond => kiloHz
origNyqFreq     = 1/(2*origTempRes) * 1000; % Division by millisecond => kiloHz

interpFreqAx    = linspace(-interpNyqFreq, interpNyqFreq, length(interpolated)) - interpNyqFreq/length(interpolated); % Have to shift by 1 freqstep for some reason
origFreqAx      = linspace(-origNyqFreq, origNyqFreq, length(PCcurves)) - origNyqFreq/length(PCcurves); % Have to shift by 1 freqstep for some reason

% Visualize ffts
if visualSwitch == true
    figure(502);
    plot(interpFreqAx, abs(interpFFT));
    xlim([-origNyqFreq, origNyqFreq]);
    title("Frequency content of PCA after timestamp correction");
    
    figure(503);
    plot(origFreqAx, abs(origFFT));
    xlim([-origNyqFreq, origNyqFreq]);
    title("Frequency content of raw PCA");
end

disp(['Section 3: timeCorrecter - Finished in ' num2str(toc) ' seconds.' newline]);
end % End of function