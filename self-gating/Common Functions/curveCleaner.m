function [cardiacPC, normalizedCardiac, breathPC, normalizedBreath, pm] = curveCleaner(prinComps, PCweights, pm, visualSwitch)
% Among the e.g. 10 first PCs, find the principal component that best represents cardiac and breathing phases.
% Use FT and denoising to cleanup the PC curve so it is ready for peak detection.

disp('Section 4: curveCleaner');
tic;
% Parameter setup
numPrinComps        = length(pm.PCselection);
timePoints          = length(prinComps);
samplingFreq        = 1/pm.newTempRes*1000; % Convert millisecond to seconds
nyqFreq             = samplingFreq/2;
freqAxis            = linspace(-nyqFreq,nyqFreq,timePoints); % Should range from e.g 77 Hz (samplingFreq = 1/6.5 ms) to 0.0160 Hz (samplingFreq/9600).
freqStep            = freqAxis(2)-freqAxis(1); % Difference between two points in the freqAxis (in Hz)
freqAxis            = freqAxis - freqStep/2; % For reasons unknown to this elastographer we have to shift the frequency axis by half a step

% Set length of timespans to visualize breathing and beating
breathVisualTime    = 5; % Seconds to visualize breathing
cardiacVisualTime   = 3; % Seconds to visualize beating
breathVisRange      = 1:(breathVisualTime*samplingFreq);
cardiacVisRange     = 1:(cardiacVisualTime*samplingFreq);

if length(prinComps) < length(breathVisRange)
    breathVisRange  = 1:length(prinComps);
end

if length(prinComps) < length(cardiacVisRange)
    cardiacVisRange = 1:length(prinComps);
end

% Set frequency range to visualize in ffts of breathing and beating
origSamplingFreq    = 1/(pm.midlineRate * pm.TR * 1e-3);
origNyqFreq         = origSamplingFreq/2;
stepsToOrigNyqFreq  = ceil(origNyqFreq/freqStep);
newMidpoint         = floor(timePoints/2);
origRangeStart      = newMidpoint - stepsToOrigNyqFreq;
origRangeEnd        = newMidpoint + stepsToOrigNyqFreq;

% Adapt ranges to rat or mouse
if pm.ratOrMouse == "rat"
    cardiacRange    = [4, 8];        % Assumed cardiac rhytm within [] Hz range
    breathingRange  = [0.5, 1.5];    % Assumed breathing rhythm within [] Hz range
elseif pm.ratOrMouse == "mouse"
    cardiacRange    = [4, 10];      
    breathingRange  = [0.5, 2];
end

% Analyze first 10 PCs to find the best cardiac and respiratory components.
% Near-zero frequencies are temporarily nulled to prevent signal drift from
% dominating actual cardiac and breathing patterns.
% Each PC is assigned a value based on what fraction of its spectral energy
% is kept in the main frequency, and based on its weighting from the
% original PCA.
FTmat           = zeros(numPrinComps, length(prinComps));
FTmaxIndices    = zeros(1,numPrinComps);
FTmaxMag        = zeros(1,numPrinComps);
FTmaxHz         = zeros(1,numPrinComps);
FTenergy        = zeros(1,numPrinComps);
nearZeroFreq    = 0.3; % Frequencies below 0.3 Hz are considered signal drift
nearZeroCutOff  = ceil(nearZeroFreq/freqStep);

PSDmat          = zeros(numPrinComps, floor(length(prinComps)/2)+1);
PSDmaxIndices   = zeros(1,numPrinComps);
PSDmaxMag       = zeros(1,numPrinComps);
PSDmaxHz        = zeros(1,numPrinComps);
PSDenergy       = zeros(1,numPrinComps);

for PC = 1:numPrinComps
    curPC                               = prinComps(PC,:);
    FTmat(PC,:)                         = fft(curPC);
    FTmat(PC,1:nearZeroCutOff)          = 0; % FT spectrum is twosided and must be nulled on both sides
    FTmat(PC,(end-nearZeroCutOff):end)  = 0; % FT spectrum is twosided and must be nulled on both sides
    [FTmaxMag(PC), FTmaxIndices(PC)]    = max(abs(squeeze(FTmat(PC, 1:round(timePoints/2)))));
    FTmaxHz(PC)                         = (FTmaxIndices(PC)-1) * freqStep;
    FTenergy(PC)                        = FTmaxMag(PC)^2/sum(abs(squeeze(FTmat(PC, 1:round(timePoints/2)))))^2;

    PSDmat(PC,:) = periodogram(curPC, rectwin(length(curPC)), length(curPC), 1/freqStep, 'psd');
    PSDmat(PC,1:nearZeroCutOff)         = 0; % PSD spectrum is one-sided for real-valued input and only needs nulling on one side
    [PSDmaxMag(PC), PSDmaxIndices(PC)]  = max(abs(squeeze(PSDmat(PC, :))));
    PSDmaxHz(PC)                        = (PSDmaxIndices(PC)-1) * freqStep;
    PSDenergy(PC)                       = PSDmaxMag(PC)/sum(PSDmat(PC,:));
    
    % Visualize
    if pm.showPCs == true
        figure(PC);
        plot(freqAxis + freqStep/2,fftshift(abs(FTmat(PC,:))));
        title(['FT Magnitude: PC ' num2str(PC) ' - Max Hz: ' num2str(FTmaxHz(PC)) ' - Energy: ' num2str(FTenergy(PC)/FTenergy(1))]);
        xlim([0,origNyqFreq]);

        figure(PC + numPrinComps);
        plot(freqStep*(0:length(PSDmat)-1), PSDmat(PC,:));
        title(['PSD: PC ' num2str(PC) ' - Max Hz: ' num2str(PSDmaxHz(PC)) ' - Energy: ' num2str(PSDenergy(PC)/PSDenergy(1))]);
        xlim([0,origNyqFreq]);

        figure(PC + 2 * numPrinComps);
        threeSecs = ceil(3000/pm.newTempRes); % 3000 milliseconds
        plot(prinComps(PC,1:threeSecs));
        title(['Original signal: PC ' num2str(PC)]);
    end
end

% Select the PCs with peaks within cardiac and breathing ranges
inCardiacRange  = (PSDmaxHz > cardiacRange(1)) .* (PSDmaxHz < cardiacRange(2));
inBreathRange   = (PSDmaxHz > breathingRange(1)) .* (PSDmaxHz < breathingRange(2));

cardiacEnergy   = PSDenergy .* inCardiacRange;% .* PCweights(1:numPrinComps).^2;
breathingEnergy = PSDenergy .* inBreathRange;% .* PCweights(1:numPrinComps).^2;

% Error if we didn't find any PCs for one or both of the categories
if sum(breathingEnergy) == 0
    error("curveCleaner could not find a PC with its main frequency within breathing range");
elseif sum(cardiacEnergy) == 0
    error("curveCleaner could not find a PC with its main frequency within cardiac range");
end

[~,cardiacPCnum]    = max(cardiacEnergy);
[~,breathPCnum]     = max(breathingEnergy);
cardiacPC           = prinComps(cardiacPCnum,:);
breathPC            = prinComps(breathPCnum,:);

pm.cardiacFreq      = FTmaxHz(cardiacPCnum);
pm.cardiacPeriod    = 1/pm.cardiacFreq * 1000; % Period in ms
pm.breathingFreq    = FTmaxHz(breathPCnum);
pm.breathingPeriod  = 1/pm.breathingFreq * 1000; % Period in ms

disp(['PC' num2str(cardiacPCnum) ' was strongest cardiac PC at ' sprintf('%.2f', pm.cardiacFreq) ' Hz (' num2str(round(pm.cardiacPeriod)) ' ms).'])
disp(['Midlines per cardiac period: ' num2str(pm.cardiacPeriod/(pm.midlineRate*pm.TR), '%.2f')]);
disp(['PC' num2str(breathPCnum) ' was strongest breath PC at ' sprintf('%.2f', pm.breathingFreq) ' Hz (' num2str(round(pm.breathingPeriod)) ' ms).'])
disp(['Midlines per breathing period: ' num2str(pm.breathingPeriod/(pm.midlineRate*pm.TR), '%.2f')]);

% Butterworth filter setup
cardiacLowThresh    = 0.5 * pm.cardiacFreq;
cardiacHighThresh   = 2 * pm.cardiacFreq;
cardiacLowOrder     = 2;
cardiacHighOrder    = 2;

breathLowThresh     = 4 * pm.breathingFreq;
breathLowOrder      = 3;
breathHardThresh    = 0.2 * pm.breathingFreq;

% Butterworth filtering
cardiacFilterLow  = 1./(1+(freqAxis/cardiacLowThresh).^(2*cardiacLowOrder));
cardiacFilterHigh = 1-1./(1+(freqAxis/cardiacHighThresh).^(2*cardiacHighOrder));
cardiacFilter     = cardiacFilterHigh.*cardiacFilterLow;

cardiacFT         = fftshift(fft(cardiacPC,[],2));
filteredCardiacFT = cardiacFT.*cardiacFilter;
filteredCardiac   = ifft(ifftshift(filteredCardiacFT),[],2,'Symmetric');

breathFilterLow  = 1./(1+(freqAxis/breathLowThresh).^(2*breathLowOrder));
breathFilterHard = abs(freqAxis)>breathHardThresh;
breathFilter     = breathFilterLow.*breathFilterHard;

breathFT         = fftshift(fft(breathPC,[],2));
filteredBreathFT = breathFT.*breathFilter;
filteredBreath   = ifft(ifftshift(filteredBreathFT),[],2,'Symmetric');

% Normalize y-scales
normalizedCardiac   = 50*filteredCardiac/max(filteredCardiac); % Ranges [-50,50]
cardiacGraphScaler  = max(normalizedCardiac)/max(cardiacPC);  % Make graphs equalish height. Quite violent, but should work for stable graphs
cardiacPC           = cardiacPC*cardiacGraphScaler;

normalizedBreath   = 50*filteredBreath/max(filteredBreath); % Ranges [-50,50]
breathGraphScaler  = max(normalizedBreath)/max(breathPC);  % Make graphs equalish height. Quite violent, but should work for stable graphs
breathPC           = breathPC*breathGraphScaler;

if visualSwitch == true    
    % Visualization
    figure(21);
    plot(cardiacPC(cardiacVisRange));
    title('Cardiac PC - Zoomed view');
    
    figure(22);
    plot(freqAxis(origRangeStart:origRangeEnd) + freqStep/2,abs(cardiacFT(origRangeStart:origRangeEnd)));
    title('Cardiac FT curve - Full');
    
    figure(23);
    plot(freqAxis(origRangeStart:origRangeEnd) + freqStep/2,abs(filteredCardiacFT(origRangeStart:origRangeEnd)));
    title('Cardiac FT curve - Filtered view');
    
    figure(24);
    plot(normalizedCardiac(cardiacVisRange));
    title('Filtered Cardiac PC - Zoomed view');
    
    figure(25);
    plot(freqAxis(origRangeStart:origRangeEnd),cardiacFilter(origRangeStart:origRangeEnd));
    title('Cardiac ButterWorth filter');

    figure(26);
    plot(breathPC(breathVisRange));
    title('Breath PC - Zoomed view');
    
    figure(27);
    plot(freqAxis(origRangeStart:origRangeEnd) + freqStep/2,abs(breathFT(origRangeStart:origRangeEnd)));
    title('Breath FT curve - Full');
    
    figure(28);
    plot(freqAxis(origRangeStart:origRangeEnd) + freqStep/2,abs(filteredBreathFT(origRangeStart:origRangeEnd)));
    title('Breath FT curve - Filtered view');
    
    figure(29);
    plot(normalizedBreath(breathVisRange));
    title('Filtered Breath PC - Zoomed view');
    
    figure(30);
    plot(freqAxis(origRangeStart:origRangeEnd),breathFilter(origRangeStart:origRangeEnd));
    title('Breath ButterWorth filter');
end

disp(['Section 4: curveCleaner - Finished in ' num2str(toc) ' seconds.' newline]);
end % Of function