function [reconkSpace, combReconkSpace, combRealSpace] = CSreconstructor(inputkSpace, pm, visualSwitch)
% Input
    % inputkSpace [coils x y slices frames MEGs] (once MRE is implemented)
    % rawObj = metadata for the scan
% Returns
    % reconData = reconstructed CS data

disp('Section 8: CS reconstruction');
tic;

pm.MEGdirections = 1; % Temp solution before future implementation of MRE
pm.numSlices     = 1;
figNumStart      = 400;
timestep         = 7.5/pm.newFrameNum;

if visualSwitch == true
    % Setup figures for later usage
    figureSetup(figNumStart);

    combGIFname     = [pm.animationFolder pm.expName '/' num2str(pm.scanNumber) ' scanTime ' num2str(round(pm.scanTime)) ' ' char(pm.absOrRel)];
    combGIFname     = [combGIFname ' breath[' num2str(pm.breathRange(1)) ',' num2str(pm.breathRange(2)) ']'];
    combGIFname     = [combGIFname ' frames ' num2str(pm.newFrameNum) ' CSrecon'];
    reconGIFname    = [combGIFname ' process'];
    phaseGIFname    = [combGIFname ' phase'];
    combGIFname     = [combGIFname '.gif'];
    reconGIFname    = [reconGIFname '.gif'];
    phaseGIFname    = [phaseGIFname '.gif'];

    if exist(combGIFname, 'file')
        delete(combGIFname);
        disp([combGIFname ' deleted for new gif.'])
    end
    if exist(reconGIFname, 'file')
        delete(reconGIFname);
        disp([reconGIFname ' deleted for new gif.'])
    end
    if exist(phaseGIFname, 'file')
        delete(phaseGIFname);
        disp([phaseGIFname ' deleted for new gif.'])
    end
end

% Basic matrices
inputkSpace         = reshape(inputkSpace,[pm.coilNum, pm.actualMatrix(1), pm.actualMatrix(2), pm.numSlices, pm.newFrameNum, pm.MechanicalPhases, pm.MEGdirections]);
hasDataMask         = inputkSpace ~= 0;          % Undersampling mask
noDataMask          = inputkSpace == 0;
tempRealImages      = zeros(size(inputkSpace));  % Compressed sensing operation matrix
reconkSpace         = zeros(size(inputkSpace));  % Compressed sensed output kspace
combReconkSpace     = zeros(pm.actualMatrix(1), pm.actualMatrix(2), pm.numSlices, pm.newFrameNum, pm.MechanicalPhases, pm.MEGdirections);
combRealSpace       = zeros(pm.actualMatrix(1), pm.actualMatrix(2), pm.numSlices, pm.newFrameNum, pm.MechanicalPhases, pm.MEGdirections);

% Loop through all data. CS must be performed on an image time series.
% Images series are unique based on their slice, coil, mechPhase, and MEG.
for slice = 1:pm.numSlices
    for coil = 1:pm.coilNum
        for mechPhase = 1:pm.MechanicalPhases
            % Transform all images at all timepoints to real space to facilitate temporal filtering
            for frame = 1:pm.newFrameNum
                for MEG = 1:pm.MEGdirections
                    tempRealImages(coil,:,:,slice,frame,mechPhase,MEG) = ifft2c(squeeze(inputkSpace(coil,:,:,slice,frame,mechPhase,MEG)));
                end
            end
    
            if visualSwitch == true && pm.showCSprocess == true
                % Visualize original data
                origkSpace = squeeze(inputkSpace(coil,:,:,slice,:,mechPhase,MEG));
                origRealIm = squeeze(tempRealImages(coil,:,:,slice,:,mechPhase,MEG));
                CSvisualizer(figNumStart,coil-1,squeeze(origkSpace(:,:,pm.selectFrame)),"kspace");
                CSvisualizer(figNumStart,coil-1+4,squeeze(origRealIm(:,:,pm.selectFrame)),"realspace");
            end
    
            % Matrix for monitoring the differences between iterations
            diffRMS = zeros(pm.MEGdirections, pm.maxIter);
    
            % Iterate to fill kspace
            for iter = 1:pm.maxIter
                % Perform temporal fft (after squeeze time is 3rd dimension)
                kspaceTemporal      = fft(squeeze(tempRealImages(coil,:,:,slice,:,mechPhase,:)), [], 3);
    
                % On the first iteration the thresholdvalue is set based on the temporal kspace data
                if iter==1
                    absData             = abs(kspaceTemporal);
                    threshVal           = prctile(absData(:), pm.CSprctThresh);
                end
                
                % Temporal kspace is thresholded according to the threshold value
                kspaceTemporalThresh = softThresh(kspaceTemporal, threshVal);
                
                % ifft back into real space
                tempRealImages(coil,:,:,slice,:,mechPhase,:) = ifft(kspaceTemporalThresh,[],3);
    
                % Compare reconstructed data to originally acquired data.
                % Stop CS when good enough.
                for MEG = 1:pm.MEGdirections
                    
                    % Transform back to 2D kspace and remove lines that were
                    % initially missing from the data to compare with original
                    % input.
                    kspaceAfterThresh = fft2c(squeeze(tempRealImages(coil,:,:,slice,:,mechPhase,MEG)));
                    maskedAfterThresh = kspaceAfterThresh .* squeeze(hasDataMask(coil,:,:,slice,:,mechPhase,MEG));
                    diffData          = squeeze(inputkSpace(coil,:,:,slice,:,mechPhase,MEG)) - maskedAfterThresh;
    
                    % Check diffRMS using mask
                    curDataMask       = squeeze(hasDataMask(coil,:,:,slice,:,mechPhase,MEG)); % [FIKS] Should this be mechPhase instead of MEG?
                    selectedDiffData  = diffData(curDataMask == true);
                    diffRMS(MEG,iter) = rms(selectedDiffData(:)); 
    
                    % Take original data, add new data for the blank areas
                    adjustedkSpace                       = squeeze(inputkSpace(coil,:,:,slice,:,mechPhase,MEG)) + kspaceAfterThresh .* squeeze(noDataMask(coil,:,:,slice,:,mechPhase,MEG));
                    tempRealImages(coil,:,:,slice,:,mechPhase,MEG) = squeeze(ifft2c(adjustedkSpace));
    
                    if visualSwitch == true && pm.showCSprocess == true
                        % Visualize
                        if mechPhase == 1
                            if coil == 1 % Arbitrarily selected the best quality coil
                                GIFmaker(figNumStart + coil-1+4, reconGIFname, timestep*3);
                            end
                            CSvisualizer(figNumStart,coil-1,squeeze(adjustedkSpace(:,:,pm.selectFrame)),"kspace");
                            CSvisualizer(figNumStart,coil-1+4,squeeze(tempRealImages(coil,:,:,slice,pm.selectFrame,mechPhase,MEG)),"realspace");
                        end
                        pause(timestep/5);
                    end
                end
                
                % Iteration stops when the RMS difference between an iteration 
                % and the original data is below 1% of the difference in the 
                % first iteration.
                if iter > 1
                    diffCurrent = diffRMS(MEG,iter) - diffRMS(MEG,iter-1);
                    diffFirst   = diffRMS(MEG,2) - diffRMS(MEG,1);
                    % if iter == 2 || mod(iter,10) == 0
                    %     disp(['Coil: ' num2str(coil) ' mechPhase: ' num2str(mechPhase) ' MEG: ' num2str(MEG) ' - Iter:' num2str(iter) ' - Current Diff = ' num2str(diffCurrent)  ' - % of first iteration: ' num2str(100*diffCurrent/diffFirst) '%']);
                    % end
                    if abs(diffCurrent/diffFirst) < pm.CSdiffThresh
                        disp(['Coil: ' num2str(coil) ' mechPhase: ' num2str(mechPhase) ' MEG: ' num2str(MEG) ' - Iter:' num2str(iter) ' - Current Diff = ' num2str(diffCurrent)  ' - % of first iteration: ' num2str(100*diffCurrent/diffFirst) '%']);
                        if visualSwitch == true
                            CSvisualizer(figNumStart,coil-1,squeeze(adjustedkSpace(:,:,pm.selectFrame)),"kspace");
                            CSvisualizer(figNumStart,coil-1+4,squeeze(tempRealImages(coil,:,:,slice,pm.selectFrame,mechPhase,MEG)),"realspace");
                        end
                        break
                    end
                end
            end % Iter loop

            % Perform lateral shift to compensate for off-center phase
            % shift (K: Somewhat arbitrary that it is exactly here, but 
            % was most convenient)
            pixelsToShift           = round(pm.phaseOffset/pm.resolution);
            shiftableTempRealImages = squeeze(tempRealImages(coil,:,:,slice,:,mechPhase,MEG));
            tempRealImages(coil,:,:,slice,:,mechPhase,MEG) = circshift(shiftableTempRealImages,pixelsToShift,2);

            % Add the images to the final matrix
            for MEG = 1:pm.MEGdirections
    
                combRealSpace(:,:,slice,:,mechPhase,MEG) = squeeze(combRealSpace(:,:,slice,:,mechPhase,MEG)) + squeeze(tempRealImages(coil,:,:,slice,:,mechPhase,MEG)).^2;
    
                % Set coil-wise reconkSpace
                reconkSpace(coil,:,:,slice,:,mechPhase,MEG) = fft2c(squeeze(tempRealImages(coil,:,:,slice,:,mechPhase,MEG)));
            end
        end % MechanicalPhases loop
    end % Coil loop
end % Slice loop
combRealSpace   = sqrt(combRealSpace);

% Calculate combReconkSpace and visualize results
upperBright = prctile(abs(combRealSpace), 99.75, "all");
for frame = 1:pm.newFrameNum
    for MEG = 1:pm.MEGdirections
        for slice = 1:pm.numSlices
            for mechPhase = 1:pm.MechanicalPhases
                % Perform the FFT on the x-y dimension of each slice
                combReconkSpace(:,:,slice,frame,mechPhase,MEG) = fftshift(fftn(ifftshift(combRealSpace(:,:,slice,frame,mechPhase,MEG))));
                
                % Visualize and record
                if visualSwitch == true
                    % All mechPhases (should) have the same magnitude
                    if mechPhase == 1
                        figure(figNumStart + 8);
                        imagesc(abs(squeeze(combRealSpace(:,:,slice,frame,mechPhase,MEG))));
                        clim([0,upperBright]);
                        title(sprintf('Magnitude Slice %d, Frame %d, mechPhase %d, MEG %d', slice, frame, mechPhase, MEG));
                        colormap gray;
                        colorbar;
                        GIFmaker(figure(figNumStart+8), combGIFname, timestep);
                    end

                    figure(figNumStart + 9);
                    complexImage = fftshift(ifftn(ifftshift(ifftshift(squeeze(reconkSpace(1,:,:,slice,frame,mechPhase,MEG)),1),2)));
                    imagesc(angle(complexImage));
                    title(sprintf('Phase Coil 1 Slice %d, Frame %d, mechPhase %d, MEG %d', slice, frame, mechPhase, MEG));
                    colormap parula;
                    colorbar;

                    % We only record one arbitrary frame
                    if frame == round(pm.newFrameNum/2)
                        GIFmaker(figure(figNumStart+9), phaseGIFname, timestep);
                    end
    
                    figure(figNumStart + 10);
                    imagesc(abs(squeeze(combReconkSpace(:,:,slice,frame,MEG)))');
                    title(sprintf('k-Space Slice %d, Frame %d, mechPhase %d, MEG %d', slice, frame, mechPhase, MEG));
                    colormap parula;
                    pause(timestep/pm.MechanicalPhases);
                end
            end
        end
    end
end

% Remove complex component from combRealSpace
combRealSpace = abs(combRealSpace);

disp(['Section 8: CS reconstruction completed in ' num2str(toc) ' seconds.' newline])
end % End of reconstructCS

%% Support functions

function figureSetup(figNumStart)
% Sets up figures that can be called later
    handle1 = figure(figNumStart);
    set(handle1, 'Position', [0,575,475,400]);
    title('Coil 1:');
    
    handle2 = figure(figNumStart + 1);
    set(handle2, 'Position', [475,575,475,400]);
    title('Coil 2:');
    
    handle3 = figure(figNumStart + 2);
    set(handle3, 'Position', [950,575,475,400]);
    title('Coil 3:');
    
    handle4 = figure(figNumStart + 3);
    set(handle4, 'Position', [1425,575,475,400]);
    title('Coil 4:');
    
    handle5 = figure(figNumStart + 4);
    set(handle5, 'Position', [0,100,475,400]);
    title('Coil 1:');
    
    handle6 = figure(figNumStart + 5);
    set(handle6, 'Position', [475,100,475,400]);
    title('Coil 2:');
    
    handle7 = figure(figNumStart + 6);
    set(handle7, 'Position', [950,100,475,400]);
    title('Coil 3:');
    
    handle8 = figure(figNumStart + 7);
    set(handle8, 'Position', [1425,100,475,400]);
    title('Coil 4:');

    handle9 = figure(figNumStart + 8);
    set(handle9, 'Position', [475,1200,475,400]);
    title('Combined magnitude:');

    handle10 = figure(figNumStart + 9);
    set(handle10, 'Position', [950,1200,475,400]);
    title('Combined phase:');

    handle11 = figure(figNumStart + 10);
    set(handle11, 'Position', [1425,1200,475,400]);
    title('Combined kspace:');
end

% Visualize CS images
function CSvisualizer(figNumStart, coil, data, kspaceOrNot)
    figure(figNumStart+coil);
    if kspaceOrNot == "kspace"
        imagesc(abs(data)');
    else
        imagesc(abs(data));
        colormap gray;
    end
    title(['Coil ' num2str(coil+1) ':']);
end

% 3D fft
function res = fft2c(x)
    fctr = size(x,1)*size(x,2);
    for n=1:size(x,3)
	    res(:,:,n) = 1/sqrt(fctr)*fftshift(fft2(ifftshift(x(:,:,n))));
    end
end

% 3D ifft
function res = ifft2c(x)
    fctr = size(x,1)*size(x,2);
    for n=1:size(x,3)
        res(:,:,n) = sqrt(fctr)*fftshift(ifft2(ifftshift(x(:,:,n))));
    end
end

function [y] = softThresh(x,lambda)
    % K: Cuts off anything below the threshold and shortens anything above by
    % the threshhold value. 
    % E.g. x = 4,   lambda = 4 => y = 0.
    %      x = 4.1, lambda = 4 => y = 0.1
    %      x = 10,  lambda = 4 => y = 6
    
    % Avoid division by zero
    x(x==0) = eps;
    
    % Perform thresholding
    y = (abs(x) > lambda).*(x.*(abs(x)-lambda)./(abs(x)));
end