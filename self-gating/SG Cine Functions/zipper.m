function [zipped, realZipped] = zipper(magnitude, pm, visualSwitch)
% Performs zero-filling interpolation to "increase" image quality for
% visual purposes

tic;
% Setup input data in the middle of expanded zeros matrix
dims        = size(magnitude);
zipped      = zeros([pm.actualMatrix*2, dims(3:4)]);
startX      = pm.actualMatrix(1)/2 + 1;
endX        = pm.actualMatrix(1)*3/2;
startY      = pm.actualMatrix(2)/2 + 1;
endY        = pm.actualMatrix(2)*3/2;

rawkSpace   = fftshift(fftshift(fft(fft(magnitude,[],1),[],2),1),2);
zipped(startX:endX,startY:endY,:,:) = rawkSpace;

% Transform to realspace
realZipped = zeros(size(zipped));
for MEG = 1:dims(3)
    for frame = 1:dims(4)
        curkSpace                   = squeeze(zipped(:,:,MEG,frame));
        realZipped(:,:,MEG,frame)   = abs(ifftn(ifftshift(curkSpace)));
    end
end

if visualSwitch == true
    % Setup GIFname and delete existing GIF if there
    GIFname = [pm.animationFolder pm.expName '/' num2str(pm.scanNumber) ' scanTime ' num2str(round(pm.scanTime)) ' ' char(pm.absOrRel)];
    GIFname = [GIFname ' breath[' num2str(pm.breathRange(1)) ',' num2str(pm.breathRange(2)) ']'];
    GIFname = [GIFname ' frames ' num2str(pm.newFrameNum)];
    GIFname = [GIFname ' CSrecon zipped.gif'];
    
    if exist(GIFname, 'file')
        delete(GIFname);
        disp([GIFname ' deleted for new gif.'])
    end
    
    % Set total time to 7.5 seconds
    timestep = 4/pm.newFrameNum;
    
    % Set upper brightness limit
    upperBright = prctile(realZipped, 99.75, "all");
    
    % Visualize and create GIF
    figHandle = figure(600);
    for frame = 1:dims(4)
        imagesc(squeeze(realZipped(:,:,1,frame)));
        clim([0,upperBright]);
        title(['Zipped real space Frame ' num2str(frame)])
        colormap gray;
        colorbar;
        GIFmaker(figHandle,GIFname,timestep);
    end
end

disp(['Section 9: Zero-filling interpolation completed in ' num2str(toc) ' seconds.'])

end % of function