function cineDataSaver(magnitude, saveName)
% Saves reconstructed magnitude and phase images together as complex valued
% images. Also saves resolution, which is needed in MRE-pipeline.

tic;

% Save it
fileName = [saveFolder '/imageData.mat'];
save(fileName, 'magnitude');

disp(['Section 10: Save data completed in ' num2str(toc) ' seconds.' newline])