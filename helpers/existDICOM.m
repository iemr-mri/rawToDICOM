function existD = existDICOM(destination, scanName)
%% existDICOM - checks if a DICOM file exists at designated destination and returns corresponding boolean
% Input
% - destination: path of DICOM save destination to check for existing file
% Output
% - existD: boolean to indicate wether or not the DICOM file already exist
    
    % Initial check to look for subject folder and make one if it doesn't exist already
    [dirPath]   = fileparts(destination);
    % "7" specifically checks if dirPath is a folder
    if exist(dirPath, 'dir') ~= 7
        mkdir(dirPath) %
    end

    % Check if DICOM file already exist and return boolean
    if exist([destination,'.dcm'], 'file')
        disp('--------')
        disp(['DICOM file ', scanName, '.dcm already exist.'])
        existD = 1;
    else
        existD = 0;
    end
end