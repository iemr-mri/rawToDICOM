function recon = reconCheck()
    disp('Self-gated DICOM files already exist for this subject.')
    recon = input('Repeat all SG reconstruction for this subject? (SAX slices must be batch processed) [yes = 1, no = 0]');
    while recon ~= 0 && recon ~= 1
        recon = input('Input not recognized. Repeat all SG reconstruction for this subject? [yes = 1, no = 0]');
    end

    if recon == 1
        disp('Repeating SG reconstruction.')
    end

    if recon == 0
        disp('Skipping SG reconstruction for this subject.')
end

