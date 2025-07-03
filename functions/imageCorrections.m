function final_im = imageCorrections(imageData, rawObj, visuParam)    
    
    %% Fixing offset in phase direction
    resolution                          = rawObj.Method.PVM_FovCm(1)/size(imageData,1); % in cm/pixels
    offset_mm                           = rawObj.Method.PVM_Phase1Offset(1); % in mm
    offset_pixels                       = (offset_mm/10)/resolution;

    im_shifted                          = circshift(imageData, -round(offset_pixels), 2); 
    imageData                           = im_shifted;

    %% Normalize image
    % Calculate normalization factor to scale maximum intensity to 30,000
    normFactor                          = 30000/max(imageData,[],'all');
    % Normalize image data
    Inorm                               = normFactor .* imageData;
    final_im                            = int16(Inorm);
    
    if contains(visuParam.VisuAcquisitionProtocol, 'LAX')
        final_im                           = rot90(final_im,2);
    end


end