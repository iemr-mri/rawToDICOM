function rotIm = orientRotation(imageData,rawObj, visuParam)
    % Rotates and flips the image matrix according to geometry data from the method file

    % view: 
    %   - LAX also contain 'LAX' in their protocol name
    %   - 'CS_191021' is the unique string in the SAX protocol
    view                                = visuParam.VisuAcquisitionProtocol;

    % Should be one of H_F, A_P or L_R
    readOrient                          = rawObj.Method.PVM_SPackArrReadOrient;

    k                                   = 0;     % constant for how many 90 degree rotations to do
    f                                   = false; % boolean for image flipping

    if contains(view, 'LAX')
        if contains(readOrient, 'H_F')
            k = 2;
        end
        if contains(readOrient, 'A_P')
            k = 1;
        end
        if contains(readOrient, 'L_R') % have not seen this case, but if it appears for coronal LAX images hopefully it works similarly to A_P
            k = 1;
        end
    end
    if contains(view, 'CS_191021')
        if contains(readOrient, 'A_P')
            k = 2;
        end
        if contains(readOrient, 'L_R')
            k = 1;
            f = true;
        end
    end
    
    imageData                           = rot90(imageData,k);
    
    if f == true
        imageData                       = flip(imageData,2);
    end

    rotIm = imageData;
end