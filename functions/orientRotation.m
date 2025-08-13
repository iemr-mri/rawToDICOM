function rotIm = orientRotation(imageData,rawObj, visuParam)
    % Rotates and flips the image matrix according to geometry data from the method file

    % Should be one of H_F, A_P or L_R
    readOrient                          = rawObj.Method.PVM_SPackArrReadOrient;
    % Should be one of sagittal, coronal and axial
    sliceOrient                         = rawObj.Method.PVM_SPackArrSliceOrient;

    k                                   = 0;    % constant for how many 90 degree rotations to do
    f                                   = 0;    % constant for direction of flip, f=0 means no flip happens

    if contains(sliceOrient, 'sagittal')
        if contains(readOrient, 'H_F')
            k = 2;
        end

        if contains(readOrient, 'A_P')
            k = 1;
            f = 2;
        end
    end

    if contains(sliceOrient, 'coronal')
        if contains(readOrient, 'H_F')
            k = 2;
        end
        if contains(readOrient, 'L_R')
            k = 1;
            f = 2;
        end
    end

    if contains(sliceOrient, 'axial')
        if contains(readOrient, 'A_P')
            k = 2;
        end
        if contains(readOrient, 'L_R')
            k = 1;
            f = 2;
        end
    end
    
    imageData                           = rot90(imageData,k);
    
    if f ~= 0
        imageData                       = flip(imageData,f);
    end

    rotIm = imageData;
end