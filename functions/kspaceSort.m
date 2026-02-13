function kspaceSorted = kspaceSort(rawObj)
    % Input
        % rawObj = image and meta data for the scan
    % Returns
        % kspaceSorted = kspace data matrix sorted according to [x, y, slices, movieFrames, flowEncDir, coils]

    %% 1) Load raw files and find k-space + relevant parameters
    kspaceRaw   = squeeze(rawObj.data{1});

    movieFrames = rawObj.Method.PVM_NMovieFrames;
    xData       = rawObj.Method.PVM_EncMatrix(1);
    yData       = rawObj.Method.PVM_EncMatrix(2);
    flowEncDir  = 1; % number of flow encoding directions
    slices      = rawObj.Method.PVM_SPackArrNSlices;
    coils       = rawObj.Method.PVM_EncNReceivers;

    kspace      = reshape(kspaceRaw, [coils, xData, movieFrames, yData, slices, flowEncDir]);
    kspace      = permute(kspace,[2 4 5 3 6 1]); % [x, y, slices, movieFrames, flowEncDir, coils]
    
    %% 2) Rearrange if kspace is undersampled
    if isfield(rawObj.Method, "CSPhaseEncList")
        CSPhaseEncList  = (rawObj.Method.CSPhaseEncList+4)*16;
        kspace_us       = zeros(xData, xData, slices, movieFrames, flowEncDir, coils);
        count           = 1;
        count_withCoil  = 1;
    
        for k = 1:slices
            for y = 1:size(kspace,2)
                for v = 1:flowEncDir
                    for t = 1:movieFrames
                        for c=1:coils
                            kspace_us(:,CSPhaseEncList(count), k, t, v, c)  = kspaceRaw(c,:,count); % [x, y, slices, movieFrames, flowEncDir, coils]
                        end
                        count = count+1;
                    end
                end
            end
        end
        kspaceSorted    = kspace_us; %the new, with zeros at the correct lines and data in the correct positions.
    
    % if kspace is not undersampled it should just return the reordered kspace
    else
        kspaceSorted    = kspace;
    end

    %% 3) Zero-fill partial echo
    if rawObj.Method.PVM_EncPft(1) > 1
        partialStart                           = round(xData*(rawObj.Method.PVM_EncPft(1) - 1)) + 1; % the starting index of the echo
        kspaceZero                             = zeros(round(xData*rawObj.Method.PVM_EncPft(1)), yData, slices, movieFrames, flowEncDir, coils);

        kspaceZero(partialStart:end,:,:,:,:,:) = kspaceSorted; % Put our kspace from step 2 into the zero-filled matrix
        
        kspaceSorted                           = kspaceZero;   % redefine kspaceSorted for function return
    end

end