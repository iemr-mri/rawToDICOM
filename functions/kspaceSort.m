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
    slices      = rawObj.Method.PVM_SPackArrNSlices; % number of slices (1 for SegFLASH)
    coils       = rawObj.Method.PVM_EncNReceivers;

    

    kspace      = reshape(kspaceRaw, [xData , coils, movieFrames, yData, slices, flowEncDir]);
    kspace      = permute(kspace,[1 4 5 3 6 2]); % [x, y, slices, movieFrames, flowEncDir, coils]
    
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
                            kspace_us(:,CSPhaseEncList(count), k, t, v, c)  = kspaceRaw(:,count_withCoil); % [x, y, slices, movieFrames, flowEncDir, coils]
                            count_withCoil                                   = count_withCoil + 1;
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
end