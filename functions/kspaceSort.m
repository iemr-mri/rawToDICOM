function kspaceSorted = kspaceSort(rawObj, visuParam)
    % Input
        % rawObj = image and meta data for the scan
    % Returns
        % kspaceSorted = kspace data matrix sorted according to [x, y, slices, movieFrames, flowEncDir, coils]
    method = rawObj.Method;
    %% 1) Load raw files and find k-space + relevant parameters
    kspaceRaw   = squeeze(rawObj.data{1});

    movieFrames = method.PVM_NMovieFrames;

    
    xData       = method.PVM_EncMatrix(1);
    yData       = method.PVM_EncMatrix(2);
    flowEncDir  = 1; % number of flow encoding directions
    slices      = method.PVM_SPackArrNSlices; % number of slices (1 for SegFLASH)
    coils       = method.PVM_EncNReceivers;

    

    kspace      = reshape(kspaceRaw, [xData , coils, movieFrames, yData, slices, flowEncDir]);
    kspace      = permute(kspace,[1 4 5 3 6 2]); % [x, y, slices, movieFrames, flowEncDir, coils]
    
    %% 2) Rearrange if kspace is undersampled by CS
    if isfield(method, "CSPhaseEncList")
        CSPhaseEncList  = (method.CSPhaseEncList+4)*16;
        kspace_us       = zeros(max(method.PVM_Matrix), max(method.PVM_Matrix), slices, movieFrames, flowEncDir, coils);
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
        kspace    = kspace_us; %the new, with zeros at the correct lines and data in the correct positions.
    end

    %% 3) Zero-fill if kspace is undersampled by partial fourier
    if max(visuParam.VisuAcqPartialFourier) ~= 1
        xP          = visuParam.VisuCoreSize(1); % x pixels in final image
        kspace_us   = zeros(xP, xY, slices, movieFrames, flowEncDir, coils);
        
        startX      = floor((xP - xData) / 2) + 1;

        kspace_us(:, :, :, :, :, :)   = kspace(:, :, :, :, :, :);

        kspace    = kspace_us; %the new, with zeros at the correct lines and data in the correct positions.
    end
        

    %% 4) If kspace is not undersampled it should just return the reordered kspace from 1)
    kspaceSorted    = kspace;
end