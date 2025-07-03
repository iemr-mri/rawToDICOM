function kspaceCS = reconstructCS(kspace)
    % Input
        % kspace = sorted kspace data matrix [x, y, slices, movieFrames, flowEncDir, coils]
    % Returns
        % reconData = reconstructed CS data
        
    %% 1) 
    usMask              = kspace~=0; % undersampling mask

    xData               = size(kspace,1); yData         = size(kspace,2);
    slices              = size(kspace,3); movieFrames   = size(kspace,4);
    MEG                 = size(kspace,5); coils         = size(kspace,6);

    maxIter             = 50;
    CS_percentThresh    = 50;

    imTemp              = zeros(xData, yData, movieFrames, MEG);  %Compressed sensing operation matrix
    kspaceCS            = zeros(size(kspace));  %Compressed sensed output kspace

    %% CS recon loop
    % Top layer - Slices
    for slice = 1:slices
        % CS reconstruction is then performed per coil
        for coil = 1:coils
            % All frames in the slice are ifft into im space
            for v = 1:MEG
                imTemp(:,:,:,v) = ifft2c(squeeze(kspace(:,:, slice, :, v, coil)));
            end
            
            % Reconstruction iteration
            for iter = 1:maxIter

                % The image data is fft along the temporal dimension (movieFrames, dim=3)
                kspaceTemporal      = fft(imTemp(:,:,:,v), [], 3);
                
                % On the first iteration the thresholdvalue is set based on the temporal kspace data
                if iter==1
                    absData             = abs(kspaceTemporal);
                    threshVal           = prctile(absData(:), CS_percentThresh);
                end
                
                % Temporal kspace is thresholded according to the threshold value
                kspaceTemporalThresh = SoftThresh(kspaceTemporal, threshVal);
                
                % ifft back into im space
                imTemp(:,:,:,v) = ifft(kspaceTemporalThresh,[],3);

                %%
                %This loop compares reconstructed data to originally acquired
                %data, which is used to determine CS endpoint
                for v = 1:MEG
                    kspaceAfterThresh   = fft2c(imTemp(:,:,:,v)).*(squeeze(usMask(:,:,slice,:,v,coil)));  
                
                    diff_data_mat       = kspaceAfterThresh - squeeze(kspace(:,:,slice,:,v,coil)); %difference between acquired data and corresponding data post thresholding
                
                    %creates diff data matrix
                    for t = 1:movieFrames
                        for y = 1:yData
                            if max(squeeze(kspace(:, y, slice, t, v, coil))) ~= 0
                                diff_data_mat_no0(:,y,t)    = diff_data_mat(:,y,t);
                            end
                        end
                    end
                
                    %rms diff value, so we can plot average differences
                    diffData(:,iter,v) = rms(diff_data_mat_no0(:));
                    
                    %inverse fourier transform to image domain
                    adjustedKspace = squeeze(kspace(:,:,slice,:,v,coil)) + fft2c(squeeze(imTemp(:,:,:,v))).*(1-squeeze(usMask(:,:,slice,:,v,coil)));
                    imTemp(:,:,:,v) = ifft2c(squeeze(adjustedKspace));
                    
                    % figure(1)
                    % imagesc(abs(squeeze(imTemp(:,:,1,v))))
                    % colormap('gray')
                    % 
                    % figure(2)
                    % imagesc(abs(squeeze(adjustedKspace(:,:,1)))')
                end
            
            
                %This section stops iterative reconstruction based on the gradient of the
                %difference between acquired and reconstructed vs iterations.
                if iter>1
                    diffIter = diffData(:,iter,v)-diffData(:,iter-1,v);
                    diffIter1st = diffData(:,2,v)-diffData(:,1,v);
                    if abs(diffIter/diffIter1st)<0.01
                        break
                    end
                end
            end
            disp(['Iterations: ', num2str(iter), ', with a final change of : ', num2str((diffIter/diffIter1st)*100), '% from initial iteration.'])
            for v = 1:MEG
                kspaceCS(:,:,slice,:,v,coil) = fft2c(squeeze(imTemp(:,:,:,v)));
            end
        end
    end
end