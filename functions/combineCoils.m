function [imData] = combineCoils(kspace)
% Function to combine kspace data from multiple coils and MEGS using the Sum-of-Squares method
% Input:
    % kspace = [xData, yData, slices, movieFrames, MEG, coils]
% Output
    % imData = [xData, yData, slices, movieFrames, MEG]
    xData               = size(kspace,1); yData         = size(kspace,2);
    slices              = size(kspace,3); movieFrames   = size(kspace,4);
    MEG                 = size(kspace,5); coils         = size(kspace,6);

    imData = zeros(xData, yData, slices, movieFrames, MEG);

    for slice = 1:slices
        for frame = 1:movieFrames
            for v = 1:MEG
                combined_image = zeros(xData, yData);
                for coil = 1:coils
                    kspace_matrix = squeeze(kspace(:,:,slice,frame,v,coil));
                    im_matrix = ifftshift(fft2(fftshift(kspace_matrix)));
                    combined_image = combined_image + abs(im_matrix).^2;
                end
                imData(:,:,slice,frame,v) = sqrt(combined_image);
            end
        end
    end
end

    % eldre måte med mean istedenfor sum of squares
    % for k=1:K
    %     for t=1:T
    %         %disp(t)
    %         StudyData.Magn(:,:,k,t) = mean(mean(abs(im(:,:,k,t,:,:)),5),6);
    %     end
    % end