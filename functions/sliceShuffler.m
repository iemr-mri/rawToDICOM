function [shuffledData] = sliceShuffler(imageData, method)
% Shuffles dataset containing multiple slices so they are ordered correctly
% One would think there is a slice order list in the metadata, but I have not found it
    sliceOrder = method.PVM_ObjOrderList + 1;
    
    shuffledData = zeros(size(imageData));

    for slice=1:length(sliceOrder)
        shuffledData(:,:,slice,:) = imageData(:,:,sliceOrder(slice),:);
    end
end