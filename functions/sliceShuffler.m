function [shuffledData] = sliceShuffler(imageData, method)
% Shuffles dataset containing multiple slices so they are ordered correctly
% One would think there is a slice order list in the metadata, but I have not found it
    
    sliceOrder   = zeros(1,size(imageData,3));
    sliceHalf    = round(length(sliceOrder)/2);
    
    % Create a vector of slice counters from 1 to N
    sliceCounter = (1:size(imageData, 3));

    % Assign values based on odd/even index
    sliceOrder(1:2:end) = sliceCounter(1:sliceHalf);      % odd indices
    sliceOrder(2:2:end) = sliceCounter(sliceHalf+1:end);   % even indices
    
    shuffledData = zeros(size(imageData), "uint16");
    for slice=1:length(sliceOrder)
        shuffledData(:,:,slice,:) = imageData(:,:,sliceOrder(slice),:);
    end

end