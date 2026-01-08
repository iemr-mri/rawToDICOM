function sortedStruct = slicePositionSort(unsortedStruct)
% Takes a struct of scans and sorts them according to the slice position found in the meta data files.
    positionList = zeros(length(unsortedStruct(2)));
    for scan=1:length(unsortedStruct)
        imagePath          = fullfile(unsortedStruct(scan).folder,unsortedStruct(scan).name);
        rawObj             = RawDataObject(imagePath, 'dataPrecision', 'double');
        position           = rawObj.Method.PVM_EffSliceOffset;

        positionList(scan) = position;
    end
    [~,index] = sort(positionList);
    sortedStruct = unsortedStruct(index);
end