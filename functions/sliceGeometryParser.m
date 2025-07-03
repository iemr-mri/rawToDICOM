function [position, orientation] = sliceGeometryParser(sliceGeo)
    
    PositionParser                      = sliceGeo{2,1};
    slicePos                            = str2num(PositionParser);
    OrientationParser                   = sliceGeo{1,1};
    orientVector                        = str2num(OrientationParser);

    orientMatrix                        = reshape(orientVector, 3,3)';
    affine                              = from_matvec(orientMatrix, slicePos);
    
    % RAS (ParaVision) to LPS (DICOM) transformation
    diag_matrix                         = diag([-1, -1, 1, 1]);
    %affine                              = diag_matrix * affine;
    
    [mat, vec]                          = to_matvec(affine);
    position                            = vec;
    orientation                         = reshape(mat',1,[]);
    %orientation                         = orientVector;
end