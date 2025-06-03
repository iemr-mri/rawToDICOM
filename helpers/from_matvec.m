function affine = from_matvec(mat, vec)
%% from_matvec
% Makes an affine transform from a matrix and a vector
% as described in https://github.com/BrkRaw/brkraw/blob/main/brkraw/api/helper/orientation.py
% Input:
%   - mat: a 3x3 matrix equivalent to the reshaped VisuCoreOrientation
%   - vec: a 1x3 vector equivalent to VisuCorePosition
% Output: 
%   - affine, 4x4 affine transformation with rmat in the top left corner and pos in the last column

    affine          = eye(4);
    affine(1:3,1:3) = mat;
    affine(1:3,4)   = vec;

end