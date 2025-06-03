function [mat, vec] = to_matvec(affine)
% Decompose a 4x4 affine matrix into a 3x3 matrix and a 1x3 vector.
% Input:
%   - affine: A 4x4 affine transformation matrix.
% Output:
%   - mat: a 3x3 matrix equivalent to the reshaped VisuCoreOrientation
%   - vec: a 1x3 vector equivalent to VisuCorePosition
    mat = affine(1:3, 1:3);
    vec = affine(1:3, 4);
end