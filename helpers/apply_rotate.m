function affine = apply_rotate(matrix, rad)
    % If angles are not specified, default to 0
    if nargin < 2
        rad = [0, 0, 0];
    end
    rad_x = rad(1); rad_y = rad(2); rad_z = rad(3);
    
    % Set up struct for rotation
    rot.x = [1, 0, 0; 
          0, cos(rad_x), -sin(rad_x); 
          0, sin(rad_x), cos(rad_x)];

    rot.y = [cos(rad_y), 0, sin(rad_y); 
          0, 1, 0; 
          -sin(rad_y), 0, cos(rad_y)];

    rot.z = [cos(rad_z), -sin(rad_z), 0; 
          sin(rad_z), cos(rad_z), 0; 
          0, 0, 1];
    
    % Find orientation matrix and vector from original affine transformation matrix
    [af_mat, af_vec] = to_matvec(matrix);

    % Rotate matrix and vector according to rotation struct
    rotated_mat = rot.z * (rot.y * (rot.x * af_mat));
    rotated_vec = rot.z * (rot.y * (rot.x * af_vec));
    
    % Convert back into affine transformation matrix
    affine =  from_matvec(rotated_mat, rotated_vec);
end