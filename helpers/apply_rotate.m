function affine = apply_rotate(mat, rad_x=0, rad_y=0, rad_z=0)
    rot.x = [1, 0, 0; 
          0, cos(rad_x), -sin(rad_x); 
          0, sin(rad_x), cos(rad_x)];

    rot.y = [cos(rad_y), 0, sin(rad_y); 
          0, 1, 0; 
          -sin(rad_y), 0, cos(rad_y)];

    rot.z = [cos(rad_z), -sin(rad_z), 0; 
          sin(rad_z), cos(rad_z), 0; 
          0, 0, 1];

    [af_mat, af_vec] = to_matvec(matrix);
    rotated_mat = rot.z * (rot.y * (rot.x * af_mat));
    rotated_vec = rot.z * (rot.y * (rot.x * af_vec));
    
    affine =  from_matvec(rotated_mat, rotated_vec);
end