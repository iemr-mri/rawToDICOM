function corrected_pose = reversed_pose_correction(pose, rmat, distance)
    reversed_pose = rmat * pose';
    reversed_pose(end) =  reversed_pose(end) + distance;
    corrected_pose = rmat' * reversed_pose;
end