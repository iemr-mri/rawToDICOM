function affine = build_affine(visuParam, method, resol)
% Builds the affine transformation matrix based on different ParaVision parameters.
% Based on method in https://github.com/BrkRaw/brkraw/blob/main/brkraw/lib/orient.py

% Input:
%   - visuParam
%   - method
% Output:
%   - affine: The transformation matrix containing the corrected VisuCoreOrientation and VisuCorePosition
    
    subj_pose                    = visuParam.VisuSubjectPosition;
    
    slicePos                     = visuParam.VisuCorePosition;
    orientVector                 = visuParam.VisuCoreOrientation;
    orientMatrix                 = reshape(orientVector, 3,3);

    %Find slice_orient, the direction where the z-value is the most prevalent
    slice_orient                 = method.PVM_SPackArrSliceOrient;
    % Convert resolution to diagonal matrix based on slice orientation
    if ismember(slice_orient, {'axial', 'sagittal'})
        resol                    = diag(resol(1)*[1; 1; 1]);
    else
        resol                    = diag(resol(1)*[1; 1; -1]);  % Negate the third diagonal element
    end

    % Combine rotation matrix and resolution
    orientMatrix = orientMatrix * resol;  % Transpose and multiply

    % Create the affine transformation matrix
    affine = from_matvec(orientMatrix, slicePos);

     % Adjusting for subject posture
    if ~isempty(subj_pose)
        switch subj_pose
            case 'Head_Supine'
                affine = apply_rotate(affine, [0, 0, pi]);  % Rotate around Z by 180 degrees
            case 'Head_Prone'
                % No adjustment needed for 'Head_Prone'
            case 'Head_Left'
                affine = apply_rotate(affine, [0, 0, pi/2]);  % Rotate around Z by 90 degrees
            case 'Head_Right'
                affine = apply_rotate(affine, [0, 0, -pi/2]);  % Rotate around Z by -90 degrees
            case {'Foot_Supine', 'Tail_Supine'}
                affine = apply_rotate(affine, [pi, 0, 0]);  % Rotate around X by 180 degrees
            case {'Foot_Prone', 'Tail_Prone'}
                affine = apply_rotate(affine, [0, pi, 0]);  % Rotate around Y by 180 degrees
            case {'Foot_Left', 'Tail_Left'}
                affine = apply_rotate(affine, [0, 0, pi/2]);  % Rotate around Z by 90 degrees
            case {'Foot_Right', 'Tail_Right'}
                affine = apply_rotate(affine, [0, 0, -pi/2]);  % Rotate around Z by -90 degrees
            otherwise
                error('NotIntegrated: Unknown subject pose.');
        end
    end

    % Base adjustment for all quadripeds
    affine = apply_rotate(affine, [-pi/2, pi, 0]);
    
    % RAS (ParaVision) to LPS (DICOM) transformation
    diag_matrix = diag([-1, -1, 1, 1]);
    affine = diag_matrix * affine;  
end