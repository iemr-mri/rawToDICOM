function orient_order = get_orient_order(orient_matrix)
    % Returns the indices of the maximum arguments of the orientation matrix
    orient_order = zeros(1,3);
    [~, orient_order(1)] = max(abs(orient_matrix(:,1)));
    [~, orient_order(2)] = max(abs(orient_matrix(:,2)));
    [~, orient_order(3)] = max(abs(orient_matrix(:,3)));
end