function orient_order = get_orient_order(orient_matrix)
    % Returns the indices of the maximum arguments of the orientation matrix
    orient_order = zeros(1,3);
    for col = 1:3
    [~, orient_order(col)] = max(abs(orient_matrix(:,col)));
    end
end