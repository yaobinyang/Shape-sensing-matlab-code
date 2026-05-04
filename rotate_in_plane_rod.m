function [r_all_rotated, D_ik_all_rotated] = rotate_in_plane_rod(r_all, D_ik_all, rod_axis_rs)
    % 1. Get original vector and its length
    vec_orig = rod_axis_rs(:,end) - rod_axis_rs(:,1);
    L_orig = norm(vec_orig);
    axis_vec_orig = vec_orig / L_orig;
    
    % 2. Get deformed vector and its length
    vec_deformed = r_all(:,end) - r_all(:,1);
    L_def = norm(vec_deformed);
    axis_vec_deformed = vec_deformed / L_def;
    
    % 3. Calculate stretch factor
    scale_factor = L_orig / L_def;
    
    % Anonymous functions for rotation
    GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0; norm(cross(A,B)) dot(A,B) 0; 0 0 1];
    FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
    UU = @(Fi,G) Fi*G*inv(Fi);
    
    % Safety check: if vectors are already perfectly parallel, avoid NaN in rotation
    if norm(cross(axis_vec_deformed, axis_vec_orig)) < 1e-10
        U = eye(3);
    else
        U = UU(FFi(axis_vec_deformed, axis_vec_orig), GG(axis_vec_deformed, axis_vec_orig));
    end
    
    r_all_rotated = zeros(size(r_all));
    D_ik_all_rotated = zeros(size(D_ik_all));
    
    for i = 1:size(r_all, 2)
        % Apply rotation AND scaling (stretching) to the positions
        r_all_rotated(:, i) = scale_factor * U * (r_all(:, i) - r_all(:, 1)) + r_all(:, 1);
        
        % Apply ONLY rotation to the directors (directors must remain unit vectors)
        D_ik_all_rotated(:, :, i) = U * D_ik_all(:, :, i);
    end
end