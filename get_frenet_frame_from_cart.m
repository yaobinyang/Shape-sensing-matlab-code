% function D_kls= get_frenet_frame_from_cart(s,rod_axis_rs,isassigneddir,assigned_dir)
% n_points=length(s);
% gradient_s_ds=gradient(s);
% partial_r_partial_s=gradient(rod_axis_rs)./(repmat(gradient_s_ds,3,1));
% parfor n=1:n_points
%     ref_rod_direct_d3o(:,n)=partial_r_partial_s(:,n)/norm(partial_r_partial_s(:,n));
% end
% partial_d3o_partial_s=gradient(ref_rod_direct_d3o)./(repmat(gradient_s_ds,3,1));
% parfor n=1:n_points
%     if isassigneddir
%         ref_rod_direct_d1o(:,n)=assigned_dir;
%     else
%     % ref_rod_direct_d1o(:,n)=[0 0 1];
%     % ref_rod_direct_d1o(:,n)=partial_d3o_partial_s(:,n)/norm(partial_d3o_partial_s(:,n));
%     ref_rod_direct_d1o(:,n)=(partial_d3o_partial_s(:,n)/norm(partial_d3o_partial_s(:,n)))/sign(partial_d3o_partial_s(2,n));
%     end
% end
% parfor n=1:n_points
%     ref_rod_direct_d2o(:,n)=cross(ref_rod_direct_d3o(:,n),ref_rod_direct_d1o(:,n));
% end
% 
% figure(9001)
% hold off
% plot3(rod_axis_rs(1,:), rod_axis_rs(2,:), rod_axis_rs(3,:),'m');
% grid on
% figure(9002)
% hold off
% plot3(rod_axis_rs(1,:), rod_axis_rs(2,:), rod_axis_rs(3,:),'m');
% grid on
% plot_index=1:round(n_points/50):n_points;
% axis equal
% hold on
% quiver3(rod_axis_rs(1,plot_index),rod_axis_rs(2,plot_index),rod_axis_rs(3,plot_index),ref_rod_direct_d1o(1,(plot_index)),ref_rod_direct_d1o(2,(plot_index)),ref_rod_direct_d1o(3,(plot_index)),'r')
% quiver3(rod_axis_rs(1,plot_index),rod_axis_rs(2,plot_index),rod_axis_rs(3,plot_index),ref_rod_direct_d2o(1,(plot_index)),ref_rod_direct_d2o(2,(plot_index)),ref_rod_direct_d2o(3,(plot_index)),'g')
% quiver3(rod_axis_rs(1,plot_index),rod_axis_rs(2,plot_index),rod_axis_rs(3,plot_index),ref_rod_direct_d3o(1,(plot_index)),ref_rod_direct_d3o(2,(plot_index)),ref_rod_direct_d3o(3,(plot_index)),'b')
% 
% D_kls(:,1,:)=ref_rod_direct_d1o;
% D_kls(:,2,:)=ref_rod_direct_d2o;
% D_kls(:,3,:)=ref_rod_direct_d3o;
% end
function D_kls = get_frenet_frame_from_cart(s, rod_axis_rs, isassigneddir, assigned_dir)

    n_points = length(s);
    D_kls = zeros(3, 3, n_points);

    % 1. Calculate Tangents (d3o) Robustly
    % ---------------------------------------
    % Calculate derivative for each row separately
    dx_ds = gradient(rod_axis_rs(1, :), s);
    dy_ds = gradient(rod_axis_rs(2, :), s);
    dz_ds = gradient(rod_axis_rs(3, :), s);

    tangent_raw = [dx_ds; dy_ds; dz_ds];
    tangent_norms = vecnorm(tangent_raw, 2, 1);
    T = tangent_raw ./ tangent_norms; % This is d3o (3xN)

    % 2. Initialize the First Frame (Index 1)
    % ---------------------------------------
    % We need ONE valid starting normal vector.

    if isassigneddir
        % Use user-provided direction, projected to be perp to Tangent
        d3_1 = T(:,1);
        % Make sure assigned_dir is a column vector
        if size(assigned_dir, 2) > 1, assigned_dir = assigned_dir'; end

        % Gram-Schmidt projection: v_perp = v - (v . t) * t
        d1_1 = assigned_dir - (dot(assigned_dir, d3_1) * d3_1);
        d1_1 = d1_1 / norm(d1_1);
    else
        % Auto-generate a valid normal for the first point
        d3_1 = T(:,1);
        % Try global 'up' vector
        guess = [0; 0; 1];
        % If tangent is too close to 'up', use 'right'
        if abs(dot(d3_1, guess)) > 0.99
            guess = [1; 0; 0];
        end
        d1_1 = cross(d3_1, guess);
        d1_1 = d1_1 / norm(d1_1);
    end

    % Calculate initial Binormal
    d2_1 = cross(T(:,1), d1_1);

    % Store first frame
    D_kls(:, 1, 1) = d1_1;
    D_kls(:, 2, 1) = d2_1;
    D_kls(:, 3, 1) = T(:,1);

    % 3. Propagate Frame (Parallel Transport)
    % ---------------------------------------
    % Instead of recalculating from scratch (which causes flips),
    % we rotate the PREVIOUS frame to align with the NEW tangent.

    for i = 1 : n_points - 1
        % Current and Next Tangents
        t_curr = T(:, i);
        t_next = T(:, i+1);

        % Current Normal/Binormal
        n_curr = D_kls(:, 1, i);
        b_curr = D_kls(:, 2, i);

        % Calculate rotation axis (cross product of tangents)
        % This is the axis we need to rotate around to get from t_curr to t_next
        v = cross(t_curr, t_next);
        c = dot(t_curr, t_next);

        % If vectors are nearly identical (straight line), just copy the frame
        if norm(v) < 1e-10 
             n_next = n_curr;
             b_next = b_curr;
        else
            % Rodrigues' Rotation Formula (Simplified for aligning two vectors)
            % This rotates n_curr by the exact angle required to align t_curr to t_next
            % without adding any extra "twist".

            % Skew-symmetric matrix of v
            vx = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];

            % Rotation Matrix R
            R = eye(3) + vx + (vx^2) * (1 / (1 + c));

            n_next = R * n_curr;
            b_next = R * b_curr;
        end

        % Normalize to prevent error accumulation
        n_next = n_next / norm(n_next);
        b_next = b_next / norm(b_next);

        % Store
        D_kls(:, 1, i+1) = n_next;
        D_kls(:, 2, i+1) = b_next;
        D_kls(:, 3, i+1) = t_next;
    end

% 4. Visualization (Verification)
    % ---------------------------------------
    figure(9002); clf; hold on;

    % Plot the curve path
    plot3(rod_axis_rs(1,:), rod_axis_rs(2,:), rod_axis_rs(3,:), 'k-', 'LineWidth', 2);

    % Downsample for cleaner arrows (e.g., every 20th point)
    idx = 1:ceil(n_points/20):n_points; 

    % --- EXTRACT DATA SAFELY ---
    % Get Positions (1 x M vectors)
    px = rod_axis_rs(1,idx);
    py = rod_axis_rs(2,idx);
    pz = rod_axis_rs(3,idx);

    % Get Normal Vector d1 (1 x M vectors)
    % We use 'reshape' to force them into Row Vectors [1, M]
    d1x = reshape(D_kls(1, 1, idx), 1, []);
    d1y = reshape(D_kls(2, 1, idx), 1, []);
    d1z = reshape(D_kls(3, 1, idx), 1, []);

    % Get Binormal Vector d2 (1 x M vectors)
    d2x = reshape(D_kls(1, 2, idx), 1, []);
    d2y = reshape(D_kls(2, 2, idx), 1, []);
    d2z = reshape(D_kls(3, 2, idx), 1, []);

    % --- PLOT ARROWS ---
    % Scale factor for visibility (optional)
    scale = 0.5; 

    quiver3(px, py, pz, d1x, d1y, d1z, scale, 'r', 'LineWidth', 2, 'DisplayName', 'd1 (Normal)');
    quiver3(px, py, pz, d2x, d2y, d2z, scale, 'g', 'LineWidth', 2, 'DisplayName', 'd2 (Binormal)');

    grid on; axis equal; view(3);
    title('Parallel Transport Frame');
    legend show;
    hold off;
end