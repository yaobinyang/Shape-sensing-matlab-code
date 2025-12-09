function D_kls= get_frenet_frame_from_cart(s,rod_axis_rs,isassigneddir,assigned_dir)
n_points=length(s);
gradient_s_ds=gradient(s);
partial_r_partial_s=gradient(rod_axis_rs)./(repmat(gradient_s_ds,3,1));
parfor n=1:n_points
    ref_rod_direct_d3o(:,n)=partial_r_partial_s(:,n)/norm(partial_r_partial_s(:,n));
end
partial_d3o_partial_s=gradient(ref_rod_direct_d3o)./(repmat(gradient_s_ds,3,1));
parfor n=1:n_points
    if isassigneddir
        ref_rod_direct_d1o(:,n)=assigned_dir;
    else
    % ref_rod_direct_d1o(:,n)=[0 0 1];
    % ref_rod_direct_d1o(:,n)=partial_d3o_partial_s(:,n)/norm(partial_d3o_partial_s(:,n));
    ref_rod_direct_d1o(:,n)=(partial_d3o_partial_s(:,n)/norm(partial_d3o_partial_s(:,n)))/sign(partial_d3o_partial_s(2,n));
    end
end
parfor n=1:n_points
    ref_rod_direct_d2o(:,n)=cross(ref_rod_direct_d3o(:,n),ref_rod_direct_d1o(:,n));
end

figure(9001)
hold off
plot3(rod_axis_rs(1,:), rod_axis_rs(2,:), rod_axis_rs(3,:),'m');
grid on
figure(9002)
hold off
plot3(rod_axis_rs(1,:), rod_axis_rs(2,:), rod_axis_rs(3,:),'m');
grid on
plot_index=1:round(n_points/50):n_points;
axis equal
hold on
quiver3(rod_axis_rs(1,plot_index),rod_axis_rs(2,plot_index),rod_axis_rs(3,plot_index),ref_rod_direct_d1o(1,(plot_index)),ref_rod_direct_d1o(2,(plot_index)),ref_rod_direct_d1o(3,(plot_index)),'r')
quiver3(rod_axis_rs(1,plot_index),rod_axis_rs(2,plot_index),rod_axis_rs(3,plot_index),ref_rod_direct_d2o(1,(plot_index)),ref_rod_direct_d2o(2,(plot_index)),ref_rod_direct_d2o(3,(plot_index)),'g')
quiver3(rod_axis_rs(1,plot_index),rod_axis_rs(2,plot_index),rod_axis_rs(3,plot_index),ref_rod_direct_d3o(1,(plot_index)),ref_rod_direct_d3o(2,(plot_index)),ref_rod_direct_d3o(3,(plot_index)),'b')

D_kls(:,1,:)=ref_rod_direct_d1o;
D_kls(:,2,:)=ref_rod_direct_d2o;
D_kls(:,3,:)=ref_rod_direct_d3o;
end