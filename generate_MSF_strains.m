function MSF_shape=generate_MSF_kinematics(R_helix,R_straight,omega,theta_helix_0,theta_straight_0,n_helix,n_straight,s,centerline_coordinates,isassigneddir,assigned_dir,uk,vk)
n_points=size(centerline_coordinates,2);
s_max=max(s);
s_min=min(s);
s_all=linspace(s_min,s_max,n_points);
x_k_helix_all=zeros(3,n_helix,n_points);
on_x_tau_k_helix_all=zeros(3,n_helix,n_points);
on_z_tau_k_helix_all=zeros(3,n_helix,n_points);
x_k_straight_all=zeros(3,n_helix,n_points);
on_x_tau_k_straight_all=zeros(3,n_helix,n_points);
on_z_tau_k_straight_all=zeros(3,n_helix,n_points);
D0_ik_all=get_frenet_frame_from_cart(s,centerline_coordinates,isassigneddir,assigned_dir);
pD_ik_ps_all=get_partial_D_ik_partial_s(D0_ik_all,s_all);
u_0k_all=zeros(3,n_points);
v_0k_all=zeros(3,n_points);
r_all=zeros(3,n_points);
D_ik_all=zeros(3,n_points);

for point=1:n_points
    s=s_all(point);
    [u_0k_at_s,v_0k_at_s]=get_ini_darboux_strain_u0_k_and_v0_k(squeeze(D0_ik_all(:,:,point)),squeeze(pD_ik_ps_all(:,:,point)));
    u_0k_all(:,point)=u_0k_at_s;
    v_0k_all(:,point)=v_0k_at_s;
    for helix=1:n_helix
        theta=theta_helix_0+2*pi*(helix-1)'/n_helix;
        x_k_helix_all(:,helix,point)=get_cable_location_x_k(R_helix,omega,theta,s);
        on_x_tau_k_helix_all(:,helix,point)=get_cable_direction_on_x_tau_k(R_helix,omega,theta,s);
        on_z_tau_k_helix_all(:,helix,point)=get_cable_direction_on_z_tau_k(squeeze(on_x_tau_k_helix_all(:,helix,point)),squeeze(D0_ik_all(:,:,point)));
        
    end
    for straight=1:n_straight
        theta=theta_straight_0+2*pi*(straight-1)'/n_straight;
        x_k_straight_all(:,straight,point)=get_cable_location_x_k(R_straight,0,theta,s);
        on_x_tau_k_straight_all(:,straight,point)=get_cable_direction_on_x_tau_k(R_straight,0,theta,s);
        on_z_tau_k_straight_all(:,straight,point)=get_cable_direction_on_z_tau_k(squeeze(on_x_tau_k_straight_all(:,helix,point)),squeeze(D0_ik_all(:,:,point)));
    end
end


plotinterval=round(n_points/100);
for helix=1:n_helix
    figure(1)
    plot3(squeeze(x_k_helix_all(1,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(2,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(3,helix,1:plotinterval:n_points)),'m');
    grid on
    hold on
    figure(2)
    plot3(squeeze(x_k_helix_all(1,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(2,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(3,helix,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    quiver3(squeeze(x_k_helix_all(1,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(2,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(3,helix,1:plotinterval:n_points)),squeeze(on_x_tau_k_helix_all(1,helix,1:plotinterval:n_points)), squeeze(on_x_tau_k_helix_all(2,helix,1:plotinterval:n_points)), squeeze(on_x_tau_k_helix_all(3,helix,1:plotinterval:n_points)),'r')
end
for straight=1:n_straight
    figure(1)
    plot3(squeeze(x_k_straight_all(1,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(2,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(3,straight,1:plotinterval:n_points)),'m');
    grid on
    hold on
    figure(2)
    plot3(squeeze(x_k_straight_all(1,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(2,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(3,straight,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    quiver3(squeeze(x_k_straight_all(1,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(2,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(3,straight,1:plotinterval:n_points)),squeeze(on_x_tau_k_straight_all(1,straight,1:plotinterval:n_points)), squeeze(on_x_tau_k_straight_all(2,straight,1:plotinterval:n_points)), squeeze(on_x_tau_k_straight_all(3,straight,1:plotinterval:n_points)),'b')

end


% 
% axis equal
% hold on
% quiver3(cart_coord_x(plot_index),cart_coord_y(plot_index),cart_coord_z(plot_index),ref_rod_direct_d1o((plot_index),1),ref_rod_direct_d1o((plot_index),2),ref_rod_direct_d1o((plot_index),3),'r')
% quiver3(cart_coord_x(plot_index),cart_coord_y(plot_index),cart_coord_z(plot_index),ref_rod_direct_d2o((plot_index),1),ref_rod_direct_d2o((plot_index),2),ref_rod_direct_d2o((plot_index),3),'g')
% quiver3(cart_coord_x(plot_index),cart_coord_y(plot_index),cart_coord_z(plot_index),ref_rod_direct_d3o((plot_index),1),ref_rod_direct_d3o((plot_index),2),ref_rod_direct_d3o((plot_index),3),'b')

end