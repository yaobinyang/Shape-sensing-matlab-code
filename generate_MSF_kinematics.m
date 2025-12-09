function [r_all,D_ik_all,helix_cable_strain_all,straight_cable_strain_all,u_0k_all,v_0k_all]=generate_MSF_kinematics(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,centerline_coordinates,isassigneddir,assigned_dir,uk_all,vk_all)
n_points=size(centerline_coordinates,2);
s_max=max(s_all);
s_min=min(s_all);
s_all=linspace(s_min,s_max,n_points);
x_k_helix_all=zeros(3,n_helix,n_points);
on_x_tau_k_helix_all=zeros(3,n_helix,n_points);
on_z_tau_k_helix_all=zeros(3,n_helix,n_points);
x_k_straight_all=zeros(3,n_helix,n_points);
on_x_tau_k_straight_all=zeros(3,n_helix,n_points);
on_z_tau_k_straight_all=zeros(3,n_helix,n_points);
D0_ik_all=get_frenet_frame_from_cart(s_all,centerline_coordinates,isassigneddir,assigned_dir);
pD_ik_ps_all=get_partial_D_ik_partial_s(D0_ik_all,s_all);
u_0k_all=zeros(3,n_points);
v_0k_all=zeros(3,n_points);
LS_u_all=zeros(3,n_points);
LS_v_all=zeros(3,n_points);
r_all=centerline_coordinates;
D_ik_all=D0_ik_all;
helix_cable_strain_all=zeros(n_helix,n_points);
straight_cable_strain_all=zeros(n_straight,n_points);
grad_s_all=gradient(s_all);


for point=1:n_points
    s=s_all(point);
    ds=grad_s_all(point);
    uk=uk_all(:,point);
    vk=vk_all(:,point);
    if point>1
        Ucross_helpermat=zeros(3,3);
        for i=1:3
            for k=1:3
                Ucross_helpermat(k,i)=0;
                for j=1:3
                    Ucross_helpermat(k,i)=Ucross_helpermat(k,i)+Tensorhelpervarepsi(j,k,i)*uk(j);
                end
            end
        end
        Ucross_helpermat_trans=Ucross_helpermat';
        [u_0k_at_s,v_0k_at_s]=get_ini_darboux_strain_u0_k_and_v0_k(squeeze(D0_ik_all(:,:,point)),squeeze(pD_ik_ps_all(:,:,point)));
        u_0k_all(:,point)=u_0k_at_s;
        v_0k_all(:,point)=v_0k_at_s;
        D_ik_s1=D_ik_all(:,:,point-1);
        norm_uk=max(norm(uk),1e-10);
        helpermat_exp=eye(3)+sin(norm_uk*ds)*Ucross_helpermat_trans/norm_uk+(1-cos(norm_uk*ds))*(Ucross_helpermat_trans*Ucross_helpermat_trans)/(norm_uk^2);
        helpermat_int=(eye(3)+Ucross_helpermat_trans*Ucross_helpermat_trans/(norm_uk^2))*ds...
            -cos(norm_uk*ds)*Ucross_helpermat_trans/(norm_uk^2)...
            -sin(norm_uk*ds)*(Ucross_helpermat_trans*Ucross_helpermat_trans)/(norm_uk^3)...
            +Ucross_helpermat_trans/norm_uk^2;
        D_ik_s=D_ik_s1*helpermat_exp;
        D_ik_all(:,:,point)=D_ik_s;
        r_s=D_ik_s1*helpermat_int*vk+r_all(:,point-1);
        r_all(:,point)=r_s;
    end
end
for point=1:n_points
    on_z_D0_ik=D0_ik_all(:,:,point);
    s=s_all(point);
    ds=grad_s_all(point);
    uk=uk_all(:,point);
    vk=vk_all(:,point);
    for helix=1:n_helix
        theta=theta_helix_0+2*pi*(helix-1)'/n_helix;
        x_k_helix=get_cable_location_x_k(R_helix,omega_helix,theta,s);
        x_k_helix_all(:,helix,point)=x_k_helix;
        on_x_tau_k_helix=get_cable_direction_on_x_tau_k(R_helix,omega_helix,theta,s);
        on_x_tau_k_helix_all(:,helix,point)=on_x_tau_k_helix;
        on_z_tau_k_helix=get_cable_direction_on_z_tau_k(on_x_tau_k_helix,on_z_D0_ik);
        on_z_tau_k_helix_all(:,helix,point)=on_z_tau_k_helix;
        J=get_jacobain_J(x_k_helix,u_0k_at_s);
        on_z_D0_j3=on_z_D0_ik(:,3);
        on_z_G_alphabeta=get_contravar_on_z_G_ij(x_k_helix,u_0k_at_s,on_z_D0_ik,J);
        on_x_F_i3=get_green_deform_grad_thirdcomp_on_x_F_i3(x_k_helix,uk,vk,J);
        on_z_E_ij=get_on_z_E_ij_nonlinear(on_z_G_alphabeta,on_x_F_i3,on_z_D0_j3);
        helix_cable_strain=on_z_tau_k_helix'*on_z_E_ij*on_z_tau_k_helix;
        helix_cable_strain_all(helix,point)=helix_cable_strain;
        E_tau_tau_helix=direct_compute_cable_strain_E_tau_tau(R_helix,omega_helix,theta,s,u_0k_at_s,v_0k_at_s,uk,vk);
        helix_cable_strain_all(helix,point)=E_tau_tau_helix;
    end
    for straight=1:n_straight
        theta=theta_straight_0+2*pi*(straight-1)'/n_straight;
        x_k_straight=get_cable_location_x_k(R_straight,omega_straight,theta,s);
        x_k_straight_all(:,straight,point)=x_k_straight;
        on_x_tau_k_straight=get_cable_direction_on_x_tau_k(R_straight,omega_straight,theta,s);
        on_x_tau_k_straight_all(:,straight,point)=on_x_tau_k_straight;
        on_z_tau_k_straight=get_cable_direction_on_z_tau_k(on_x_tau_k_straight,on_z_D0_ik);
        on_z_tau_k_straight_all(:,straight,point)=on_z_tau_k_straight;
        J=get_jacobain_J(x_k_straight,u_0k_at_s);
        on_z_D0_j3=on_z_D0_ik(:,3);
        on_z_G_alphabeta=get_contravar_on_z_G_ij(x_k_straight,u_0k_at_s,on_z_D0_ik,J);
        on_x_F_i3=get_green_deform_grad_thirdcomp_on_x_F_i3(x_k_straight,uk,vk,J);
        on_z_E_ij=get_on_z_E_ij_nonlinear(on_z_G_alphabeta,on_x_F_i3,on_z_D0_j3);
        straight_cable_strain=on_z_tau_k_straight'*on_z_E_ij*on_z_tau_k_straight;
        straight_cable_strain_all(straight,point)=straight_cable_strain;
        E_tau_tau_straight=direct_compute_cable_strain_E_tau_tau(R_straight,omega_straight,theta,s,u_0k_at_s,v_0k_at_s,uk,vk);
        straight_cable_strain_all(straight,point)=E_tau_tau_straight;
    end

end
%
% parfor point=1:n_points
%     point
%     s=s_all(point);
%     [u,v]=nonlinear_least_square_solve_darboux(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,squeeze(u_0k_all(:,point)),squeeze(v_0k_all(:,point)),squeeze(helix_cable_strain_all(:,point)),squeeze(straight_cable_strain_all(:,point)));
%     LS_u_all(:,point)=u;
%     LS_v_all(:,point)=v;
% end

v_0k_all(3,1)=1;

plotinterval=1;
for helix=1:n_helix
    figure(900001)
    plot3(squeeze(x_k_helix_all(1,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(2,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(3,helix,1:plotinterval:n_points)),'r');
    grid on
    hold on
    figure(900002)
    plot3(squeeze(x_k_helix_all(1,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(2,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(3,helix,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    quiver3(squeeze(x_k_helix_all(1,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(2,helix,1:plotinterval:n_points)), squeeze(x_k_helix_all(3,helix,1:plotinterval:n_points)),squeeze(on_x_tau_k_helix_all(1,helix,1:plotinterval:n_points)), squeeze(on_x_tau_k_helix_all(2,helix,1:plotinterval:n_points)), squeeze(on_x_tau_k_helix_all(3,helix,1:plotinterval:n_points)),R_helix,'Color','r',LineWidth=R_helix)
end
for straight=1:n_straight
    figure(900001)
    plot3(squeeze(x_k_straight_all(1,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(2,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(3,straight,1:plotinterval:n_points)),'g');
    grid on
    hold on
    figure(900002)
    plot3(squeeze(x_k_straight_all(1,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(2,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(3,straight,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    quiver3(squeeze(x_k_straight_all(1,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(2,straight,1:plotinterval:n_points)), squeeze(x_k_straight_all(3,straight,1:plotinterval:n_points)),squeeze(on_x_tau_k_straight_all(1,straight,1:plotinterval:n_points)), squeeze(on_x_tau_k_straight_all(2,straight,1:plotinterval:n_points)), squeeze(on_x_tau_k_straight_all(3,straight,1:plotinterval:n_points)),R_helix,'Color','b',LineWidth=R_helix)
end
plotinterval=ceil(n_points/100);
figure(900003)
plot3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),'m');
grid on
hold on
figure(900004)
plot3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),'m');
grid on
axis equal
hold on
quiver3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),squeeze(D_ik_all(1,1,1:plotinterval:n_points))', squeeze(D_ik_all(2,1,1:plotinterval:n_points))', squeeze(D_ik_all(3,1,1:plotinterval:n_points))','r')
quiver3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),squeeze(D_ik_all(1,2,1:plotinterval:n_points))', squeeze(D_ik_all(2,2,1:plotinterval:n_points))', squeeze(D_ik_all(3,2,1:plotinterval:n_points))','g')
quiver3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),squeeze(D_ik_all(1,3,1:plotinterval:n_points))', squeeze(D_ik_all(2,3,1:plotinterval:n_points))', squeeze(D_ik_all(3,3,1:plotinterval:n_points))','b')



%
% axis equal
% hold on
% quiver3(cart_coord_x(plot_index),cart_coord_y(plot_index),cart_coord_z(plot_index),ref_rod_direct_d1o((plot_index),1),ref_rod_direct_d1o((plot_index),2),ref_rod_direct_d1o((plot_index),3),'r')
% quiver3(cart_coord_x(plot_index),cart_coord_y(plot_index),cart_coord_z(plot_index),ref_rod_direct_d2o((plot_index),1),ref_rod_direct_d2o((plot_index),2),ref_rod_direct_d2o((plot_index),3),'g')
% quiver3(cart_coord_x(plot_index),cart_coord_y(plot_index),cart_coord_z(plot_index),ref_rod_direct_d3o((plot_index),1),ref_rod_direct_d3o((plot_index),2),ref_rod_direct_d3o((plot_index),3),'b')

end