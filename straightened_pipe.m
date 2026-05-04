close all
clear
clc
load pipe_clip.mat
%% settings
%%
n_points=1953;
sample_spacing=0.1;
plot_number=200;
measure_error_level=20e-6;
n_samples=100;
n_substep=1;
strain_N=readmatrix("data\ClaremontNorth36_N.csv");
Times = readcell('data\ClaremontNorth36_N.csv','range',"1:1");
strain_N=strain_N(2:end,:);
strain_S=readmatrix("data\ClaremontNorth36_S.csv");
strain_S=strain_S(2:end,:);
strain_Z=readmatrix("data\ClaremontNorth36_Z.csv");
strain_Z=strain_Z(2:end,:);
plot_index=1:ceil(n_points/plot_number):n_points;
s = linspace(0, sample_spacing*(n_points-1), n_points)';
cart_coords=cart_coords(30:end-30,:);
cart_coord_x = cart_coords(:,1)*1e5;
cart_coord_y = cart_coords(:,2)*1e5;
cart_coord_z = cart_coords(:,3);
p=0.1;
cart_coord_x_sp = csaps(s,cart_coord_x,p,s);
cart_coord_y_sp = csaps(s,cart_coord_y,p,s);
cart_coord_z_sp = csaps(s,cart_coord_z,p,s);
%%
% s_all=[0:1:100];
s_all=s';
% rod_axis_rs=[cart_coord_x_sp';cart_coord_y_sp';cart_coord_z_sp'];
rod_axis_rs=[s_all;zeros(size(s_all));zeros(size(s_all))];
isassigneddir=true;
assigned_dir=[0,1,0];
get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
R_cables=[0.45,0.45,0.45];
theta_cables=[pi/4,pi/2,3*pi/4];
omega_cables=[0,0,0];
n_cables=3;
%%
% R_helix=R_cable;R_straight=R_cable;omega_helix=0;omega_straight=0;theta_helix_0=0;theta_straight_0=theta_cable;n_helix=0;n_straight=3;
%%
% R_helix=R_rope/sqrt(3);R_straight=R_rope*(1+2/sqrt(3));omega_helix=1;omega_straight=0;theta_helix_0=-pi/3;theta_straight_0=0;n_helix=6;n_straight=4;
% R_helix=0.01;R_straight=0.005;omega_helix=10.0;omega_straight=0;theta_helix_0=0;theta_straight_0=0;n_helix=6;n_straight=3;
% design_uk=[1;0.5;0.1];
% % design_vk=[-0.1;0.3;1.1];
% design_uk_all=0.1*[0.5;0;0.1]*ones(size(s_all))+0.2*[0;0.5;0]*sin(2*pi*s_all/s_all(end));
% design_vk_all=[0;0;1]*ones(size(s_all))+1e-4*[0;0;1]*ones(size(s_all))+1e-5*[1;0;0]*ones(size(s_all))+0*0.01*[-0.05;0.4;0.1]*ones(size(s_all))+0*[0.1;0;0.01]*sin(2*pi*s_all/s_all(end));

%% generate strains
 % [r_all,D_ik_all,helix_cable_strain_all,straight_cable_strain_all,u_0k_all,v_0k_all]=generate_MSF_kinematics(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,design_uk_all,design_vk_all);
% figure(9)
% plot(s_all,helix_cable_strain_all')
% hold on
% plot(s_all,straight_cable_strain_all')
 % check_null_space_of_cable(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all(1),u_0k_all(:,1),v_0k_all(:,1),u_0k_all(:,1),v_0k_all(:,1))

%
% readneversensvithvideo;
% helix_cable_strain_all=[];
%% 

vw_v3 = VideoWriter("straightened_Axial elongation","MPEG-4");
vw_v3.FrameRate=1/0.6;
open(vw_v3)
vw = VideoWriter("straightened_pipe","MPEG-4");
vw.FrameRate=1/0.6;
open(vw)
k_frame_inter=1;
n_points=length(s_all);
trajectoryendpoint=zeros(3,length(Times));
for k = 1:k_frame_inter:length(Times)
    slicetimeid=k
    cables_strain_all=1*[strain_N(:,slicetimeid)';strain_S(:,slicetimeid)';strain_Z(:,slicetimeid)']*1e-6;
    plot_figure_id=0;
    [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep);
    figure(1001)
    plotinterval=ceil(n_points/100);
    plot3(squeeze(r_all_NFCA(1,1:plotinterval:n_points)), squeeze(r_all_NFCA(2,1:plotinterval:n_points)), squeeze(r_all_NFCA(3,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    quiver3(squeeze(r_all_NFCA(1,1:plotinterval:n_points)), squeeze(r_all_NFCA(2,1:plotinterval:n_points)), squeeze(r_all_NFCA(3,1:plotinterval:n_points)),squeeze(D_ik_all_NFCA(1,1,1:plotinterval:n_points))', squeeze(D_ik_all_NFCA(2,1,1:plotinterval:n_points))', squeeze(D_ik_all_NFCA(3,1,1:plotinterval:n_points))','r')
    quiver3(squeeze(r_all_NFCA(1,1:plotinterval:n_points)), squeeze(r_all_NFCA(2,1:plotinterval:n_points)), squeeze(r_all_NFCA(3,1:plotinterval:n_points)),squeeze(D_ik_all_NFCA(1,2,1:plotinterval:n_points))', squeeze(D_ik_all_NFCA(2,2,1:plotinterval:n_points))', squeeze(D_ik_all_NFCA(3,2,1:plotinterval:n_points))','g')
    quiver3(squeeze(r_all_NFCA(1,1:plotinterval:n_points)), squeeze(r_all_NFCA(2,1:plotinterval:n_points)), squeeze(r_all_NFCA(3,1:plotinterval:n_points)),squeeze(D_ik_all_NFCA(1,3,1:plotinterval:n_points))', squeeze(D_ik_all_NFCA(2,3,1:plotinterval:n_points))', squeeze(D_ik_all_NFCA(3,3,1:plotinterval:n_points))','b')
    hold off
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    xlim([0,220])
    ylim([-10,10])
    zlim([-10,10])

    view(488.2878,37)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
	title(caption, 'FontSize', 15);
    frame = getframe(gcf);
    writeVideo(vw,frame)
    trajectoryendpoint(:,k)=r_all_NFCA(:,end);

    figure(1004)
    col=v_k_all_NFCA(3,1:plotinterval:n_points)-1;
    clf(1004)
    surface([squeeze(r_all_NFCA(1,1:plotinterval:n_points));squeeze(r_all_NFCA(1,1:plotinterval:n_points))],...
        [squeeze(r_all_NFCA(2,1:plotinterval:n_points));squeeze(r_all_NFCA(2,1:plotinterval:n_points))],...
        [squeeze(r_all_NFCA(3,1:plotinterval:n_points));squeeze(r_all_NFCA(3,1:plotinterval:n_points))],...
        [col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5);
    colorbar;
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    xlim([0,220])
    ylim([-10,10])
    zlim([-10,10])
    clim([-2e-4,2e-4])
    view(519.2878,61)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title({"Axial elongation [1]";caption}, 'FontSize', 15);
        frame = getframe(gcf);
    writeVideo(vw_v3,frame)
end
close(vw_v3)
close(vw)
trajectoryendpointini=trajectoryendpoint(:,1);
for i=1:3
    for t=1:length(Times)
    trajectoryendpoint(i,t)=trajectoryendpoint(i,t)-trajectoryendpointini(i);
    end
end
vw = VideoWriter("straightened_trajectory","MPEG-4");
vw.FrameRate=1/0.6;
open(vw)
k_frame_inter=1;
n_points=length(s_all);
for k = 1:k_frame_inter:length(Times)
    slicetimeid=k
    figure(1002)
    plot3(trajectoryendpoint(1,1:slicetimeid),trajectoryendpoint(2,1:slicetimeid),trajectoryendpoint(3,1:slicetimeid),LineWidth=3.0,Color='r')
    grid on
    axis equal
    hold on
    hold off
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    xlim([-2,2])
    ylim([-2,2])
    zlim([-2,2])
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title({"Trajectory of the Selected End Point";caption}, 'FontSize', 15);
    view(488.2878-5*slicetimeid,37+0*slicetimeid)
    frame = getframe(gcf);
    writeVideo(vw,frame)

end
close(vw)



%% compute the noise-free Euler–Bernoulli approximation
% plot_figure_id=1001;
% [r_all_NEBA,D_ik_all_NEBA,D0_ik_all_NEBA,u0_k_all_NEBA,v0_k_all_NEBA,u_k_all_NEBA,v_k_all_NEBA]=EBbeam_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%% compute the Deterministic Cosserat analytical solution
% plot_figure_id=1011;
% [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%% compute the Deterministic Cosserat simulation
% plot_figure_id=1021;
% [r_all_NFCS,D_ik_all_NFCS,D0_ik_all_NFCS,u0_k_all_NFCS,v0_k_all_NFCS,u_k_all_NFCS,v_k_all_NFCS]=Cosserat_simulation(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%% compute the stochastic Cosserat analytical solution
% plot_figure_id=1031;
% [r_all_SACS,D_ik_all_SACS,D0_ik_all_SACS,u0_k_all_SACS,v0_k_all_SACS,u_k_all_SACS,v_k_all_SACS,Cu_k_all_SACS,Cv_k_all_SACS,MSE_all_SACS]=Cosserat_analytical_SDE(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%% Numerical Stochastic Cosserat theory with Milstein scheme
% plot_figure_id=1041;
% [r_all_NSCM,D_ik_all_NSCM,D0_ik_all_NSCM,u0_k_all_NSCM,v0_k_all_NSCM,u_k_all_NSCM,v_k_all_NSCM,Cu_k_all_NSCM,Cv_k_all_NSCM,MSE_all_NSCM]...
%     =Cosserat_Milstein_Samples(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
%% Numerical Stochastic Cosserat theory with RK scheme
% plot_figure_id=1051;
% [r_all_NSCRK,D_ik_all_NNSCRK,D0_ik_all_NSCRK,u0_k_all_NSCRK,v0_k_all_NSCRK,u_k_all_NSCRK,v_k_all_NSCRK,Cu_k_all_NSCRK,Cv_k_all_NSCRK,MSE_all_NSCRK]...
%     =Cosserat_WeakRK_Samples(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
%% Monte-carlo simulation EB
% plot_figure_id=1061;
% [r_all_MCEBS,D_ik_all_MCEBS,D0_ik_all_MCEBS,u0_k_all_MCEBS,v0_k_all_MCEBS,u_k_all_MCEBS,v_k_all_MCEBS,Cu_k_all_MCEBS,Cv_k_all_MCEBS,MSE_all_MCEBS]...
%     =EB_analytical_Montecarlo(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
%% Monte-carlo simulation Cosserat
% plot_figure_id=1071;
% [r_all_MCCS,D_ik_all_MCCS,D0_ik_all_MCCS,u0_k_all_MCCS,v0_k_all_MCCS,u_k_all_MCCS,v_k_all_MCCS,Cu_k_all_MCCS,Cv_k_all_MCCS,MSE_all_MCCS]...
%     =Cosserat_Simulation_Montecarlo(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
% %% Monte-carlo simulation fully nonlinear Cosserat
% plot_figure_id=1081;
% [r_all_MCFNCS,D_ik_all_MCFNCS,D0_ik_all_MCFNCS,u0_k_all_MCFNCS,v0_k_all_MCFNCS,u_k_all_MCFNCS,v_k_all_MCFNCS,Cu_k_all_MCFNCS,Cv_k_all_MCFNCS,MSE_all_MCFNCS]...
%     =Fully_nonlin_Cosserat_Simulation_Montecarlo(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
