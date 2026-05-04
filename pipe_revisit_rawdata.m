close all
clear
clc
% load pipe_clip.mat
%% settings
%%
scaleeff=1;
n_points=1953;
sample_spacing=0.1;
plot_number=200;
measure_error_level=scaleeff*20e-6;
scaleeff=1;
n_samples=100;
n_substep=1;
cart_coords=readmatrix("claremont_data/pipe_coords_xyz_36pipe.csv");
% strain_N=readmatrix("data/newdata/ClaremontNorth36_N_tc.csv");
strain_N=readmatrix("claremont_data/ClaremontNorth36_N_measured_notc.csv");
Times = readcell('claremont_data/ClaremontNorth36_N_measured_notc.csv','range',"1:1");
strain_N=strain_N(2:end,:);
strain_S=readmatrix("claremont_data/ClaremontNorth36_S_measured_notc.csv");
strain_S=strain_S(2:end,:);
strain_Z=readmatrix("claremont_data/ClaremontNorth36_Z_measured_notc.csv");
strain_Z=strain_Z(2:end,:);
plot_index=1:ceil(n_points/plot_number):n_points;
s = linspace(0, sample_spacing*(n_points-1), n_points)';
% cart_coords=cart_coords(30:end-30,:);

cart_coord_x = cart_coords(:,1)-cart_coords(1,1);
cart_coord_y = cart_coords(:,2)-cart_coords(1,2);
cart_coord_z = cart_coords(:,3)-cart_coords(1,3);
% cart_coord_x=flip(cart_coord_x);
% % cart_coord_x=-cart_coord_x;
% cart_coord_y=flip(cart_coord_y);
% cart_coord_z=flip(cart_coord_z);
% strain_N=flip(strain_N);
% strain_S=flip(strain_S);
% strain_Z=flip(strain_Z);
p=0.1;
cart_coord_x_sp = csaps(s,cart_coord_x,p,s);
cart_coord_y_sp = csaps(s,cart_coord_y,p,s);
cart_coord_z_sp = csaps(s,cart_coord_z,p,s);
%%
% s_all=[0:1:100];
s_all=s';
rod_axis_rs=[cart_coord_x_sp';cart_coord_y_sp';cart_coord_z_sp'];
isassigneddir=false;
assigned_dir=[1,0,0];
get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
R_cables=[0.45,0.45,0.45];
%NSZ
theta_cables=-[pi/4,pi/2,3*pi/4];
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
    figure(1004)
    hold on
       figure(1001)
    hold on 
        figure(1005)
    hold on
    figure(1056)
    hold on
vw_v3 = VideoWriter("Axial elongation","Uncompressed AVI");
vw_v3.FrameRate=1/0.6;
open(vw_v3)
vw = VideoWriter("pipe","Uncompressed AVI");
vw.FrameRate=1/0.6;
open(vw)
k_frame_inter=1;
n_points=length(s_all);
trajectoryendpoint=zeros(3,length(Times));
trajectorypoint1=zeros(3,length(Times));
trajectorypoint2=zeros(3,length(Times));
for k = 1:k_frame_inter:length(Times)
    slicetimeid=k
    cables_strain_all=scaleeff*[strain_N(:,slicetimeid)';strain_S(:,slicetimeid)';strain_Z(:,slicetimeid)']*1e-6;
    plot_figure_id=0;
    % [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep);
    % [r_all_NFCA,D_ik_all_NFCA]=rotate_in_plane_rod(r_all_NFCA,D_ik_all_NFCA,rod_axis_rs);
    [r_all_NSCRK,D_ik_all_NSCRK,D0_ik_all_NSCRK,u0_k_all_NSCRK,v0_k_all_NSCRK,u_k_all_NSCRK,v_k_all_NSCRK,Cu_k_all_NSCRK,Cv_k_all_NSCRK,MSE_all_NSCRK,s_r_all_NSCRK,s_D_ik_all_NSCRK]...
        =Cosserat_WeakRK_Samples_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
    [r_all_NSCRK,D_ik_all_NSCRK]=rotate_in_plane_rod(r_all_NSCRK,D_ik_all_NSCRK,rod_axis_rs);
    figure(1001)
    plotinterval=ceil(n_points/100);
    plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,1,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,1,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,1,1:plotinterval:n_points))','r')
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,2,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,2,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,2,1:plotinterval:n_points))','g')
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,3,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,3,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,3,1:plotinterval:n_points))','b')
    hold off
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    % xlim([-200,0])
    % ylim([-20,60])
    % zlim([95,140])
    % view(-12.283639988412114,39.810373502729341)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
	title(caption, 'FontSize', 15);
    print(['claremontquiver' caption '.png'],'-dpng')
    frame = getframe(gcf);
    writeVideo(vw,frame)
    trajectoryendpoint(:,k)=r_all_NSCRK(:,end);
    trajectorypoint1(:,k)=r_all_NSCRK(:,467);
    trajectorypoint2(:,k)=r_all_NSCRK(:,1334);
    figure(1004)
    col=v_k_all_NSCRK(3,1:plotinterval:n_points)-1;
    clf(1004)
    surface([squeeze(r_all_NSCRK(1,1:plotinterval:n_points));squeeze(r_all_NSCRK(1,1:plotinterval:n_points))],...
        [squeeze(r_all_NSCRK(2,1:plotinterval:n_points));squeeze(r_all_NSCRK(2,1:plotinterval:n_points))],...
        [squeeze(r_all_NSCRK(3,1:plotinterval:n_points));squeeze(r_all_NSCRK(3,1:plotinterval:n_points))],...
        [col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5);
    colorbar;
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    % xlim([-200,0])
    % ylim([-20,60])
    % zlim([95,140])
    % clim([-2e-4,2e-4])
    % view(540,-90)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title({"Axial elongation [1]";caption}, 'FontSize', 15);
        frame = getframe(gcf);
    writeVideo(vw_v3,frame)
    %%
    if mod(k,8)==1
    figure(1005)
    hold on
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    currenttime.Format = 'yyyy-MM';
    caption = sprintf('Time = %s', string(currenttime));
    captionshort=sprintf('%s', string(currenttime));
    thisplot=plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'linewidth',3, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    figure(1056)
    hold on
    thisplot=plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'Color',thiscolor,'linewidth',3, 'DisplayName', captionshort);
    % thiscolor=get(thisplot, 'Color');
    for isample=1:n_samples
        r_all_isample=squeeze(s_r_all_NSCRK(isample,:,:));
        D_ik_all_isample=squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
        [r_all_isample,D_ik_all_isample]=rotate_in_plane_rod(r_all_isample,D_ik_all_isample,rod_axis_rs);
        plot3(squeeze(r_all_isample(1,1:plotinterval:n_points)), squeeze(r_all_isample(2,1:plotinterval:n_points)), squeeze(r_all_isample(3,1:plotinterval:n_points)),'Color',[thiscolor,0.01],'linewidth',1,'HandleVisibility','off');
    end
    end
end
close(vw_v3)
close(vw)
trajectoryendpointini=trajectoryendpoint(:,1);
for i=1:3
    trajectoryendpoint(i,:)=trajectoryendpoint(i,:)-trajectoryendpoint(i,1);
    trajectorypoint1(i,:)=trajectorypoint1(i,:)-trajectorypoint1(i,1);
    trajectorypoint2(i,:)=trajectorypoint2(i,:)-trajectorypoint2(i,1);
end
%%
vw = VideoWriter("trajectory","Uncompressed AVI");
vw.FrameRate=1/0.6;
open(vw)
    figure(1002)
    hold on
k_frame_inter=1;
n_points=length(s_all);
for k = 1:k_frame_inter:length(Times)
    slicetimeid=k
    figure(1002)
    hold on
    plot3(trajectoryendpoint(1,1:slicetimeid),trajectoryendpoint(2,1:slicetimeid),trajectoryendpoint(3,1:slicetimeid),LineWidth=3.0,Color='r')
    plot3(trajectorypoint1(1,1:slicetimeid),trajectorypoint1(2,1:slicetimeid),trajectorypoint1(3,1:slicetimeid),LineWidth=3.0,Color='g')
    plot3(trajectorypoint2(1,1:slicetimeid),trajectorypoint2(2,1:slicetimeid),trajectorypoint2(3,1:slicetimeid),LineWidth=3.0,Color='b')
    grid on
    axis equal
    hold on
    hold off
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    % xlim([-2,2])
    % ylim([-2,2])
    % zlim([-2,2])
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title({"Trajectory of the Selected Points";caption}, 'FontSize', 15);
    legend("end","point 1","point 2")
    % view(488.2878-5*slicetimeid,37+0*slicetimeid)
    frame = getframe(gcf);
    writeVideo(vw,frame)

end
close(vw)
%%
vw = VideoWriter("relativetrajectory","Uncompressed AVI");
vw.FrameRate=1/0.6;
open(vw)
    figure(1012)
    hold on
k_frame_inter=1;
n_points=length(s_all);
for k = 1:k_frame_inter:length(Times)
    slicetimeid=k
    figure(1012)
    hold on
    relativetrajectory=trajectorypoint1-trajectorypoint2;
    plot3(relativetrajectory(1,1:slicetimeid),relativetrajectory(2,1:slicetimeid),relativetrajectory(3,1:slicetimeid),LineWidth=3.0,Color='r')
    grid on
    axis equal
    hold on
    hold off
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    % xlim([-2,2])
    % ylim([-2,2])
    % zlim([-2,2])
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title({"Trajectory of the relative displacment 1-2";caption}, 'FontSize', 15);
    view(488.2878-5*slicetimeid,37+0*slicetimeid)
    frame = getframe(gcf);
    writeVideo(vw,frame)

end
close(vw)
%%
positionvec=[10         478        1917         518];
figure(1005)
set(gcf,"Position",positionvec)
hold on
xlabel("$x$ [mm]",'interpreter',"latex",'fontsize',14)
ylabel("$y$ [mm]",'interpreter',"latex",'fontsize',14)
zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
legend('Location','eastoutside')

% axis equal
    % xlim([-200,0])
    % ylim([-20,60])
    % zlim([95,140])
    % view(-12.283639988412114,39.810373502729341)
grid on
title({"3D deformation"}, 'FontSize', 15);
set(findall(gcf,'-property','FontSize'),'FontSize',15)
print('claremont3Dstaticplot.png','-dpng')
view(0,0)
title({"2D deformation plane xz"}, 'FontSize', 15);
print('claremont2Dstaticplot xz.png','-dpng')
% ylim(1e3*[-3e-4,3e-4])
print('claremont2Dstaticplot yz.png','-dpng')
view(0,90)
title({"2D deformation plane xy"}, 'FontSize', 15);
% set(gcf,"Position",[965,703,630,472])
print('claremont2Dstaticplot xy.png','-dpng')
view(-90,0)
set(gcf,"Position",[255,604,744,472])
title({"2D deformation-MC-plane yz"}, 'FontSize', 15);
%%
figure(1056)
set(gcf,"Position",positionvec)
hold on
xlabel("$x$ [mm]",'interpreter',"latex",'fontsize',14)
ylabel("$y$ [mm]",'interpreter',"latex",'fontsize',14)
zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
legend('Location','eastoutside')

    % xlim([-200,0])
    % ylim([-20,60])
    % zlim([95,140])
    % view(-12.283639988412114,39.810373502729341)
grid on
title({"3D deformation-MC"}, 'FontSize', 15);
set(findall(gcf,'-property','FontSize'),'FontSize',15)
print('claremont3DstaticplotMC.png','-dpng')
view(0,0)
title({"2D deformation-MC-plane xz"}, 'FontSize', 15);
print('claremont2DstaticplotMC xz.png','-dpng')

% ylim(1e3*[-3e-4,3e-4])
print('claremont2DstaticplotMC yz.png','-dpng')
view(0,90)
title({"2D deformation-MC-plane xy"}, 'FontSize', 15);
% 
print('claremont2DstaticplotMC xy.png','-dpng')
view(-90,0)
set(gcf,"Position",[255,604,744,472])
title({"2D deformation-MC-plane yz"}, 'FontSize', 15);
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
