close all
clear
clc
%% settings
%%
s_all=[0:0.0026:4];
%%
% s_all=[0:1:100];
rod_axis_rs=[zeros(size(s_all));zeros(size(s_all));s_all];
isassigneddir=true;
assigned_dir=[1,0,0];
get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
R_cable=sqrt((5.031/2)^2+(4.745/2)^2)*1e-3;
theta_cable=atan(5.031/4.745);
%%
R_helix=R_cable;R_straight=R_cable;omega_helix=0;omega_straight=0;theta_helix_0=-theta_cable-pi/2;theta_straight_0=theta_cable-3*pi/2;n_helix=2;n_straight=2;
%%
% R_helix=R_rope/sqrt(3);R_straight=R_rope*(1+2/sqrt(3));omega_helix=1;omega_straight=0;theta_helix_0=-pi/3;theta_straight_0=0;n_helix=6;n_straight=4;
% R_helix=0.01;R_straight=0.005;omega_helix=10.0;omega_straight=0;theta_helix_0=0;theta_straight_0=0;n_helix=6;n_straight=3;
% design_uk=[1;0.5;0.1];
% % design_vk=[-0.1;0.3;1.1];
design_uk_all=0.1*[0.5;0;0.1]*ones(size(s_all))+0.2*[0;0.5;0]*sin(2*pi*s_all/s_all(end));
design_vk_all=[0;0;1]*ones(size(s_all))+1e-4*[0;0;1]*ones(size(s_all))+1e-5*[1;0;0]*ones(size(s_all))+0*0.01*[-0.05;0.4;0.1]*ones(size(s_all))+0*[0.1;0;0.01]*sin(2*pi*s_all/s_all(end));
measure_error_level=20e-6;
n_samples=100;
n_substep=1;
%% generate strains
[r_all,D_ik_all,helix_cable_strain_all,straight_cable_strain_all,u_0k_all,v_0k_all]=generate_MSF_kinematics(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,design_uk_all,design_vk_all);
% figure(9)
% plot(s_all,helix_cable_strain_all')
% hold on
% plot(s_all,straight_cable_strain_all')
check_null_space_of_cable(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all(1),u_0k_all(:,1),v_0k_all(:,1),u_0k_all(:,1),v_0k_all(:,1))

%
readneversens;
%% 
vw = VideoWriter("neversens.avi");
open(vw)
k_frame_inter=1;
n_points=length(s_all);
for k = 1:k_frame_inter:length(time4)
    slicetimeid=k
    dat1=dat1-dat1(2,:);
    dat2=dat2-dat2(2,:);
    dat3=dat3-dat3(2,:);
    dat4=dat4-dat4(2,:);
    reading_strainu1=dat1(slicetimeid,:);
    loc_strainu1=loc1-loc1(1);
    reading_strainu2=dat2(slicetimeid,:);
    loc_strainu2=loc2-loc2(1);
    reading_strainl1=dat3(slicetimeid,:);
    loc_strainl1=loc3-loc3(1);
    reading_strainl2=dat4(slicetimeid,:);
    loc_strainl2=loc4-loc4(1);
    strainu1=interp1(loc_strainu1,reading_strainu1,s_all*loc_strainu1(end)/s_all(end));
    strainu2=interp1(loc_strainu2,reading_strainu2,s_all*loc_strainu2(end)/s_all(end));
    strainl1=interp1(loc_strainl1,reading_strainl1,s_all*loc_strainl1(end)/s_all(end));
    strainl2=interp1(loc_strainl2,reading_strainl2,s_all*loc_strainl2(end)/s_all(end));
    nans=isnan(strainl2);
    strainl2(nans) = interp1(s_all(~nans), strainl2(~nans), s_all(nans));
    nans=isnan(strainu1);
    strainu1(nans) = interp1(s_all(~nans), strainu1(~nans), s_all(nans));
    nans=isnan(strainu2);
    strainu2(nans) = interp1(s_all(~nans), strainu2(~nans), s_all(nans));
    nans=isnan(strainl1);
    strainl1(nans) = interp1(s_all(~nans), strainl1(~nans), s_all(nans));
    helix_cable_strain_all=[strainu1;strainl2]*1e-6;
    straight_cable_strain_all=[strainl1;strainu2]*1e-6;
    plot_figure_id=0;
    [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
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
    xlim([-1,1])
    ylim([-1,1])
    zlim([0,4.5])
    view(1.642537145473441e+02,10.25614035087716)
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
