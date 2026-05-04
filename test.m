close all
clear
clc
%% settings
s_all=[0:0.1:10];
rod_axis_rs=[s_all;zeros(size(s_all));zeros(size(s_all))];
isassigneddir=true;
assigned_dir=[0,1,0];
get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
R_rope=0.01;
% R_helix=R_rope/sqrt(3);R_straight=R_rope*(1+2/sqrt(3));omega_helix=1.0;omega_straight=0;theta_helix_0=0;theta_straight_0=pi/2;n_helix=3;n_straight=3;
% R_helix=R_rope/sqrt(3);R_straight=R_rope*(1+2/sqrt(3));omega_helix=1.0;omega_straight=0;theta_helix_0=0;theta_straight_0=pi/2;n_helix=6;n_straight=6;
R_helix=R_rope;R_straight=R_rope;omega_helix=7.5157;omega_straight=0;theta_helix_0=0;theta_straight_0=pi/2;n_helix=3;n_straight=4;
% R_helix=0.05;R_straight=0.05;omega_helix=5.0;omega_straight=0;theta_helix_0=0;theta_straight_0=pi/2;n_helix=3;n_straight=2;
% design_uk=[1;0.5;0.1];
% design_vk=[-0.1;0.3;1.1];
design_uk_all=[1;2;0.1]*ones(size(s_all))+0*[0.1;-0.2;0.1]*sin(2*pi*s_all/10)+[0;0;0]*(s_all/max(s_all)<0.5)+[0;0;0]*(s_all/max(s_all)>0.5);
design_vk_all=[0;0;1]*ones(size(s_all))+0*[-0;0.5;1.0]*ones(size(s_all))+0*[-0.05;0.4;0.1]*ones(size(s_all))+0*[0.1;0;0.01]*sin(2*pi*s_all/10);
measure_error_level=10e-6;
n_samples=100;
%% generate strains
[r_all,D_ik_all,helix_cable_strain_all,straight_cable_strain_all,u_0k_all,v_0k_all]=generate_MSF_kinematics(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,design_uk_all,design_vk_all);
figure(9)
% plot(s_all,helix_cable_strain_all')
hold on
plot(s_all,straight_cable_strain_all')
% for point=1:1
%     sp=s_all(point);
%     u0_k=u_0k_all(:,point);
%     v0_k=v_0k_all(:,point);
%     epsilon_p_helix=helix_cable_strain_all(:,point);
%     epsilon_p_straight=straight_cable_strain_all(:,point);
%     [A_pq,B_pq_l]=nonlin_rod_assymbling_matrixs(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,sp,u0_k,v0_k,epsilon_p_helix,epsilon_p_straight);
% end


% % %%
% plot_figure_id=1051;
% measure_error_level=10e-6;
% n_substep=1;
% n_samples=100;
% % [r_all_calc,D_ik_all_calc,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all]=Cosserat_Milstein(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%  % [r_all_calc,D_ik_all_calc,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all]=Cosserat_WeakRK(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
% % [r_all_calc,D_ik_all_calc,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all,MSE_all]=Cosserat_analytical(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
% [r_all_calc,D_ik_all_calc,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all,MSE_all]=Cosserat_noisefree(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
% % [r_all_calc,D_ik_all_calc,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all,MSE_all]=EBbeam_noisefree(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
% [s_r_all,s_D_ik_all,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all,MSE_all]...
%     =Cosserat_Milstein_Samples(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);

% 
% 
% %%
% 
% 
% 
% figure(101)
% plot(s_all,straight_cable_strain_all')
% hold on
% plot(s_all,helix_cable_strain_all')
% figure(102)
% plot(s_all,MSE_all)
% % %%
% syms R omega theta s
% syms u0_k v0_k uk vk [3,1]
% E_tau_tau=direct_compute_cable_strain_E_tau_tau(R,omega,theta,s,u0_k,v0_k,uk,vk);

%
%% compute the noise-free Euler–Bernoulli approximation
plot_figure_id=1011;
n_substep=1;
[r_all_NEBA,D_ik_all_NEBA,D0_ik_all_NEBA,u0_k_all_NEBA,v0_k_all_NEBA,u_k_all_NEBA,v_k_all_NEBA]=EBbeam_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%% compute the noise-free Cosserat approximation
plot_figure_id=1001;
n_substep=1;
[r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%% compute the Deterministic Cosserat analytical solution
plot_figure_id=1021;
n_substep=1;
[r_all_DACS,D_ik_all_DACS,D0_ik_all_DACS,u0_k_all_DACS,v0_k_all_DACS,u_k_all_DACS,v_k_all_DACS]=Cosserat_analytical(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%% compute the stochastic Cosserat analytical solution
plot_figure_id=1031;
n_substep=1;
[r_all_SACS,D_ik_all_SACS,D0_ik_all_SACS,u0_k_all_SACS,v0_k_all_SACS,u_k_all_SACS,v_k_all_SACS,Cu_k_all_SACS,Cv_k_all_SACS,MSE_all_SACS]=Cosserat_analytical_SDE(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%% Numerical Stochastic Cosserat theory with Milstein scheme
plot_figure_id=1041;
n_substep=1;
[r_all_NSCM,D_ik_all_NSCM,D0_ik_all_NSCM,u0_k_all_NSCM,v0_k_all_NSCM,u_k_all_NSCM,v_k_all_NSCM,Cu_k_all_NSCM,Cv_k_all_NSCM,MSE_all_NSCM]...
    =Cosserat_Milstein_Samples(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
%% Numerical Stochastic Cosserat theory with RK scheme
plot_figure_id=1051;
n_substep=1;
[r_all_NSCRK,D_ik_all_NNSCRK,D0_ik_all_NSCRK,u0_k_all_NSCRK,v0_k_all_NSCRK,u_k_all_NSCRK,v_k_all_NSCRK,Cu_k_all_NSCRK,Cv_k_all_NSCRK,MSE_all_NSCRK]...
    =Cosserat_WeakRK_Samples(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
%% Monte-carlo simulation EB approximation
plot_figure_id=1061;
n_substep=1;
[r_all_MCEBA,D_ik_all_MCEBA,D0_ik_all_MCEBA,u0_k_all_MCEBA,v0_k_all_MCEBA,u_k_all_MCEBA,v_k_all_MCEBA,Cu_k_all_MCEBA,Cv_k_all_MCEBA,MSE_all_MCEBA]...
    =EB_analytical_Montecarlo(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
%% Monte-carlo simulation Cosserat analytical
plot_figure_id=1071;
n_substep=1;
[r_all_MCCS,D_ik_all_MCCS,D0_ik_all_MCCS,u0_k_all_MCCS,v0_k_all_MCCS,u_k_all_MCCS,v_k_all_MCCS,Cu_k_all_MCCS,Cv_k_all_MCCS,MSE_all_MCCS]...
    =Cosserat_analytical_Montecarlo(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
%% Save
name='data'
name1 = [strrep(datestr(now), ':', '-')];
mkdir(name1);
save([name1 '\' name])
figureids=[1011:10:1071];
figureids=[figureids figureids+1];
for fid=figureids
    figure(fid);
    savefig([name1 '\' num2str(fid)]);
end