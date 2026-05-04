close all
clear
clc
%% settings
%%
s_all=[0:0.0026:4.5682];
%%
% s_all=[0:1:100];
rod_axis_rs=[zeros(size(s_all));zeros(size(s_all));s_all];
isassigneddir=true;
assigned_dir=[1,0,0];
get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
R_rope=8.84e-3;
%%
R_helix=R_rope/sqrt(3);R_straight=R_rope*(1+2/sqrt(3));omega_helix=7.5157;omega_straight=0.073;theta_helix_0=-pi/3;theta_straight_0=0;n_helix=3;n_straight=1;
%%
% R_helix=R_rope/sqrt(3);R_straight=R_rope*(1+2/sqrt(3));omega_helix=1;omega_straight=0;theta_helix_0=-pi/3;theta_straight_0=0;n_helix=6;n_straight=4;
% R_helix=0.01;R_straight=0.005;omega_helix=10.0;omega_straight=0;theta_helix_0=0;theta_straight_0=0;n_helix=6;n_straight=3;
% design_uk=[1;0.5;0.1];
% % design_vk=[-0.1;0.3;1.1];
design_uk_all=0.1*[0.5;0;0.1]*ones(size(s_all))+0.2*[0;0.5;0]*sin(2*pi*s_all/s_all(end));
design_vk_all=[0;0;1]*ones(size(s_all))+1e-4*[0;0;1]*ones(size(s_all))+1e-5*[1;0;0]*ones(size(s_all))+0*0.01*[-0.05;0.4;0.1]*ones(size(s_all))+0*[0.1;0;0.01]*sin(2*pi*s_all/s_all(end));
measure_error_level=20e-6;
n_samples=1000;
n_substep=1;
%% generate strains
[r_all,D_ik_all,helix_cable_strain_all,straight_cable_strain_all,u_0k_all,v_0k_all]=generate_MSF_kinematics(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,design_uk_all,design_vk_all);
% figure(9)
% plot(s_all,helix_cable_strain_all')
% hold on
% plot(s_all,straight_cable_strain_all')
check_null_space_of_cable(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all(1),u_0k_all(:,1),v_0k_all(:,1),u_0k_all(:,1),v_0k_all(:,1))

%
% readlunafile;
% interpstrain;
load('ropephotodata/shapedata.mat')
shapei=3;
strainhelix1=strainhelix1s(shapei,:);
strainhelix2=strainhelix2s(shapei,:);
strainhelix3=strainhelix1s(shapei,:);
strainstraight0=strainhelix0s(shapei,:);
helix_cable_strain_all=1*[strainhelix1;strainhelix2;strainhelix3]*1e-6;
straight_cable_strain_all=1*strainstraight0*1e-6;
% filtercoeff=1e-3;
% straight_cable_strain_all=lowpass(straight_cable_strain_all',filtercoeff)';
% helix_cable_strain_all=lowpass(helix_cable_strain_all',filtercoeff)';
% figure(10)
% plot(s_all,helix_cable_strain_all')
% hold on
% plot(s_all,straight_cable_strain_all')
%% compute the stochastic Cosserat analytical solution
% plot_figure_id=1031;
% n_substep=1;
% [r_all_SACS,D_ik_all_SACS,D0_ik_all_SACS,u0_k_all_SACS,v0_k_all_SACS,u_k_all_SACS,v_k_all_SACS,Cu_k_all_SACS,Cv_k_all_SACS,MSE_all_SACS]=Cosserat_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
% % % figure(14)
% plot(u_k_all_SACS')
% % filtercoeff=1.0;
% % design_uk_all=lowpass(u_k_all_SACS',filtercoeff)';
% design_uk_all=u_k_all_SACS;
% figure(13)
% plot(design_uk_all')
% design_vk_all=v_k_all_SACS;
% [r_all,D_ik_all,helix_cable_strain_all_filtered,straight_cable_strain_all_filtered,u_0k_all,v_0k_all]=generate_MSF_kinematics(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,design_uk_all,design_vk_all);
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

% [r_all,D_ik_all,helix_cable_strain_all_NFCA,straight_cable_strain_all_NFCA,u_0k_all,v_0k_all]=generate_MSF_kinematics(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,u_k_all_NFCA,v_k_all_NFCA);
% figure(11)
% plot(helix_cable_strain_all_NFCA')
% hold on
% plot(straight_cable_strain_all_NFCA')


% %% compute the noise-free Cosserat approximationMCEBA
% plot_figure_id=1001;
% n_substep=1;
% [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
% % compute the noise-free Euler–Bernoulli approximation
% plot_figure_id=1001;
% [r_all_NEBA,D_ik_all_NEBA,D0_ik_all_NEBA,u0_k_all_NEBA,v0_k_all_NEBA,u_k_all_NEBA,v_k_all_NEBA]=EBbeam_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
% %% compute the Deterministic Cosserat analytical solution
% plot_figure_id=1011;
% [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
% %% compute the Deterministic Cosserat simulation
% plot_figure_id=1021;
% [r_all_NFCS,D_ik_all_NFCS,D0_ik_all_NFCS,u0_k_all_NFCS,v0_k_all_NFCS,u_k_all_NFCS,v_k_all_NFCS]=Cosserat_simulation(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
%% compute the stochastic Cosserat analytical solution
plot_figure_id=1031;
[r_all_SACS,D_ik_all_SACS,D0_ik_all_SACS,u0_k_all_SACS,v0_k_all_SACS,u_k_all_SACS,v_k_all_SACS,Cu_k_all_SACS,Cv_k_all_SACS,MSE_all_SACS]=Cosserat_analytical_SDE(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
% %% Numerical Stochastic Cosserat theory with Milstein scheme
% plot_figure_id=1041;
% [r_all_NSCM,D_ik_all_NSCM,D0_ik_all_NSCM,u0_k_all_NSCM,v0_k_all_NSCM,u_k_all_NSCM,v_k_all_NSCM,Cu_k_all_NSCM,Cv_k_all_NSCM,MSE_all_NSCM]...
%     =Cosserat_Milstein_Samples(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
% %% Numerical Stochastic Cosserat theory with RK scheme
% plot_figure_id=1051;
% [r_all_NSCRK,D_ik_all_NNSCRK,D0_ik_all_NSCRK,u0_k_all_NSCRK,v0_k_all_NSCRK,u_k_all_NSCRK,v_k_all_NSCRK,Cu_k_all_NSCRK,Cv_k_all_NSCRK,MSE_all_NSCRK]...
%     =Cosserat_WeakRK_Samples(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
% %% Monte-carlo simulation EB
% plot_figure_id=1061;
% [r_all_MCEBS,D_ik_all_MCEBS,D0_ik_all_MCEBS,u0_k_all_MCEBS,v0_k_all_MCEBS,u_k_all_MCEBS,v_k_all_MCEBS,Cu_k_all_MCEBS,Cv_k_all_MCEBS,MSE_all_MCEBS]...
%     =EB_analytical_Montecarlo(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
% %% Monte-carlo simulation Cosserat
% plot_figure_id=1071;
% [r_all_MCCS,D_ik_all_MCCS,D0_ik_all_MCCS,u0_k_all_MCCS,v0_k_all_MCCS,u_k_all_MCCS,v_k_all_MCCS,Cu_k_all_MCCS,Cv_k_all_MCCS,MSE_all_MCCS]...
%     =Cosserat_Simulation_Montecarlo(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
% % %% Monte-carlo simulation fully nonlinear Cosserat
% % plot_figure_id=1081;
% % [r_all_MCFNCS,D_ik_all_MCFNCS,D0_ik_all_MCFNCS,u0_k_all_MCFNCS,v0_k_all_MCFNCS,u_k_all_MCFNCS,v_k_all_MCFNCS,Cu_k_all_MCFNCS,Cv_k_all_MCFNCS,MSE_all_MCFNCS]...
% %     =Fully_nonlin_Cosserat_Simulation_Montecarlo(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
%%
n_cables=4;
R_cables=[R_straight,R_helix,R_helix,R_helix];
cables_strain_all=1.0*[straight_cable_strain_all;helix_cable_strain_all];
 omega_cables=[omega_straight,omega_helix,omega_helix,omega_helix];
 theta_cables=[theta_straight_0+2*pi*([1:n_straight]-1)/n_straight,theta_helix_0+2*pi*([1:n_helix]-1)/n_helix];
plot_figure_id=0;
 [r_all_NSCRK,D_ik_all_NSCRK,D0_ik_all_NSCRK,u0_k_all_NSCRK,v0_k_all_NSCRK,u_k_all_NSCRK,v_k_all_NSCRK,Cu_k_all_NSCRK,Cv_k_all_NSCRK,MSE_all_NSCRK,s_r_all_NSCRK,s_D_ik_all_NSCRK]...
        =Cosserat_WeakRK_Samples_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);

%% Save
% name='data'
% name1 = ['rope' strrep(datestr(now), ':', '-')];
% mkdir(name1);
% save([name1 '\' name])
%%
figure(900001)
hline = findobj(gcf, 'type', 'line')
set(hline(:),'Linewidth',2)
set(hline(1),'Color','b')
xlabel('x [m]',Interpreter='latex')
ylabel('y [m]',Interpreter='latex')
zlabel('z [m]',Interpreter='latex')
axis equal
zlim([0,0.2])
xlim([-0.007,0.02])
ylim([-0.007,0.01])
print(gcf,'foo.pdf','-dpdf','-r300');
%%
figure(700002)
ax1 = axes();
Imm = imread('ropephotodata/shape1.JPG');
uistack(ax1,'bottom');
set(ax1,'handlevisibility','off', ...
            'visible','off')
image(Imm)
ax2 = axes('Color', 'none');

hold on
plotinterval=10;
n_points=length(s_all);
plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'m',LineWidth=0.5,DisplayName='\mbox{\boldmath$r$}');
hold on
for isample=1:n_samples
    r_all_isample=squeeze(s_r_all_NSCRK(isample,:,:));
    D_ik_all_isample=squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
    % [r_all_isample,D_ik_all_isample]=rotate_in_plane_rod(r_all_isample,D_ik_all_isample,rod_axis_rs);
    plot3(squeeze(r_all_isample(1,1:plotinterval:n_points)), squeeze(r_all_isample(2,1:plotinterval:n_points)), squeeze(r_all_isample(3,1:plotinterval:n_points)),'Color',[1 0 1,0.02],'linewidth',0.001,'HandleVisibility','off');
end
plotinterval=200;
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,1,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,1,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,1,1:plotinterval:n_points))',0.2,'LineWidth',0.5,'Color',[1 0 0 1],DisplayName='$\mbox{\boldmath$d$}_1$')
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,2,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,2,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,2,1:plotinterval:n_points))',0.2,'LineWidth',0.5,'Color',[0 1 0 1],DisplayName='$\mbox{\boldmath$d$}_2$')
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,3,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,3,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,3,1:plotinterval:n_points))',0.2,'LineWidth',0.5,'Color',[0 0 1 1],DisplayName='$\mbox{\boldmath$d$}_3$')
legend(Interpreter='latex')
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
h.ZAxis.Visible = 'off';
hold on
% xlabel('x [m]',Interpreter='latex')
% ylabel('y [m]',Interpreter='latex')
% zlabel('z [m]',Interpreter='latex')
axis equal
% xlim([-0.4,0.2])
% ylim([-0.2,0.8])
% zlim([0,4.2])
% view(gca,[-53.6002789366093 21.9491441384439]);
% grid(gca,'on');
% axis(gca,'tight');
% hold(gca,'off');
% Set the remaining axes properties
% set(gca,'DataAspectRatio',[1 1 1]);
% Create legend
%%shape1camera
set(ax2,'CameraPosition',...
    [-21.7442856942752 -19.4303916009191 4.5702702348023],'CameraTarget',...
    [-0.0242841606193389 0.122408282422188 2.51487214345462],'CameraUpVector',...
    [0.671000923345908 -0.73558006009508 0.0931651010811266],'CameraViewAngle',...
    9.5816282408368,'Color','none','DataAspectRatio',[1 1 1],...
    'PlotBoxAspectRatio',[1.26923676822294 1 1.93629793019323],'Projection',...
    'perspective');
% Create legend
legend1 = legend(ax2,'show');
set(legend1,...
    'Position',[0.150891740196385 0.482171312461336 0.12344106878154 0.14003984291715],...
    'Interpreter','latex');
% %%shape3camera
% set(ax2,'CameraPosition',...
%     [-19.7070262946338 -18.5910447497699 14.6646041739029],'CameraTarget',...
%     [0.117645490544377 -0.224836563857804 2.52248682802173],'CameraUpVector',...
%     [0.65177378739663 -0.754486777099267 -0.0770755035339639],'CameraViewAngle',...
%     8.77467640354203,'Color','none','DataAspectRatio',[1 1 1],...
%     'PlotBoxAspectRatio',[1.26923676822294 1 1.93629793019323]);
% % Create legend
% legend1 = legend(ax2,'show');
% set(legend1,...
%     'Position',[0.150891740196385 0.477720893724562 0.12344106878154 0.148940680390697],...
%     'Interpreter','latex');
%%shape2camera
% set(ax2,'CameraPosition',...
%     [30.407137662931 7.66000760741913 3.94982444251594],'CameraTarget',...
%     [-0.0650409972530577 0.205771757147424 2.60552347413288],'CameraUpVector',...
%     [-0.234427085810119 0.969979486717434 -0.0646818118636014],...
%     'CameraViewAngle',9.81437938309367,'Color','none','DataAspectRatio',[1 1 1],...
%     'PlotBoxAspectRatio',[1.26923676822294 1 1.93629793019323]);
% % Create legend
% legend1 = legend(ax2,'show');
% set(legend1,...
%     'Position',[0.151856123546241 0.477720893724562 0.121512302081828 0.148940680390697],...
%     'Interpreter','latex');

print(gcf,'ropeshape.png','-dpng','-r500');
%%
name='data'
name1 = ['rope_largesample' strrep(datestr(now), ':', '-')];
mkdir(name1);
save([name1 '\' name])
figureids=[1001:10:1081];
figureids=[figureids figureids+1];
for fid=figureids
    figure(fid);
    savefig([name1 '\' num2str(fid)]);
end