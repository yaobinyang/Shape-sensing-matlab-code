%% WeakRK
function [r_all,D_ik_all,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all]...
    =Cosserat_analytical(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep)
n_points=length(s_all);
n_measure=n_helix+n_straight;
D0_ik_all=get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
pD_ik_ps_all=get_partial_D_ik_partial_s(D0_ik_all,s_all);
grad_s=gradient(s_all);
mean_value_m_of_Y_all=zeros(12,n_points);
% covariance_P_of_Y_all=zeros(12,12,n_points);
% vecP_all=zeros(12*12,n_points);
mean_value_m_of_Y_all(:,1)=[rod_axis_rs(:,1);D0_ik_all(:,1,1);D0_ik_all(:,2,1);D0_ik_all(:,3,1)];
covariance_P_of_Y_all(:,:,1)=mean_value_m_of_Y_all(:,1)*mean_value_m_of_Y_all(:,1)';
vecP_all(:,1)=reshape(covariance_P_of_Y_all(:,:,1),[12*12,1]);
u0_k_all=zeros(3,n_points);
v0_k_all=zeros(3,n_points);
u_k_all=zeros(3,n_points);
v_k_all=zeros(3,n_points);
% Cu_k_all=zeros(3,n_measure,n_points);
% Cv_k_all=zeros(3,n_measure,n_points);
A_all=zeros(12,12,n_points);
B_all=zeros(n_measure,12,12,n_points);
% MSE_all=zeros(1,n_points);
r_all=rod_axis_rs;
D_ik_all=D0_ik_all;
parfor point=1:n_points
    point
    s=s_all(point);
    epsilon_p_helix=helix_cable_strain_all(:,point);
    epsilon_p_straight=straight_cable_strain_all(:,point);
    [u0_k,v0_k]=get_ini_darboux_strain_u0_k_and_v0_k(squeeze(D0_ik_all(:,:,point)),squeeze(pD_ik_ps_all(:,:,point)));
    u0_k_all(:,point)=u0_k;
    v0_k_all(:,point)=v0_k;
    [A,B,u_k,v_k,Cu_k,Cv_k]=nonlin_rod_assymbling_matrixs(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,epsilon_p_helix,epsilon_p_straight,measure_error_level,grad_s(point));
    u_k_all(:,point)=u_k;
    v_k_all(:,point)=v_k;
    % Cu_k_all(:,:,point)=Cu_k;
    % Cv_k_all(:,:,point)=Cv_k;
    A_all(:,:,point)=A
    % B_all(:,:,:,point)=B
end
permutedA=permute(A_all,[3,1,2]);
% permutedB=permute(B_all,[4,1,2,3]);
helperfunA =@(t) squeeze(interp1(s_all',permutedA,t,"spline"));
% helperfunB =@(t) squeeze(interp1(s_all',permutedB,t,"spline"));
opts = odeset('RelTol',1e-9,'AbsTol',1e-12,'Refine',n_substep,'MaxStep',min(grad_s)/n_substep);
[t,solm]=ode45(@(s,m)helperfunA(s)*m,s_all,mean_value_m_of_Y_all(:,1),opts);
% [t,solvecp]=ode45(@(s,vecP)get_rhs_for_vecP(vecP,helperfunA(s),helperfunB(s)),s_all,vecP_all(:,1),opts);

for point=2:n_points
    % point
    % % helperintmatA=integral(helperfunA,s_all(1),s_all(point),ArrayValued=true,AbsTol=1e-12);
    % % helperintmatB=integral(helperfunB,s_all(1),s_all(point),ArrayValued=true,AbsTol=1e-12);
    % mean_value_m_of_Y_all(:,point)=expm(helperintmatA)*mean_value_m_of_Y_all(:,1);
    mean_value_m_of_Y_all(:,point)=solm(point,:);
    % covariance_P_of_Y_all(:,:,point)=reshape(solvecp(point,:)',[12,12]);
    r_all(:,point)=mean_value_m_of_Y_all(1:3,point);
    D_ik_all(:,1,point)=mean_value_m_of_Y_all(4:6,point);
    D_ik_all(:,2,point)=mean_value_m_of_Y_all(7:9,point);
    D_ik_all(:,3,point)=mean_value_m_of_Y_all(10:12,point);
    % MSE_all(point)=sqrt(covariance_P_of_Y_all(1,1,point)+covariance_P_of_Y_all(2,2,point)+covariance_P_of_Y_all(3,3,point)-norm(mean_value_m_of_Y_all(1:3,point))^2);
end
if plot_figure_id
    plotinterval=round(n_points/100);
    figure(plot_figure_id)
    plot3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),'m');
    grid on
    hold on
    figure(plot_figure_id+1)
    plot3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    quiver3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),squeeze(D_ik_all(1,1,1:plotinterval:n_points))', squeeze(D_ik_all(2,1,1:plotinterval:n_points))', squeeze(D_ik_all(3,1,1:plotinterval:n_points))','r')
    quiver3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),squeeze(D_ik_all(1,2,1:plotinterval:n_points))', squeeze(D_ik_all(2,2,1:plotinterval:n_points))', squeeze(D_ik_all(3,2,1:plotinterval:n_points))','g')
    quiver3(squeeze(r_all(1,1:plotinterval:n_points)), squeeze(r_all(2,1:plotinterval:n_points)), squeeze(r_all(3,1:plotinterval:n_points)),squeeze(D_ik_all(1,3,1:plotinterval:n_points))', squeeze(D_ik_all(2,3,1:plotinterval:n_points))', squeeze(D_ik_all(3,3,1:plotinterval:n_points))','b')
end