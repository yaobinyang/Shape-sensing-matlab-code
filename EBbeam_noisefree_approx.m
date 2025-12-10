%% WeakRK
function [r_all,D_ik_all,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all,MSE_all]...
    =EBbeam_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep)
n_points=length(s_all);
n_measure=n_helix+n_straight;
D0_ik_all=get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
pD_ik_ps_all=get_partial_D_ik_partial_s(D0_ik_all,s_all);
grad_s=gradient(s_all);
mean_value_m_of_Y_all=zeros(12,n_points);
covariance_P_of_Y_all=zeros(12,12,n_points);
vecP_all=zeros(12*12,n_points);
mean_value_m_of_Y_all(:,1)=[rod_axis_rs(:,1);D0_ik_all(:,1,1);D0_ik_all(:,2,1);D0_ik_all(:,3,1)];
covariance_P_of_Y_all(:,:,1)=mean_value_m_of_Y_all(:,1)*mean_value_m_of_Y_all(:,1)';
vecP_all(:,1)=reshape(covariance_P_of_Y_all(:,:,1),[12*12,1]);
u0_k_all=zeros(3,n_points);
v0_k_all=zeros(3,n_points);
u_k_all=zeros(3,n_points);
v_k_all=zeros(3,n_points);
Cu_k_all=zeros(3,n_measure,n_points);
Cv_k_all=zeros(3,n_measure,n_points);
A_all=zeros(12,12,n_points);
B_all=zeros(n_measure,12,12,n_points);
MSE_all=zeros(1,n_points);
r_all=rod_axis_rs;
D_ik_all=D0_ik_all;
parfor point=1:n_points
    s=s_all(point);
    epsilon_p_helix=helix_cable_strain_all(:,point);
    epsilon_p_straight=straight_cable_strain_all(:,point);
    [u0_k,v0_k]=get_ini_darboux_strain_u0_k_and_v0_k(squeeze(D0_ik_all(:,:,point)),squeeze(pD_ik_ps_all(:,:,point)));
    u0_k_all(:,point)=u0_k;
    v0_k_all(:,point)=v0_k;
    [u_k,v_k]=nonlinear_least_square_solve_darboux(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,epsilon_p_helix,epsilon_p_straight);
    u_k_all(:,point)=u_k;
    v_k_all(:,point)=v_k;
end

for point=1:n_points
    s=s_all(point);
    ds=grad_s(point);
    uk=u_k_all(:,point);
    vk=v_k_all(:,point);
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
        u0_k_all(:,point)=u_0k_at_s;
        v0_k_all(:,point)=v_0k_at_s;
        D_ik_s1=D_ik_all(:,:,point-1);
        norm_uk=max(norm(uk),1e-10);
        helpermat_exp=eye(3)+Ucross_helpermat_trans*ds;
        helpermat_int=eye(3)*ds+0.5*Ucross_helpermat_trans*ds^2;
        D_ik_s=D_ik_s1*helpermat_exp;
        D_ik_all(:,:,point)=D_ik_s;
        r_s=D_ik_s1*helpermat_int*[0;0;vk(3)]+r_all(:,point-1);
        r_all(:,point)=r_s;
    end
end

if plot_figure_id
    plotinterval=ceil(n_points/100);
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