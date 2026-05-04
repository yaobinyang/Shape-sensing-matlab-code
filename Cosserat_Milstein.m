%% Milstein
function [r_all,D_ik_all,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all]...
    =Cosserat_Milstein(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep)
n_points=length(s_all);
n_measure=n_helix+n_straight;
D0_ik_all=get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
pD_ik_ps_all=get_partial_D_ik_partial_s(D0_ik_all,s_all);
grad_s=gradient(s_all);
Ys_Milstein=zeros(12,n_points);
Ys_Milstein(:,1)=[rod_axis_rs(:,1);D0_ik_all(:,1,1);D0_ik_all(:,2,1);D0_ik_all(:,3,1)];
u0_k_all=zeros(3,n_points);
v0_k_all=zeros(3,n_points);
u_k_all=zeros(3,n_points);
v_k_all=zeros(3,n_points);
Cu_k_all=zeros(3,n_measure,n_points);
Cv_k_all=zeros(3,n_measure,n_points);
A_all=zeros(12,12,n_points);
B_all=zeros(n_measure,12,12,n_points);
r_all=rod_axis_rs;
D_ik_all=D0_ik_all;
parfor point=1:n_points
    s=s_all(point);
    epsilon_p_helix=helix_cable_strain_all(:,point);
    epsilon_p_straight=straight_cable_strain_all(:,point);
    [u0_k,v0_k]=get_ini_darboux_strain_u0_k_and_v0_k(squeeze(D0_ik_all(:,:,point)),squeeze(pD_ik_ps_all(:,:,point)));
    u0_k_all(:,point)=u0_k;
    v0_k_all(:,point)=v0_k;
    [A,B,u_k,v_k,Cu_k,Cv_k]=nonlin_rod_assymbling_matrixs(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,epsilon_p_helix,epsilon_p_straight,measure_error_level);
    u_k_all(:,point)=u_k;
    v_k_all(:,point)=v_k;
    Cu_k_all(:,:,point)=Cu_k;
    Cv_k_all(:,:,point)=Cv_k;
    A_all(:,:,point)=A
    B_all(:,:,:,point)=B
end
for point=2:n_points
    point
    ds=grad_s(point)/n_substep;
    Y_update=Ys_Milstein(:,point-1);
    for step=1:n_substep
        A=((n_substep-step)/n_substep)*squeeze(A_all(:,:,point-1))+(step/n_substep)*squeeze(A_all(:,:,point));
        B=((n_substep-step)/n_substep)*squeeze(B_all(:,:,:,point-1))+(step/n_substep)*squeeze(B_all(:,:,:,point));
        Y=Y_update;
        DeltaW_j=sqrt(ds)*randn([n_measure 1]);
        for p=1:12
            for q=1:12
                Y_update(p)=Y_update(p)+A(p,q)*Y(q)*ds;
                for j=1:n_measure
                    Y_update(p)=Y_update(p)+B(j,p,q)*Y(q)*DeltaW_j(j);
                    for l=1:12
                        Y_update(p)=Y_update(p)-0.5*B(j,l,q)*Y(q)*B(j,p,l)*DeltaW_j(j);
                        for r=1:n_measure
                            Y_update(p)=Y_update(p)+0.5*B(j,l,q)*Y(q)*B(r,p,l)*DeltaW_j(j)*DeltaW_j(r);
                        end
                    end
                end
            end
        end
    end
    Ys_Milstein(:,point)=Y_update;
    r_all(:,point)=Ys_Milstein(1:3,point);
    D_ik_all(:,1,point)=Ys_Milstein(4:6,point);
    D_ik_all(:,2,point)=Ys_Milstein(7:9,point);
    D_ik_all(:,3,point)=Ys_Milstein(10:12,point);
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
end