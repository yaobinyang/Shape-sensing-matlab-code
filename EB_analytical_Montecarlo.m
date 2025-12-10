function [mean_r_all,mean_D_ik_all,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all,MSE_all]...
    =EB_analytical_Montecarlo(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
n_points=length(s_all);
n_measure=n_helix+n_straight;
MSE_all=zeros(1,n_points);
Cu_k_all=zeros(3,n_measure,n_points);
Cv_k_all=zeros(3,n_measure,n_points);
s_r_all=zeros(n_samples,3,n_points);
s_D_ik_all=zeros(n_samples,3,3,n_points);
D0_ik_all=[];
u0_k_all=[];
v0_k_all=[];
u_k_all=[];
v_k_all=[];
parfor sample=1:n_samples
    disp("now on sample = ")
    disp(sample)
    errored_helix_cable_strain_all=helix_cable_strain_all+measure_error_level*randn(size(helix_cable_strain_all));
    errored_straight_cable_strain_all=straight_cable_strain_all+measure_error_level*randn(size(straight_cable_strain_all));
    [r_all,D_ik_all,D0_ik_all_,u0_k_all_,v0_k_all_,u_k_all_,v_k_all_]...
        =EBbeam_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,errored_helix_cable_strain_all,errored_straight_cable_strain_all,measure_error_level,0,n_substep);
    s_r_all(sample,:,:)=r_all;
    s_D_ik_all(sample,:,:,:)=D_ik_all;
end
mean_r_all=squeeze(mean(s_r_all,1));
std_r_all=squeeze(std(s_r_all,0,1));
mean_D_ik_all=squeeze(mean(s_D_ik_all,1));
MSE_all(:)=sqrt(std_r_all(1,:).^2+std_r_all(2,:).^2+std_r_all(3,:).^2);
plotinterval=ceil(n_points/100);
transpancy=0.2;
if plot_figure_id
    figure(plot_figure_id)
    plot3(squeeze(mean_r_all(1,1:plotinterval:n_points)), squeeze(mean_r_all(2,1:plotinterval:n_points)), squeeze(mean_r_all(3,1:plotinterval:n_points)), 'Color',[0 0 0 1]);
    grid on
    hold on
    figure(plot_figure_id+1)
    plot3(squeeze(mean_r_all(1,1:plotinterval:n_points)), squeeze(mean_r_all(2,1:plotinterval:n_points)), squeeze(mean_r_all(3,1:plotinterval:n_points)), 'Color',[0 0 0 1]);
    grid on
    axis equal
    hold on
    quiver3(squeeze(mean_r_all(1,1:plotinterval:n_points)), squeeze(mean_r_all(2,1:plotinterval:n_points)), squeeze(mean_r_all(3,1:plotinterval:n_points)),squeeze(mean_D_ik_all(1,1,1:plotinterval:n_points))', squeeze(mean_D_ik_all(2,1,1:plotinterval:n_points))', squeeze(mean_D_ik_all(3,1,1:plotinterval:n_points))','Color',[1 0 0 1])
    quiver3(squeeze(mean_r_all(1,1:plotinterval:n_points)), squeeze(mean_r_all(2,1:plotinterval:n_points)), squeeze(mean_r_all(3,1:plotinterval:n_points)),squeeze(mean_D_ik_all(1,2,1:plotinterval:n_points))', squeeze(mean_D_ik_all(2,2,1:plotinterval:n_points))', squeeze(mean_D_ik_all(3,2,1:plotinterval:n_points))','Color',[0 1 0 1])
    quiver3(squeeze(mean_r_all(1,1:plotinterval:n_points)), squeeze(mean_r_all(2,1:plotinterval:n_points)), squeeze(mean_r_all(3,1:plotinterval:n_points)),squeeze(mean_D_ik_all(1,3,1:plotinterval:n_points))', squeeze(mean_D_ik_all(2,3,1:plotinterval:n_points))', squeeze(mean_D_ik_all(3,3,1:plotinterval:n_points))','Color',[0 0 1 1])

    for sample=1:n_samples

        figure(plot_figure_id)
        plot3(squeeze(s_r_all(sample,1,1:plotinterval:n_points)), squeeze(s_r_all(sample,2,1:plotinterval:n_points)), squeeze(s_r_all(sample,3,1:plotinterval:n_points)), 'Color',[0 0 0 transpancy]);
        grid on
        hold on
        figure(plot_figure_id+1)
        plot3(squeeze(s_r_all(sample,1,1:plotinterval:n_points)), squeeze(s_r_all(sample,2,1:plotinterval:n_points)), squeeze(s_r_all(sample,3,1:plotinterval:n_points)),'Color',[0 0 0 transpancy]);
        grid on
        axis equal
        hold on
        quiver3(squeeze(s_r_all(sample,1,1:plotinterval:n_points)), squeeze(s_r_all(sample,2,1:plotinterval:n_points)), squeeze(s_r_all(sample,3,1:plotinterval:n_points)),squeeze(s_D_ik_all(sample,1,1,1:plotinterval:n_points)), squeeze(s_D_ik_all(sample,2,1,1:plotinterval:n_points)), squeeze(s_D_ik_all(sample,3,1,1:plotinterval:n_points)),'Color',[1 0 0 transpancy])
        quiver3(squeeze(s_r_all(sample,1,1:plotinterval:n_points)), squeeze(s_r_all(sample,2,1:plotinterval:n_points)), squeeze(s_r_all(sample,3,1:plotinterval:n_points)),squeeze(s_D_ik_all(sample,1,2,1:plotinterval:n_points)), squeeze(s_D_ik_all(sample,2,2,1:plotinterval:n_points)), squeeze(s_D_ik_all(sample,3,2,1:plotinterval:n_points)),'Color',[0 1 0 transpancy])
        quiver3(squeeze(s_r_all(sample,1,1:plotinterval:n_points)), squeeze(s_r_all(sample,2,1:plotinterval:n_points)), squeeze(s_r_all(sample,3,1:plotinterval:n_points)),squeeze(s_D_ik_all(sample,1,3,1:plotinterval:n_points)), squeeze(s_D_ik_all(sample,2,3,1:plotinterval:n_points)), squeeze(s_D_ik_all(sample,3,3,1:plotinterval:n_points)),'Color',[0 0 1 transpancy])
    end
end