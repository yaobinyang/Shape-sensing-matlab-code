%% Milstein
function [mean_r_all,mean_D_ik_all,D0_ik_all,u0_k_all,v0_k_all,u_k_all,v_k_all,Cu_k_all,Cv_k_all,MSE_all]...
    =Cosserat_Milstein_Samples(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples)
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
MSE_all=zeros(1,n_points);
Cu_k_all=zeros(3,n_measure,n_points);
Cv_k_all=zeros(3,n_measure,n_points);
A_all=zeros(12,12,n_points);
B_all=zeros(n_measure,12,12,n_points);
s_r_all=zeros(n_samples,3,n_points);
s_D_ik_all=zeros(n_samples,3,3,n_points);
s_r_all=repmat(rod_axis_rs,[1,1,n_samples]);
s_r_all=permute(s_r_all,[3,1,2]);
s_D_ik_all=repmat(D0_ik_all,[1,1,1,n_samples]);
s_D_ik_all=permute(s_D_ik_all,[4,1,2,3]);
parfor point=1:n_points
    s=s_all(point);
    ds=grad_s(point);
    epsilon_p_helix=helix_cable_strain_all(:,point);
    epsilon_p_straight=straight_cable_strain_all(:,point);
    [u0_k,v0_k]=get_ini_darboux_strain_u0_k_and_v0_k(squeeze(D0_ik_all(:,:,point)),squeeze(pD_ik_ps_all(:,:,point)));
    u0_k_all(:,point)=u0_k;
    v0_k_all(:,point)=v0_k;
    [A,B,u_k,v_k,Cu_k,Cv_k]=nonlin_rod_assymbling_matrixs(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,epsilon_p_helix,epsilon_p_straight,measure_error_level,ds);
    u_k_all(:,point)=u_k;
    v_k_all(:,point)=v_k;
    Cu_k_all(:,:,point)=Cu_k;
    Cv_k_all(:,:,point)=Cv_k;
    A_all(:,:,point)=A
    B_all(:,:,:,point)=B
end
permutedA=permute(A_all,[1,2,3]);
permutedB=permute(B_all,[1,2,3,4]);
spline_A=spline(s_all',permutedA);
spline_B=spline(s_all',permutedB);
helperfunA =@(t) ppval(spline_A,t);
helperfunB =@(t) ppval(spline_B,t);
parfor sample=1:n_samples
    disp("now on sample = ")
    disp(sample)
    Ys_Milstein=zeros(12,n_points);
    Ys_Milstein(:,1)=[rod_axis_rs(:,1);D0_ik_all(:,1,1);D0_ik_all(:,2,1);D0_ik_all(:,3,1)];
    for point=2:n_points
        ds=grad_s(point-1)/n_substep;
        s=s_all(point-1);
        Y_update=Ys_Milstein(:,point-1);
        for step=1:n_substep
            A=helperfunA(s+step*ds);
            B=helperfunB(s+step*ds);
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
        s_r_all(sample,:,point)=Ys_Milstein(1:3,point);
        s_D_ik_all(sample,:,:,point)=[Ys_Milstein(4:6,point) Ys_Milstein(7:9,point) Ys_Milstein(10:12,point)];
    end
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
end