close all
clear
clc
%% settings
%%
s_all=[0:0.0026:4]; % 0:your gauge pitch:total length
%%
% s_all=[0:1:100];
rod_axis_rs=[zeros(size(s_all));zeros(size(s_all));s_all];
isassigneddir=true;
assigned_dir=[1,0,0];
get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
R_cable=sqrt((5.031/2)^2+(4.745/2)^2)*1e-3;
theta_cable=atan(5.031/4.745);
%%Stochastic 
R_helix=R_cable;R_straight=R_cable;omega_helix=0;omega_straight=0;theta_helix_0=-theta_cable-pi/2;theta_straight_0=theta_cable-1*pi/2;n_helix=2;n_straight=1;
%%
measure_error_level=20e-6;
n_samples=100;
n_substep=1;
%% generate strains
% [r_all,D_ik_all,helix_cable_strain_all,straight_cable_strain_all,u_0k_all,v_0k_all]=generate_MSF_kinematics(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,design_uk_all,design_vk_all);
% figure(9)
% plot(s_all,helix_cable_strain_all')
% hold on
% plot(s_all,straight_cable_strain_all')
% check_null_space_of_cable(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all(1),u_0k_all(:,1),v_0k_all(:,1),u_0k_all(:,1),v_0k_all(:,1))

%
readnervesensvithvideo;
%% 
vw = VideoWriter("nervesens.avi");
vw.FrameRate=10/0.6;
vw_together = VideoWriter("nervesens_together.avi");
vw_together.FrameRate=10/0.6;
open(vw)
open(vw_together)
k_frame_inter=10;
n_points=length(s_all);

vidObj = VideoReader("nervesenswithvideo/converted_00000.mp4");
ratio=87/540;
width=1158-702;
lengthim=width/ratio;
movpt=[702 1048; 1158 1053; 1057 133; 871 131];
fixpt=[702 1048; 702+width 1048; 702+width 1048-lengthim; 702 1048-lengthim];
tform = fitgeotform2d(movpt,fixpt,"projective");
invtform = invert(tform);
framerate=100;
frameid=0;

% for k = 1:k_frame_inter:length(time4)
for k = [1362 2481 4122 4841 5432 5982]
    slicetimeid=k
    dat1=dat1-dat1(2,:);
    dat2=dat2-dat2(2,:);
    % dat3=dat3-dat3(2,:);
    dat4=dat4-dat4(2,:);
    reading_strainu1=dat1(slicetimeid,:);
    loc_strainu1=loc1-loc1(1);
    reading_strainu2=dat2(slicetimeid,:);
    loc_strainu2=loc2-loc2(1);
    % reading_strainl1=dat3(slicetimeid,:);
    % loc_strainl1=loc3-loc3(1);
    reading_strainl2=dat4(slicetimeid,:);
    loc_strainl2=loc4-loc4(1);
    strainu1=interp1(loc_strainu1,reading_strainu1,s_all*loc_strainu1(end)/s_all(end));
    strainu2=interp1(loc_strainu2,reading_strainu2,s_all*loc_strainu2(end)/s_all(end));
    % strainl1=interp1(loc_strainl1,reading_strainl1,s_all*loc_strainl1(end)/s_all(end));
    strainl2=interp1(loc_strainl2,reading_strainl2,s_all*loc_strainl2(end)/s_all(end));
    nans=isnan(strainl2);
    strainl2(nans) = interp1(s_all(~nans), strainl2(~nans), s_all(nans));
    nans=isnan(strainu1);
    strainu1(nans) = interp1(s_all(~nans), strainu1(~nans), s_all(nans));
    nans=isnan(strainu2);
    strainu2(nans) = interp1(s_all(~nans), strainu2(~nans), s_all(nans));
    % nans=isnan(strainl1);
    % strainl1(nans) = interp1(s_all(~nans), strainl1(~nans), s_all(nans));
    helix_cable_strain_all=[strainu1;strainl2]*1e-6;
    straight_cable_strain_all=[strainu2]*1e-6;
    plot_figure_id=0;
    [r_all_SACS,D_ik_all_SACS,D0_ik_all_SACS,u0_k_all_SACS,v0_k_all_SACS,u_k_all_SACS,v_k_all_SACS,Cu_k_all_SACS,Cv_k_all_SACS,MSE_all_SACS]=Cosserat_analytical_SDE(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
    % [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,helix_cable_strain_all,straight_cable_strain_all,measure_error_level,plot_figure_id,n_substep);
    figure(1001)
    plotinterval=ceil(n_points/100);
    plot3(squeeze(r_all_SACS(1,1:plotinterval:n_points)), squeeze(r_all_SACS(2,1:plotinterval:n_points)), squeeze(r_all_SACS(3,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    quiver3(squeeze(r_all_SACS(1,1:plotinterval:n_points)), squeeze(r_all_SACS(2,1:plotinterval:n_points)), squeeze(r_all_SACS(3,1:plotinterval:n_points)),squeeze(D_ik_all_SACS(1,1,1:plotinterval:n_points))', squeeze(D_ik_all_SACS(2,1,1:plotinterval:n_points))', squeeze(D_ik_all_SACS(3,1,1:plotinterval:n_points))','r')
    quiver3(squeeze(r_all_SACS(1,1:plotinterval:n_points)), squeeze(r_all_SACS(2,1:plotinterval:n_points)), squeeze(r_all_SACS(3,1:plotinterval:n_points)),squeeze(D_ik_all_SACS(1,2,1:plotinterval:n_points))', squeeze(D_ik_all_SACS(2,2,1:plotinterval:n_points))', squeeze(D_ik_all_SACS(3,2,1:plotinterval:n_points))','g')
    quiver3(squeeze(r_all_SACS(1,1:plotinterval:n_points)), squeeze(r_all_SACS(2,1:plotinterval:n_points)), squeeze(r_all_SACS(3,1:plotinterval:n_points)),squeeze(D_ik_all_SACS(1,3,1:plotinterval:n_points))', squeeze(D_ik_all_SACS(2,3,1:plotinterval:n_points))', squeeze(D_ik_all_SACS(3,3,1:plotinterval:n_points))','b')
    hold off
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    xlim([-1,1])
    ylim([-1,1])
    zlim([0,4.5])
    view(462.4878,10.2476)
    frame = getframe(gcf);
    writeVideo(vw,frame)
    videotime=seconds(time4(slicetimeid)-time4(1))+1;
    title(sprintf("Current Time = %.3f sec",videotime))
    % videotime=vidObj.CurrentTime
    while(hasFrame(vidObj) && vidObj.CurrentTime<videotime)
        frameid=frameid+1;
        frame = readFrame(vidObj);
        if (vidObj.CurrentTime>videotime)
            warppedframe = imwarp(frame,tform);
            framerim=200;
            croppedframe = imcrop(warppedframe,[2710-framerim 1802-framerim 3161-2710+2*framerim 4075-1802+2*framerim]);
            % binedfame=imbinarize(im2gray(croppedframe));
            % edgedframe = edge(binedfame,'Canny');
            figure(1003)
            croppedframe = flip(croppedframe ,1);
            image([-0.8124,0.8054],[-0.3533,4.694],croppedframe)
            ax = gca;
            ax.YDir = 'normal';
            % set(gca, 'YDir','reverse')
            axis equal
            hold on
            % [z_frame,y_frame]=find(edgedframe);
            % scalefac=(2270-232)/4;
            % z_frame=2270-z_frame;
            % z_frame=z_frame/scalefac;
            % y_frame=y_frame/scalefac;$\boldmath $r$ $
            % y_frame=y_frame-0.4867;
            % [z_frame,sort_idenx]=sort(z_frame);
            % y_frame=y_frame(sort_idenx);
            % plot(y_frame,z_frame)
            ylim([-0.5,4.5]);
            xlim([-1,1]);
            meanline_x=squeeze(r_all_SACS(2,1:plotinterval:n_points));
            meanline_y=squeeze(r_all_SACS(3,1:plotinterval:n_points));
            errdir_x=squeeze(D_ik_all_SACS(2,2,1:plotinterval:n_points))';
            errdir_y=squeeze(D_ik_all_SACS(3,2,1:plotinterval:n_points))';
            stdp_x=meanline_x+errdir_x.*MSE_all_SACS(1:plotinterval:n_points);
            stdp_y=meanline_y+errdir_y.*MSE_all_SACS(1:plotinterval:n_points);
            stdm_x=meanline_x-errdir_x.*MSE_all_SACS(1:plotinterval:n_points);
            stdm_y=meanline_y-errdir_y.*MSE_all_SACS(1:plotinterval:n_points);
            plot(meanline_x, meanline_y,'-m',LineWidth=2);
            plot(stdp_x, stdp_y,'--k',LineWidth=1);
            plot(stdm_x, stdm_y,'--k',LineWidth=1);
            legend('$E[${\boldmath$r$}$]$','$E[${\boldmath$r$}$] \pm Std[${\boldmath$r$}$]$','interpreter',"latex")
            hold on
            title(sprintf("Current Time = %.3f sec",vidObj.CurrentTime))
            set(gcf,"Position",[675,50.33333333333333,321.3333333333333,670.6666666666666])
            axis equal
            xlim([-1,1])
            
            % quiver( squeeze(r_all_SACS(2,1:plotinterval:n_points)), squeeze(r_all_SACS(3,1:plotinterval:n_points)), squeeze(D_ik_all_SACS(2,2,1:plotinterval:n_points))', squeeze(D_ik_all_SACS(3,2,1:plotinterval:n_points))','g')
            % quiver( squeeze(r_all_SACS(2,1:plotinterval:n_points)), squeeze(r_all_SACS(3,1:plotinterval:n_points)), squeeze(D_ik_all_SACS(2,3,1:plotinterval:n_points))', squeeze(D_ik_all_SACS(3,3,1:plotinterval:n_points))','b')
            hold off
            frame = getframe(gcf);
            writeVideo(vw_together,frame)
            xlabel("$x$ [m]",'interpreter',"latex")
            ylabel("$y$ [m]",'interpreter',"latex")
            saveas(gcf,sprintf("Current Time = %.3f sec.pdf",vidObj.CurrentTime))
            % pause(0.01/vidObj.FrameRate)
            break
        end
    end
end
    
close(vw)
close(vw_together)
clear vidObj
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
