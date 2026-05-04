close all
clear
clc
% load andrewdatarevisit/aydata.mat
Tnumber=4;
filename="Strain_T"+num2str(Tnumber)+".xlsx";
s1a_strain=readmatrix("andrewdatarevisit/"+filename,Range="B2",Sheet="S1A");
s1b_strain=readmatrix("andrewdatarevisit/"+filename,Range="B2",Sheet="S1B");
s2a_strain=readmatrix("andrewdatarevisit/"+filename,Range="B2",Sheet="S2A");
s2b_strain=readmatrix("andrewdatarevisit/"+filename,Range="B2",Sheet="S2B");
offsets=readmatrix("andrewdatarevisit/"+filename,Range=[2 1 size(s1a_strain,1)+1 1],Sheet="S1A");
loadstep=[1:size(s1a_strain,2);1:size(s1a_strain,2)];
plotframe=find(loadstep(2,:)~=0);
load("andrewdatarevisit/loadSteps_T"+num2str(Tnumber)+".mat")
% loadstep(2,:)=(2*loadstep(2,:)-1)*100*4.44822;
loadstep(2,:)=loadSteps;
%% settings
%%
n_points=size(s1a_strain,1);
hh=mean(gradient(offsets));%0.2552
% sample_spacing=0.1;
plot_number=200;
measure_error_level=20e-6*sqrt(0.2552/hh);
n_samples=1000;
n_substep=10;
plot_index=1:ceil(n_points/plot_number):n_points;
s = linspace(0, max(offsets), n_points)';
s1a_interpstrain=interp1(offsets,s1a_strain,s,"linear","extrap")*-1;
s1b_interpstrain=interp1(offsets,s1b_strain,s,"linear","extrap")*-1;
s2a_interpstrain=interp1(offsets,s2a_strain,s,"linear","extrap")*-1;
s2b_interpstrain=interp1(offsets,s2b_strain,s,"linear","extrap")*-1;
FrameTimes=[1:1:size(s1a_strain,2)];
cart_coord_x = zeros(size(s));
cart_coord_y = zeros(size(s));
cart_coord_z = -flip(s);

%%
% s_all=[0:1:100];
s_all=s';
rod_axis_rs=[cart_coord_x';cart_coord_y';cart_coord_z'];
isassigneddir=true;
assigned_dir=[1,0,0];
get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);
R_cables=[0.4572,0.4572,0.4572,0.4572];
%NSZ
if Tnumber==1
    theta_cables=-[pi,0,pi/2,3*pi/2];
elseif Tnumber ==2
    theta_cables=-[pi/2,pi*3/2,0,pi];
elseif Tnumber ==3
    theta_cables=-[0,pi,pi*3/2,pi/2];
elseif Tnumber ==4
    theta_cables=-[pi/2,3*pi/2,pi,0];
end
%T1
 % theta_cables=-[pi,0,pi/2,3*pi/2];
%T2
% theta_cables=-[pi/2,pi*3/2,0,pi];
%T3
% theta_cables=-[0,pi,pi*3/2,pi/2];
%T4
% theta_cables=-[pi/2,3*pi/2,pi,0];
omega_cables=[0,0,0,0];
n_cables=4;

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
frameratethisone=1.0;
vw_v2 = VideoWriter("AYPILE_Axialxz elongation","Motion JPEG AVI");
vw_v2.FrameRate=frameratethisone;
open(vw_v2)
vw_v3 = VideoWriter("AYPILE_Axialyz elongation","Motion JPEG AVI");
vw_v3.FrameRate=frameratethisone;
open(vw_v3)
vw = VideoWriter("AYPILE","Motion JPEG AVI");
vw.FrameRate=frameratethisone;
open(vw)
vw_v4 = VideoWriter("AYPILE_strain","Motion JPEG AVI");
vw_v4.FrameRate=frameratethisone;
open(vw_v4)
k_frame_inter=1;
positionvec=[618,47,640,1292];
n_points=length(s_all);
trajectoryendpoint=zeros(3,length(FrameTimes));
figure(1001)
hold on
figure(1004)
hold on
figure(2004)
hold on
figure(70001)
set(gcf,"Position",positionvec)
hold on
r_all_NSCRK_allframe=zeros(3,n_points,size(loadstep,2));
xlsfilename="andrewdatarevisit/"+filename+"deformedpositions.xlsx";

for k = plotframe(1:end)
% for k = 1:size(loadstep,2)-1
    slicetimeid=k;
    cables_strain_all=1*[s1a_interpstrain(:,slicetimeid)';s1b_interpstrain(:,slicetimeid)';s2a_interpstrain(:,slicetimeid)';s2b_interpstrain(:,slicetimeid)'];
    cables_strain_all=flip(cables_strain_all,2);
    %%
    figure(70001)
        clf(70001)
    captionshort=sprintf("%s kN", string(loadstep(2,slicetimeid)));
    plot(cables_strain_all,cart_coord_z,'linewidth',3)
    xlim([-4.5e-4,0])
    legend("s1a","s1b","s2a","s2b")
    xlabel("DFOS strain [$\mu \varepsilon$]",'interpreter',"latex",'fontsize',14)
    ylabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    caption = sprintf('Frame #%d of %d, load = %s kN', slicetimeid, length(FrameTimes), string(loadstep(2,slicetimeid)));
    title(captionshort, 'FontSize', 15);
    frame = getframe(gcf);
    writeVideo(vw_v4,frame)

    %%
    plot_figure_id=0;
    % [r_all_SACS,D_ik_all_SACS,D0_ik_all_SACS,u0_k_all_SACS,v0_k_all_SACS,u_k_all_SACS,v_k_all_SACS,Cu_k_all_SACS,Cv_k_all_SACS,MSE_all_SACS]=Cosserat_analytical_SDE_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep);
    % [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep);

    [r_all_NSCRK,D_ik_all_NSCRK,D0_ik_all_NSCRK,u0_k_all_NSCRK,v0_k_all_NSCRK,u_k_all_NSCRK,v_k_all_NSCRK,Cu_k_all_NSCRK,Cv_k_all_NSCRK,MSE_all_NSCRK,s_r_all_NSCRK,s_D_ik_all_NSCRK]...
        =Cosserat_WeakRK_Samples_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);

    [r_all_NSCRK,D_ik_all_NSCRK]=rotate_in_plane_rod(r_all_NSCRK,D_ik_all_NSCRK,rod_axis_rs);
    r_all_NSCRK_allframe(:,:,slicetimeid)=r_all_NSCRK;
    writecell({'s [m]','x [m]','y [m]','z [m]'},xlsfilename,"Sheet",slicetimeid,'Range','A1')
    writematrix(r_all_NSCRK',xlsfilename,"Sheet",slicetimeid,'Range','B2')
    writematrix(s,xlsfilename,"Sheet",slicetimeid,'Range','A2')
    figure(1001)
    plotinterval=ceil(n_points/100);
    plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    % xlim([-3e-3,3e-3])
    % ylim([-3e-3,3e-3])
    % zlim([-11,2])
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,1,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,1,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,1,1:plotinterval:n_points))','r')
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,2,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,2,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,2,1:plotinterval:n_points))','g')
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,3,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,3,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,3,1:plotinterval:n_points))','b')
    hold off
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)

    % xlim([0,200])
    % ylim([-20,60])
    % zlim([95,140])
    % view(540,-90)
    caption = sprintf('Frame #%d of %d, load = %s kN', slicetimeid, length(FrameTimes), string(loadstep(2,slicetimeid)));
    title(caption, 'FontSize', 15);
    frame = getframe(gcf);
    writeVideo(vw,frame)
    trajectoryendpoint(:,k)=r_all_NSCRK(:,end);

    figure(1004)
    hold on
    col=v_k_all_NSCRK(3,1:plotinterval:n_points)-1;
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
    % axis equal
    xlim([-3e-3,3e-3])
    ylim([-3e-3,3e-3])
    zlim([-11,2])
    clim([-2e-4,0])
    view(0,0)
    % caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(FrameTimes), string(FrameTimes(slicetimeid)));
    title({"Axial elongation - plane xz [1]";caption}, 'FontSize', 15);
    frame = getframe(gcf);
    writeVideo(vw_v2,frame)
    clf(1004)
    figure(2004)
    hold on

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
    % axis equal
    xlim([-3e-3,3e-3])
    ylim([-3e-3,3e-3])
    zlim([-11,2])
    clim([-2e-4,0])
    view(90,0)
    title({"Axial elongation - plane yz [1]";caption}, 'FontSize', 15);
    frame = getframe(gcf);
    writeVideo(vw_v3,frame)
    clf(2004)
    figure(1005)
    hold on
    captionshort=sprintf("%s kN", string(round(loadstep(2,slicetimeid))));
    thisplot=plot3(1e3*squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), 1e3*squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'linewidth',3, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    %%
    figure(1056)
    hold on
    thisplot=plot3(1e3*squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), 1e3*squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'Color',thiscolor,'linewidth',3, 'DisplayName', captionshort);
    % thiscolor=get(thisplot, 'Color');
    for isample=1:n_samples
        r_all_isample=squeeze(s_r_all_NSCRK(isample,:,:));
        D_ik_all_isample=squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
        [r_all_isample,D_ik_all_isample]=rotate_in_plane_rod(r_all_isample,D_ik_all_isample,rod_axis_rs);
        plot3(squeeze(1e3*r_all_isample(1,1:plotinterval:n_points)), 1e3*squeeze(r_all_isample(2,1:plotinterval:n_points)), squeeze(r_all_isample(3,1:plotinterval:n_points)),'Color',[thiscolor,0.01],'linewidth',1,'HandleVisibility','off');
    end
end
close(vw_v3)
close(vw_v2)
close(vw)
close(vw_v4)
%%
positionvec=[618,47,640,1292];
figure(1005)
set(gcf,"Position",positionvec)
hold on
xlabel("$x$ [mm]",'interpreter',"latex",'fontsize',14)
ylabel("$y$ [mm]",'interpreter',"latex",'fontsize',14)
zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
legend
% axis equal
xlim(1e3*[-4.2e-3,4.2e-3])
ylim(1e3*[-4.2e-3,4.2e-3])
zlim([-12,0])
grid on
view(-51.9,10.8108)
title({"T"+num2str(Tnumber)+" 3D deformation"}, 'FontSize', 15);
legend('Location', 'Best')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
saveas(gcf,"andrewdatarevisit/"+ "T"+num2str(Tnumber)+"AY3Dstaticplot.png")
view(0,0)
title({"T"+num2str(Tnumber)+" 2D deformation plane xz"}, 'FontSize', 15);
legend('Location', 'Best')
saveas(gcf,"andrewdatarevisit/"+ "T"+num2str(Tnumber)+"AY2Dstaticplot xz.png")
view(90,0)
title({"T"+num2str(Tnumber)+" 2D deformation plane yz"}, 'FontSize', 15);
ylim(1e3*[-4.2e-3,4.2e-3])
legend('Location', 'Best')
saveas(gcf,"andrewdatarevisit/"+ "T"+num2str(Tnumber)+"AY2Dstaticplot yz.png")
view(0,90)
xlim(1e3*[-4.2e-3,4.2e-3])
ylim(1e3*[-4.2e-3,4.2e-3])
zlim([-12,0])
title({"T"+num2str(Tnumber)+" 2D deformation plane xy"}, 'FontSize', 15);
legend('Location', 'Best')
set(gcf,"Position",[965,703,630,472])
saveas(gcf,"andrewdatarevisit/"+ "T"+num2str(Tnumber)+"AY2Dstaticplot xy.png")
%%
axislim=2e-3;
figure(1056)
set(gcf,"Position",positionvec)
hold on
xlabel("$x$ [mm]",'interpreter',"latex",'fontsize',14)
ylabel("$y$ [mm]",'interpreter',"latex",'fontsize',14)
zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
legend
% axis equal
xlim(1e3*[-axislim,axislim])
ylim(1e3*[-axislim,axislim])
zlim([-12,0])
grid on
view(-51.9,10.8108)
% view(-17.0754,10.3311)
title({"T"+num2str(Tnumber)+" 3D deformation-MC"}, 'FontSize', 15);
legend('Location', 'Best')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
saveas(gcf,"andrewdatarevisit/"+ "T"+num2str(Tnumber)+"AY3DstaticplotMC.png")
view(0,0)
title({"T"+num2str(Tnumber)+" 2D deformation-MC-plane xz"}, 'FontSize', 15);
legend('Location', 'Best')
saveas(gcf,"andrewdatarevisit/"+ "T"+num2str(Tnumber)+"AY2DstaticplotMC xz.png")
view(90,0)
title({"T"+num2str(Tnumber)+" 2D deformation-MC-plane yz"}, 'FontSize', 15);
legend('Location', 'Best')
saveas(gcf,'AY2DstaticplotMC yz.png')
view(0,90)
title({"T"+num2str(Tnumber)+" 2D deformation-MC-plane xy"}, 'FontSize', 15);
legend('Location', 'northeast')
set(gcf,"Position",[965,703,630,472])
saveas(gcf,"andrewdatarevisit/"+ "T"+num2str(Tnumber)+"AY2DstaticplotMC xy.png")
%%

trajectoryendpointini=trajectoryendpoint(:,1);
for i=1:3
    for t=1:length(FrameTimes)
        trajectoryendpoint(i,t)=trajectoryendpoint(i,t)-trajectoryendpointini(i);
    end
end
vw = VideoWriter("AYPILE_trajectory","Motion JPEG AVI");
vw.FrameRate=frameratethisone;
open(vw)
k_frame_inter=1;
n_points=length(s_all);
figure(1002)
hold on
for k = 1:k_frame_inter:length(FrameTimes)
    slicetimeid=k
    figure(1002)
    plot3(trajectoryendpoint(1,1:slicetimeid),trajectoryendpoint(2,1:slicetimeid),trajectoryendpoint(3,1:slicetimeid),LineWidth=3.0,Color='r')
    grid on
    axis equal
    hold on
    hold off
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    xlim([-2.5e-3,2.5e-3])
    ylim([-2.5e-3,2.5e-3])
    zlim([-2.5e-3,2.5e-3])
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(FrameTimes), string(FrameTimes(slicetimeid)));
    title({"Trajectory of the Selected End Point";caption}, 'FontSize', 15);
    % view(488.2878-5*slicetimeid,37+0*slicetimeid)
    frame = getframe(gcf);
    writeVideo(vw,frame)

end
close(vw)
figure(1022)
plot(loadstep(2,plotframe(1:end)),trajectoryendpoint(3,plotframe(1:end)),'-or',LineWidth=1.5,MarkerSize=10.0)
xlabel("Pile load $F$ [kN]",'interpreter',"latex",'fontsize',14)
ylabel({"Vertical displacement $\Delta r_z$ [m]"},'interpreter',"latex",'fontsize',14)
title("Vertical displacement - load", 'FontSize', 15)
saveas(gcf,'VerticaldisplacementloadAYpile')
%%
figure(1007)
plotinterval=10;
plot(s,real(MSE_all_NSCRK),'-or',LineWidth=1.5,MarkerSize=10.0,MarkerIndices=1:plotinterval:n_points)
xlabel("Rod Coordinate $s$ [m]",'interpreter',"latex",'fontsize',14)
ylabel({"Mean Squared Deviation of ";"the Rod Axis Estimation $MSD$ [m]"},'interpreter',"latex",'fontsize',14)

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
