close all
clear
clc
% load pipe_clip.mat
%% settings
%%
scaleeff=1.0;
n_gauge_points=1983;
gaugepitch=0.1;
oversampling=100;
sample_spacing=gaugepitch/oversampling;
n_points=1+(n_gauge_points-1)*oversampling;
plot_number=200;
measure_error_level=1.0*scaleeff*0.001e-6*sqrt(oversampling);
n_samples=1;
n_substep=1;
strain_N_coarse=readmatrix("/home/srgnsfdealii/Desktop/claremont_shape_sensing/Pipe_Strain_N_pipe_temp_comp.csv")*1e-6;
Times = readcell("/home/srgnsfdealii/Desktop/claremont_shape_sensing/Pipe_Strain_N_pipe_temp_comp.csv",'range',"1:1");
strain_N_coarse=strain_N_coarse(2:end,:);
strain_S_coarse=readmatrix("/home/srgnsfdealii/Desktop/claremont_shape_sensing/Pipe_Strain_S_pipe_temp_comp.csv")*1e-6;
strain_S_coarse=strain_S_coarse(2:end,:);
strain_Z_coarse=readmatrix("/home/srgnsfdealii/Desktop/claremont_shape_sensing/Pipe_Strain_Z_pipe_temp_comp.csv")*1e-6;
strain_Z_coarse=strain_Z_coarse(2:end,:);
s_coarse=linspace(0, gaugepitch*(n_gauge_points-1), n_gauge_points)';

% strain_N=0.0*(strain_N);
% strain_S=0.0*(strain_S);
% strain_Z=0.0*(strain_Z);
% strain_N_coarse(550:850,2:end)=16.66*strain_N_coarse(550:850,2:end);
% strain_S_coarse(550:850,2:end)=16.66*strain_N_coarse(550:850,2:end);
% strain_Z_coarse(550:850,2:end)=16.66*strain_N_coarse(550:850,2:end);
% strain_N_coarse(:,2:end)=0.01;
% strain_S_coarse(:,2:end)=0.01;
% strain_Z_coarse(:,2:end)=0.01;
% strain_N(:,1)=0;
% strain_S(:,1)=0;
% strain_Z(:,1)=0;
plot_index=1:ceil(n_points/plot_number):n_points;
s = linspace(0, sample_spacing*(n_points-1), n_points)';
cart_coords_coarse=readmatrix("/home/srgnsfdealii/Desktop/claremont_shape_sensing/pipe_coords_xyz_36pipe.csv");
strain_N=interp1(s_coarse,strain_N_coarse,s,"spline");
strain_S=interp1(s_coarse,strain_S_coarse,s,"spline");
strain_Z=interp1(s_coarse,strain_Z_coarse,s,"spline");
cart_coords=interp1(s_coarse,cart_coords_coarse,s,"spline");

cart_coords=cart_coords(:,:);
cart_coord_x = cart_coords(:,1);
cart_coord_y = cart_coords(:,2);
cart_coord_z = cart_coords(:,3);

figure
plot3(cart_coord_x,cart_coord_y,cart_coord_z)
cart_coord_x = cart_coord_x-cart_coord_x(1);
cart_coord_y = cart_coord_y-cart_coord_y(1);
cart_coord_z = cart_coord_z-cart_coord_z(1);
% cart_coord_x=flip(cart_coord_x);
% % cart_coord_x=-cart_coord_x;
% cart_coord_y=flip(cart_coord_y);
% cart_coord_z=flip(cart_coord_z);
% strain_N=flip(strain_N);
% strain_S=flip(strain_S);
% strain_Z=flip(strain_Z);
p=0.1;
cart_coord_x_sp = csaps(s,cart_coord_x,p,s);
cart_coord_y_sp = csaps(s,cart_coord_y,p,s);
cart_coord_z_sp = csaps(s,cart_coord_z,p,s);
cart_coord_x_sp = cart_coord_x_sp-cart_coord_x_sp(1);
cart_coord_y_sp = cart_coord_y_sp-cart_coord_y_sp(1);
cart_coord_z_sp = cart_coord_z_sp-cart_coord_z_sp(1);
% cart_coord_x_sp=cos(2.0*pi*s/max(s));
% cart_coord_y_sp=sin(2.0*pi*s/max(s));
% cart_coord_z_sp=0*s;
figure
plot3(cart_coord_x_sp,cart_coord_y_sp,cart_coord_z_sp)
%%
% s_all=[0:1:100];
d_cart_coord_x_sp=gradient(cart_coord_x_sp);
d_cart_coord_y_sp=gradient(cart_coord_y_sp);
d_cart_coord_z_sp=gradient(cart_coord_z_sp);
segment_lengths = sqrt(d_cart_coord_x_sp.^2 + d_cart_coord_y_sp.^2 + d_cart_coord_z_sp.^2);
s = cumsum(segment_lengths);
s=s-s(1);
s_all=s';
rod_axis_rs=[cart_coord_x_sp';cart_coord_y_sp';cart_coord_z_sp'];
isassigneddir=false;
assigned_dir=[1,0,0];
get_frenet_frame_from_cart(s_all,rod_axis_rs,isassigneddir,assigned_dir);

%%
R_cables=[0.45,0.45,0.45];
%NSZ theta should be minus
theta_cables=[1*pi/4,pi*3/4,pi/2];
omega_cables=[0,0,0];
n_cables=3;
[~, fault_index] = min(abs(s-70));
strike=130*pi/180;
dip=0;
%%
% reconstruct to keep consistency
cables_strain_all=scaleeff*[strain_N(:,1)';strain_S(:,1)';strain_Z(:,1)'];
plot_figure_id=0;
% [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep);
% [r_all_NFCA,D_ik_all_NFCA]=rotate_in_plane_rod(r_all_NFCA,D_ik_all_NFCA,rod_axis_rs);
% [r_all_ANASDE0,D_ik_all_ANASDE0,D0_ik_all_ANASDE0,u0_k_all_ANASDE0,v0_k_all_ANASDE0,u_k_all_ANASDE0,v_k_all_ANASDE0,Cu_k_all_ANASDE0,Cv_k_all_ANASDE0,MSE_all_ANASDE0]...
%     =Cosserat_analytical_SDE_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,0*cables_strain_all,measure_error_level,plot_figure_id,n_substep);
[r_all_ANASDE0,D_ik_all_ANASDE0,D0_ik_all_ANASDE0,u0_k_all_ANASDE0,v0_k_all_ANASDE0,u_k_all_ANASDE0,v_k_all_ANASDE0,Cu_k_all_ANASDE0,Cv_k_all_ANASDE0,MSE_all_ANASDE0]...
    =Cosserat_WeakRK_Samples_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep,1);
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


figure(1004)

hold on
figure(1001)

hold on
figure(1005)

figure(10051)
plot(squeeze(rod_axis_rs(1,1:1:n_points)), squeeze(rod_axis_rs(2,1:1:n_points)),'linewidth',3, 'DisplayName', "initial");
hold on
plot(squeeze(r_all_ANASDE0(1,1:1:n_points)), squeeze(r_all_ANASDE0(2,1:1:n_points)),'linewidth',3, 'DisplayName', "reconstruct");
figure(10052)
axis equal
plot3(squeeze(rod_axis_rs(1,1:1:n_points)), squeeze(rod_axis_rs(2,1:1:n_points)), squeeze(rod_axis_rs(3,1:1:n_points)),'linewidth',3, 'DisplayName', "initial");
hold on
plot3(squeeze(r_all_ANASDE0(1,1:1:n_points)), squeeze(r_all_ANASDE0(2,1:1:n_points)), squeeze(r_all_ANASDE0(3,1:1:n_points)),'linewidth',3, 'DisplayName', "reconstruct");
axis equal
figure(1056)

hold on
vw_v3 = VideoWriter("Axial elongation","Uncompressed AVI");
vw_v3.FrameRate=1/0.6;
open(vw_v3)
vw = VideoWriter("pipe","Uncompressed AVI");
vw.FrameRate=1/0.6;
open(vw)
k_frame_inter=1;
n_points=length(s_all);
all_points_history=zeros(3,n_points,length(Times));
all_axis_history=zeros(3,3,n_points,length(Times));
for k = [1:1:length(Times)-1 length(Times)]
    slicetimeid=k
    cables_strain_all=scaleeff*[strain_N(:,slicetimeid)';strain_S(:,slicetimeid)';strain_Z(:,slicetimeid)'];
    plot_figure_id=0;
    % [r_all_NFCA,D_ik_all_NFCA,D0_ik_all_NFCA,u0_k_all_NFCA,v0_k_all_NFCA,u_k_all_NFCA,v_k_all_NFCA]=Cosserat_noisefree_approx_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep);
    % [r_all_NFCA,D_ik_all_NFCA]=rotate_in_plane_rod(r_all_NFCA,D_ik_all_NFCA,rod_axis_rs);
    [r_all_NSCRK,D_ik_all_NSCRK,D0_ik_all_NSCRK,u0_k_all_NSCRK,v0_k_all_NSCRK,u_k_all_NSCRK,v_k_all_NSCRK,Cu_k_all_NSCRK,Cv_k_all_NSCRK,MSE_all_NSCRK,s_r_all_NSCRK,s_D_ik_all_NSCRK]...
        =Cosserat_WeakRK_Samples_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
    [r_all_NSCRK,D_ik_all_NSCRK]=rotate_at_point(r_all_NSCRK,D_ik_all_NSCRK,r_all_ANASDE0,D_ik_all_ANASDE0,fault_index);
    all_points_history(:,:,slicetimeid)=r_all_NSCRK;
    all_axis_history(:,:,:,slicetimeid)=D_ik_all_NSCRK;
    %%
    figure(1001)
    plotinterval=ceil(n_points/100);
    plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'m');
    grid on
    axis equal
    hold on
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,1,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,1,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,1,1:plotinterval:n_points))','r')
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,2,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,2,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,2,1:plotinterval:n_points))','g')
    quiver3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),squeeze(D_ik_all_NSCRK(1,3,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(2,3,1:plotinterval:n_points))', squeeze(D_ik_all_NSCRK(3,3,1:plotinterval:n_points))','b')
    hold off
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    xlim([-183.2869283785117 6.341750763666314])
    ylim([-14.757046048603588 53.488444303517625])
    zlim([-23.955800875800776 0.990397380988718])
    view(-12.283639988412114,39.810373502729341)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(['claremontquiver' caption '.png'],'-dpng')
    frame = getframe(gcf);
    % FORCE RESIZE: Ensure the captured image is exactly 960x540
    if size(frame.cdata, 1) ~= 960 || size(frame.cdata, 2) ~= 540
        frame.cdata = imresize(frame.cdata, [960, 540]);
    end
    writeVideo(vw,frame)
    %%
    figure(1004)
    col=v_k_all_NSCRK(3,1:plotinterval:n_points)-1;
    clf(1004)
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
    axis equal
    xlim([-183.2869283785117 6.341750763666314])
    ylim([-14.757046048603588 53.488444303517625])
    zlim([-23.955800875800776 0.990397380988718])
    clim([-2e-4,2e-4])
    view(0,90)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title({"Axial elongation [1]";caption}, 'FontSize', 15);
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(['claremontv3' caption '.png'],'-dpng')
    frame = getframe(gcf);
        % FORCE RESIZE: Ensure the captured image is exactly 960x540
    if size(frame.cdata, 1) ~= 960 || size(frame.cdata, 2) ~= 540
        frame.cdata = imresize(frame.cdata, [960, 540]);
    end
    writeVideo(vw,frame)
    writeVideo(vw_v3,frame)
    %%
    figure(1005)
    
    hold on
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    % currenttime.Format = 'yyyy-MM';
    caption = sprintf('Time = %s', string(currenttime));
    captionshort=sprintf('%s', string(currenttime));
    thisplot=plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'linewidth',3, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    % xlim([-183.2869283785117 6.341750763666314])
    % ylim([-14.757046048603588 53.488444303517625])
    % zlim([-23.955800875800776 0.990397380988718])
    view(-12.283639988412114,39.810373502729341)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    legend
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(['claremontcloud3D' caption '.png'],'-dpng')
        %%
    figure(10051)

    hold on
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    % currenttime.Format = 'yyyy-MM';
    caption = sprintf('Time = %s', string(currenttime));
    captionshort=sprintf('%s', string(currenttime));
    thisplot=plot(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)),'linewidth',3, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    % xlim([-183.2869283785117 6.341750763666314])
    % ylim([-14.757046048603588 53.488444303517625])
    % zlim([-23.955800875800776 0.990397380988718])
    % view(-12.283639988412114,39.810373502729341)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    legend
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(['claremontcloud2D' caption '.png'],'-dpng')
            %%
    figure(10055)

    hold on
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    % currenttime.Format = 'yyyy-MM';
    caption = sprintf('Time = %s', string(currenttime));
    captionshort=sprintf('%s', string(currenttime));
    thisplot=plot(s,v_k_all_NSCRK(3,:),'linewidth',3, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    % xlim([-183.2869283785117 6.341750763666314])
    % ylim([-14.757046048603588 53.488444303517625])
    % zlim([-23.955800875800776 0.990397380988718])
    % view(-12.283639988412114,39.810373502729341)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    xlabel("Points",'interpreter',"latex",'fontsize',14)
    ylabel("$v_3$ [1]",'interpreter',"latex",'fontsize',14)
    legend
    set(gca, 'XDir', 'reverse');
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(['elongation' caption '.png'],'-dpng')
    %%
                %%
    figure(10056)

    hold on
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    % currenttime.Format = 'yyyy-MM';
    caption = sprintf('Time = %s', string(currenttime));
    captionshort=sprintf('%s', string(currenttime));
    thisplot=plot(s,1.0+mean(cables_strain_all,1),'linewidth',3, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    % xlim([-183.2869283785117 6.341750763666314])
    % ylim([-14.757046048603588 53.488444303517625])
    % zlim([-23.955800875800776 0.990397380988718])
    % view(-12.283639988412114,39.810373502729341)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    xlabel("Points",'interpreter',"latex",'fontsize',14)
    ylabel("$avgstrain+1$ [1]",'interpreter',"latex",'fontsize',14)
    legend
    set(gca, 'XDir', 'reverse');
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(['strain' caption '.png'],'-dpng')
    %%
    figure(1056)
    hold on
    thisplot=plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'Color',thiscolor,'linewidth',1, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    for isample=1:n_samples
        r_all_isample=squeeze(s_r_all_NSCRK(isample,:,:));
        D_ik_all_isample=squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
        [r_all_isample,D_ik_all_isample]=rotate_at_point(r_all_isample,D_ik_all_isample,r_all_ANASDE0,D_ik_all_ANASDE0,fault_index);
        plot3(squeeze(r_all_isample(1,1:plotinterval:n_points)), squeeze(r_all_isample(2,1:plotinterval:n_points)), squeeze(r_all_isample(3,1:plotinterval:n_points)),'Color',[thiscolor,0.01],'linewidth',1,'HandleVisibility','off');
    end
    % plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'linewidth',3, 'DisplayName', captionshort);
    xlim([-183.2869283785117 6.341750763666314])
    ylim([-14.757046048603588 53.488444303517625])
    zlim([-23.955800875800776 0.990397380988718])
    view(-12.283639988412114,39.810373502729341)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    legend
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(['claremontcloud3DwithMC' caption '.png'],'-dpng')
    %%
    faultvec_para=[cos(strike),sin(strike),0];
    faultvec_perp = [-sin(strike), cos(strike), 0];
    %%
    figure(1006)
    hold on
    projectiondisplacement_parallel=faultvec_para*(r_all_NSCRK-r_all_ANASDE0);
    thisplot=plot(s,projectiondisplacement_parallel,'Color',thiscolor,'linewidth',1, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    for isample=1:n_samples
        r_all_isample=squeeze(s_r_all_NSCRK(isample,:,:));
        D_ik_all_isample=squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
        [r_all_isample,D_ik_all_isample]=rotate_at_point(r_all_isample,D_ik_all_isample,r_all_ANASDE0,D_ik_all_ANASDE0,fault_index);

        projectiondisplacement_parallel=faultvec_para*(r_all_isample-r_all_ANASDE0);
        plot(s,projectiondisplacement_parallel,'Color',[thiscolor,0.01],'linewidth',1,'HandleVisibility','off');
    end
    
    xlim([0 200])
    % ylim([-1 1])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Displacement Parallel to the Fault $u_{\parallel}$ [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(['claremontparallel' caption '.png'],'-dpng')
    %%
    figure(1007)
    hold on
    projectiondisplacement_perp=faultvec_perp*(r_all_NSCRK-r_all_ANASDE0);
    thisplot=plot(s,projectiondisplacement_perp,'Color',thiscolor,'linewidth',1, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    for isample=1:n_samples
        r_all_isample=squeeze(s_r_all_NSCRK(isample,:,:));
        D_ik_all_isample=squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
        [r_all_isample,D_ik_all_isample]=rotate_at_point(r_all_isample,D_ik_all_isample,r_all_ANASDE0,D_ik_all_ANASDE0,fault_index);

        projectiondisplacement_perp=faultvec_perp*(r_all_isample-r_all_ANASDE0);
        plot(s,projectiondisplacement_perp,'Color',[thiscolor,0.01],'linewidth',1,'HandleVisibility','off');
    end
    xlim([0 200])
    % ylim([-1 1])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Displacement Perpendicular to the Fault $u_{\perp}$ [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(['claremontperpendicular' caption '.png'],'-dpng')

%% --- Optimized Saving Section ---
    folder_name = "results_excel_" + num2str(scaleeff) + "_scale";
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);
    end
    
    % 1. Pre-allocate a BIG matrix to hold all data for this frame in RAM
    %    Columns: [x, y, z, SampleID]
    %    Rows: n_points * (n_samples + 1 base)
    % total_rows = n_points * (n_samples + 1);
    total_rows = n_points;
    big_data_matrix = zeros(total_rows, 7);
    
    % 2. Fill the Base Curve (ID = 0)
    %    Transpose to (n_points x 3)
    idx_start = 1;
    idx_end = n_points;
    big_data_matrix(idx_start:idx_end, 1:3) = r_all_NSCRK.'; 
    big_data_matrix(idx_start:idx_end, 4) = 0; % ID 0 marks the Base
    big_data_matrix(idx_start:idx_end, 5)=projectiondisplacement_parallel';
    big_data_matrix(idx_start:idx_end, 6)=projectiondisplacement_perp';
    big_data_matrix(idx_start:idx_end, 7)=s;
    
    % % 3. Loop for Calculation Only (Fast in Memory)
    % for isample = 1:n_samples
    %     % Retrieve and Rotate
    %     r_temp = squeeze(s_r_all_NSCRK(isample,:,:));
    %     D_temp = squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
    %     [r_temp, ~] = rotate_at_point(r_temp, D_temp, r_all_ANASDE0, D_ik_all_ANASDE0, fault_index);
    % 
    %     % Calculate indices for the big matrix
    %     idx_start = idx_end + 1;
    %     idx_end = idx_start + n_points - 1;
    % 
    %     % Store in Big Matrix
    %     big_data_matrix(idx_start:idx_end, 1:3) = r_temp.';
    %     big_data_matrix(idx_start:idx_end, 4) = isample; % Store Sample ID
    % end

    % 4. Generate Filename
    currenttime = Times{slicetimeid};
    if iscell(currenttime), currenttime = currenttime{1}; end
    
    % --- OPTION A: Save as CSV (Fastest, ~0.1 seconds) ---
    % Excel opens .csv files automatically. This is much faster than .xlsx
    filename_csv = sprintf('%s/claremont_frame_%03d_%s.csv', folder_name, slicetimeid, string(currenttime));
    
    % Create a table with headers so it looks nice in Excel
    T = array2table(big_data_matrix, 'VariableNames', {'x', 'y', 'z', 'SampleID','displacement_parallel','displacement_perpendicular','curve_lin_coord'});
    writetable(T, filename_csv);
    
    disp(['Saved Fast CSV: ', filename_csv]);
end
close(vw_v3)
close(vw)
