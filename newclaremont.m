close all
clear
clc
% load pipe_clip.mat
%% settings
%%
scaleeff=1.0;
n_gauge_points=1973;
gaugepitch=0.1;
oversampling=10;
sample_spacing=gaugepitch/oversampling;
n_points=1+(n_gauge_points-1)*oversampling;
plot_number=200;
measure_error_level=1.0*scaleeff*40e-6*sqrt(oversampling);
n_samples=100;
n_substep=1;
strain_N_coarse=readmatrix("/home/srgnsfdealii/Desktop/claremont_shape_sensing/Pipe_Strain_N_pipe_temp_comp.csv")*1e-6;
Times = readcell("/home/srgnsfdealii/Desktop/claremont_shape_sensing/Pipe_Strain_N_pipe_temp_comp.csv",'range',"1:1");
strain_N_coarse=strain_N_coarse(2:end,:);
strain_S_coarse=readmatrix("/home/srgnsfdealii/Desktop/claremont_shape_sensing/Pipe_Strain_S_pipe_temp_comp.csv")*1e-6;
strain_S_coarse=strain_S_coarse(2:end,:);
strain_Z_coarse=readmatrix("/home/srgnsfdealii/Desktop/claremont_shape_sensing/Pipe_Strain_Z_pipe_temp_comp.csv")*1e-6;
strain_Z_coarse=strain_Z_coarse(2:end,:);
s_coarse=linspace(0, gaugepitch*(n_gauge_points-1), n_gauge_points)';
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
p=0.1;
cart_coord_x_sp = csaps(s,cart_coord_x,p,s);
cart_coord_y_sp = csaps(s,cart_coord_y,p,s);
cart_coord_z_sp = csaps(s,cart_coord_z,p,s);
cart_coord_x_sp = cart_coord_x_sp-cart_coord_x_sp(1);
cart_coord_y_sp = cart_coord_y_sp-cart_coord_y_sp(1);
cart_coord_z_sp = cart_coord_z_sp-cart_coord_z_sp(1);
figure
plot3(cart_coord_x_sp,cart_coord_y_sp,cart_coord_z_sp)
%%
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
[~, fault_index] = min(abs(s-66));
strike=130*pi/180;
dip=0;
%%
% reconstruct to keep consistency
cables_strain_all=scaleeff*[strain_N(:,1)';strain_S(:,1)';strain_Z(:,1)'];
plot_figure_id=0;
[r_all_ANASDE0,D_ik_all_ANASDE0,D0_ik_all_ANASDE0,u0_k_all_ANASDE0,v0_k_all_ANASDE0,u_k_all_ANASDE0,v_k_all_ANASDE0,Cu_k_all_ANASDE0,Cv_k_all_ANASDE0,MSE_all_ANASDE0]...
    =Cosserat_WeakRK_Samples_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep,1);
%% Initialize Media Output Directory
out_dir = 'outputmedia';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
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
% Define video writers for outputmedia folder
vw_v3 = VideoWriter(fullfile(out_dir, "Axial_elongation.avi"), "Uncompressed AVI");
vw_v3.FrameRate=1/0.6; open(vw_v3);
vw = VideoWriter(fullfile(out_dir, "pipe.avi"), "Uncompressed AVI");
vw.FrameRate=1/0.6; open(vw);
vw_para = VideoWriter(fullfile(out_dir, "Displacement_Parallel.avi"), "Uncompressed AVI");
vw_para.FrameRate=1/0.6; open(vw_para);
vw_perp = VideoWriter(fullfile(out_dir, "Displacement_Perpendicular.avi"), "Uncompressed AVI");
vw_perp.FrameRate=1/0.6; open(vw_perp);

% NEW: Video writers for Pipe-local displacements
vw_pipe_para = VideoWriter(fullfile(out_dir, "Displacement_Pipe_Parallel_d3.avi"), "Uncompressed AVI");
vw_pipe_para.FrameRate=1/0.6; open(vw_pipe_para);
vw_pipe_perp = VideoWriter(fullfile(out_dir, "Displacement_Pipe_Perpendicular_d1.avi"), "Uncompressed AVI");
vw_pipe_perp.FrameRate=1/0.6; open(vw_pipe_perp);

k_frame_inter=1;
n_points=length(s_all);
all_points_history=zeros(3,n_points,length(Times));
all_axis_history=zeros(3,3,n_points,length(Times));
for k = [1:1:length(Times)-8]
    slicetimeid=k
    cables_strain_all=scaleeff*[strain_N(:,slicetimeid)';strain_S(:,slicetimeid)';strain_Z(:,slicetimeid)'];
    plot_figure_id=0;
    
    [r_all_NSCRK,D_ik_all_NSCRK,D0_ik_all_NSCRK,u0_k_all_NSCRK,v0_k_all_NSCRK,u_k_all_NSCRK,v_k_all_NSCRK,Cu_k_all_NSCRK,Cv_k_all_NSCRK,MSE_all_NSCRK,s_r_all_NSCRK,s_D_ik_all_NSCRK]...
        =Cosserat_WeakRK_Samples_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s_all,rod_axis_rs,isassigneddir,assigned_dir,cables_strain_all,measure_error_level,plot_figure_id,n_substep,n_samples);
    [r_all_NSCRK,D_ik_all_NSCRK]=rotate_at_point(r_all_NSCRK,D_ik_all_NSCRK,r_all_ANASDE0,D_ik_all_ANASDE0,fault_index);
    % [r_all_NSCRK,D_ik_all_NSCRK]=rotate_in_plane_rod(r_all_NSCRK,D_ik_all_NSCRK,r_all_ANASDE0);
    all_points_history(:,:,slicetimeid)=r_all_NSCRK;
    all_axis_history(:,:,:,slicetimeid)=D_ik_all_NSCRK;
    
    %% Figure 1001
    figure(1001)
    clf(1001) % Clear figure to only keep current frame
    hold on
    plotinterval=ceil(n_points/100);
    plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'m');
    grid on
    axis equal
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
    print(fullfile(out_dir, ['claremontquiver' caption '.png']),'-dpng')
    
    frame = getframe(gcf);
    if size(frame.cdata, 1) ~= 540 || size(frame.cdata, 2) ~= 960
        frame.cdata = imresize(frame.cdata, [540, 960]); 
    end
    writeVideo(vw,frame)
    
    %% Figure 1004
    figure(1004)
    clf(1004) % Clear figure
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
    axis equal
    xlim([-183.2869283785117 6.341750763666314])
    ylim([-14.757046048603588 53.488444303517625])
    zlim([-23.955800875800776 0.990397380988718])
    clim([-2e-4,2e-4])
    view(0,90)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title({"Axial elongation [1]";caption}, 'FontSize', 15);
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(fullfile(out_dir, ['claremontv3' caption '.png']),'-dpng')
    
    frame = getframe(gcf);
    if size(frame.cdata, 1) ~= 540 || size(frame.cdata, 2) ~= 960
        frame.cdata = imresize(frame.cdata, [540, 960]);
    end
    writeVideo(vw,frame)
    writeVideo(vw_v3,frame)
    
    %% Figure 1005
    figure(1005)
    clf(1005)
    hold on
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    caption = sprintf('Time = %s', string(currenttime));
    captionshort=sprintf('%s', string(currenttime));
    thisplot=plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'linewidth',3, 'DisplayName', captionshort);
    view(-12.283639988412114,39.810373502729341)
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    legend
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(fullfile(out_dir, ['claremontcloud3D' caption '.png']),'-dpng')
    
    %% Figure 10051
    figure(10051)
    clf(10051)
    hold on
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    caption = sprintf('Time = %s', string(currenttime));
    captionshort=sprintf('%s', string(currenttime));
    thisplot=plot(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)),'linewidth',3, 'DisplayName', captionshort);
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
    zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
    axis equal
    legend
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(fullfile(out_dir, ['claremontcloud2D' caption '.png']),'-dpng')
            
    %% Figure 10055
    figure(10055)
    clf(10055)
    hold on
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    caption = sprintf('Time = %s', string(currenttime));
    captionshort=sprintf('%s', string(currenttime));
    thisplot=plot(s,v_k_all_NSCRK(3,:),'linewidth',3, 'DisplayName', captionshort);
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    xlabel("Points",'interpreter',"latex",'fontsize',14)
    ylabel("$v_3$ [1]",'interpreter',"latex",'fontsize',14)
    legend
    set(gca, 'XDir', 'reverse');
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(fullfile(out_dir, ['elongation' caption '.png']),'-dpng')
    
    %% Figure 10056
    figure(10056)
    clf(10056)
    hold on
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    caption = sprintf('Time = %s', string(currenttime));
    captionshort=sprintf('%s', string(currenttime));
    thisplot=plot(s,1.0+mean(cables_strain_all,1),'linewidth',3, 'DisplayName', captionshort);
    caption = sprintf('Frame #%d of %d, t = %s', slicetimeid, length(Times), string(Times(slicetimeid)));
    title(caption, 'FontSize', 15);
    xlabel("Points",'interpreter',"latex",'fontsize',14)
    ylabel("$avgstrain+1$ [1]",'interpreter',"latex",'fontsize',14)
    legend
    set(gca, 'XDir', 'reverse');
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(fullfile(out_dir, ['strain' caption '.png']),'-dpng')
    
    %% Figure 1056
    figure(1056)
    clf(1056)
    hold on
    thisplot=plot3(squeeze(r_all_NSCRK(1,1:plotinterval:n_points)), squeeze(r_all_NSCRK(2,1:plotinterval:n_points)), squeeze(r_all_NSCRK(3,1:plotinterval:n_points)),'linewidth',1, 'DisplayName', captionshort);
    thiscolor=get(thisplot, 'Color');
    for isample=1:n_samples
        r_all_isample=squeeze(s_r_all_NSCRK(isample,:,:));
        D_ik_all_isample=squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
        [r_all_isample,D_ik_all_isample]=rotate_at_point(r_all_isample,D_ik_all_isample,r_all_ANASDE0,D_ik_all_ANASDE0,fault_index);
        % [r_all_isample,D_ik_all_isample]=rotate_in_plane_rod(r_all_isample,D_ik_all_isample,r_all_ANASDE0);
        plot3(squeeze(r_all_isample(1,1:plotinterval:n_points)), squeeze(r_all_isample(2,1:plotinterval:n_points)), squeeze(r_all_isample(3,1:plotinterval:n_points)),'Color',[thiscolor,0.01],'linewidth',1,'HandleVisibility','off');
    end
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
    print(fullfile(out_dir, ['claremontcloud3DwithMC' caption '.png']),'-dpng')
    
    %% Fault vectors
    faultvec_para=[cos(strike),sin(strike),0];
    faultvec_perp = [-sin(strike), cos(strike), 0];
    
    % Initial Directors for Pipe Projections
    % D_ik_all_ANASDE0 dimensions: [3, 3, n_points] -> d1 is index 1, d3 is index 3
    d1_init = squeeze(D_ik_all_ANASDE0(:, 1, :)); 
    d3_init = squeeze(D_ik_all_ANASDE0(:, 3, :)); 
    
    %% Calculate Sample Statistics for Shading
    base_projection_para = faultvec_para*(r_all_NSCRK-r_all_ANASDE0);
    base_projection_perp = faultvec_perp*(r_all_NSCRK-r_all_ANASDE0);
    
    % Dot products for pipe-local projections (summing across the x,y,z components)
    base_projection_pipe_para = sum((r_all_NSCRK-r_all_ANASDE0) .* d3_init, 1);
    base_projection_pipe_perp = sum((r_all_NSCRK-r_all_ANASDE0) .* d1_init, 1);
    
    all_samples_para = zeros(n_samples, n_points);
    all_samples_perp = zeros(n_samples, n_points);
    all_samples_pipe_para = zeros(n_samples, n_points);
    all_samples_pipe_perp = zeros(n_samples, n_points);
    
    for isample=1:n_samples
        r_all_isample=squeeze(s_r_all_NSCRK(isample,:,:));
        D_ik_all_isample=squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
        [r_all_isample,~]=rotate_at_point(r_all_isample,D_ik_all_isample,r_all_ANASDE0,D_ik_all_ANASDE0,fault_index);
        % [r_all_isample,~]=rotate_in_plane_rod(r_all_isample,D_ik_all_isample,r_all_ANASDE0);
        all_samples_para(isample,:) = faultvec_para*(r_all_isample-r_all_ANASDE0);
        all_samples_perp(isample,:) = faultvec_perp*(r_all_isample-r_all_ANASDE0);
        
        % Sample dot products for pipe-local projections
        all_samples_pipe_para(isample,:) = sum((r_all_isample-r_all_ANASDE0) .* d3_init, 1);
        all_samples_pipe_perp(isample,:) = sum((r_all_isample-r_all_ANASDE0) .* d1_init, 1);
    end
    
    % Compute standard deviation across all samples (dimension 1)
    std_para = std(all_samples_para, 0, 1);
    std_perp = std(all_samples_perp, 0, 1);
    std_pipe_para = std(all_samples_pipe_para, 0, 1);
    std_pipe_perp = std(all_samples_pipe_perp, 0, 1);
    
    % Prepare X-axis fill vector (forward then backward to create a polygon)
    s_fill = [s; flipud(s)];
    
  %% --- Displacement Parallel to the Fault ---
    
    % 1. Accumulated Figure for saving the final image (1006)
    figure(1006)
    hold on
    thisplot = plot(s, base_projection_para, 'linewidth', 1.5, 'DisplayName', captionshort);
    thiscolor = get(thisplot, 'Color'); 
    
    upper_para = base_projection_para + std_para;
    lower_para = base_projection_para - std_para;
    para_fill = [upper_para, fliplr(lower_para)]';
    fill(s_fill, para_fill, thiscolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    xlim([0 200])
    ylim([-8e-3,8e-3])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    title("Displacement Parallel to the Fault", 'FontSize', 15); % <-- NEW TITLE
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Displacement Parallel to the Fault $u_{\parallel}$ [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(fullfile(out_dir, ['claremontparallel_accumulated' caption '.png']),'-dpng')

    % 2. Single-frame Figure for capturing the video (2006)
    figure(2006)
    clf(2006) 
    hold on
    plot(s, base_projection_para, 'Color', thiscolor, 'linewidth', 1.5, 'DisplayName', captionshort);
    fill(s_fill, para_fill, thiscolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    xlim([0 200])
    ylim([-8e-3,8e-3])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    title({"Displacement Parallel to the Fault"}, 'FontSize', 15); % <-- UPDATED TITLE
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Displacement Parallel to the Fault $u_{\parallel}$ [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    
    frame_para = getframe(gcf);
    if size(frame_para.cdata, 1) ~= 540 || size(frame_para.cdata, 2) ~= 960
        frame_para.cdata = imresize(frame_para.cdata, [540, 960]);
    end
    writeVideo(vw_para, frame_para);
    
    %% --- Displacement Perpendicular to the Fault ---
    % 1. Accumulated Figure for saving the final image (1007)
    figure(1007)
    hold on
    thisplot = plot(s, base_projection_perp, 'linewidth', 1.5, 'DisplayName', captionshort);
    thiscolor = get(thisplot, 'Color'); 
    
    upper_perp = base_projection_perp + std_perp;
    lower_perp = base_projection_perp - std_perp;
    perp_fill = [upper_perp, fliplr(lower_perp)]';
    fill(s_fill, perp_fill, thiscolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    xlim([0 200])
    ylim([-8e-3,8e-3])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    title("Displacement Perpendicular to the Fault", 'FontSize', 15); % <-- NEW TITLE
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Displacement Perpendicular to the Fault $u_{\perp}$ [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(fullfile(out_dir, ['claremontperpendicular_accumulated' caption '.png']),'-dpng')

    % 2. Single-frame Figure for capturing the video (2007)
    figure(2007)
    clf(2007) 
    hold on
    plot(s, base_projection_perp, 'Color', thiscolor, 'linewidth', 1.5, 'DisplayName', captionshort);
    fill(s_fill, perp_fill, thiscolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    xlim([0 200])
    ylim([-8e-3,8e-3])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    title({"Displacement Perpendicular to the Fault"}, 'FontSize', 15); % <-- UPDATED TITLE
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Displacement Perpendicular to the Fault $u_{\perp}$ [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    
    frame_perp = getframe(gcf);
    if size(frame_perp.cdata, 1) ~= 540 || size(frame_perp.cdata, 2) ~= 960
        frame_perp.cdata = imresize(frame_perp.cdata, [540, 960]);
    end
    writeVideo(vw_perp, frame_perp);
    %% --- NEW: Displacement Parallel to the Pipe (d3) ---
    % 1. Accumulated Figure (1008)
    figure(1008)
    hold on
    thisplot = plot(s, base_projection_pipe_para, 'linewidth', 1.5, 'DisplayName', captionshort);
    thiscolor = get(thisplot, 'Color'); 
    
    upper_pipe_para = base_projection_pipe_para + std_pipe_para;
    lower_pipe_para = base_projection_pipe_para - std_pipe_para;
    pipe_para_fill = [upper_pipe_para, fliplr(lower_pipe_para)]';
    fill(s_fill, pipe_para_fill, thiscolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    xlim([0 200])
    ylim([-8e-3,8e-3])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Axial Displacement (d3) [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(fullfile(out_dir, ['claremont_pipe_parallel_accumulated' caption '.png']),'-dpng')

    % 2. Video Figure (2008)
    figure(2008)
    clf(2008) 
    hold on
    plot(s, base_projection_pipe_para, 'Color', thiscolor, 'linewidth', 1.5, 'DisplayName', captionshort);
    fill(s_fill, pipe_para_fill, thiscolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    xlim([0 200])
    ylim([-8e-3,8e-3])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    title({"Displacement Parallel to Pipe (d3)"}, 'FontSize', 15);
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Axial Displacement [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    
    frame_pipe_para = getframe(gcf);
    if size(frame_pipe_para.cdata, 1) ~= 540 || size(frame_pipe_para.cdata, 2) ~= 960
        frame_pipe_para.cdata = imresize(frame_pipe_para.cdata, [540, 960]);
    end
    writeVideo(vw_pipe_para, frame_pipe_para);

    %% --- NEW: Displacement Perpendicular to the Pipe (d1) ---
    % 1. Accumulated Figure (1009)
    figure(1009)
    hold on
    thisplot = plot(s, base_projection_pipe_perp, 'linewidth', 1.5, 'DisplayName', captionshort);
    thiscolor = get(thisplot, 'Color'); 
    
    upper_pipe_perp = base_projection_pipe_perp + std_pipe_perp;
    lower_pipe_perp = base_projection_pipe_perp - std_pipe_perp;
    pipe_perp_fill = [upper_pipe_perp, fliplr(lower_pipe_perp)]';
    fill(s_fill, pipe_perp_fill, thiscolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    xlim([0 200])
    ylim([-8e-3,8e-3])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Transverse Displacement (d1) [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    print(fullfile(out_dir, ['claremont_pipe_perpendicular_accumulated' caption '.png']),'-dpng')

    % 2. Video Figure (2009)
    figure(2009)
    clf(2009) 
    hold on
    plot(s, base_projection_pipe_perp, 'Color', thiscolor, 'linewidth', 1.5, 'DisplayName', captionshort);
    fill(s_fill, pipe_perp_fill, thiscolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    xlim([0 200])
    ylim([-8e-3,8e-3])
    set(gca, 'XDir', 'reverse');
    legend('Location', 'eastoutside')
    title({"Displacement Perpendicular to Pipe (d1)"}, 'FontSize', 15);
    xlabel("Offset from Start point: $s$ [m]",'interpreter',"latex",'fontsize',14)
    ylabel("Transverse Displacement [m]",'interpreter',"latex",'fontsize',14)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 960, 540]);
    
    frame_pipe_perp = getframe(gcf);
    if size(frame_pipe_perp.cdata, 1) ~= 540 || size(frame_pipe_perp.cdata, 2) ~= 960
        frame_pipe_perp.cdata = imresize(frame_pipe_perp.cdata, [540, 960]);
    end
    writeVideo(vw_pipe_perp, frame_pipe_perp);
    
%% --- Optimized Saving Section ---
    folder_name = fullfile(out_dir, "results_excel_" + num2str(scaleeff) + "_scale");
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);
    end
    
    total_rows = n_points * (n_samples + 1);
    big_data_matrix = zeros(total_rows, 7);
    
    idx_start = 1;
    idx_end = n_points;
    big_data_matrix(idx_start:idx_end, 1:3) = r_all_NSCRK.'; 
    big_data_matrix(idx_start:idx_end, 4) = 0; 
    big_data_matrix(idx_start:idx_end, 5)=base_projection_para';
    big_data_matrix(idx_start:idx_end, 6)=base_projection_perp';
    big_data_matrix(idx_start:idx_end, 7)=s;
    
    for isample = 1:n_samples
        r_temp = squeeze(s_r_all_NSCRK(isample,:,:));
        D_temp = squeeze(s_D_ik_all_NSCRK(isample,:,:,:));
        [r_temp, ~] = rotate_at_point(r_temp, D_temp, r_all_ANASDE0, D_ik_all_ANASDE0, fault_index);
        % [r_temp, ~] = rotate_in_plane_rod(r_temp, D_temp, r_all_ANASDE0);
        idx_start = idx_end + 1;
        idx_end = idx_start + n_points - 1;
        big_data_matrix(idx_start:idx_end, 1:3) = r_temp.';
        big_data_matrix(idx_start:idx_end, 4) = isample; 
        
        % Reuse pre-calculated projections
        big_data_matrix(idx_start:idx_end, 5)=all_samples_para(isample,:)';
        big_data_matrix(idx_start:idx_end, 6)=all_samples_perp(isample,:)';
        big_data_matrix(idx_start:idx_end, 7)=s;
    end
    
    currenttime = Times{slicetimeid};
    if iscell(currenttime), currenttime = currenttime{1}; end
    
    filename_csv = sprintf('%s/claremont_frame_%03d_%s.csv', folder_name, slicetimeid, string(currenttime));
    
    T = array2table(big_data_matrix, 'VariableNames', {'x', 'y', 'z', 'SampleID','displacement_parallel','displacement_perpendicular','curve_lin_coord'});
    writetable(T, filename_csv);
    
    disp(['Saved Fast CSV: ', filename_csv]);
end

% Close all VideoWriters
close(vw_v3)
close(vw)
close(vw_para)
close(vw_perp)
close(vw_pipe_para)
close(vw_pipe_perp)