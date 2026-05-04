close all
clear
clc
% n_points=2012;
% n_points=1953;
sample_spacing=0.1;
plot_number=200;
measure_error_level=20e-6;
n_samples=10;
n_substep=1;
cart_coords=readmatrix("claremont_data/22pipe/pipe_coords_xyz_22pipe.csv");
% strain_N=readmatrix("data/newdata/ClaremontNorth36_N_tc.csv");
filestrain="tc_pipe";
strain_N=readmatrix("claremont_data/22pipe/ClaremontSouth22_N_"+filestrain+".csv");
Times = readcell("claremont_data/22pipe/ClaremontSouth22_N_"+filestrain+".csv",'range',"1:1");
strain_N=strain_N(2:end,:);
strain_S=readmatrix("claremont_data/22pipe/ClaremontSouth22_S_"+filestrain+".csv");
strain_S=strain_S(2:end,:);
strain_Z=readmatrix("claremont_data/22pipe/ClaremontSouth22_Z_"+filestrain+".csv");
strain_Z=strain_Z(2:end,:);
n_points=size(strain_Z,1);
% strain_Z=flip(strain_Z,1);
plot_index=1:ceil(n_points/plot_number):n_points;
s = linspace(0, sample_spacing*(n_points-1), n_points)';
slicetimeid=42
%%
strain_N=detrend(strain_N);
strain_S=detrend(strain_S);
strain_Z=detrend(strain_Z);
%%
strain_diff=strain_N(:,slicetimeid)*1e-6-strain_S(:,slicetimeid)*1e-6;
% s=[1:length(strain_N(:,slicetimeid))]*0.1;
curveture=strain_diff/0.64;
latdisp=cumtrapz(s,cumtrapz(s,curveture));
figure
plot(s,latdisp)
nsample=1000;
errorlevel=50*1e-6;
ftsz=20;
z=0.64;
L=max(s);
colors=colororder;
simpsupportlegendline={};
axiallegendline={};
legs={};
figure(10001)
% for slicetimeid=[1:8:33 34:6:42]
colorid=0;
for slicetimeid=[1:7:42 42]
    colorid=colorid+1;
    noiseN=errorlevel*randn(length(strain_N(:,slicetimeid)),nsample);
    noiseS=errorlevel*randn(length(strain_N(:,slicetimeid)),nsample);
    noiseZ=errorlevel*randn(length(strain_N(:,slicetimeid)),nsample);
    plotcolor=colors(mod(colorid,7)+1,:);
    strain_diff=strain_N(:,slicetimeid)*ones([1,nsample])*1e-6+noiseN-strain_S(:,slicetimeid)*ones([1,nsample])*1e-6-noiseS;
    % strain_avg=(strain_N(:,slicetimeid)*ones([1,nsample])*1e-6+noiseN+strain_S(:,slicetimeid)*ones([1,nsample])*1e-6+noiseS+strain_Z(:,slicetimeid)*ones([1,nsample])*1e-6+noiseZ)/3;
    strain_avg=(strain_N(:,slicetimeid)*ones([1,nsample])*1e-6+noiseN+strain_S(:,slicetimeid)*ones([1,nsample])*1e-6+noiseS)/2;
    curveture=strain_diff/z;
    latdisp=cumtrapz(s,cumtrapz(s,curveture));
    axial_disp=cumtrapz(s,strain_avg);
    %% bring to zero
    simpsupportdisp=latdisp-0*s/L*(latdisp(end,:)+0e-3);
    % avgdisp=mean(simpsupportdisp,"all");
    [~,indexfixed]=min(abs(s-69));
    %% rotate to flat
    [~,simpsupportdisp_tangent]=gradient(simpsupportdisp,1:nsample,s);
    simpsupportdisp=simpsupportdisp-((s-s(indexfixed))).*simpsupportdisp_tangent(indexfixed,:);
    %% rotate to negatively equal
    % simpsupportdisp=simpsupportdisp-((s-s(indexfixed))).*(simpsupportdisp(end,:)+simpsupportdisp(1,:))/(s(end)+s(1)-2*s(indexfixed));
    avgdisp=simpsupportdisp(indexfixed,:);
    simpsupportdisp=simpsupportdisp-(1+0*s)*avgdisp;
    %%
    figure(10001)

    plot(s,simpsupportdisp,color=[plotcolor,0.01]);
    hold on
    simpsupportmean=mean(simpsupportdisp,2);
    simpsupportcvar=var(simpsupportdisp,0,2);
    simpsupportcstd=sqrt(simpsupportcvar);

    plot(s,simpsupportmean,color=[plotcolor,1],LineWidth=1.0);
    hold on
    plot(s,simpsupportmean+simpsupportcstd,color=[plotcolor,1],LineWidth=0.5,LineStyle="--");
    plot(s,simpsupportmean-simpsupportcstd,color=[plotcolor,1],LineWidth=0.5,LineStyle="--");
    simpsupportlegendline= [simpsupportlegendline plot(nan,color=[plotcolor,1],LineWidth=2.0)] ;
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    currenttime.Format = 'yyyy_MM';
    caption = sprintf('%s', string(currenttime));
    legs= [legs caption];
    set(gca, 'XDir','reverse')
    % set(gca, 'YDir','reverse')
    %%
    figure(20001)
    plot(s,axial_disp,color=[plotcolor,0.01]);
    hold on
    axialmean=mean(axial_disp,2);
    axialvar=var(axial_disp,0,2);
    axialstd=sqrt(axialvar);
    plot(s,axialmean,color=[plotcolor,1],LineWidth=1.0);
    hold on
    plot(s,axialmean+axialstd,color=[plotcolor,1],LineWidth=0.5,LineStyle="--");
    plot(s,axialmean-axialstd,color=[plotcolor,1],LineWidth=0.5,LineStyle="--");
    axiallegendline= [axiallegendline plot(nan,color=[plotcolor,1],LineWidth=2.0)] ;
end
%%
figure(10001)
simpsupportlegendline= [simpsupportlegendline plot(nan,color='k',LineWidth=2.0)] ;
legs= [legs '$E[w(s)]$'];
simpsupportlegendline= [simpsupportlegendline plot(nan,color='k',LineWidth=0.5,LineStyle="--")] ;
legs= [legs '$E[w(s)] \pm Std[w(s)]$'];
% legend(legendline, legs, 'location', 'best', 'Interpreter','latex')
legend(simpsupportlegendline, legs, 'location', 'eastoutside', 'Interpreter','latex')
% set(gcf, 'Position', get(0, 'Screensize'));
ylabel("Lateral displacement $w(s)$ [m]",'interpreter',"latex",'fontsize',14)
xlabel("Distance $s$ [m]",'interpreter',"latex",'fontsize',14)
title('Laterial Displacement')
ax = gca;
ax.FontSize = ftsz;
set(gcf,"Position",[10         478        1917         518])
grid on
xlim([0,max(s)]);
saveas(gcf,'claremont_data/22pipe/lateraldisp_tcpipe_faultzero.png')
figure(10002)
coeff=2*errorlevel^2/z^2*sample_spacing;
% fplot(@(s)sqrt(coeff*(s^4/3/L-2*s^3/3+L*s^2/3)),LineWidth=2)
fplot(@(ss)abs(sqrt(coeff*(ss-s(indexfixed)).^3/3)),[0,300],LineWidth=2)
hold on
plot(s,simpsupportcstd,"o",LineWidth=2,MarkerIndices=1:20:length(s))
xlim([0,max(s)]);
% ylim([-0.5,0.5])
ylabel("Deviation $Std[w(s)]$ [m]",'interpreter',"latex",'fontsize',14)
xlabel("Distance $s$ [m]",'interpreter',"latex",'fontsize',14)
legend('$Std[w(s)]$ - Analytical','$Std[w(s)]$ - Monte Carlo','interpreter',"latex")
title('Propogation of Error')
grid on
ax = gca;
ax.FontSize = ftsz;
set(gcf,"Position",[10         478        1917         518])
saveas(gcf,'claremont_data/22pipe/lateraldisperror_tcpipe_faultzero.png')
%%
figure(20001)
axiallegendline= [axiallegendline plot(nan,color='k',LineWidth=2.0)] ;
legs= [legs '$E[w(s)]$'];
axiallegendline= [axiallegendline plot(nan,color='k',LineWidth=0.5,LineStyle="--")] ;
legs= [legs '$E[w(s)] \pm Std[w(s)]$'];
% legend(legendline, legs, 'location', 'best', 'Interpreter','latex')
legend(axiallegendline, legs, 'location', 'eastoutside', 'Interpreter','latex')
% set(gcf, 'Position', get(0, 'Screensize'));
ylabel("Axial displacement $w(s)$ [m]",'interpreter',"latex",'fontsize',14)
xlabel("Distance $s$ [m]",'interpreter',"latex",'fontsize',14)
title('Axial Displacement')
xlim([0,max(s)]);
ax = gca;
ax.FontSize = ftsz;
set(gcf,"Position",[10         478        1917         518])
grid on
saveas(gcf,'claremont_data/22pipe/axialdisp_tcpipe_faultzero.png')
figure(20002)
coeff=errorlevel^2*sample_spacing/2;
fplot(@(s)sqrt(coeff*(s)),LineWidth=2)
hold on
plot(s,axialstd,"o",LineWidth=2,MarkerIndices=1:20:length(s))
xlim([0,max(s)]);
% ylim([0,6e-3])
ylabel("Deviation $Std[w(s)]$ [m]",'interpreter',"latex",'fontsize',14)
xlabel("Distance $s$ [m]",'interpreter',"latex",'fontsize',14)
legend('$Std[w(s)]$ - Analytical','$Std[w(s)]$ - Monte Carlo','interpreter',"latex")
title('Propogation of Error')
grid on
ax = gca;
ax.FontSize = ftsz;
set(gcf,"Position",[10         478        1917         518])
saveas(gcf,'claremont_data/22pipe/axialdisperror_tcpipe_faultzero.png')
%%

figure(30001)



vw = VideoWriter("claremont_data/22pipe/simpsupportvideo_tcpipe_faultzero","Motion JPEG AVI");
vw.FrameRate=1/0.6;
open(vw)
for slicetimeid=1:1:42
    figure(30001)
    simpsupportlegendline={};
    axiallegendline={};
    legs={};
    hold off
    noiseN=errorlevel*randn(length(strain_N(:,slicetimeid)),nsample);
    noiseS=errorlevel*randn(length(strain_N(:,slicetimeid)),nsample);
    noiseZ=errorlevel*randn(length(strain_N(:,slicetimeid)),nsample);
    plotcolor=colors(mod(slicetimeid,7)+1,:);
    strain_diff=strain_N(:,slicetimeid)*ones([1,nsample])*1e-6+noiseN-strain_S(:,slicetimeid)*ones([1,nsample])*1e-6-noiseS;
    strain_avg=(strain_N(:,slicetimeid)*ones([1,nsample])*1e-6+noiseN+strain_S(:,slicetimeid)*ones([1,nsample])*1e-6+noiseS+strain_Z(:,slicetimeid)*ones([1,nsample])*1e-6+noiseZ)/3;
    curveture=strain_diff/z;
    latdisp=cumtrapz(s,cumtrapz(s,curveture));
    axial_disp=cumtrapz(s,strain_avg);
    simpsupportdisp=latdisp-s/L*latdisp(end,:);
    plot(s,simpsupportdisp,color=[plotcolor,0.01]);
    hold on
    simpsupportmean=mean(simpsupportdisp,2);
    simpsupportcvar=var(simpsupportdisp,0,2);
    simpsupportcstd=sqrt(simpsupportcvar);

    plot(s,simpsupportmean,color=[plotcolor,1],LineWidth=1.0);
    hold on
    plot(s,simpsupportmean+simpsupportcstd,color=[plotcolor,1],LineWidth=0.5,LineStyle="--");
    plot(s,simpsupportmean-simpsupportcstd,color=[plotcolor,1],LineWidth=0.5,LineStyle="--");
    simpsupportlegendline= [simpsupportlegendline plot(nan,color=[plotcolor,1],LineWidth=2.0)] ;
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    currenttime.Format = 'yyyy_MM';
    caption = sprintf('%s', string(currenttime));
    legs= [legs caption];
    simpsupportlegendline= [simpsupportlegendline plot(nan,color='k',LineWidth=2.0)] ;
    legs= [legs '$E[w(s)]$'];
    simpsupportlegendline= [simpsupportlegendline plot(nan,color='k',LineWidth=0.5,LineStyle="--")] ;
    legs= [legs '$E[w(s)] \pm Std[w(s)]$'];
    % legend(legendline, legs, 'location', 'best', 'Interpreter','latex')
    ylim([-0.18,0.18])
    % set(gcf, 'Position', get(0, 'Screensize'));
    ylabel("Lateral displacement $w(s)$ [m]",'interpreter',"latex",'fontsize',14)
    xlabel("Distance $s$ [m]",'interpreter',"latex",'fontsize',14)
    title('Laterial Displacement')
    legend(simpsupportlegendline, legs, 'location', 'eastoutside', 'Interpreter','latex')
    xlim([0,max(s)]);
    ax = gca;
    ax.FontSize = ftsz;
    set(gcf,"Position",[10         478        1917         518])
    grid on
    frame = getframe(gcf);
    writeVideo(vw,frame)
    hold off
end
close(vw)
%%

figure(40001)



vw = VideoWriter("claremont_data/22pipe/axialdispvideo_tcpipe_faultzero","Motion JPEG AVI");
vw.FrameRate=1/0.6;
open(vw)
for slicetimeid=1:1:42
    figure(40001)
    simpsupportlegendline={};
    axiallegendline={};
    legs={};
    hold off
    noiseN=errorlevel*randn(length(strain_N(:,slicetimeid)),nsample);
    noiseS=errorlevel*randn(length(strain_N(:,slicetimeid)),nsample);
    noiseZ=errorlevel*randn(length(strain_N(:,slicetimeid)),nsample);
    plotcolor=colors(mod(slicetimeid,7)+1,:);
    strain_diff=strain_N(:,slicetimeid)*ones([1,nsample])*1e-6+noiseN-strain_S(:,slicetimeid)*ones([1,nsample])*1e-6-noiseS;
    strain_avg=(strain_N(:,slicetimeid)*ones([1,nsample])*1e-6+noiseN+strain_S(:,slicetimeid)*ones([1,nsample])*1e-6+noiseS+strain_Z(:,slicetimeid)*ones([1,nsample])*1e-6+noiseZ)/3;
    curveture=strain_diff/z;
    latdisp=cumtrapz(s,cumtrapz(s,curveture));
    axial_disp=cumtrapz(s,strain_avg);
    simpsupportdisp=latdisp-s/L*latdisp(end,:);
    plot(s,axial_disp,color=[plotcolor,0.01]);
    hold on
    simpsupportmean=mean(simpsupportdisp,2);
    simpsupportcvar=var(simpsupportdisp,0,2);
    simpsupportcstd=sqrt(simpsupportcvar);
    axialmean=mean(axial_disp,2);
    axialvar=var(axial_disp,0,2);
    axialstd=sqrt(axialvar);
    plot(s,axialmean,color=[plotcolor,1],LineWidth=1.0);
    hold on
    plot(s,axialmean+axialstd,color=[plotcolor,1],LineWidth=0.5,LineStyle="--");
    plot(s,axialmean-axialstd,color=[plotcolor,1],LineWidth=0.5,LineStyle="--");
    simpsupportlegendline= [simpsupportlegendline plot(nan,color=[plotcolor,1],LineWidth=2.0)] ;
    currenttime=Times(slicetimeid);
    currenttime=currenttime{1};
    currenttime.Format = 'yyyy_MM';
    caption = sprintf('%s', string(currenttime));
    legs= [legs caption];
    simpsupportlegendline= [simpsupportlegendline plot(nan,color='k',LineWidth=2.0)] ;
    legs= [legs '$E[w(s)]$'];
    simpsupportlegendline= [simpsupportlegendline plot(nan,color='k',LineWidth=0.5,LineStyle="--")] ;
    legs= [legs '$E[w(s)] \pm Std[w(s)]$'];
    % legend(legendline, legs, 'location', 'best', 'Interpreter','latex')
    ylim([-0.02,0.02])
    % set(gcf, 'Position', get(0, 'Screensize'));
    ylabel("Axial displacement $w(s)$ [m]",'interpreter',"latex",'fontsize',14)
    xlabel("Distance $s$ [m]",'interpreter',"latex",'fontsize',14)
    title('Axial Displacement')
    legend(simpsupportlegendline, legs, 'location', 'eastoutside', 'Interpreter','latex')
    xlim([0,max(s)]);
    ax = gca;
    ax.FontSize = ftsz;
    set(gcf,"Position",[10         478        1917         518])
    grid on
    frame = getframe(gcf);
    writeVideo(vw,frame)
    hold off
end
close(vw)