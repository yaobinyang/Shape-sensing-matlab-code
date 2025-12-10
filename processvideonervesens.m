clear
clc
close all
vidObj = VideoReader("neversenswithvideo\00000.MTS");
frame = read(vidObj,20);
figure
imshow(frame)
axis on
% fixpt=100*0.0127*[0 10;10 10;10 0;0 0];
% fixpt=[1028 849;1086 849;1086 733;1028 733];
% movpt=[1028 849; 1086 849; 1080 776; 1026 773];
ratio=87/540;
width=1158-702;
length=width/ratio;
movpt=[702 1048; 1158 1053; 1057 133; 871 131];
fixpt=[702 1048; 702+width 1048; 702+width 1048-length; 702 1048-length];

tform = fitgeotform2d(movpt,fixpt,"projective");
invtform = invert(tform);
warppedframe = imwarp(frame,tform);
figure
imshow(warppedframe)
axis on
figure
    croppedframe = imcrop(warppedframe,[2710 1802 3161-2710 4075-1802]);
    imshow(croppedframe)
figure
binedfame=imbinarize(im2gray(croppedframe));
imshow(binedfame)
edgedframe = edge(binedfame,'Canny');
imshow(edgedframe)
figure
[z_frame,y_frame]=find(edgedframe);
scalefac=(2270-232)/4;
z_frame=2270-z_frame;
z_frame=z_frame/scalefac;
y_frame=y_frame/scalefac;
[z_frame,sort_idenx]=sort(z_frame);
y_frame=y_frame(sort_idenx);
plot(y_frame,z_frame)
axis equal
figure
framerate=100;
frameid=0;
while(hasFrame(vidObj))
    frameid=frameid+1;
    frame = readFrame(vidObj);
    if(mod(frameid,framerate)==0)
        
        warppedframe = imwarp(frame,tform);
        framerim=200;
        croppedframe = imcrop(warppedframe,[2710-framerim 1802-framerim 3161-2710+2*framerim 4075-1802+2*framerim]);
        % binedfame=imbinarize(im2gray(croppedframe));
        % edgedframe = edge(binedfame,'Canny');
        figure(701)
        imshow(croppedframe)
        % [z_frame,y_frame]=find(edgedframe);
        % scalefac=(2270-232)/4;
        % z_frame=2270-z_frame;
        % z_frame=z_frame/scalefac;
        % y_frame=y_frame/scalefac;
        % y_frame=y_frame-0.4867;
        % [z_frame,sort_idenx]=sort(z_frame);
        % y_frame=y_frame(sort_idenx);
        % plot(y_frame,z_frame)
        % ylim([0,4.5]);
        % xlim([-1,1]);
        axis equal
        title(sprintf("Current Time = %.3f sec",vidObj.CurrentTime))
        pause(0.01/vidObj.FrameRate)
    end
end

clear vidObj
