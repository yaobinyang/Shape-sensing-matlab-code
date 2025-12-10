
% datinfo1=readcell("DFOSshapesensing_nervesens\ShapeSensing_NerveSensor_3DSensor_2024-10-04_21-38-07_ch1_gages.tsv",FileType="text",Delimiter='\t',Range=[1 1 200 200]);
dat1=readmatrix("nervesenswithvideo/shapesensing_video_2024-11-19_02-59-48_ch2_gages.tsv",FileType="text",Delimiter='\t',Range=[35 6]);
loc1=readmatrix("nervesenswithvideo/shapesensing_video_2024-11-19_02-59-48_ch2_gages.tsv",FileType="text",Delimiter='\t',Range=[34 6 34 6+size(dat1,2)]);
time1=readtable("nervesenswithvideo/shapesensing_video_2024-11-19_02-59-48_ch2_gages.tsv",FileType="text",Delimiter='\t',Range=[35 1 35+size(dat1,1) 1]);
time1=time1.Var1;
% datinfo2=readcell("DFOSshapesensing_nervesens\ShapeSensing_NerveSensor_3DSensor_2024-10-04_21-38-07_ch2_gages.tsv",FileType="text",Delimiter='\t',Range=[1 1 200 200]);
dat2=readmatrix("nervesenswithvideo/shapesensing_video_2024-11-19_02-59-48_ch1_gages.tsv",FileType="text",Delimiter='\t',Range=[35 6]);
loc2=readmatrix("nervesenswithvideo/shapesensing_video_2024-11-19_02-59-48_ch1_gages.tsv",FileType="text",Delimiter='\t',Range=[34 6 34 6+size(dat2,2)]);
time2=readtable("nervesenswithvideo/shapesensing_video_2024-11-19_02-59-48_ch1_gages.tsv",FileType="text",Delimiter='\t',Range=[35 1 35+size(dat2,1) 1]);
time2=time2.Var1;
% % datinfo3=readcell("DFOSshapesensing_nervesens\ShapeSensing_NerveSensor_3DSensor_2024-10-04_21-38-07_ch3_gages.tsv",FileType="text",Delimiter='\t',Range=[1 1 200 200]);
% dat3=readmatrix("DFOSshapesensing_nervesens\ShapeSensing_NerveSensor_3DSensor_2024-10-04_21-38-07_ch3_gages.tsv",FileType="text",Delimiter='\t',Range=[35 6]);
% loc3=readmatrix("DFOSshapesensing_nervesens\ShapeSensing_NerveSensor_3DSensor_2024-10-04_21-38-07_ch3_gages.tsv",FileType="text",Delimiter='\t',Range=[34 6 34 6+size(dat3,2)]);
% time3=readtable("DFOSshapesensing_nervesens\ShapeSensing_NerveSensor_3DSensor_2024-10-04_21-38-07_ch3_gages.tsv",FileType="text",Delimiter='\t',Range=[35 1 35+size(dat3,1) 1]);
% time3=time3.Var1;
% datinfo4=readcell("DFOSshapesensing_nervesens\ShapeSensing_NerveSensor_3DSensor_2024-10-04_21-38-07_ch4_gages.tsv",FileType="text",Delimiter='\t',Range=[1 1 200 200]);
dat4=readmatrix("nervesenswithvideo/shapesensing_video_2024-11-19_02-59-48_ch4_gages.tsv",FileType="text",Delimiter='\t',Range=[35 6]);
loc4=readmatrix("nervesenswithvideo/shapesensing_video_2024-11-19_02-59-48_ch4_gages.tsv",FileType="text",Delimiter='\t',Range=[34 6 34 6+size(dat4,2)]);
time4=readtable("nervesenswithvideo/shapesensing_video_2024-11-19_02-59-48_ch4_gages.tsv",FileType="text",Delimiter='\t',Range=[35 1 35+size(dat4,1) 1]);
time4=time4.Var1;

% loc1=[loc1;loc1s];
% loc2=[loc2;loc2s];
% loc3=[loc3;loc3s];
% loc4=[loc4;loc4s];
% time1=[time1;time1s];
% time2=[time2;time2s];
% % time3=[time3;time3s];
% time4=[time4;time4s];

startid1=469;
endid1=1997;
startid2=477;
endid2=2005;
startid3=478;
endid3=2008;
startid4=475;
endid4=2006;

ds_factor_loc=1;
ds_factor_time=1;
dat1=dat1(1:ds_factor_time:length(time1),startid1:ds_factor_loc:endid1);
loc1=loc1(startid1:ds_factor_loc:endid1);
time1=time1(1:ds_factor_time:length(time1));
dat2=dat2(1:ds_factor_time:length(time2),startid2:ds_factor_loc:endid2);
loc2=loc2(startid2:ds_factor_loc:endid2);
time2=time2(1:ds_factor_time:length(time2));
% dat3=dat3(1:ds_factor_time:length(time3),startid3:ds_factor_loc:endid3);
% loc3=loc3(startid3:ds_factor_loc:endid3);
% time3=time3(1:ds_factor_time:length(time3));
dat4=dat4(1:ds_factor_time:length(time4),startid4:ds_factor_loc:endid4);
loc4=loc4(startid4:ds_factor_loc:endid4);
time4=time4(1:ds_factor_time:length(time4));

