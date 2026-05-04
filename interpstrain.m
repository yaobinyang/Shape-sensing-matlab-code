strainhelix1s=[]
strainhelix2s=[]
strainhelix3s=[]
strainhelix0s=[]
for slicetimeid=[40 120 175]
% slicetimeid=175;
%40 120 175
dat0=dat0-dat0(2,:);
dat1=dat1-dat1(2,:);
dat2=dat2-dat2(2,:);
dat3=dat3-dat3(2,:);
reading_strainhelix1=dat1(slicetimeid,:);
loc_strainhelix1=loc1-loc1(1);
reading_strainhelix2=dat3(slicetimeid,:);
loc_strainhelix2=loc3-loc3(1);
reading_strainhelix3=dat2(slicetimeid,:);
loc_strainhelix3=loc2-loc2(1);
reading_strainstraight0=dat0(slicetimeid,:);
loc_strainstraight0=loc0-loc0(1);
strainhelix1i=interp1(loc_strainhelix1,reading_strainhelix1,s_all*loc_strainhelix1(end)/s_all(end));
strainhelix2i=interp1(loc_strainhelix2,reading_strainhelix2,s_all*loc_strainhelix2(end)/s_all(end));
strainhelix3i=interp1(loc_strainhelix3,reading_strainhelix3,s_all*loc_strainhelix3(end)/s_all(end));
strainstraight0i=interp1(loc_strainstraight0,reading_strainstraight0,s_all*loc_strainstraight0(end)/s_all(end));
nans=isnan(strainstraight0i);
strainstraight0i(nans) = interp1(s_all(~nans), strainstraight0i(~nans), s_all(nans));
nans=isnan(strainhelix1i);
strainhelix1i(nans) = interp1(s_all(~nans), strainhelix1i(~nans), s_all(nans));
nans=isnan(strainhelix2i);
strainhelix2i(nans) = interp1(s_all(~nans), strainhelix2i(~nans), s_all(nans));
nans=isnan(strainhelix3i);
strainhelix3i(nans) = interp1(s_all(~nans), strainhelix3i(~nans), s_all(nans));
strainhelix1s=[strainhelix1s;strainhelix1i];
strainhelix2s=[strainhelix2s;strainhelix2i];
strainhelix3s=[strainhelix3s;strainhelix3i];
strainhelix0s=[strainhelix0s;strainstraight0i]
% strainstraight0=lowpass(strainstraight0',0.001)';
end
figure(77771)
plot(s_all,[strainhelix1i;strainhelix2i;strainhelix3i;strainstraight0i]',LineWidth=2)
legend('Helix 1','Helix 2','Helix 3','Straight')
xlabel('Offset [m]','Interpreter','latex')
ylabel('Strain Measurement [$\mu\varepsilon$]','Interpreter','latex')
title('Shape - Bending Arch')
set(gca,"FontSize",14)
