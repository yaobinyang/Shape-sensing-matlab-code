slicetimeid=113;
%40 120 175
dat1=dat1-dat1(2,:);
dat2=dat2-dat2(2,:);
dat3=dat3-dat3(2,:);
dat4=dat4-dat4(2,:);
reading_strainu1=dat1(slicetimeid,:);
loc_strainu1=loc1-loc1(1);
reading_strainu2=dat2(slicetimeid,:);
loc_strainu2=loc2-loc2(1);
reading_strainl1=dat3(slicetimeid,:);
loc_strainl1=loc3-loc3(1);
reading_strainl2=dat4(slicetimeid,:);
loc_strainl2=loc4-loc4(1);
strainu1=interp1(loc_strainu1,reading_strainu1,s_all*loc_strainu1(end)/s_all(end));
strainu2=interp1(loc_strainu2,reading_strainu2,s_all*loc_strainu2(end)/s_all(end));
strainl1=interp1(loc_strainl1,reading_strainl1,s_all*loc_strainl1(end)/s_all(end));
strainl2=interp1(loc_strainl2,reading_strainl2,s_all*loc_strainl2(end)/s_all(end));
nans=isnan(strainl2);
strainl2(nans) = interp1(s_all(~nans), strainl2(~nans), s_all(nans));
nans=isnan(strainu1);
strainu1(nans) = interp1(s_all(~nans), strainu1(~nans), s_all(nans));
nans=isnan(strainu2);
strainu2(nans) = interp1(s_all(~nans), strainu2(~nans), s_all(nans));
nans=isnan(strainl1);
strainl1(nans) = interp1(s_all(~nans), strainl1(~nans), s_all(nans));

% figure
% plot([strainu1;strainu2;strainl1;strainl2]')
