figure(101)
hold off
plot(rod_axis_rs(3,:),rod_axis_rs(1,:),'k',LineWidth=2)
hold on
plot(r_all_NFCA(3,:),r_all_NFCA(1,:),'g',LineWidth=2)
plot(r_all_NFCA(3,:),r_all_NFCA(2,:),'r',LineWidth=2)
plot(r_all_NFCA(3,:),sqrt(r_all_NFCA(2,:).^2+r_all_NFCA(1,:).^2+(rod_axis_rs(3,:)-r_all_NFCA(3,:)).^2),'b',LineWidth=2)
xlabel('Rope Direction [m]')
ylabel('Estimated Deformation [m]')
legend('Orig','Deformed Shape Projection on Y','Deformed Shape Projection on X','Deformed Shape 2-Norm')
p=r_all_NFCA(:,end);
% figure(102)
% for point=1:length(s_all)
%     dist(point)=norm()
% end

% deform=0;
% a=fsolve(@(x) x^2+((deform-r_all_NFCA(1,end)*x)/r_all_NFCA(2,end))^2-1,1);
% b=(deform-r_all_NFCA(1,end)*a)/r_all_NFCA(2,end);
% 
% figure(104)
% hold off
% plot(rod_axis_rs(3,:),rod_axis_rs(1,:),'k',LineWidth=2)
% hold on
% plot(r_all_NFCA(3,:),r_all_NFCA(1,:),'g',LineWidth=2)
% plot(r_all_NFCA(3,:),r_all_NFCA(2,:),'r',LineWidth=2)
% plot(r_all_NFCA(3,:),r_all_NFCA(1,:)*a+r_all_NFCA(2,:)*b,'b',LineWidth=2)
% xlabel('Rope Direction [m]')
% ylabel('Estimated Deformation [m]')
% legend('Orig','Deformed Shape Projection on Y','Deformed Shape Projection on X','Deformed Shape 2-Norm')

% p=r_all_NFCA(:,end);
% options = optimoptions('fmincon','Display','iter',...
%     "Algorithm","interior-point",...
%     "EnableFeasibilityMode",true,...
%     "SubproblemAlgorithm","cg");
% fun=@(abc) helperfun(abc,r_all_NFCA);
% nonlincon=@(abc) helpcon(abc,p);
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% abc0 = [0;0;0];
% lb=-0.1*[pi;pi;pi];
% ub=0.1*[pi;pi;pi];
% abc = fmincon(fun,abc0,A,b,Aeq,beq,lb,ub,nonlincon,options);
% Rox=rotx(abc(1));
% Roy=roty(abc(2));
% Roz=rotz(abc(3));
% R=Rox*Roy*Roz;
% figure(105)
% rotated_r_all_NFCA=R*r_all_NFCA;
% hold off
% plot(rod_axis_rs(3,:),rod_axis_rs(1,:),'k',LineWidth=2)
% hold on
% plot(rotated_r_all_NFCA(3,:),rotated_r_all_NFCA(1,:),'g',LineWidth=2)
% plot(rotated_r_all_NFCA(3,:),rotated_r_all_NFCA(2,:),'r',LineWidth=2)
% plot(rotated_r_all_NFCA(3,:),sqrt(rotated_r_all_NFCA(2,:).^2+rotated_r_all_NFCA(1,:).^2+(rod_axis_rs(3,:)-rotated_r_all_NFCA(3,:)).^2),'b',LineWidth=2)
% xlabel('Rope Direction [m]')
% ylabel('Estimated Deformation [m]')
% legend('Orig','Deformed Shape Projection on Y','Deformed Shape Projection on X','Deformed Shape 2-Norm')
% 
% % p=r_all_NFCA(:,end);
% % options = optimoptions('ga','Display','iter');
% % fun=@(abc) helperfun(abc,r_all_NFCA);
% % nonlincon=@(abc) helpcon(abc,p);
% % A = [];
% % b = [];
% % Aeq = [];
% % beq = [];
% % abc0 = [0;0;0];
% % lb=-[pi;pi;pi];
% % ub=[pi;pi;pi];
% % abc = ga(fun,3,A,b,Aeq,beq,lb,ub,[],options);
% % Rox=rotx(abc(1));
% % Roy=roty(abc(2));
% % Roz=rotz(abc(3));
% % R=Rox*Roy*Roz;
% % figure(115)
% % rotated_r_all_NFCA=R*r_all_NFCA;
% % hold off
% % plot(rod_axis_rs(3,:),rod_axis_rs(1,:),'k',LineWidth=2)
% % hold on
% % plot(rotated_r_all_NFCA(3,:),rotated_r_all_NFCA(1,:),'g',LineWidth=2)
% % plot(rotated_r_all_NFCA(3,:),rotated_r_all_NFCA(2,:),'r',LineWidth=2)
% % plot(rotated_r_all_NFCA(3,:),sqrt(rotated_r_all_NFCA(2,:).^2+rotated_r_all_NFCA(1,:).^2+(rod_axis_rs(3,:)-rotated_r_all_NFCA(3,:)).^2),'b',LineWidth=2)
% % xlabel('Rope Direction [m]')
% % ylabel('Estimated Deformation [m]')
% % legend('Orig','Deformed Shape Projection on Y','Deformed Shape Projection on X','Deformed Shape 2-Norm')
% 
% 
figure(106)
plot(strainstraight0)

filtercoeff=1e-3;
design_uk_all=lowpass(u_k_all_NFCA',filtercoeff)';
figure(107)
plot(design_uk_all')
figure(117)
plot(u_k_all_NFCA')
% design_vk_all=v_k_all_NFCA;
% [r_all_lowpass,D_ik_all_lowpass,helix_cable_strain_all_filtered,straight_cable_strain_all_filtered,u_0k_all,v_0k_all]=generate_MSF_kinematics(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s_all,rod_axis_rs,isassigneddir,assigned_dir,design_uk_all,design_vk_all);
% 
% p=r_all_lowpass(:,end);
% options = optimoptions('ga','Display','iter');
% fun=@(abc) helperfun(abc,r_all_lowpass);
% nonlincon=@(abc) helpcon(abc,p);
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% abc0 = [0;0;0];
% lb=-[pi;pi;pi];
% ub=[pi;pi;pi];
% abc = ga(fun,3,A,b,Aeq,beq,lb,ub,nonlincon,options);
% Rox=rotx(abc(1));
% Roy=roty(abc(2));
% Roz=rotz(abc(3));
% R=Rox*Roy*Roz;
% figure(125)
% rotated_r_all_NFCA=R*r_all_lowpass;
% hold off
% plot(rod_axis_rs(3,:),rod_axis_rs(1,:),'k',LineWidth=2)
% hold on
% plot(rotated_r_all_NFCA(3,:),rotated_r_all_NFCA(1,:),'g',LineWidth=2)
% plot(rotated_r_all_NFCA(3,:),rotated_r_all_NFCA(2,:),'r',LineWidth=2)
% plot(rotated_r_all_NFCA(3,:),sqrt(rotated_r_all_NFCA(2,:).^2+rotated_r_all_NFCA(1,:).^2+(rod_axis_rs(3,:)-rotated_r_all_NFCA(3,:)).^2),'b',LineWidth=2)
% xlabel('Rope Direction [m]')
% ylabel('Estimated Deformation [m]')
% legend('Orig','Deformed Shape Projection on Y','Deformed Shape Projection on X','Deformed Shape 2-Norm')
% figure(116)
% plot(strainstraight0)

function outofplannorm=helperfun(abc,r_all_NFCA)
Rox=rotx(abc(1));
Roy=roty(abc(2));
Roz=rotz(abc(3));
R=Rox*Roy*Roz;
mappedmeasures=R*r_all_NFCA;
outofplannorm=mappedmeasures(2,:)*mappedmeasures(2,:)';
end

function outofplannorm=helperfun1(abc,r_all_NFCA)
Rox=rotx(abc(1));
Roy=roty(abc(2));
Roz=rotz(abc(3));
R=Rox*Roy*Roz;
mappedmeasures=R*r_all_NFCA;
outofplannorm=mappedmeasures(2,:)*mappedmeasures(2,:)'+mappedmeasures(3,:)*mappedmeasures(3,:)';
end

function e=consfun(abc,p)
Rox=rotx(abc(1));
Roy=roty(abc(2));
Roz=rotz(abc(3));
R=Rox*Roy*Roz;
e=R*p-norm(p)*[0;0;1];
end

function [c,ceq] = helpcon(abc,p)
c=[];
ceq = norm(consfun(abc,p));
end