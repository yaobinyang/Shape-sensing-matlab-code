hasanalytical=true;
MSE_all_proposed_2=zeros(1,length(s_all));
MSE_all_proposed_3=zeros(1,length(s_all));
MSE_all_proposed_1=zeros(1,length(s_all));
MSE_all_proposed_4=zeros(1,length(s_all));
MSE_all_Nishiguchi=zeros(1,length(s_all));
% MSE_all_MCCS=zeros(size(s_all));
% MSE_all_MCEBA=zeros(size(s_all));
% MSE_all_MCEBS=zeros(size(s_all));
grads=gradient(s_all);
normvest=max(vecnorm(v_k_all_SACS,2,1));
elongationest=mean(v_k_all_SACS,2);
elongationest=elongationest(3);
for point=1:length(s_all)
    C_A=zeros(3,3);
    for k=1:3
        for i=1:3
            for j=1:3
                for p=1:n_straight+n_helix
                    for q=1:n_straight+n_helix
                        for l=1:3
                            for m=1:3
                            C_A(k,l)=C_A(k,l)+Tensorhelpervarepsi(k,i,j)*Cu_k_all_SACS(j,p,point)*(p==q) ...
                                *Tensorhelpervarepsi(l,i,m)*Cu_k_all_SACS(m,q,point);
                            end
                        end
                    end
                end
            end
        end
    end
    for p=1:n_straight+n_helix
        MSE_all_proposed_1(point)=MSE_all_proposed_1(point)+(Cu_k_all_SACS(1,p,point)^2+Cu_k_all_SACS(2,p,point)^2+0*Cu_k_all_SACS(3,p,point)^2)*s_all(point)^3/3;
        MSE_all_proposed_1(point)=MSE_all_proposed_1(point)+(Cv_k_all_SACS(1,p,point)^2+Cv_k_all_SACS(2,p,point)^2+Cv_k_all_SACS(3,p,point)^2)*s_all(point);
        MSE_all_proposed_2(point)=MSE_all_proposed_2(point)+(Cu_k_all_SACS(1,p,point)^2+Cu_k_all_SACS(2,p,point)^2+Cu_k_all_SACS(3,p,point)^2)*s_all(point)^3/3;
        MSE_all_proposed_2(point)=MSE_all_proposed_2(point)+(Cv_k_all_SACS(1,p,point)^2+Cv_k_all_SACS(2,p,point)^2+Cv_k_all_SACS(3,p,point)^2)*s_all(point);
        MSE_all_proposed_3(point)=MSE_all_proposed_3(point)+(Cu_k_all_SACS(1,p,point)^2+Cu_k_all_SACS(2,p,point)^2+0*Cu_k_all_SACS(3,p,point)^2)*elongationest^2*s_all(point)^3/3;
        MSE_all_proposed_3(point)=MSE_all_proposed_3(point)+(Cv_k_all_SACS(1,p,point)^2+Cv_k_all_SACS(2,p,point)^2+Cv_k_all_SACS(3,p,point)^2)*s_all(point);
        MSE_all_proposed_4(point)=MSE_all_proposed_4(point)+(Cu_k_all_SACS(1,p,point)^2+Cu_k_all_SACS(2,p,point)^2+Cu_k_all_SACS(3,p,point)^2)*normvest^2*s_all(point)^3/3;
        MSE_all_proposed_4(point)=MSE_all_proposed_4(point)+(Cv_k_all_SACS(1,p,point)^2+Cv_k_all_SACS(2,p,point)^2+Cv_k_all_SACS(3,p,point)^2)*s_all(point);
    end
    MSE_all_proposed_1(point)=sqrt(MSE_all_proposed_1(point))*measure_error_level*sqrt(grads(point));
    MSE_all_proposed_2(point)=sqrt(MSE_all_proposed_2(point))*measure_error_level*sqrt(grads(point));
    MSE_all_proposed_3(point)=sqrt(MSE_all_proposed_3(point))*measure_error_level*sqrt(grads(point));
    MSE_all_proposed_4(point)=sqrt(MSE_all_proposed_4(point))*measure_error_level*sqrt(grads(point));
    lams=eig(C_A);
    lambda_max=max(lams);
    MSE_all_Nishiguchi(point)=sqrt(lambda_max*s_all(point)^3/3)*measure_error_level*sqrt(grads(point));
end
% plot(s_all,[MSE_all_SACS;MSE_all_NSCRK;MSE_all_NSCM;MSE_all_MCCS;MSE_all_MCEBA;MSE_all_proposed_2;MSE_all_proposed_1])
% legend({'Analytical SDE Solution';'Numerical SDE Solution - Weak RK';'Numerical SDE Solution - Milstein'; ...
%     'Monte-carlo simulation Cosserat';'Monte-carlo simulation EB approximation';'Extension of Nishiguchi, 2018';'Proposed Simple Solution'})
close all

figure(10001)
plot(s_all,[MSE_all_SACS;MSE_all_NSCRK;MSE_all_NSCM;MSE_all_proposed_1;MSE_all_proposed_2;MSE_all_proposed_3;MSE_all_proposed_4;MSE_all_Nishiguchi;MSE_all_MCCS;MSE_all_MCEBS])
legend({'Proposed SDE Solution';'Numerical SDE Solution - Weak RK';'Numerical SDE Solution - Milstein'; ...
    'Proposed Estimation 1';'Proposed Estimation 2';'Proposed Estimation 3';'Proposed Estimation 4';'Generalization of Nishiguchi, 2018';'Monte-carlo Simulation Cosserat';'Monte-carlo Simulation EB Approximation'})

plotindex=1:10:length(s_all);
legs={'Proposed SDE Solution'; ...
    'Proposed Estimation 1';'Proposed Estimation 2';'Proposed Estimation 3';'Proposed Estimation 4';'Generalization of Nishiguchi, 2018';'Monte-carlo Simulation Cosserat'}
linestyles={'-or';'--*g';'--*b';'--*c';'--*m';':diamondb';'-k'};
plotdata=[MSE_all_SACS;MSE_all_proposed_1;MSE_all_proposed_2;MSE_all_proposed_3;MSE_all_proposed_4;MSE_all_Nishiguchi;MSE_all_MCCS];
figure(10002)
hold on
for i=1:length(linestyles)
    plot(s_all(plotindex),plotdata(i,plotindex),linestyles{i},LineWidth=1.5,MarkerSize=10.0)
end
legend(legs)
xlabel("Rod Coordinate $s$ [m]",'interpreter',"latex",'fontsize',14)
ylabel({"Mean Squared Deviation of ";"the Rod Axis Estimation $MSD$ [m]"},'interpreter',"latex",'fontsize',14)

figure(10003)

for i=1:length(linestyles)
    semilogy(s_all(plotindex),plotdata(i,plotindex),linestyles{i},LineWidth=1.5,MarkerSize=10.0)
    hold on
end
legend(legs)
figure(10004)

for i=1:length(linestyles)
    loglog(s_all(plotindex),plotdata(i,plotindex),linestyles{i},LineWidth=1.5,MarkerSize=10.0)
    hold on
end
legend(legs)
xlabel("Rod Coordinate $s$ [m]",'interpreter',"latex",'fontsize',14)
ylabel({"Mean Squared Deviation of ";"the Rod Axis Estimation $MSD$ [m]"},'interpreter',"latex",'fontsize',14)
n_points=size(r_all,2);
r_all_analytical=rod_axis_rs;
D_analytical=D_ik_all;
D0=D_analytical(:,:,1);
for point=1:n_points
    s=s_all(point);
    uk=design_uk_all(:,1);
    vk=design_vk_all(:,1);
    norm_uk=max(norm(uk),1e-10);
    Ucross_helpermat=zeros(3,3);
    for i=1:3
        for k=1:3
            Ucross_helpermat(k,i)=0;
            for j=1:3
                Ucross_helpermat(k,i)=Ucross_helpermat(k,i)+Tensorhelpervarepsi(j,k,i)*uk(j);
            end
        end
    end
    Ucross_helpermat_trans=Ucross_helpermat';
    helpermat_exp=eye(3)+sin(norm_uk*s)*Ucross_helpermat_trans/norm_uk+(1-cos(norm_uk*s))*(Ucross_helpermat_trans*Ucross_helpermat_trans)/(norm_uk^2);
    D_analytical(:,:,point)=D0*helpermat_exp;
    helpermat_int=(eye(3)+Ucross_helpermat_trans*Ucross_helpermat_trans/(norm_uk^2))*s...
        -cos(norm_uk*s)*Ucross_helpermat_trans/(norm_uk^2)...
        -sin(norm_uk*s)*(Ucross_helpermat_trans*Ucross_helpermat_trans)/(norm_uk^3)...
        +Ucross_helpermat_trans/norm_uk^2;
    r_s=D0*helpermat_int*vk+r_all_analytical(:,1);
    r_all_analytical(:,point)=r_s;
end
if hasanalytical
    r_all_analytical=r_all_SACS;
end
plotindex=1:10:length(s_all);
figure(10005)
linestyles={'-^r',;'-vr';'-or';'--*g';'--*b';'--*c';'--*m';'-k'};
grid on
hold on
plot3(squeeze(r_all_NEBA(1,plotindex)), squeeze(r_all_NEBA(2,plotindex)), squeeze(r_all_NEBA(3,plotindex)),linestyles{1},LineWidth=1.5,MarkerSize=10.0);
plot3(squeeze(r_all_NFCA(1,plotindex)), squeeze(r_all_NFCA(2,plotindex)), squeeze(r_all_NFCA(3,plotindex)),linestyles{2},LineWidth=1.5,MarkerSize=10.0);
plot3(squeeze(r_all_SACS(1,plotindex)), squeeze(r_all_SACS(2,plotindex)), squeeze(r_all_SACS(3,plotindex)),linestyles{3},LineWidth=1.5,MarkerSize=10.0);
plot3(squeeze(r_all_NSCRK(1,plotindex)), squeeze(r_all_NSCRK(2,plotindex)), squeeze(r_all_NSCRK(3,plotindex)),linestyles{4},LineWidth=1.5,MarkerSize=10.0);
plot3(squeeze(r_all_NSCM(1,plotindex)), squeeze(r_all_NSCM(2,plotindex)), squeeze(r_all_NSCM(3,plotindex)),linestyles{5},LineWidth=1.5,MarkerSize=10.0);
plot3(squeeze(r_all_MCCS(1,plotindex)), squeeze(r_all_MCCS(2,plotindex)), squeeze(r_all_MCCS(3,plotindex)),linestyles{6},LineWidth=1.5,MarkerSize=10.0);
plot3(squeeze(r_all_MCEBS(1,plotindex)), squeeze(r_all_MCEBS(2,plotindex)), squeeze(r_all_MCEBS(3,plotindex)),linestyles{7},LineWidth=1.5,MarkerSize=10.0);
plot3(squeeze(r_all_analytical(1,plotindex)), squeeze(r_all_analytical(2,plotindex)), squeeze(r_all_analytical(3,plotindex)),linestyles{8},LineWidth=1.5,MarkerSize=10.0);
legend({'EB Piecewise';'Cosserat Piecewise';'Proposed SDE Solution';'Numerical SDE Solution - Weak RK';'Numerical SDE Solution - Milstein'; ...
    'Monte-carlo Simulation Cosserat';'Monte-carlo Simulation EB Approximation';'Analytical Solution'})
xlabel("$x$ [m]",'interpreter',"latex",'fontsize',14)
ylabel("$y$ [m]",'interpreter',"latex",'fontsize',14)
zlabel("$z$ [m]",'interpreter',"latex",'fontsize',14)
%axis equal
plotindex=1:5:length(s_all);
figure(10006)
linestyles={'-or';'--*g';'--*b';'--*c';'--*m';'-^r';'-vr'};
hold on
plot(s_all(plotindex),vecnorm(r_all_SACS(:,plotindex)-r_all_analytical(:,plotindex)),linestyles{1},LineWidth=1.5,MarkerSize=10.0)
plot(s_all(plotindex),vecnorm(r_all_NSCRK(:,plotindex)-r_all_analytical(:,plotindex)),linestyles{2},LineWidth=1.5,MarkerSize=10.0)
plot(s_all(plotindex),vecnorm(r_all_NSCM(:,plotindex)-r_all_analytical(:,plotindex)),linestyles{3},LineWidth=1.5,MarkerSize=10.0)
plot(s_all(plotindex),vecnorm(r_all_MCCS(:,plotindex)-r_all_analytical(:,plotindex)),linestyles{4},LineWidth=1.5,MarkerSize=10.0)
plot(s_all(plotindex),vecnorm(r_all_MCEBS(:,plotindex)-r_all_analytical(:,plotindex)),linestyles{5},LineWidth=1.5,MarkerSize=10.0)
plot(s_all(plotindex),vecnorm(r_all_NEBA(:,plotindex)-r_all_analytical(:,plotindex)),linestyles{6},LineWidth=1.5,MarkerSize=10.0)
plot(s_all(plotindex),vecnorm(r_all_NFCA(:,plotindex)-r_all_analytical(:,plotindex)),linestyles{7},LineWidth=1.5,MarkerSize=10.0)
legend({'Proposed SDE Solution';'Numerical SDE Solution - Weak RK';'Numerical SDE Solution - Milstein'; ...
    'Monte-carlo Simulation Cosserat';'Monte-carlo Simulation EB Approximation';'EB Piecewise';'Cosserat Piecewise'})
grid on
hold on
xlabel("Rod Coordinate $s$ [m]",'interpreter',"latex",'fontsize',14)
ylabel("Bias of the Rod Axis Estimation $Bias$ [m]",'interpreter',"latex",'fontsize',14)