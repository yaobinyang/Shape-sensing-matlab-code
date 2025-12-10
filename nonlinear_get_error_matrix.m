function [Cu_k,Cv_k]=nonlinear_get_error_matrix(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,u_k,v_k,measure_error_level,spatial_resolution)
helpmat=zeros([n_helix+n_straight,6]);
unknownerrorlevel=0;
if measure_error_level>1e-12
    unknownerrorlevel=1e-2/(measure_error_level*sqrt(spatial_resolution));
end
for helix=1:n_helix
    theta=theta_helix_0+2*pi*(helix-1)'/n_helix;
    [pe_pu, pe_pv]=partial_cablestrain_partial_darboux(R_helix,omega_helix,theta,s,u0_k, v0_k, u_k, v_k);
    helpmat(helix,:)=[pe_pu;pe_pv];
end
for straight=1:n_straight
    theta=theta_straight_0+2*pi*(straight-1)'/n_straight;
    [pe_pu, pe_pv]=partial_cablestrain_partial_darboux(R_straight,omega_straight,theta,s,u0_k, v0_k, u_k, v_k);
    helpmat(n_helix+straight,:)=[pe_pu;pe_pv];
end
if n_straight+n_helix>6
nullspace=null(helpmat'*helpmat);
if size(nullspace,2)==0
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=helpmat1(1:3,:);
    Cv_k=helpmat1(4:6,:);
elseif size(nullspace,2)>0 && size(nullspace,2)<3
    helpmat(:,[4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=helpmat1(1:3,:);
    Cv_k=[unknownerrorlevel*ones(2,n_straight+n_helix);helpmat1(3,:)];
elseif size(nullspace,2)>0 && size(nullspace,2)==3
    helpmat(:,[3,4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,n_straight+n_helix)];
    Cv_k=[unknownerrorlevel*ones(2,n_straight+n_helix);helpmat1(3,:)];
elseif size(nullspace,2)>0 && size(nullspace,2)==4
    helpmat(:,[3,4,5,6])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,n_straight+n_helix)];
    Cv_k=unknownerrorlevel*ones(3,n_straight+n_helix);
end
elseif n_straight+n_helix==6
    nullspace=null(helpmat'*helpmat);
if size(nullspace,2)==0
    helpmat1=inv(helpmat);
    Cu_k=helpmat1(1:3,:);
    Cv_k=helpmat1(4:6,:);
elseif size(nullspace,2)>0 && size(nullspace,2)<3
    helpmat(:,[4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=helpmat1(1:3,:);
    Cv_k=[unknownerrorlevel*ones(2,n_straight+n_helix);helpmat1(3,:)];
elseif size(nullspace,2)>0 && size(nullspace,2)==3
    helpmat(:,[3,4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,n_straight+n_helix)];
    Cv_k=[unknownerrorlevel*ones(2,n_straight+n_helix);helpmat1(3,:)];
elseif size(nullspace,2)>0 && size(nullspace,2)==4
    helpmat(:,[3,4,5,6])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,n_straight+n_helix)];
    Cv_k=unknownerrorlevel*ones(3,n_straight+n_helix);
end
elseif n_straight+n_helix==5
helpmat(:,3)=[];
helpmat1=inv(helpmat);
Cu_k=[helpmat1(1:2,:);zeros(1,5)];
Cv_k=helpmat1(3:5,:);
elseif n_straight+n_helix==4
% helpmat(:,4:5)=[];
% helpmat1=inv(helpmat);
% Cu_k=helpmat1(1:3,:);
% Cv_k=[zeros(2,4);helpmat1(4,:)];
    helpmat(:,[3,4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,n_straight+n_helix)];
    Cv_k=[unknownerrorlevel*ones(2,n_straight+n_helix);helpmat1(3,:)];
elseif n_straight+n_helix==3
    helpmat(:,[3,4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,n_straight+n_helix)];
    Cv_k=[unknownerrorlevel*ones(2,n_straight+n_helix);helpmat1(3,:)];
elseif n_straight+n_helix==2
helpmat(:,3:6)=[];
helpmat1=inv(helpmat);
Cu_k=[helpmat1(1:2,:);zeros(1,2)];
Cv_k=[zeros(3,2)];
elseif n_straight+n_helix==1
helpmat(:,2:6)=[];
helpmat1=inv(helpmat);
Cu_k=[helpmat1(1,:);zeros(2,1)];
Cv_k=[zeros(3,1)];
end
Cu_k=Cu_k*measure_error_level*sqrt(spatial_resolution);
Cv_k=Cv_k*measure_error_level*sqrt(spatial_resolution);
end