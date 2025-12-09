function [Cu_k,Cv_k]=nonlinear_get_error_matrix_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s,u0_k,v0_k,u_k,v_k,measure_error_level,spatial_resolution)
helpmat=zeros([n_cables,6]);
unknownerrorlevel=0;
helpmat=zeros([n_cables,6]);
for cable=1:n_cables
    omega=omega_cables(cable);
    theta=theta_cables(cable);
    R=R_cables(cable);
    [pe_pu, pe_pv]=partial_cablestrain_partial_darboux(R,omega,theta,s,u0_k, v0_k, u_k, v_k);
    helpmat(cable,:)=[pe_pu;pe_pv]/(measure_error_level*sqrt(spatial_resolution));
    % helpmat(cable,:)=[pe_pu;pe_pv];
end
if n_cables>1
nullspace=null(helpmat'*helpmat);
if size(nullspace,2)==0
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=helpmat1(1:3,:);
    Cv_k=helpmat1(4:6,:);
elseif size(nullspace,2)>0 && size(nullspace,2)<3
    helpmat(:,[4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=helpmat1(1:3,:);
    Cv_k=[unknownerrorlevel*[eye(2),ones(2,n_cables-2)];helpmat1(3,:)];
elseif size(nullspace,2)>0 && size(nullspace,2)==3
    helpmat(:,[3,4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*[eye(1),ones(1,n_cables-1)]];
    Cv_k=[unknownerrorlevel*[eye(2),ones(2,n_cables-2)];helpmat1(3,:)];
elseif size(nullspace,2)>0 && size(nullspace,2)==4
    helpmat(:,[3,4,5,6])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,n_cables)];
    Cv_k=unknownerrorlevel*ones(3,n_cables);
end
elseif n_cables==6
    nullspace=null(helpmat'*helpmat);
if size(nullspace,2)==0
    helpmat1=inv(helpmat);
    Cu_k=helpmat1(1:3,:);
    Cv_k=helpmat1(4:6,:);
elseif size(nullspace,2)>0 && size(nullspace,2)<3
    helpmat(:,[4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=helpmat1(1:3,:);
    Cv_k=[unknownerrorlevel*ones(2,n_cables);helpmat1(3,:)];
elseif size(nullspace,2)>0 && size(nullspace,2)==3
    helpmat(:,[3,4,5])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,n_cables)];
    Cv_k=[unknownerrorlevel*ones(2,n_cables);helpmat1(3,:)];
elseif size(nullspace,2)>0 && size(nullspace,2)==4
    helpmat(:,[3,4,5,6])=[];
    helpmat1=(helpmat'*helpmat)\helpmat';
    Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,n_cables)];
    Cv_k=unknownerrorlevel*ones(3,n_cables);
end
elseif n_cables==5
helpmat(:,3)=[];
helpmat1=inv(helpmat);
Cu_k=[helpmat1(1:2,:);zeros(1,5)];
Cv_k=helpmat1(3:5,:);
elseif n_cables==4
% helpmat(:,4:5)=[];
% helpmat1=inv(helpmat);
% Cu_k=helpmat1(1:3,:);
% Cv_k=[zeros(2,4);helpmat1(4,:)];
helpmat(:,[3,4,5])=[];
helpmat1=(helpmat'*helpmat)\helpmat';
Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,4)];
Cv_k=[unknownerrorlevel*ones(2,4);helpmat1(3,:)];
elseif n_cables==3
helpmat(:,[3,4,5])=[];
helpmat1=(helpmat'*helpmat)\helpmat';
Cu_k=[helpmat1(1:2,:);unknownerrorlevel*ones(1,3)];
Cv_k=[unknownerrorlevel*ones(2,3);helpmat1(3,:)];
elseif n_cables==2
helpmat(:,3:6)=[];
helpmat1=inv(helpmat);
Cu_k=[helpmat1(1:2,:);zeros(1,2)];
Cv_k=[zeros(3,2)];
elseif n_cables==1
helpmat(:,2:6)=[];
helpmat1=inv(helpmat);
Cu_k=[helpmat1(1,:);zeros(2,1)];
Cv_k=[zeros(3,1)];
end
end