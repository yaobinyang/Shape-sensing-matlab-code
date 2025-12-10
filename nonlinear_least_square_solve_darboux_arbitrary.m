function [u_k,v_k]=nonlinear_least_square_solve_darboux_arbitrary(R_cables,omega_cables,theta_cables,n_cables,s,u0_k,v0_k,epsilon_p_cables)

helpmat=zeros([n_cables,6]);
for cable=1:n_cables
    omega=omega_cables(cable);
    theta=theta_cables(cable);
    R=R_cables(cable);
    [pe_pu, pe_pv]=partial_cablestrain_partial_darboux(R,omega,theta,s,u0_k, v0_k, u0_k, v0_k);
    helpmat(cable,:)=[pe_pu;pe_pv];
end

if n_cables>6
    nullspace=null(helpmat'*helpmat);
    if size(nullspace,2)==0
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_cables];
        u_k=u0_k+[Delta_uv(1:3)];
        v_k=v0_k+[Delta_uv(4:6)];
    elseif size(nullspace,2)>0 && size(nullspace,2)<3
        helpmat(:,[4,5])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_cables];
        u_k=u0_k+[Delta_uv(1:3)];
        v_k=v0_k+[0;0;Delta_uv(4)];
    elseif size(nullspace,2)>0 && size(nullspace,2)==3
        helpmat(:,[3,4,5])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_cables];
        u_k=u0_k+[Delta_uv(1:2);0];
        v_k=v0_k+[0;0;Delta_uv(3)];
    elseif size(nullspace,2)>0 && size(nullspace,2)==4
        helpmat(:,[3,4,5,6])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_cables];
        u_k=u0_k+[Delta_uv(1:2);0];
        v_k=v0_k+[0;0;0];
    end
elseif n_cables==6
    nullspace=null(helpmat'*helpmat);
    if size(nullspace,2)==0
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_cables];
        u_k=u0_k+[Delta_uv(1:3)];
        v_k=v0_k+[Delta_uv(4:6)];
    elseif size(nullspace,2)>0 && size(nullspace,2)<3
        helpmat(:,[4,5])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_cables];
        u_k=u0_k+[Delta_uv(1:3)];
        v_k=v0_k+[0;0;Delta_uv(4)];
    elseif size(nullspace,2)>0 && size(nullspace,2)==3
        helpmat(:,[3,4,5])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_cables];
        u_k=u0_k+[Delta_uv(1:2);0];
        v_k=v0_k+[0;0;Delta_uv(3)];
    elseif size(nullspace,2)>0 && size(nullspace,2)==4
        helpmat(:,[3,4,5,6])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_cables];
        u_k=u0_k+[Delta_uv(1:2);0];
        v_k=v0_k+[0;0;0];
    end
else
    helpmat(:,[3,4,5])=[];
    Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_cables];
    u_k=u0_k+[Delta_uv(1:2);0];
    v_k=v0_k+[0;0;Delta_uv(3)];
end


