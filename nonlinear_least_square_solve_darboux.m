function [u_k,v_k]=nonlinear_least_square_solve_darboux(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,epsilon_p_helix,epsilon_p_straight)

helpmat=zeros([n_helix+n_straight,6]);
for helix=1:n_helix
    theta=theta_helix_0+2*pi*(helix-1)'/n_helix;
    [pe_pu, pe_pv]=partial_cablestrain_partial_darboux(R_helix,omega_helix,theta,s,u0_k, v0_k, u0_k, v0_k);
    helpmat(helix,:)=[pe_pu;pe_pv];
end
for straight=1:n_straight
    theta=theta_straight_0+2*pi*(straight-1)'/n_straight;
    [pe_pu, pe_pv]=partial_cablestrain_partial_darboux(R_straight,omega_straight,theta,s,u0_k, v0_k, u0_k, v0_k);
    helpmat(n_helix+straight,:)=[pe_pu;pe_pv];
end
if n_helix+n_straight>6
    nullspace=null(helpmat'*helpmat);
    if size(nullspace,2)==0
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_helix;epsilon_p_straight];
        u_k=u0_k+[Delta_uv(1:3)];
        v_k=v0_k+[Delta_uv(4:6)];
    elseif size(nullspace,2)>0 && size(nullspace,2)<3
        helpmat(:,[4,5])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_helix;epsilon_p_straight];
        u_k=u0_k+[Delta_uv(1:3)];
        v_k=v0_k+[0;0;Delta_uv(4)];
    elseif size(nullspace,2)>0 && size(nullspace,2)==3
        helpmat(:,[3,4,5])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_helix;epsilon_p_straight];
        u_k=u0_k+[Delta_uv(1:2);0];
        v_k=v0_k+[0;0;Delta_uv(3)];
    elseif size(nullspace,2)>0 && size(nullspace,2)==4
        helpmat(:,[3,4,5,6])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_helix;epsilon_p_straight];
        u_k=u0_k+[Delta_uv(1:2);0];
        v_k=v0_k+[0;0;0];
    end
elseif n_helix+n_straight==6
    nullspace=null(helpmat'*helpmat);
    if size(nullspace,2)==0
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_helix;epsilon_p_straight];
        u_k=u0_k+[Delta_uv(1:3)];
        v_k=v0_k+[Delta_uv(4:6)];
    elseif size(nullspace,2)>0 && size(nullspace,2)<3
        helpmat(:,[4,5])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_helix;epsilon_p_straight];
        u_k=u0_k+[Delta_uv(1:3)];
        v_k=v0_k+[0;0;Delta_uv(4)];
    elseif size(nullspace,2)>0 && size(nullspace,2)==3
        helpmat(:,[3,4,5])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_helix;epsilon_p_straight];
        u_k=u0_k+[Delta_uv(1:2);0];
        v_k=v0_k+[0;0;Delta_uv(3)];
    elseif size(nullspace,2)>0 && size(nullspace,2)==4
        helpmat(:,[3,4,5,6])=[];
        Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_helix;epsilon_p_straight];
        u_k=u0_k+[Delta_uv(1:2);0];
        v_k=v0_k+[0;0;0];
    end
else
    helpmat(:,[3,4,5])=[];
    Delta_uv=(helpmat'*helpmat)\helpmat'*[epsilon_p_helix;epsilon_p_straight];
    u_k=u0_k+[Delta_uv(1:2);0];
    v_k=v0_k+[0;0;Delta_uv(3)];
end


