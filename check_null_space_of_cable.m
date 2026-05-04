function nullspace=check_null_space_of_cable(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,u_k,v_k)
helpmat=zeros([n_helix+n_straight,6]);
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
helpmat
inv(helpmat'*helpmat)
nullspace=null(helpmat'*helpmat);
end