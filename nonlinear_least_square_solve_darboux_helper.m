function errorsquare=nonlinear_least_square_solve_darboux_helper(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,epsilon_p_helix,epsilon_p_straight,u_k,v_k)
errorsquare=0;
for helix=1:n_helix
    theta=theta_helix_0+2*pi*(helix-1)'/n_helix;
    errorsquare=errorsquare+(direct_compute_cable_strain_E_tau_tau(R_helix,omega_helix,theta,s,u0_k,v0_k,u_k,v_k)-epsilon_p_helix(helix))^2;
end
for straight=1:n_straight
    theta=theta_straight_0+2*pi*(straight-1)'/n_straight;
    errorsquare=errorsquare+(direct_compute_cable_strain_E_tau_tau(R_straight,omega_straight,theta,s,u0_k,v0_k,u_k,v_k)-epsilon_p_straight(straight))^2;
end