function E_tautau=on_tau_strain(on_z_E_ij,on_z_tau_k)
E_tautau=on_z_tau_k'*on_z_E_ij*on_z_tau_k;
end