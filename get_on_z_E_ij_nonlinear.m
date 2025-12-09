function on_z_E_ij_nonlinear=get_on_z_E_ij_nonlinear(on_z_G_alphabeta,on_x_F_i3,on_z_D0_j3)
on_z_E_ij_nonlinear=zeros([3,3]);
for alpha=1:2
    on_z_E_ij_nonlinear=on_z_E_ij_nonlinear+on_z_G_alphabeta(:,alpha)*on_z_G_alphabeta(:,alpha)'...
        +on_x_F_i3(alpha)*(on_z_G_alphabeta(:,alpha)*on_z_D0_j3'+on_z_D0_j3*on_z_G_alphabeta(:,alpha)');
end
for i=1:3
    on_z_E_ij_nonlinear=on_z_E_ij_nonlinear+(on_x_F_i3(i))^2*(on_z_D0_j3*on_z_D0_j3');
end
on_z_E_ij_nonlinear=on_z_E_ij_nonlinear-eye(3);
on_z_E_ij_nonlinear=on_z_E_ij_nonlinear/2;
end