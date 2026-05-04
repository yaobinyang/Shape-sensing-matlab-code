function on_z_G_ij=get_contravar_on_z_G_ij(x,u0k,on_z_D0_ik,J)
on_z_G_ij=on_z_D0_ik;
for i=1:3
for j=1:3
    on_z_G_ij(:,i)=on_z_G_ij(:,i)+Tensorhelpervarepsi(i,j,3)*x(j)*u0k(3)*on_z_D0_ik(:,3)/J;
end
on_z_G_ij(:,i)=on_z_G_ij(:,i)+Tensorhelperdelta(i,3)*(on_z_D0_ik(:,3)/J-on_z_D0_ik(:,3));
end
% on_z_g_i
end