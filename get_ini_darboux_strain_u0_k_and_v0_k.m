function [u0_k,v0_k]=  get_ini_darboux_strain_u0_k_and_v0_k(D0_ik,pD_ik_ps)
v0_k=[0;0;1];
u0_k=[0;0;0];
for j=1:3
    u0_k(j)=0;
    for k=1:3
        for l=1:3
            u0_k(j)=u0_k(j)+0.5*Tensorhelpervarepsi(j,k,l)*(D0_ik(:,l)'*pD_ik_ps(:,k));
        end
    end
end
end
