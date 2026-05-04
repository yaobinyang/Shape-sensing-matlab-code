function [u0_k,v0_k]=  get_ini_darboux_strain_u0_k_and_v0_k(D0_ik,pD_ik_ps)
v0_k=[0;0;1];
u0_k=[0;0;0];
for m=1:3
    u0_k(m)=0;
    for l=1:3
        for k=1:3
            for i=1:3
                for j=1:3
                    u0_k(m)=u0_k(m)+0.5*Tensorhelpervarepsi(i,j,l)*(D0_ik(i,k)*pD_ik_ps(j,k)*D0_ik(l,m));
                end
            end
        end
    end
end
end
