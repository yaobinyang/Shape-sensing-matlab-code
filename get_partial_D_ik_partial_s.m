function pD_ik_ps_all=get_partial_D_ik_partial_s(D_ik_all,s_all)

grads=gradient(s_all);

gradD=zeros(size(D_ik_all));
gradS=zeros(size(D_ik_all));
for i=1:3
    for j=1:3
       gradD(i,j,:)=gradient(squeeze(D_ik_all(i,j,:)));
       gradS(i,j,:)=grads;
    end
end
pD_ik_ps_all=gradD./gradS;
end