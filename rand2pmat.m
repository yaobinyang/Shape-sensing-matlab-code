function R_2p_mat=rand2pmat(noisenumber, Delta_s)
R=rand([noisenumber noisenumber]);
R_2p_mat=zeros([noisenumber noisenumber]);
tmp=sign(R-0.5);
for i=1:noisenumber
    for j=1:noisenumber
        tmp(i,j)=tmp(i,j)*(i<j)*Delta_s;
        tmp(i,j)=tmp(i,j)-(tmp(i,j)+Delta_s)*(i==j);
        tmp(i,j)=tmp(i,j)-(tmp(i,j)+tmp(j,i))*(i>j);
    end
end
R_2p_mat=tmp;
end