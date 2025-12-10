function F=get_rhs_for_vecP(vecP,A,B)
vecP=reshape(vecP,[12,12]);
F=A*vecP+vecP*A';
for l=1:size(B,1)
    F=F+squeeze(B(l,:,:))*vecP*squeeze(B(l,:,:))';
end
F=reshape(F,[12*12,1]);
end