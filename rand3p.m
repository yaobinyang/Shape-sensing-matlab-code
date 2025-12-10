function R_3p=rand3p(sizeofvec,Delta_s)
R=rand(sizeofvec);
R_3p=zeros(sizeofvec);
R_3p(R<=1/6)=-sqrt(3*Delta_s);
R_3p(R>5/6)=sqrt(3*Delta_s);
end