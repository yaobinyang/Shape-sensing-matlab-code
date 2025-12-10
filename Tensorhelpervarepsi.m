function varepsi=Tensorhelpervarepsi(i,j,k)
E=eye(3);
varepsi=cross(E(i,:),E(j,:))*E(k,:)';
end