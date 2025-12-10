function J=get_jacobain_J(x,u0_k)
J=1;
for alpha=1:2
    for beta=1:2
        J=J+Tensorhelpervarepsi(alpha,beta,3)*x(beta)*u0_k(alpha);
    end
end
end