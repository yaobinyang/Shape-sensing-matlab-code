function [pe_pu, pe_pv]=partial_cablestrain_partial_darboux(R,omega,theta,s,u0_k, v0_k, uk, vk)
alpha=(2*pi*R*omega)/sqrt(1+(2*pi*R*omega)^2);
beta=sqrt(1-alpha^2);
alpha_divid_R=(2*pi*omega)/sqrt(1+(2*pi*R*omega)^2);
x_k=[R*cos(2*pi*omega*s+theta),R*sin(2*pi*omega*s+theta),s];
J=get_jacobain_J(x_k,u0_k);
pe_pu1=x_k(2)*beta^2*(vk(3)+uk(1)*x_k(2)-uk(2)*x_k(1))/J^2;
pe_pu2=-x_k(1)*beta^2*(vk(3)+uk(1)*x_k(2)-uk(2)*x_k(1))/J^2;
pe_pu3=beta^2*(x_k(1)*vk(2)-x_k(2)*vk(1)+uk(3)*R^2)/J^2+alpha*R/J-u0_k(3)*beta*R^2/J;
pe_pv1=u0_k(3)*x_k(2)*beta/J^2-alpha_divid_R*x_k(2)/(J)+beta*(vk(1)-uk(3)*x_k(2));
pe_pv2=-u0_k(3)*x_k(1)*beta/J^2+alpha_divid_R*x_k(1)/(J)+beta*(vk(2)+uk(3)*x_k(1));
pe_pv3=beta^2*(vk(3)+uk(1)*x_k(2)-uk(2)*x_k(1))/J^2;
pe_pu=[pe_pu1;pe_pu2;pe_pu3];
pe_pv=[pe_pv1;pe_pv2;pe_pv3];
end