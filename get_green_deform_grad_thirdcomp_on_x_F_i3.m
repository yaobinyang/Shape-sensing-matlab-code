function on_x_F_i3=get_green_deform_grad_thirdcomp_on_x_F_i3(x,uk,vk,J)
on_x_F_i3=vk;
for i=1:3
for j=1:3
    on_x_F_i3(i)=on_x_F_i3(i)-(Tensorhelpervarepsi(i,j,3)*x(j)*uk(3));
    for k=1:3
        on_x_F_i3(i)=on_x_F_i3(i)-Tensorhelperdelta(i,3)*Tensorhelpervarepsi(3,j,k)*x(j)*uk(k);

    end
end
on_x_F_i3(i)=on_x_F_i3(i)/J;
end
end