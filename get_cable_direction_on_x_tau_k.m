function tau_k=get_cable_direction_on_x_tau_k(R,omega,theta,s)
alpha=(2*pi*R*omega)/sqrt(1+(2*pi*R*omega)^2);
beta=sqrt(1-alpha^2);
tau_k=[-alpha*sin(2*pi*omega*s-theta);alpha*cos(2*pi*omega*s-theta);beta];
end