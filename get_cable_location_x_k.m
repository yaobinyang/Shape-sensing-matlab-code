function x_k=get_cable_location_x_k(R,omega,theta,s)
x_k=[R*cos(2*pi*omega*s-theta),R*sin(2*pi*omega*s-theta),s];
end