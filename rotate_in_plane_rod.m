function [r_all_rotated,D_ik_all_rotated]=rotate_in_plane_rod(r_all,D_ik_all,rod_axis_rs)
axis_vec_orig=rod_axis_rs(:,end)-rod_axis_rs(:,1);
axis_vec_orig=axis_vec_orig/norm(axis_vec_orig);
axis_vec_deformed=r_all(:,end)-r_all(:,1);
axis_vec_deformed=axis_vec_deformed/norm(axis_vec_deformed);
GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;norm(cross(A,B)) dot(A,B)  0;0 0 1];
FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
UU = @(Fi,G) Fi*G*inv(Fi);
U = UU(FFi(axis_vec_deformed,axis_vec_orig), GG(axis_vec_deformed,axis_vec_orig));
r_all_rotated=0*r_all;
D_ik_all_rotated=0*D_ik_all;
for i=1:size(r_all,2)
r_all_rotated(:,i)=U*(r_all(:,i)-r_all(:,1))+r_all(:,1);
D_ik_all_rotated(:,:,i)=U*D_ik_all(:,:,i);
end
end