function [r_all_rotated,D_ik_all_rotated]=rotate_at_point(r_all,D_ik_all,rod_axis_rs,D0_ik_all, ...
    rotate_index)
r_all_rotated=0*r_all;
D_ik_all_rotated=0*D_ik_all;
for i=1:size(r_all,2)
    U= D0_ik_all(:,:,rotate_index)*D_ik_all(:,:,rotate_index)';
    r_all_rotated(:,i)=U*(r_all(:,i)-r_all(:,rotate_index))+rod_axis_rs(:,rotate_index);
    D_ik_all_rotated(:,:,i)=U*D_ik_all(:,:,i);

end
% for i = 1:size(r_all,2)
%     now_rotating_percentage=i/size(r_all,2)*100
%     parfor isample = 1:size(s_r_all,1)
%         U= D0_ik_all(:,:,rotate_index)*squeeze(s_D_ik_all(isample,:,:,rotate_index))';
%         s_r_all_rotated(isample,:,i)=U*(squeeze(s_r_all(isample,:,i)-s_r_all(isample,:,rotate_index))')...
%             +rod_axis_rs(:,rotate_index);
%         s_D_ik_all_rotated(isample,:,:,i)=U*squeeze(s_D_ik_all(isample,:,:,i));
%     end
% end
end