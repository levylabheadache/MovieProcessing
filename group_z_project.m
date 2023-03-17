function [mov_zproj] = group_z_project(mov,group_z_value);

sz = size(mov);
% groupZ_value = sz(3)/slice_number;
slice_number = floor(sz(3)/group_z_value);
if slice_number*group_z_value ~= sz(3)
    display(['Doesn"t equally divide so dropping ' num2str(sz(3)-slice_number*group_z_value) ' frames from end'])
    mov = mov(:,:,1:slice_number*group_z_value);
end
mov2 = reshape(mov,sz(1),sz(2),group_z_value,slice_number);
mov_zproj = mean(mov2,3);
% std_zproj = std(mov2,[],3);
% SE_zproj = (std(mov2,[],3))./slice_number;
mov_zproj = squeeze(mov_zproj); % Group Z project by appropriate number
% std_zproj = squeeze(std_zproj);
% SE_zproj = squeeze(SE_zproj);
