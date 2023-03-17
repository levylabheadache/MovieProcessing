
Nz = size(ref_vol,3);
z_ref = round(Nz/2);

metric = registration.metric.MeanSquares;
optimizer = registration.optimizer.RegularStepGradientDescent;

imView = imref2d(size(ref_vol,[1,2]));

%temp_ref = ref_vol(:,:,z); figure;
temp_trim = [0,0]; % y, x
dX = zeros(Nz,1); dY = zeros(Nz,1); 

blank_tform = imregtform( ref_vol(:,:,z), ref_vol(:,:,z), 'translation', optimizer, metric );

final_tform = repmat(blank_tform, 1, Nz);
% middle to bottom
for z = z_ref:Nz-1 
    z_neigh = z+1;
    temp_ref = ref_vol(:,:,z);
    %temp_ref = temp_ref(1:end-temp_trim(1), temp_trim(2)+1:end);
    temp_target = ref_vol(:,:,z_neigh);
    %temp_target = temp_target(1:end-temp_trim(1), temp_trim(2)+1:end);

    temp_tform(z_neigh) = imregtform( temp_ref, temp_target, 'translation', optimizer, metric );
    dX(z_neigh) = temp_tform(z_neigh).T(3);
    dY(z_neigh) = temp_tform(z_neigh).T(6);

    %final_tform(z_neigh).T = [0,0,0; 0,0,0; sum(dX(z_ref:z_neigh)), sum(dY(z_ref:z_neigh)), 0];
    final_tform(z_neigh).T(3) = sum(dY(z_ref:z_neigh));
    final_tform(z_neigh).T(6) = sum(dX(z_ref:z_neigh));
    %{
    %temp_trim = temp_trim + abs(round(temp_tform(z_neigh).T([3,6])));
    temp_frame = imwarp( ref_vol(:,:,z_neigh), temp_tform(z_neigh), 'OutputView',imView); % apply optotune dewarping
    sp(1) = subplot(3,1,1); imshow(temp_ref, []); title(sprintf('z = %i', z));
    sp(2) = subplot(3,1,2); imshow(temp_target, []); title(sprintf('z = %i', z_neigh)); 
    sp(3) = subplot(3,1,3); imshow(temp_frame, []); title(sprintf('z = %i, aligned', z_neigh)); 
    %imshow(temp_frame(1:end-temp_trim(1), abs(temp_trim(2))+1:end), []); title(sprintf('z = %i, aligned', z_neigh)); 
    linkaxes(sp,'xy');
    impixelinfo;
    pause;
    %}
end
% middle to top
for z = z_ref:-1:2
    z_neigh = z-1;
    temp_ref = ref_vol(:,:,z);
    %temp_ref = temp_ref(1:end-temp_trim(1), temp_trim(2)+1:end);
    temp_target = ref_vol(:,:,z_neigh);
    %temp_target = temp_target(1:end-temp_trim(1), temp_trim(2)+1:end);

    temp_tform(z_neigh) = imregtform( temp_ref, temp_target, 'translation', optimizer, metric );
    dX(z_neigh) = temp_tform(z_neigh).T(3);
    dY(z_neigh) = temp_tform(z_neigh).T(6);

    final_tform(z_neigh).T(3) = sum(dY(z_neigh:z_ref));
    final_tform(z_neigh).T(6) = sum(dX(z_neigh:z_ref));
    %{
    %temp_trim = temp_trim + abs(round(temp_tform(z_neigh).T([3,6])));
    sp(1) = subplot(3,1,1); imshow(temp_ref, []); title(sprintf('z = %i', z));
    sp(2) = subplot(3,1,2); imshow(temp_target, []); title(sprintf('z = %i', z_neigh)); 
    temp_frame = imwarp( ref_vol(:,:,z_neigh), temp_tform(z_neigh), 'OutputView',imView); % apply optotune dewarping    
    sp(3) = subplot(3,1,3); imshow(temp_frame, []); title(sprintf('z = %i, aligned', z_neigh)); 
    %imshow(temp_frame(1:end-temp_trim(1), abs(temp_trim(2))+1:end), []); title(sprintf('z = %i, aligned', z_neigh)); 
    linkaxes(sp,'xy');
    impixelinfo;
    pause;
    %}
end

align_vol = nan(size(ref_vol));
for z = 1:Nz
    align_vol(:,:,z) = imwarp( ref_vol(:,:,z), final_tform(z), 'OutputView',imView); % apply optotune dewarping    
end

saveastiff(uint16(align_vol), 'D:\2photon\AB40\220921_FOV2\AB40_220921_FOV2_run1\opto_test.tif')



%%
Yprofile = sum(optProjMean,[1,3]);

Yprofile = sum(optProjMean(:,:,:,1),[2,3]);
Xprofile = squeeze(sum(optProjMean(:,:,:,1),[1,3]));