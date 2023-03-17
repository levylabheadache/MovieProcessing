function [R,C,reg] = RectifyStackDFT(vol, start, upscale)
[~,~,Nz] = size(vol); % Nx,Ny
R = zeros(Nz,1);
C = zeros(Nz,1);    
reg = zeros(size(vol));

%figure('WindowState','maximized');

target = vol(:,:,start); %imgaussfilt(vol(:,:,start), 3); % 
%forward
for i = start:Nz
    source = vol(:,:,i); %imgaussfilt(vol(:,:,i), 3); % 
    S = dftregistrationAlex(fft2(target),fft2(source),upscale);
    R(i) = S(1);
    C(i) = S(2);
    %subplot(1,2,1); imshow(target, []); title(sprintf('z = %i: target',i));
    %subplot(1,2,2); imshow(source, []); title(sprintf('source. [Row, Col] = [%2.1f, %2.1f]',R(i), C(i))); pause;
    target = imtranslate(source,[S(2),S(1)]);
    reg(:,:,i) = target;
end

%backwards
target = vol(:,:,start);
for i = flip(1:start)
    source = vol(:,:,i); %imgaussfilt(vol(:,:,i), 3); %
    S = dftregistrationAlex(fft2(target),fft2(source),upscale);
    R(i) = S(1);
    C(i) = S(2);
    %subplot(1,2,1); imshow(target, []); title(sprintf('z = %i: target',i));
    %subplot(1,2,2); imshow(source, []); title(sprintf('source. [Row, Col] = [%2.1f, %2.1f]',R(i), C(i))); pause;
    target = imtranslate(source,[S(2),S(1)]);
    reg(:,:,i) = target;
end
end