function [ShowEdgeFig, sp] = ShowEdges(edges, showIm, varargin) %
if nargin > 2 
    binaryIm = varargin{1}; 
    Nsp = 3;
else
    binaryIm = []; 
    Nsp = 2;
end
imSize = size(showIm);
opt = {[0.01,0.03], [0.07,0.04], [0.02,0.02]};  % {[vert, horz], [bottom, top], [left, right] }

ShowEdgeFig = figure('WindowState','maximized');
sp(1) = subtightplot(1,Nsp,1,opt{:});
imshow( showIm, [] );  
hold on;
line( edges(1)*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
line( (imSize(2)-edges(2))*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
line( [1, imSize(2)], edges(3)*[1,1], 'color','r', 'LineStyle','--');
line( [1, imSize(2)], (imSize(1)-edges(4))*[1,1], 'color','r', 'LineStyle','--');

sp(2) = subtightplot(1,Nsp,2,opt{:});
imshow( imadjust(showIm), [] ); 
hold on;
line( edges(1)*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
line( (imSize(2)-edges(2))*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
line( [1, imSize(2)], edges(3)*[1,1], 'color','r', 'LineStyle','--');
line( [1, imSize(2)], (imSize(1)-edges(4))*[1,1], 'color','r', 'LineStyle','--');
title( sprintf('edges: [L, R, T, B] = [%i, %i, %i, %i]', edges ) );
%title( 'Contrast Enhanced' );

if Nsp > 2
    sp(3) = subtightplot(1,Nsp,3,opt{:});
    imshow( binaryIm ); hold on;
    line( edges(1)*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
    line( (imSize(2)-edges(2))*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
    line( [1, imSize(2)], edges(3)*[1,1], 'color','r', 'LineStyle','--');
    line( [1, imSize(2)], (imSize(1)-edges(4))*[1,1], 'color','r', 'LineStyle','--');
    title( sprintf('edges: [L, R, T, B] = [%i, %i, %i, %i]', edges ) );
end
linkaxes(sp,'xy');
impixelinfo;
end