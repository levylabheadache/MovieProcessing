function edges = GetEdges(inputIm, varargin) 
% edges  =[ leftEdge, rightEdge, topEdge, bottomEdge ]

IP = inputParser;
addRequired( IP, 'inputIm', @isnumeric )
addParameter( IP, 'max', [300, 300, 150, 150], @isnumeric )
addParameter( IP, 'minInt', 2000, @isnumeric ) %
addParameter( IP, 'show', false, @islogical)
parse( IP, inputIm, varargin{:} );  %  exptDir, exptName
edgeMax = IP.Results.max;
minInt = IP.Results.minInt;
show = IP.Results.show;

% Determine edges from an image stack bounded by 0s
imSize = size(inputIm);
binIm = true( imSize );
binIm(inputIm < minInt) = false;

minFreq = 3;
% Left
[~, leftInd] = max( binIm, [], 2 );
[leftVals, ~, iLeft] = unique( leftInd );
leftValFreq = accumarray(iLeft, 1); 
leftVals( leftValFreq < minFreq ) = NaN;
leftVals( leftVals > edgeMax(1) ) = NaN;
leftEdge = max( leftVals );

% Right
[~, rightInd] = max( flip(binIm,2), [], 2 );
[rightVals, ~, iRight] = unique( rightInd );
rightValFreq = accumarray(iRight, 1); 
rightVals( rightValFreq < minFreq ) = NaN;
rightVals( rightVals > edgeMax(2) ) = NaN;
rightEdge = max( rightVals );

% Top
[~, topInd] = max( binIm, [], 1 );
[topVals, ~, iTop] = unique( topInd );
topValFreq = accumarray(iTop, 1); 
topVals( topValFreq < minFreq ) = NaN;
topVals( topVals > edgeMax(3) ) = NaN;
topEdge = max( topVals );

% Bottom
[~, bottomInd] = max( flip(binIm,1), [], 1 );
[bottomVals, ~, iBottom] = unique( bottomInd );
bottomValFreq = accumarray(iBottom, 1); 
bottomVals( bottomValFreq < minFreq ) = NaN;
bottomVals( bottomVals > edgeMax(4) ) = NaN;
bottomEdge = max( bottomVals );

% Final edges
edges = [leftEdge, rightEdge, topEdge, bottomEdge] + 1;

% Show the results
if show
    [ShowEdgeFig, ShowEdgeAx] = ShowEdges(edges, inputIm, binIm);
    subplot(ShowEdgeAx(1)); title( sprintf('minInt = %i', minInt ) );
    %{
    opt = {[0.01,0.03], [0.07,0.04], [0.02,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
    figure('WindowState','maximized'); 
    sp(1) = subtightplot(1,3,1,opt{:});
    imshow( inputIm, [] );  % imadjust(inputIm)
    hold on;  
    line( edges(1)*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
    line( (imSize(2)-edges(2))*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
    line( [1, imSize(2)], edges(3)*[1,1], 'color','r', 'LineStyle','--');
    line( [1, imSize(2)], (imSize(1)-edges(4))*[1,1], 'color','r', 'LineStyle','--');
    title( sprintf('minInt = %i', minInt ) );
    
    sp(2) = subtightplot(1,3,2,opt{:});
    imshow( imadjust(inputIm), [] );  % imadjust(inputIm)
    hold on;  
    line( edges(1)*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
    line( (imSize(2)-edges(2))*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
    line( [1, imSize(2)], edges(3)*[1,1], 'color','r', 'LineStyle','--');
    line( [1, imSize(2)], (imSize(1)-edges(4))*[1,1], 'color','r', 'LineStyle','--');
    title( 'Contrast Enhanced' );

    sp(3) = subtightplot(1,3,3,opt{:});
    imshow( binIm ); hold on; 
    impixelinfo;
    line( edges(1)*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
    line( (imSize(2)-edges(2))*[1,1], [1, imSize(1)], 'color','r', 'LineStyle','--');
    line( [1, imSize(2)], edges(3)*[1,1], 'color','r', 'LineStyle','--');
    line( [1, imSize(2)], (imSize(1)-edges(4))*[1,1], 'color','r', 'LineStyle','--');
    title( sprintf('edges: [L, R, T, B] = [%i, %i, %i, %i]', edges ) );
    linkaxes(sp,'xy');
    %}
end
end