function Y = MovingPercentile(X, pct, win, alignment)
% X = matrix to calculate moving percentile on (along first dimension)
% pct = percentile to calculate
% window = width of window to calculate percentile of
% alignment = 'pre', 'post' or 'center' position of window relative to index
tic 
if mod(win,2) == 0
    warning('k must be an odd integer, adding one');
    win = win + 1;
end
Xsize = size(X);
Ndim = ndims(X);
colonsCell = repmat({':'},[1,Ndim-1]);
Ypar = repmat( {nan([1,Xsize(2:end)])}, Xsize(1), 1 );%cell(1,sizeX(1));
switch alignment
    case 'pre'
         parfor ii = 1:Xsize(1)
            idx = ii-win:ii-1; 
            idx = idx((idx > 0) & (idx <= Xsize(1))); 
            if numel(idx) > floor(win/2)%~isempty(idx)
               X_window = X(idx, colonsCell{:});
               Ypar{ii} = prctile(X_window, pct, 1); % perform percentile action
            end
         end
    case 'post'
        parfor ii = 1:Xsize(1)
            idx = ii+1:ii+win; 
            idx = idx((idx > 0) & (idx <= Xsize(1))); 
            if numel(idx) > floor(win/2)%~isempty(idx)
               X_window = X(idx, colonsCell{:});
               Ypar{ii} = prctile(X_window, pct, 1); 
            end
        end
    case 'center'
        parfor ii = 1:Xsize(1)
            idx = ii-floor(win/2):ii+floor(win/2); 
            idx = idx((idx > 0) & (idx <= Xsize(1))); 
            if numel(idx) > floor(win/2)%~isempty(idx)
               X_window = X(idx, colonsCell{:});
               Ypar{ii} = prctile(X_window, pct, 1); 
            end
        end
end
Y = vertcat( Ypar{:} );
toc
%
figure; plot( [X(:,1), Y(:,1)] )

end
