function [fileNames, filePaths] = FileFinder( mainDir, varargin )
%FileFinder Finds files within a directory which contain the "filetype" string. Useful for picking out files based on filetype (such as .xlsx).
% same as fileFind, except returns 2 separate cell arrays instead of one
% If filetype is not a string, FileFind will pick out folders only
% -----Inputs------
% maindir = string containing address of main directory to search
% filetype = string specifying the type of files to find (optional, default returns only folders. Use '*' to retain all non-directory files of any type)
% keepExt = logical whether or not to retain the extension in 'files' ( optional, default = false. )
% criteria = function for further selection beyond file type (ex: only keep files whose names start with a number) (optional)
% -----Outputs------
% files = cell array of matching files' names (first column) and their full addresses (second column)

%parse inputs
IP = inputParser;
addRequired( IP, 'maindir', @ischar )
addParameter( IP, 'type', '*', @(x) ( ischar(x) || isnumeric(x) ) )
addParameter( IP, 'keepExt', false, @islogical )
addParameter(IP, 'contains', '', @(x) ( ischar(x) || iscell(x) || isstring(x) )) % @ischar
addParameter( IP, 'criteria', [] )
parse( IP, mainDir, varargin{:} );
fileType = IP.Results.type;
keepExt = IP.Results.keepExt;
criteria = IP.Results.criteria;
containStr = IP.Results.contains;
if ~isempty(containStr) && ischar(containStr), containStr = {containStr}; end
if ispc, fileSep = '\'; else, fileSep = '/'; end
if ~strcmpi(mainDir(end), fileSep), mainDir = [mainDir,fileSep]; end 

% Find files of the desired type within the main dir
mainfiles = dir(mainDir); 
mainfiles(1:2) = []; % dir returns '.' and '..' even for an empty folder
files = cell( size(mainfiles,1) , 2 ); % files = [{filename} , {full path}];
if strcmpi(fileType,'*') % get all non-directories
    for i = find(~[mainfiles.isdir]) %1:length(mainfiles) 
        files{i,1} = mainfiles(i).name;
        files{i,2} = [mainDir, mainfiles(i).name]; % 
    end
    %{
    for i = 1:length(mainfiles) 
        dotPos = strfind(lower(mainfiles(i).name), '.'); % find the place within the filename at which the file keepExt begins
        if ~isempty(dotPos) %contains( lower(mainfiles(i).name), lower(filetype) )  && isempty( strfind( mainfiles(i).name, '$' ) )
            if keepExt
                files{i,1} = mainfiles(i).name;
            else
                files{i,1} = mainfiles(i).name(1:dotPos-1);
            end
            files{i,2} = [mainDir, mainfiles(i).name];
        end 
    end
    %}
elseif ischar(fileType) % get all files of with a specific extension
    for i = 1:length(mainfiles) 
        tempName = lower(mainfiles(i).name);
        extStart = strfind(tempName , lower(fileType) ); % find the place within the filename at which the file extension begins
        if ~isempty(extStart)
            if strcmpi(tempName(extStart(end):end), fileType) && isempty( strfind( mainfiles(i).name, '$' ) )  % contains( lower(mainfiles(i).name), lower(fileType) )
                if keepExt
                    files{i,1} = mainfiles(i).name;
                else
                    files{i,1} = mainfiles(i).name(1:extStart(end)-2);
                end
                files{i,2} = [mainDir, mainfiles(i).name];
            end
        end
    end
else
    for i = 1:length(mainfiles)
        if mainfiles(i).isdir && ~strcmpi(mainfiles(i).name , '.') && ~strcmpi(mainfiles(i).name , '..') 
            files{i,1} = mainfiles(i).name;
            files{i,2} = [mainDir, mainfiles(i).name,'\'];
        end
    end
end   
files(cellfun(@isempty,files(:,1)),:) = [];

% Apply contains criteria (optional)
%files = files(cellfun(@contains, files(:,1), repmat({containStr}, size(files,1), 1 ) ), :);
for cont = 1:numel(containStr)
    files = files(cellfun(@contains, files(:,1), repmat({containStr{cont}}, size(files,1), 1 ) ), :);
end

% Apply any additional criteria beyond file-type and contains
if ~isempty( criteria )
    fileind = cell2mat( cellfun( criteria, files(:,1), 'UniformOutput',false ) );
    files = files(fileind,:);
end
fileNames = files(:,1);
filePaths = files(:,2);
end