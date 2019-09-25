
function [file,subject] = search_db(S)

% Full path to raw files in Storm database at CFIN, Aarhus University
% FORMAT file = search_db(S)
% 
% Input :
%
% S.project  - project Code or ID
% S.modality - modality   
% S.study    - study                (optional)
% S.series   - series               (optional)
% S.status   - subject identifier   (optional)
% S.aux      - series aux files     (optional)
% S.info     - series info          (optional)
%
%
% Output :
% 
% file  - full path to raw files    (cell hierarchy)
% index - subject no. for indexing  (cell array)
%
% _____________________________________________________________
% (C) 2015-17 Martin Dietz, CFIN, Aarhus University



% mount point
% --------------------------------------------------------
 
basepath = '/projects'; 

if isunix && exist('/Volumes/projects','dir')
    basepath = '/Volumes/projects';
elseif ispc && exist('H:','dir');
    basepath = 'Z:';
end


% defaults
% --------------------------------------------------------

server   = 'http://hyades00.pet.auh.dk';
module   = '/modules/StormDb/extract/';
format   = '&format=db';
included = '&included=1';
remove   = num2str(double(numel(basepath) > 0));


% modality types
% --------------------------------------------------------

modtype = urlread(strcat(server, module,'modalitytypes'));
modtype = textscan(modtype,'%s');
modtype = modtype{:};


% project ID and Code
% ------------------------------------------------------

projCode = '';

if isnumeric(S.project) || ~isnan(str2double(S.project))
    projID = num2str(S.project);
else
    projID   = '0';
    projCode = S.project;
end

projID   = strcat('&projectId=',projID);
projCode = strcat('&projectCode=',projCode);


% login
% --------------------------------------------------------

h = dblogin(server, module);


% subjects
% ------------------------------------------------------

subject = urlread(strcat(server,module,'subjects?', ...
          h,projID,projCode,included));

subject = textscan(subject,'%s','delimiter','\n');
subject = subject{:};


% if password has changed
      
if any(cell2mat(regexp(subject,'error')))
    if isunix
        delete(fullfile(getenv('HOME'),'.stormdblogin'))
    elseif ispc
        delete(fullfile(getenv('USERPROFILE'),'stormdblogin'))
    end
    
    h = dblogin(server, module);
    
    subject = urlread(strcat(server,module,'subjects?', ...
          h,projID,projCode,included));
    
    subject = textscan(subject,'%s','delimiter','\n');
    subject = subject{:};  
end


% subject meta (group or treatment)

prop = urlread(strcat(server,module,'subjectmetas?',h,projID,projCode));
prop = unpack(textscan(prop,'%s','delimiter','\n'));

k = isfield(S,prop);

if any(k)
    s = false(size(subject));
    
    for i = 1:length(subject)
        x = urlread(strcat(server,module,'subjectmeta?', ...
            h,projID,projCode,'&subjectNo=',subject{i}, ...
            '&prop=',urlencode(prop{k})));
        
        s(i) = any(strcmpi(S.(prop{k}),x));
    end
    
    % update
    
    subject = subject(s);
end    


% study time

study = cell(size(subject));

for i = 1:length(subject)
    
    x = urlread(strcat(server,module,'studies?', ...
        h,projID,projCode,'&subjectNo=',subject{i}, ...
        format,included));
    
    if ~isempty(x)
        study(i) = textscan(x,'%s','delimiter','\n');
    else
        study{i} = char;
    end
end


% study meta

if isfield(S,'meta')
    prop = urlread(strcat(server,module,'studymetas?',...
           h,projID,projCode));    
    
    for i = 1:length(subject)
        m = false(size(study{i}));
        
        for j = 1:length(m)
            
            x = urlread(strcat(server,module,'studymeta?', ...
                h,projID,projCode,'&subjectNo=',subject{i}, ...
                '&study=',urlencode(study{i}{j}),...
                '&prop=',urlencode(prop)));
        
            x = textscan(x,'%s','delimiter','\n');
            
            m(j) = strcmpi(char(S.meta),x{:});
        end
        
        % update
        
        study{i} = study{i}(m);
    end
end


% modality

modality = cell(size(subject));
search   = isfield(S,'modality') && any(strcmpi(S.modality,modtype));
s        = false(size(subject));

for i = 1:length(subject)
    
    k = false(size(study{i}));
    m = cell(size(k));
    
    for j = 1:length(k)
        
        x = urlread(strcat(server, module,'modalities?',...
            h,projID,projCode,'&subjectNo=',subject{i},...
            '&study=',urlencode(study{i}{j})));    
        
        if ~isempty(x)
            x = textscan(x,'%s','delimiter','\n');
        else
            x = char;
        end
        
        
        % search modality
        
        if search
            k(j) = any(strcmpi(S.modality,x{:}));
            m{j} = cellstr(S.modality);
        else
            m{j} = x{:};
        end
        
    end
    
    % update
    
    if search
        modality{i} = m(k);
        study{i}    = study{i}(k);
    else
        modality{i} = m;
    end
    
    s(i) = ~isempty(study{i});
end


% update

subject  = subject(s);
study    = study(s);
modality = modality(s);


% series

serieno   = cell(size(subject)); 
search    = isfield(S,'series');

if search
    seriedesc = cellstr(S.series);
else
    seriedesc = cellstr('');
end

for i = 1:length(subject)
    for j = 1:length(study{i})
        for k = 1:length(modality{i}{j})
            for p = 1:length(seriedesc)
                
                x = urlread(strcat(server,module,'serienumbers?',...
                    h,projID,projCode,'&subjectNo=',subject{i},...
                    '&study=',urlencode(study{i}{j}),...
                    '&modality=',modality{i}{j}{k},...
                    '&seriedescr=',seriedesc{p},included));
                
                if ~isempty(x)
                    x = textscan(x,'%s','delimiter','\n');
                else
                    x = char;
                end
                
                % sort series
                
                if search
                    serieno{i}{j,1}{k,1}{p,1} = unpack(x);
                else
                    serieno{i}{j,1}{k,1} = x{:};
                end
            end
        end
    end
end


% raw files (sorted)

file = cell(size(subject));

for i = 1:length(subject)
    for j = 1:length(study{i})
        for k = 1:length(modality{i}{j})
            for p = 1:length(serieno{i}{j}{k})
                
                x = urlread(strcat(server,module,'files?',...
                    h,projID,projCode,'&subjectNo=',subject{i},...
                    '&study=',urlencode(study{i}{j}),...
                    '&modality=',modality{i}{j}{k},...
                    '&serieNo=',serieno{i}{j}{k}{p},...
                    '&removeProjects=',remove));
                
                if ~isempty(x)
                    
                    x = textscan(x,'%s','delimiter','\n');
                    f = x{:};
                    
                    for s = 1:numel(f)
                        f{s} = fullfile(basepath,f{s});
                    end
                    
                else
                    f = char;
                end
                
                
                % unpack
                
                file{i}{j,1}{k,1}{p,1} = unpack(f);
            end
            
            file{i}{j}{k} = unpack(file{i}{j}{k});
        end
        
        file{i}{j} = unpack(file{i}{j});
    end
    
    file{i} = unpack(file{i});
end



% login
% ------------------------------------------

function h = dblogin(server, module)

if isunix
    P = fullfile(getenv('HOME'),'.stormdblogin');
elseif ispc
    P = fullfile(getenv('USERPROFILE'),'stormdblogin');
end

if ~exist(P,'file')
    
    u = input('Username: ','s');
    p = input('Password: ','s');
   
    h = urlread(strcat(server,module,...
        'login/username/',u,'/password/',p));
    
    if ~any(strfind(h,'@'))
        error('Could not log in to Hyades')
    end
    
    f = fopen(P,'w');
    fprintf(f,'%s',h);
    fclose(f);
else
    f = fopen(P,'r');
    h = fgetl(f);
    fclose(f);
end

if ispc
    fileattrib(P,'+h');
end


% unpack
% -----------------------------------------

function x = unpack(x)

while iscell(x) && (length(x) == 1)
    x = x{:};
end

