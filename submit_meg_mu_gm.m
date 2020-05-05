%% PREPROCESSING 
% submitting MEG-music-autism for preprocessing on the cluster


% clusterconfig('scheduler','none'); % Run everything without any clusterizing code at all
% clusterconfig('scheduler','local'); % Run everything clusterized but only to the machine you're working on (NB! NOT Hyades01!)
% clusterconfig('scheduler','cluster'); % Run everything truly clusterized
% clusterconfig('long_running',1); % for jobs with a duration > 1 hr

% clusterconfig('long_running', 0, 'slot', 4); % the job shouldn't take more than 1 hr
clusterconfig('long_running', 1); % the job shouldn't take more than 1 hr

grouplist{1} = [1 3 5 6 7 8 11 12 15 41 44 48 49 50 52 55 62 65]; % All CTRL participants, identified by their exam nr (subj0059 is excluded from MuMUFE due to technical problems, see the database (wiki) for the full report of exclusions)
grouplist{2} = [19 21 24 27 30 37 38 39 43 45 46 53 54 56 57 60 61 63 64];                             % All ASD participants, identified by their exam nr (subj0058 is excluded due to medical history, see the database (wiki) for the full report of exclusions)

groupNames = {'CTRL','ASD'};
groupIDs = [1 2];

% input-structure specific to MEG-MuMUFE 
for h = 1:length(groupIDs)
    input{h} = groupIDs(h);
end

cell2mat(input(:))

jobid = job2cluster(@meg_mu_grandaverage, input);


% % example function
% function out = my_own_function(patient)
%   out = [];
%   filedir = fullfile('/my/file/location', patient);
%  spm_something(filedir);
 
% clusterconfig('scheduler', [])     % [] = no clusterization; 'local' = clusterized on the local machine; 'cluster' = truly clusterized
% clusterconfig('wait', 1);          % in order for results = jobresults(job_id) not to be submitted until the jub has actually finished
%  
%  
%  
% e=joberrors(job_id);
% o=joboutput(job_id);

% % kill a job
% jobdestroy(job_id);

% % who's running what?
% clusterjobs();

