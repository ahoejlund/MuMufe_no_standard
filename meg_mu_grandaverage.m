function out = meg_mu_grandaverage(input)

out = [];

%% GRAND-AVERAGE

%% Setup default function paths
addpath /projects/MINDLAB2011_14-music-autism/scripts/MEG   % search_db has been cloned from meeg-cfin github and put in our scripts-fodler 
addpath /projects/MINDLAB2011_14-music-autism/scripts/MEG/fieldtrip-20161231/
% addpath /Users/au183362/Documents/MATLAB/toolboxes/fieldtrip-20161228/

ft_defaults; 

%% Setup file paths and variables
savePath = '/projects/MINDLAB2011_14-music-autism/scratch/MEG/';

groupNames = {'CTRL','ASD'};
serieNames = {'m1','m2','ns'};
extName = 'ft';    

condNames = {
    'std', ...
    'pitch_m', ...
    'timbre_m', ...
    'loc_m', ...
    'intens_m', ...
    'slide_m', ...
    'duration_m', ...
    'rhythm_m', ...
    'pitch_ns', ...
    'timbre_ns', ...
    'loc_ns', ...
    'intens_ns', ...
    'slide_ns', ...
    'duration_ns', ...
    'rhythm_ns', ...
    'std2_ns', ...
    'std1', ...
    'std2', ...
    'std4', ...
    'std124', ...
    'std_key', ...
    'std_meter', ...
    'std_key_meter', ...
    'std1_all', ...
    'std1_ns', ...
    'std4_ns' ...
    'std124_ns', ...
    'std_key_ns', ...
    'std_meter_ns', ...
    'std_key_meter_ns', ...
    'std1_all_ns'};

%% Grand average
%% Load the data
grouplist{1} = [1 3 5 6 7 8 11 12 15 41 44 48 49 50 52 55 62 65];  % CTRL
grouplist{2} = [19 21 24 27 30 37 38 39 43 45 46 53 54 56 57 60 61 63 64];  % ASD

% grouplist{1} = [1 3 5 6 7 8 11 12 15 41 48 49 50 52 55 62 65];  % CTRL
% grouplist{2} = [19 21 24 27 30 37 38 39 43 46 53 54 56 57 60 61 63 64];  % ASD

for j = length(condNames):-1:1
    for i = length(grouplist{input}):-1:1
        load(fullfile(savePath, groupNames{input}, serieNames{3}, extName, sprintf('%04d', grouplist{input}(i)), sprintf('%04d_%s_evoked.mat', grouplist{input}(i), condNames{j})))
        data{i} = data_struct_tlck;
        clear data_struct_tlck
    end
    GM = ft_timelockgrandaverage([], data{:});
    save(fullfile(savePath, 'grandmean/ft', sprintf('%s_GM_%s.mat', groupNames{input}, condNames{j})), 'GM');
    clear data GM
end    


%% Command-line cluster submission using submit_to_cluster and the highmem.q

% for k in 1 2; do submit_to_cluster -q highmem.q -m 40G "matlab -r 'meg_mu_grandaverage($k)'"; done