function out = meg_mu_preproc(input)

%input = [1 1 1];

%% MUSIC-AUTISM project, from /projects/MINDLAB2011_14-music-autism/

% PROJECT: MINDLAB2011_14-music-autism
% DATESTAMP: 2019-03-05
% OWNER: Andreas H�jlund/Line Gebauer/Cecilie M�ller (partly based on Niels
% Chr. Hansen's preprocessing script)
% This function preprocesses already maxfiltered data
out = [];

%% Setup default function paths
addpath /projects/MINDLAB-templates/scripts/database                   % for search_database
addpath /projects/MINDLAB2011_14-music-autism/scripts/MEG   % search_db has been cloned from meeg-cfin github and put in our scripts-fodler 
% addpath /usr/local/common/matlab_toolbox/fieldtrip/latest/               % latest CFIN-version of Fieldtrip
addpath /projects/MINDLAB2011_14-music-autism/scripts/MEG/fieldtrip-20161231/

ft_defaults; 

%% Setup file paths
Pdir = '/projects/MINDLAB2011_14-music-autism/';
savePath = '/projects/MINDLAB2011_14-music-autism/scratch/MEG/';

%% This cell was run on Tuesday, Mar 4, 2014 for the present project and should NOT be run again - it copies ctc and cal-files to the project folder
% mkdir(fullfile(Pdir, 'misc/sss/'))
% cmd = sprintf('! cp /neuro/databases/sss/sss_cal_Mar11-May13.dat %smisc/sss/',Pdir);
% eval(cmd)
% 
% mkdir(fullfile(Pdir, 'misc/ctc/'))
% cmd = sprintf('! cp /neuro/databases/ctc/ct_sparse_Mar11-May13.fif %smisc/ctc/',Pdir);
% eval(cmd)

%% Setup names and query the database
% Query database for files (sorted by name)
% ------------------------------------------------
S = [];
S.project = 'MINDLAB2011_14-music-autism'; 
S.modality  = 'MEG';
S.Gruppe = sprintf('%d',input(2));  % input(3)
S.series  = {'m1', 'm2', 'ns'};

F = search_db(S);

% Query database for bad channels (sorted)
% ------------------------------------------------
S = [];
S.project = 'MINDLAB2011_14-music-autism'; 
S.modality  = 'MEG';
S.Gruppe = sprintf('%d',input(2));  % input(3)
S.info = 'badchannels';

B = search_db(S);


%% Setup subject and variable specifics
% g = input(1); % participant's sequential number
n = input(1); % participant's id number
o = input(2); % group number
first = 200;
keys = 300;
rhys = 400;
meters = 500;
% stdrds = 600; 
groupNames = {'CTRL','ASD'};
serieNames = {'m1','m2','ns'};
extName = 'ft';    
makeDir = 1;                                                                % (default = 0); option to create the extName-folders if needed (i.e. on first run of a new analysis)
dirCount = 1;                                                               % counts the number of dirs created (and cancels the if-loop when it reaches 4)
loadFT = []; %'bfmadcefspm12'; % [];                                                             % (default = empty); name of the file that should be loaded as the first step
loadExtName = []; %extName; % [];                                                          % (default = empty); MUST BE FILLED OUT IF loadFT is specified; extName for the FT-file that potentially needs to be loaded as the first step of the FT analysis (i.e. it's a shortcut if certain early processing steps can be reused in a subsequent re-run of analyses)
allConcat = 1; %1; %0;                                                              % (default = 0); option to specify that all conditions have been merged
maxTempFileName = 'tsss_mc';                                                                            % File name for output of first step of max-filtering (i.e. tsss and mc)
maxFullFileName = 'tsss_mc_trans_ds';                                                                          % File name for "fully" max-filtered data (i.e. including tsss, mc, trans_default and downsampling), subject number and serie number are added as prefixes
mainFolder = 3;                                                                 % specifies which of the 4 series folders should be used for subsequent processing/saving of the finally merged files
mainCond = 'ns';
mTrials = 3456;
nsTrials = 3427;
trialCheck = cell(length(F{1}));
serie_ids = {[1 2], 3};

% Cond labels
condList = {
    'std', ...
    'pitch-m', ...
    'timbre-m', ...
    'loc-m', ...
    'intens-m', ...
    'slide-m', ...
    'duration-m', ...
    'rhythm-m', ...
    'pitch-ns', ...
    'timbre-ns', ...
    'loc-ns', ...
    'intens-ns', ...
    'slide-ns', ...
    'duration-ns', ...
    'rhythm-ns', ...
    'std2-ns', ...
    'std1', ...
    'std2', ...
    'std4', ...
    'std124', ...
    'std-key', ...
    'std-meter', ...
    'std-key-meter', ...
    'std1-all', ...
    'std1-ns', ...
    'std4-ns' ...
    'std124-ns', ...
    'std-key-ns', ...
    'std-meter-ns', ...
    'std-key-meter-ns', ...
    'std1-all-ns'};

% Cond labels by series
condSeries = {
    {'std', 'pitch-m', 'timbre-m', 'loc-m', 'intens-m', 'slide-m', 'duration-m', 'rhythm-m', ...
    'std1', 'std2', 'std4', 'std124', 'std-key', 'std-meter', 'std-key-meter', 'std1-all'}, ...
    {'pitch-ns', 'timbre-ns', 'loc-ns', 'intens-ns', 'slide-ns', 'duration-ns', 'rhythm-ns', ...
    'std2-ns', 'std1-ns', 'std4-ns', 'std124-ns', 'std-key-ns', 'std-meter-ns', 'std-key-meter-ns', 'std1-all-ns'}};

%%
% condsid = {'aud', 'vis', 'AV'};
% 
% cond_names = {'bslH', 'bslL', 'visU', 'highU', 'highD', 'visD', ...
%     'lowU', 'lowD', 'std', 'bslCmb', 'sameCmb', 'oppCmb', 'visCmb'};
% mmn_names  = {'MMN_bslH', 'MMN_bslL', 'MMN_visU', 'MMN_highU', 'MMN_highD', ...
%     'MMN_visD', 'MMN_lowU', 'MMN_lowD', 'MMN_dummy', 'MMN_bsl', 'MMN_same', ... 
%     'MMN_opp', 'MMN_vis', 'diff_same_bsl', 'diff_opp_bsl', 'diff_same_opp'};
% std_names = {'stdbslH', 'stdbslL', 'stdvisU', 'stdhighU', 'stdhighD', ...
%     'stdvisD', 'stdlowU', 'stdlowD', 'stddummy', 'stdbslCmb', 'stdsameCmb', ...
%     'stdoppCmb', 'stdvisCmb'};
% 
% triggers = {[1 2 23], [3 6 23], [1:8 23]}; % Trigger 23 is every 3rd standard (not color) which doesn't involve a "visual event" like the deviant.
% trigs = triggers{input(2)};
% numbtrls = [400 400 1600];


%% TRIAL DEFINITION
% for j = 1:length(condList)  %length(unique(labels_concat))   %length(condList)
%     if ismember(condList{input(3)}, condSeries{1})
%         series_idx = [1 2];
%     elseif ismember(condList{input(3)}, condSeries{2})
%         series_idx = 3;
%     end
    series_idx = serie_ids{input(3)};
    for h = series_idx %1:length(serieNames)
        if ~exist(fullfile(savePath, groupNames{o}, serieNames{3}, extName, sprintf('%04d', n)), 'dir')
            mkdir(fullfile(savePath, groupNames{o}, serieNames{3}, extName, sprintf('%04d', n)))
        end
        %%% Define trial structure for ALL trials
        cfg = [];
        cfg.dataset = fullfile(savePath, groupNames{o}, serieNames{h}, ...
            'tsss', sprintf('%04d_%s_tsss_mc.fif',n,serieNames{h}));
        cfg.adjust_timeline = -0.014; % measured by AHN & CB 2013-12-18 (matches AHN's memory of 13.6 ms from spring 2011, and the 14 ms measurement from 2013-11-01)
        cfg.trialdef.prestim = 0.100+cfg.adjust_timeline; % preparing for the 14 ms offset (= 396)
        cfg.trialdef.poststim = 0.410-cfg.adjust_timeline; % SOA 205 (prep for 14 ms offset = 424)
        cfg.trialdef.eventtype = 'STI101';  % Specifying where markers are located.
        cfg.trialdef.eventvalue = 1:133; % Trigger values to include.
        cfg = ft_definetrial(cfg);
        orig_cfg = cfg; % saving the originally defined trial structure before shaping it to match our epoching criteria
        
        %%% Weed out spurious, "wrong", tranposition and rhythm trials
        % spurious
        spurious = find(diff(cfg.trl(:,1))==1);
        if ~isempty(spurious)
            cfg.trl(spurious,:) = [];
            save(fullfile(savePath, groupNames{o}, serieNames{3}, extName, sprintf('%04d', n), ...
                sprintf('index_%04d_%s',n,serieNames{h})),'spurious','cfg');
        end
        
        % "wrong"
        if strcmp(serieNames{h},'ns')
            % spuriuos 17-trial in the ns-sequence at 221; spurious 8-and-10-mix-up in lines 481 and 489 - thus deleting the 4 bars (478-493) from
            % immdediately before to immediately after (see sequence_ns_final-comments.txt in docs); as well as 1810-1817 due to spurious 6'er in 1810
            cfg.trl([221:225 478:493 1810:1817],:) = [];
        end
        
        % transposition
        stds = find(cfg.trl(:,4) < 13); % all std trl
        stds_next = stds(1:end-1) + 1; % trial after each std trl
        stds1_idx = find(diff([cfg.trl(stds(1:end-1), 4), cfg.trl(stds_next, 4)], [], 2)==7); % finding all 1st note stds (for which the following note must be +7)
        stds1 = stds(stds1_idx(2:end)); % discarding the very first trial
        stds4 = stds1 - 1; % trial before each 1st note std
        key_stds1 = find(diff([cfg.trl(stds1, 4), cfg.trl(stds4, 4)], [], 2)~=7); % finding all those 1st notes for which the preceding tone is not exactly +7 > key change
        cfg.trl(1,4) = cfg.trl(1,4) + first; % renaming the very first trial to +200
        cfg.trl(stds1(key_stds1),4) = cfg.trl(stds1(key_stds1),4) + keys; % renaming all key changes to +300
        
        % rhythm + "meter" deviants/standards
        durs = find((cfg.trl(:,4) > 114) & (cfg.trl(:,4) < 134)); % triggers from 115-133 (both incl.) are shortened notes (needed for rhythm deviants)
        rhythms = durs + 1; % the actual rhyhtm "deviants" are then the 4th note trials immediately following the shortened 3rd notes
        meterstds = rhythms + 1; % because the 4th notes follwing the shortened 3rd notes are NOT prolonged accordingly, the meter is offset from then onwards, and we therefore discard the followring std as well
        cfg.trl(rhythms,4) = cfg.trl(rhythms,4) + rhys; % renaming all actual rhythm deviants (i.e. 4th notes) to +400
        cfg.trl(meterstds,4) = cfg.trl(meterstds,4) + meters; % renaming all subsequent meter offset 1st note stds to +500 (hence, some will be +800 because of both key and meter change
        %                     % if necessary to account for those meter changes that are also key changes
        %                     cfg.trl(meterstds(ismember(meterstds,stds1(key_stds1))),4) = cfg.trl(meterstds(ismember(meterstds,stds1(key_stds1))),4) + meters - keys;
        
        % re-labeling
        labels = cell(size(cfg.trl,1),1);
        ex1_labels = cell(size(cfg.trl,1),1);
        ex2_labels = cell(size(cfg.trl,1),1);
        stds1_clean = stds1(~ismember(stds1, [stds1(key_stds1); meterstds(~ismember(meterstds,stds1(key_stds1)))])); % subtracting out both key change and meter offset stds
        stds2_clean = stds1 + 1; % except for the very first trial (see above)
        stds4_clean = stds1(~ismember(stds1 + 3, rhythms)) + 3; % also does not include the very first trial, i.e. the 4th note
        if ~strcmp(serieNames{h},'ns')
            stds3 = stds1 + 2; % except for the very first trial (see above)
            stds3_clean = stds3(cfg.trl(stds3, 4) < 20);
        end
        
        if ~strcmp(serieNames{h},'ns')
            labels(stds3_clean) = {'std'};
            labels(stds1_clean) = {'std1'};
            labels(stds2_clean) = {'std2'};
            labels(stds4_clean) = {'std4'};
            ex1_labels([stds1_clean; stds2_clean; stds4_clean]) = {'std124'}; % stds1 + stds2 + stds4
            labels((cfg.trl(:,4) > 299) & (cfg.trl(:,4) < 400)) = {'std-key'};
            labels((cfg.trl(:,4) > 499) & (cfg.trl(:,4) < 600)) = {'std-meter'};
            labels((cfg.trl(:,4) > 799) & (cfg.trl(:,4) < 900)) = {'std-key-meter'};
            ex2_labels([stds1_clean; find(ismember(cfg.trl(:,4), [300:399 500:599 800:899]))]) = {'std1-all'}; % ALL stds1 (i.e. also incl. key, meter and key+meter changes)
            labels((cfg.trl(:,4) > 19) & (cfg.trl(:,4) < 39)) = {'pitch-m'};
            labels((cfg.trl(:,4) > 38) & (cfg.trl(:,4) < 58)) = {'timbre-m'};
            labels((cfg.trl(:,4) > 57) & (cfg.trl(:,4) < 77)) = {'loc-m'};
            labels((cfg.trl(:,4) > 76) & (cfg.trl(:,4) < 96)) = {'intens-m'};
            labels((cfg.trl(:,4) > 95) & (cfg.trl(:,4) < 115)) = {'slide-m'};
            labels(durs) = {'duration-m'};
            labels(rhythms) = {'rhythm-m'};
        elseif strcmp(serieNames{h},'ns')
            %                         labels(stds3_clean) = {'std-ns'};
            labels(stds1_clean) = {'std1-ns'};
            labels(stds2_clean) = {'std2-ns'};
            labels(stds4_clean) = {'std4-ns'};
            ex1_labels([stds1_clean; stds2_clean; stds4_clean]) = {'std124-ns'}; % stds1 + stds2 + stds4
            labels((cfg.trl(:,4) > 299) & (cfg.trl(:,4) < 400)) = {'std-key-ns'};
            labels((cfg.trl(:,4) > 499) & (cfg.trl(:,4) < 600)) = {'std-meter-ns'};
            labels((cfg.trl(:,4) > 799) & (cfg.trl(:,4) < 900)) = {'std-key-meter-ns'};
            ex2_labels([stds1_clean; find(ismember(cfg.trl(:,4), [300:399 500:599 800:899]))]) = {'std1-all-ns'}; % ALL stds1 (i.e. also incl. key, meter and key+meter changes)
            labels((cfg.trl(:,4) > 19) & (cfg.trl(:,4) < 39)) = {'pitch-ns'};
            labels((cfg.trl(:,4) > 38) & (cfg.trl(:,4) < 58)) = {'timbre-ns'};
            labels((cfg.trl(:,4) > 57) & (cfg.trl(:,4) < 77)) = {'loc-ns'};
            labels((cfg.trl(:,4) > 76) & (cfg.trl(:,4) < 96)) = {'intens-ns'};
            labels((cfg.trl(:,4) > 95) & (cfg.trl(:,4) < 115)) = {'slide-ns'};
            labels(durs) = {'duration-ns'};
            labels(rhythms) = {'rhythm-ns'};
        end
        labels = labels(1:size(cfg.trl,1),1); % discarding any last trials/triggers that weren't registered properly by the MEG
        ex1_labels = ex1_labels(1:size(cfg.trl,1),1); % discarding any last trials/triggers that weren't registered properly by the MEG
        ex2_labels = ex2_labels(1:size(cfg.trl,1),1); % discarding any last trials/triggers that weren't registered properly by the MEG
        cfg.trl = cfg.trl(~cellfun(@isempty, labels),:); % removing the stds from the very first trial
        labels = labels(~cellfun(@isempty, labels),1); % same for the labels (so that they match the cfg.trl structure)
        ex1_labels = ex1_labels(~cellfun(@isempty, labels),1); % same for the labels (so that they match the cfg.trl structure)
        ex2_labels = ex2_labels(~cellfun(@isempty, labels),1); % same for the labels (so that they match the cfg.trl structure)
        full_trl = cfg.trl;
        temp_cfg = cfg;
        %     toc
        
        for j = 1:length(condSeries{input(3)})  %length(unique(labels_concat))   %length(condList)
            cfg = temp_cfg;
            % select only the trials pertaining to the relevant condition
            if ismember(condSeries{input(3)}{j}, {'std124', 'std124-ns'})
                cfg.trl = cfg.trl(match_str(ex1_labels, condSeries{input(3)}{j}),:);
            elseif ismember(condSeries{input(3)}{j}, {'std1-all', 'std1-all-ns'})
                cfg.trl = cfg.trl(match_str(ex2_labels, condSeries{input(3)}{j}),:);
            else
                cfg.trl = cfg.trl(match_str(labels, condSeries{input(3)}{j}),:);
            end
            trl = cfg.trl;
            
            % Baseline correction
            cfg.demean = 'yes';
            cfg.baselinewindow = [-Inf 0];
            
            data = ft_preprocessing(cfg);
            
            % ARTEFACT DETECTION
            cfg = [];
            %     cfg.datafile = fullfile(savePath, groupNames{o}, serieNames{h}, 'tsss', sprintf('%04d_%s_tsss_mc.fif',n,serieNames{h}));
            %     cfg.headerfile = fullfile(savePath, groupNames{o}, serieNames{h}, 'tsss', sprintf('%04d_%s_tsss_mc.fif',n,serieNames{h}));
            cfg.trl = trl;
            cfg.continuous = 'yes';
            cfg.memory = 'low';
            % channel selection, cutoff and padding
            cfg.artfctdef.zvalue.channel    = 'MEG';
            cfg.artfctdef.zvalue.cutoff     = 30;
            cfg.artfctdef.zvalue.trlpadding = 0;
            cfg.artfctdef.zvalue.artpadding = 0;
            cfg.artfctdef.zvalue.fltpadding = 0;
            % algorithmic parameters
            cfg.artfctdef.zvalue.cumulative    = 'yes';
            cfg.artfctdef.zvalue.medianfilter  = 'yes';
            cfg.artfctdef.zvalue.medianfiltord = 9;
            cfg.artfctdef.zvalue.absdiff       = 'yes';
            [cfg, artifact_jump] = ft_artifact_zvalue(cfg, data);
            
            % ARTEFACT REJECTION
            cfg.artfctdef.reject = 'complete'; % this rejects complete trials
            cfg.artfctdef.jump.artifact = artifact_jump;
            data = ft_rejectartifact(cfg, data);
            
            % FILTERING
            cfg = [];
            cfg.lpfilter = 'yes';
            cfg.lpfreq = 40;
%             cfg.hpfilter = 'yes';
%             cfg.hpsfreq = 1;
            data = ft_preprocessing(cfg, data);
            
            if ~strcmp(serieNames{h},'ns')
                conds{h} = data;
                full_trials{h} = full_trl;
                trials{h} = trl;
            else
                full_trials = full_trl;
                trials = trl;
            end
            
            % APPEND M1 & M2
            % Because location of sensors with respect to the head is lost with this
            % function, we save them and load them again.
            if strcmp(serieNames{h},'m2')
                grad = conds{1}.grad;
                data = ft_appenddata([],conds{:});
                data.grad = grad;
                data.subject = sprintf('%04d', n);
            end
            
            if ~strcmp(serieNames{h},'m1')
                % TIMELOCK ANALYSIS AND IMPLICITLY SELECTING DATA FOR SINGLE TRIAL ANALYSIS
                % NB! Use ft_timelockanalysis WITH cfg.keeptrials = 'yes' to get data into
                % 'rpt_chan_time' dimensions (i.e. a trial-field with just one big 3D
                % matrix rather than N cells with matrices of chan_time (where N is number
                % of trials)). This structure allows ft_math to be applied to the single
                % trials.
                
                % averaging (aka. timelock analysis)
                cfg = [];
                cfg.channel = 'MEG';
                cfg.keeptrials = 'yes';
                cfg.removemean = 'no';
                
                %     cfg.trials = match_str(labels_concat, condList{input(3)});
                cfg.trials = 'all'; % cuz we're looping over conditions
                
                cfg_old = cfg;
                cfg_old.keeptrials = 'no';
                data_struct = ft_timelockanalysis(cfg, data);
                data_struct_tlck = ft_timelockanalysis(cfg_old, data);
                save(fullfile(savePath, groupNames{o}, serieNames{3}, extName, sprintf('%04d', n), sprintf('%04d_%s_evoked.mat', n, strrep(condSeries{input(3)}{j},'-','_'))), 'data_struct_tlck');
                save(fullfile(savePath, groupNames{o}, serieNames{3}, extName, sprintf('%04d', n), sprintf('%04d_%s_preprocessed.mat', n, strrep(condSeries{input(3)}{j},'-','_'))), 'data_struct');
                save(fullfile(savePath, groupNames{o}, serieNames{3}, extName, sprintf('%04d', n), sprintf('%04d_%s_checks.mat', n, strrep(condSeries{input(3)}{j},'-','_'))), ...
                    'full_trials', 'trials');
                clear data_struct data_struct_tlck
            end
            clear full_trials trials
        end
        clear temp_cfg
    end
% end


%% MMNs - single trial level
% % preserving grad-information
% grad_temp = data_struct.(cond_names{input(3)}).grad;
        
% cfg = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'trial';
% % mmn_trigs = 1:length(trigs)+size(cmb_trigs,1);
% 
% % MMNs 
% % for j = length(mmn_trigs):-1:1 
% j = input(3);
% data_struct.(mmn_names{j}) = ft_math(cfg,data_struct.(cond_names{j}),data_struct.(std_names{j}));
% % data_struct_bsl.(mmn_names{j}) = ft_math(cfg,data_struct_bsl.(cond_names{j}),data_struct_bsl.(std_names{j}));
% 
% % Average (timelockanalysis) of the newly calculated MMNs at the single
% % trial level, as well as of the deviant and the standard + bsl-corr
% cfg = [];
% cfg.keeptrials = 'yes';
% cfg.removemean = 'no';
% data_struct_tlck.(cond_names{j}) = ft_timelockanalysis(cfg, data_struct.(cond_names{j}));
% data_struct_tlck.(std_names{j}) = ft_timelockanalysis(cfg, data_struct.(std_names{j}));
% data_struct_tlck.(mmn_names{j}) = ft_timelockanalysis(cfg, data_struct.(mmn_names{j}));
% 
% 
% %% Old-school MMN - avg and then combineplanar
% 
% 
% % cfg = [];
% % cfg.operation = 'subtract';
% % cfg.parameter = 'avg';
% % % mmn_trigs = 1:length(trigs)+size(cmb_trigs,1);
% % 
% % % MMNs 
% % % for j = length(mmn_trigs):-1:1 
% % j = input(3);
% % data_struct_tlck.(sprintf('%s_old',mmn_names{j})) = ft_math(cfg,data_struct_tlck.(sprintf('%s_old',cond_names{j})),data_struct_tlck.(sprintf('%s_old',std_names{j})));
% 
% % FIRST COMBINEPLANAR, THEN BASELINE CORRECTION, THEN AVERAGE...
% cfg = [];
% data_struct_tlck.(sprintf('%s_old',cond_names{j})) = ft_combineplanar(cfg, data_struct_tlck.(sprintf('%s_old',cond_names{j})));
% data_struct_tlck.(sprintf('%s_old',std_names{j})) = ft_combineplanar(cfg, data_struct_tlck.(sprintf('%s_old',std_names{j})));
% data_struct_tlck.(sprintf('%s_old',mmn_names{j})) = ft_combineplanar(cfg, data_struct_tlck.(sprintf('%s_old',mmn_names{j})));
% 
% % baseline correction
% cfg = [];
% cfg.baseline = [-0.100 0];
% cfg.parameter = {'avg'};
% data_struct_tlck.(sprintf('%s_old',cond_names{j})) = ft_timelockbaseline(cfg, data_struct_tlck.(sprintf('%s_old',cond_names{j})));
% data_struct_tlck.(sprintf('%s_old',std_names{j})) = ft_timelockbaseline(cfg, data_struct_tlck.(sprintf('%s_old',std_names{j})));
% data_struct_tlck.(sprintf('%s_old',mmn_names{j})) = ft_timelockbaseline(cfg, data_struct_tlck.(sprintf('%s_old',mmn_names{j})));
% 
% % ...AND THEN SUBTRACT?
% cfg = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'avg';
% % mmn_trigs = 1:length(trigs)+size(cmb_trigs,1);
% 
% % MMNs 
% % for j = length(mmn_trigs):-1:1 
% j = input(3);
% data_struct_tlck.(sprintf('%s_old',mmn_names{j})) = ft_math(cfg,data_struct_tlck.(sprintf('%s_old',cond_names{j})),data_struct_tlck.(sprintf('%s_old',std_names{j})));

%% DIFFS of MMNm - not relevant right now
% 
% diff_trigs = [
%     23 22; % same minus bsl
%     24 22; % opp minus bsl
%     23 24]; % same minus opp (see more below)
% 
% % diffs of MMNs
% diff_dummy = dummy_length+length(mmn_trigs);
% for k = size(diff_trigs,1):-1:1
%     data_struct.(cond_names{diff_dummy+k}) = ft_math(cfg,data_struct.(cond_names{diff_trigs(k,1)}),data_struct.(cond_names{diff_trigs(k,2)}));
%     data_struct.(cond_names{diff_dummy+k}).grad = grad_temp;
% end

% %%
% % data_struct indices:
% % 1 = aud baseline, high	tone HIGH, ball MID (no ball at all in aud block)
% % 2 = aud baseline, low  	tone LOW, ball MID (no ball at all in aud block)
% % 3 = vis only, up 		tone STD, ball UP
% % 4 = same dir. high/up		tone HIGH, ball UP
% % 5 = opp.dir.high/down		tone HIGH, ball DOWN
% % 6 = vis only, down 		tone STD, ball DOWN
% % 7 = opp.dir.low/up		tone LOW, ball UP
% % 8 = same dir. down		tone LOW, ball DOWN
% % 9 = standard
% % 10 = 1 + 2 > aud baseline, combined (both HIGH + LOW)
% % 11 = 4 + 8 > same direction,combined  
% % 12 = 5 + 7 > opposite directions, combined
% % 13 = 3 + 6 > vis only, combined
% % 14 = MMN > aud baseline, high
% % 15 = MMN > aud baseline, low
% % 16 = MMN > vis only, up 
% % 17 = MMN > same dir. high/up
% % 18 = MMN > opp.dir.high/down	
% % 19 = MMN > vis only, down
% % 20 = MMN > opp.dir.low/up
% % 21 = MMN > same dir. down
% % 22 = MMN > 1+2 aud baseline, combined
% % 23 = MMN > 4+8 > same direction,combined  
% % 24 = MMN > 5+7 > opposite directions,combined  
% % 25 = MMN > 3+6 > vis only, combined
% % 26 = diff > same minus bsl (MMNs), combined
% % 27 = diff > opp minus bsl (MMNs), combined
% % 28 = diff > same minus opp (MMNs), combined 
% 
% %% SAVE PREPROCESSED DATA to disk for further analysis
% 
% data_struct_tlck.(cond_names{j}).previous = [];
% data_struct_tlck.(std_names{j}).previous = [];
% data_struct_tlck.(mmn_names{j}).previous = [];
% 
% data_struct_tlck.(sprintf('%s_old',cond_names{j})).previous = [];
% data_struct_tlck.(sprintf('%s_old',std_names{j})).previous = [];
% data_struct_tlck.(sprintf('%s_old',mmn_names{j})).previous = [];
% 
% % checking variable size
% s = whos('data_struct_tlck');
% if (s.bytes/1024^3) >= 2
%     save(fullfile(Pdir,Sdir,subjid,sprintf('%s_%s_preprocessed.mat',condsid{input(2)},cond_names{input(3)})),'data_struct_tlck','-v7.3');
% else
%     save(fullfile(Pdir,Sdir,subjid,sprintf('%s_%s_preprocessed.mat',condsid{input(2)},cond_names{input(3)})),'data_struct_tlck');
% end

