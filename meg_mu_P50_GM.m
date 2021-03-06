function out = meg_mu_P50_GM(input)

% PROJECT: MINDLAB2011_14-music-autism
% DATESTAMP: 2020-01-28
% OWNER: Cecilie M�ller, Andreas Hoejlund, Line Gebauer
% This function identifies the channel-pair (or top 4 pairs) in each hemisphere showing max
% P50m amplitude in the std condition in the w-std paradigm (as defined as 
% the peak gradiometer-pair value betw 50-100 ms)

% input should be flagged as group (1=CTRL or 2=ASD) and 
% 1=P50m; 2=peakMMN [NB! not set up yet]; 4=peakMMN-GM [NB! not set up yet]

%% 
out = [];

% if isempty(input)
%     input = 1:2;
% end

%% Setup default function paths
addpath /projects/MINDLAB-templates/scripts/database                   % for search_database
addpath /projects/MINDLAB2011_14-music-autism/scripts/MEG   % search_db has been cloned from meeg-cfin github and put in our scripts-fodler 
% addpath /usr/local/common/matlab_toolbox/fieldtrip/latest/               % latest CFIN-version of Fieldtrip
addpath /projects/MINDLAB2011_14-music-autism/scripts/MEG/fieldtrip-20161231/

ft_defaults; 

%% Setup file paths

Pdir = '/projects/MINDLAB2011_14-music-autism/';
Savedir = '/projects/MINDLAB2011_14-music-autism/scratch/MEG/';
% NB!: Sdir further specified within some cells 
Wdir = '/projects/MINDLAB2011_14-music-autism/scratch/MEG/grandmean/ft';

conds = {'w-std', 'no-std'}; % only 'w-std' to be used at the outset
groups = {'CTRL', 'ASD'};
group = groups{input(1)};
sides = {'left', 'right'};
chan_types = {'GRADS', 'MAGS'};
layouts = {'neuromag306cmb.lay', 'neuromag306mag.lay'};
ylabels = {'RMS (fT/cm)', 'Amplitude (fT)'};

cond_name = 'std'; % only working with the std condition

ylims = [-2e-12 5e-12
    -1.5e-13 1.5e-13];

%% Load data (and combine grads)

load(fullfile(Wdir, sprintf('%s_GM_%s.mat', group, cond_name)));
cfg = [];
GM = ft_combineplanar(cfg, GM);
cfg = [];
cfg.baseline = [-0.1 0];
GM = ft_timelockbaseline(cfg, GM);

%% Overview and selection of Neuromag MEG sensors (left/right hemi)
% % selecting 48 grad-pairs in left and right hemispheres (i.e., not
% % including the 6 grad-pairs on the midline 
% cfg = [];
% cfg.layout = 'neuromag306cmb.lay'; % could also be 'neuromag306all.lay' or 'neuromag306mag.lay'
% cfg.layout = ft_prepare_layout(cfg);
% 
% %%%%% NB! Not very useful for func localizer, but indeed for GFP (where it makes sense to look across the entire hemisphere) 
% % selecting all the chan-pairs with a negative position relative to the
% % midline (= 0 in 2D space) > all left chans (and opposite for right chans,
% % i.e., >0.001); last two rows of layout.label/pos are "COMNT" and "SCALE"
% % which we're not interested in 
% hemi.(sides{1}) = cfg.layout.label(cfg.layout.pos(1:end-2,1)<0 & cfg.layout.pos(1:end-2,2)>-0.105); 
% hemi.(sides{2}) = cfg.layout.label(cfg.layout.pos(1:end-2,1)>0.001 & cfg.layout.pos(1:end-2,2)>-0.105);  
% 
% % % plot overview (not very informative)
% % figure, ft_plot_lay(cfg.layout);
% 
% % Make sure that the selected chan-pairs truly reflect all left and right
% % hemi-sensors separately
% %%%%% Try using ft_clusterplot, but it needs a bit of hacking before it'll
% %%%%% work


%% Max 40-70 ms = P50m

if ismember(1,input(2))
    
    Sdir = fullfile(Savedir, 'plots/ft/stds/P50');
    if ~exist(Sdir,'dir')
        mkdir(Sdir)
    end
    
    if ~isfield(GM,'fsample')
        fsample = 1000;
    else
        fsample = GM.fsample;
    end
    
    top = 1:4; % number of chans on each side - don't forget to update line 151 accordingly (title-line)
    
    for h = 1:length(chan_types)
        cfg = [];
        cfg.layout = layouts{h}; % could also be 'neuromag306all.lay' or 'neuromag306mag.lay'
        cfg.layout = ft_prepare_layout(cfg);
        hemi.(sides{1}) = cfg.layout.label(cfg.layout.pos(1:end-2,1)<0 & cfg.layout.pos(1:end-2,2)>-0.105); 
        hemi.(sides{2}) = cfg.layout.label(cfg.layout.pos(1:end-2,1)>0.001 & cfg.layout.pos(1:end-2,2)>-0.105);
        
        for i = 1:length(sides)
            cfg                 = [];
            cfg.latency         = [0.050 0.100];
            cfg.channel         = hemi.(sides{i});
            
            % cfg_sing               = [];
            % cfg_sing.latency       = [0.1];
            
            P50        = ft_selectdata(cfg, GM);
            for j = length(P50.label):-1:1
                if strcmp(chan_types{h}, 'MAGS') && strcmp(sides{i}, 'right')
                    % inverting the amp-values in order to detect the troughs in the right hemi
                    [yy,ii]       = findpeaks(-P50.avg(j,:)); % find peaks in each of the selected chans 
                else
                    [yy,ii]       = findpeaks(P50.avg(j,:)); % find peaks in each of the selected chans
                end
                if isempty(yy)
                    Y(j) = -Inf; I(j) = -Inf; % if no peak, leave "empty"
                else
                    [Y(j),jj]       = max(yy); % find max amp among identified peaks
                    I(j)            = ii(jj); % assing latency index of the identified max peak
                    clear jj
                end
                clear yy ii
            end
            if isempty(Y)
                if strcmp(chan_types{h}, 'MAGS') && strcmp(sides{i}, 'right')
                    [Y,I]       = max(-P50.avg,[],2); % if no peaks in any of the hemi-chans, just pick the max across all chans
                else
                    [Y,I]       = max(P50.avg,[],2); % if no peaks in any of the hemi-chans, just pick the max across all chans
                end
            end
            %         [Y,I]       = max(P50.avg,[],2);
            [X,J]       = max(Y);
            %         [XX,JJ]     = sort(Y,1,'descend');
            if strcmp(chan_types{h}, 'MAGS') && strcmp(sides{i}, 'right')
                maxP50.amp = -X;
            else
                maxP50.amp = X;
            end
            maxP50.lat = (I(J)/fsample) - 1/fsample + P50.time(1);
            maxP50.chan = P50.label{J};
            
            cfg_sing = [];
            cfg_sing.channel = maxP50.chan;
            cfg_sing.ylim = ylims(h,:);
            cfg_sing.linewidth = 3;
            figure, ft_singleplotER(cfg_sing, GM);
            hold on, plot(maxP50.lat, maxP50.amp, 'ro', 'linewidth', 2);
            text(maxP50.lat+.01, maxP50.amp+.04e-12, [num2str(maxP50.lat*1000) ' ms'])
            title(sprintf('P50m - %s (%s): %s - %s', chan_types{h}, sides{i}, maxP50.chan, group))
            
            xlabel('Time (s)')
            ylabel(ylabels{h})
            ax = gca();
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
            ax.TickDir = 'out';
            box off
            ax.XLabel.Position(2) = ylims(h,1)*0.5;
            ax.YLabel.Position = [0.01 ylims(h,2)*0.99 1];
            ax.YLabel.HorizontalAlignment = 'left';
            ax.Layer = 'top';
            
            saveas(gcf,fullfile(Sdir,sprintf('maxP50_%s_%s_%s_GM.pdf',group, chan_types{h}, sides{i})),'pdf');
            
            if strcmp(chan_types{h}, 'MAGS') && strcmp(sides{i}, 'right')
                [topY,topI]       = max(-P50.avg,[],2);
                [topXX,topJJ]     = sort(-topY,1,'ascend');
            else
                [topY,topI]       = max(P50.avg,[],2);
                [topXX,topJJ]     = sort(topY,1,'descend');
            end
            
            
            topP50.amp = topXX(top);
            topP50.lat = (topI(topJJ(top))./fsample) - 1/fsample + P50.time(1);
            topP50.chan = P50.label(topJJ(top));
            cfg_sing.channel = topP50.chan;
            figure, ft_singleplotER(cfg_sing, GM);
            hold on, plot(topP50.lat, topP50.amp, 'ro', 'linewidth', 2);
            temp_top = char(topP50.chan); temp_top = cellstr(temp_top(:,4:7));
            title(sprintf('P50m - %s (%s): %s + %s + %s + %s - %s', chan_types{h}, sides{i}, temp_top{:}, group))
            
            xlabel('Time (s)')
            ylabel('RMS (fT/cm)')
            ax = gca();
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
            ax.TickDir = 'out';
            box off
            ax.XLabel.Position(2) = ylims(h,1)*0.5;
            ax.YLabel.Position = [0.01 ylims(h,2)*0.99 1];
            ax.YLabel.HorizontalAlignment = 'left';
            ax.Layer = 'top';
            
            saveas(gcf,fullfile(Sdir,sprintf('topP50_%s_%s_%s_GM.pdf',group, chan_types{h}, sides{i})),'pdf');
            
            save(fullfile(Wdir,sprintf('P50_%s_%s_%s.mat',group, chan_types{h}, sides{i})),'maxP50','Y','I','J','P50','topP50');
            clear Y I J maxP50 P50 topP50
        end
        clear hemi
        close all
    end
end


%%
% % Selecting neighbouring channels with respect to the max P50m-chans:
% cfg.method = 'template';
% cfg.template = 'neuromag306cmb_neighb.mat';
% 
% for s = 1:length(sides)
%     load(fullfile(Wdir,sprintf('P50_%s.mat',sides{s})))
%     clear Y I J P50 topP50 gfp_top gfp_neigh gfp_hemi
%     cfg.channel = maxP50.chan;
%     neighbours = ft_prepare_neighbours(cfg, GM);
%     chans.(sides{s}) = {neighbours.label, neighbours.neighblabel{:}};
% end

%% [NB! NOT UPDATED YET!] Max 160-260 ms in diff wave (of aud) = MMNm

if ismember(2,input(2))
    
    % std = load(fullfile(Wdir, subj, 'aud_tlck.mat')); std = std.data_struct.std;
    
    if ~isfield(std,'fsample')
        fsample = 1000;
    else
        fsample = std.fsample;
    end
    
%     top = 1:5; % number of chans on each side - don't forget to update line 151 accordingly (title-line)
    
    
    for i = 1:length(sides)
        load(fullfile(Wdir,sprintf('P50_%s.mat',sides{i})),'maxP50','Y','I','J','P50','topP50');
        
        cfg                 = [];
        cfg.latency         = [0.160 0.260];
        cfg.channel         = topP50.chan;
        cfg.avgoverchan     = 'yes';
        
        % cfg_sing               = [];
        % cfg_sing.latency       = [0.1];
        
        MMN        = ft_selectdata(cfg, mmn);
        [Y,I]       = findpeaks(MMN.avg); % find peak in the average across the 4 top channels
%         for j = length(MMN.label):-1:1
%             [yy,ii]       = findpeaks(MMN.avg(j,:)); % find peaks in each of the selected chans
%             if isempty(yy)
%                 Y(j) = -Inf; I(j) = -Inf; % if no peak, leave "empty"
%             else
%                 [Y(j),jj]       = max(yy); % find max amp among identified peaks
%                 I(j)            = ii(jj); % passing latency index of the identified max peak
%                 clear jj
%             end
%             clear yy ii
%         end
        if isempty(Y)
            [Y,I]       = max(MMN.avg,[],2); % if no peaks in any of the hemi-chans, just pick the max across all chans
        end
        [X,J]       = max(Y); % find max peak among the different channels' peaks
%         [XX,JJ]     = sort(Y,1,'descend');
        maxMMN.amp = X; % assign identified max amp to the .amp-field
        maxMMN.lat = (I(J)/fsample) - 1/fsample + MMN.time(1); % assign relevant latency
        maxMMN.chan = topP50.chan; % assign relevant chan-ids from the top-selection - if going back to 1 chan only, don't forget to correct title-line accordingly (line 291)
                
        cfg_sing = [];
        cfg_sing.channel = maxMMN.chan;
        cfg_sing.ylim = ylims(1,:);
        figure, ft_singleplotER(cfg_sing, std, dev, mmn);
        hold on, plot(maxMMN.lat, maxMMN.amp, 'ro', 'linewidth', 1);
        text(maxMMN.lat+.01, maxMMN.amp+.04e-12, [num2str(maxMMN.lat*1000) ' ms'])
        title(sprintf('MMN aud - %s: %s + %s + %s + %s - GM',sides{i}, maxMMN.chan{:}))
        saveas(gcf,fullfile(Wdir,sprintf('maxMMN_%s_GM.pdf',sides{i})),'pdf');
        
        save(fullfile(Wdir,sprintf('MMN_%s.mat',sides{i})),'maxMMN','Y','I','J','MMN');
        clear Y I J MMN maxMMN
    end
end


%% NB! Prob not relevant anymore (max 200-300 ms in diff wave (of aud) = P3a)

% if ismember(3,input)
%     
%     % std = load(fullfile(Wdir, subj, 'aud_tlck.mat')); std = std.data_struct.std;
%     
%     if ~isfield(std,'fsample')
%         fsample = 1000;
%     else
%         fsample = std.fsample;
%     end
%     
% %     top = 1:5; % number of chans on each side - don't forget to update line 151 accordingly (title-line)
%     
%     
%     for i = 1:length(sides)
%         cfg                 = [];
%         cfg.latency         = [0.200 0.300];
%         cfg.channel         = chans.(sides{i});
%         
%         % cfg_sing               = [];
%         % cfg_sing.latency       = [0.1];
%         
%         P3a        = ft_selectdata(cfg, mmn);
%         for j = length(P3a.label):-1:1
%             [yy,ii]       = findpeaks(P3a.avg(j,:)); % find peaks in each of the selected chans
%             if isempty(yy)
%                 Y(j) = []; I(j) = []; % if no peak, leave empty
%             else
%                 [Y(j),jj]       = max(yy); % find max amp among identified peaks
%                 I(j)            = ii(jj); % assing latency index of the identified max peak
%                 clear jj
%             end
%             clear yy ii
%         end
%         if isempty(Y)
%             [Y(j),I(j)]       = max(P3a.avg,[],2); % if no peaks in any of the hemi-chans, just pick the max across all chans
%         end
%         [X,J]       = max(Y);
% %         [XX,JJ]     = sort(Y,1,'descend');
%         maxP3a.amp = X;
%         maxP3a.lat = (I(J)/fsample) - 1/fsample + P3a.time(1);
%         maxP3a.chan = P3a.label{J};
%         
%         cfg_sing = [];
%         cfg_sing.channel = maxP3a.chan;
%         cfg_sing.ylim = ylims(1,:);
%         figure, ft_singleplotER(cfg_sing, std, dev, mmn);
%         hold on, plot(maxP3a.lat, maxP3a.amp, 'ro', 'linewidth', 1);
%         text(maxP3a.lat+.01, maxP3a.amp+.04e-12, [num2str(maxP3a.lat*1000) ' ms'])
%         title(sprintf('P3a aud - %s: %s - GM',sides{i}, maxP3a.chan))
%         saveas(gcf,fullfile(Wdir,sprintf('maxP3a_%s_GM.pdf',sides{i})),'pdf');
%         
%         save(fullfile(Wdir,sprintf('P3a_%s.mat',sides{i})),'maxP3a','Y','I','J','P3a');
%     end
% end

%% Extract peak latency for the MMN in the greatGM (bslCmb + sameCmb + oppCmb)

if ismember(4,input(2))
    load(fullfile(Wdir, 'greatGMall_AV_all_devCmb_tlck.mat'));
    std = data.std;
    dev = data.dev;
    mmn = data.mmn;
    clear data
    
    if ~isfield(std,'fsample')
        fsample = 1000;
    else
        fsample = std.fsample;
    end
    
    top = 1:5; % number of chans on each side - don't forget to update line 151 accordingly (title-line)
    
    
    for i = 1:length(sides)
        load(fullfile(Wdir,sprintf('MMN_%s.mat',sides{i})))
        
        cfg                 = [];
        cfg.latency         = [0.160 0.260];
        cfg.channel         = maxMMN.chan;
        cfg.avgoverchan     = 'yes';

        MMN_cmb        = ft_selectdata(cfg, mmn);
        [Y,I]       = findpeaks(MMN_cmb.avg);
        if isempty(Y)
            [Y,I]       = max(MMN_cmb.avg,[],2);
        end
        [X,J]       = max(Y);
        greatMMN.amp = X;
        greatMMN.lat = (I(J)/fsample) - 1/fsample + MMN_cmb.time(1);
        greatMMN.chan = maxMMN.chan;

        cfg_sing = [];
        cfg_sing.channel = maxMMN.chan;
        cfg_sing.ylim = ylims(1,:);
        figure, ft_singleplotER(cfg_sing, std, dev, mmn);
        hold on, plot(greatMMN.lat, greatMMN.amp, 'ro', 'linewidth', 1);
        text(greatMMN.lat+.01, greatMMN.amp+.04e-12, [num2str(greatMMN.lat*1000) ' ms'])
        title(sprintf('MMN cmb - %s: %s + %s + %s + %s - GM',sides{i}, greatMMN.chan{:}))    % max chan is not based on the top4-selection - if going back to 1 chan only, don't forget to correct this line accordingly  
        legend({'std','bsl-same-oppCmb','mmnCmb'})
        saveas(gcf,fullfile(Wdir,sprintf('greatMMN_%s_GM.pdf',sides{i})),'pdf');
        
        save(fullfile(Wdir,sprintf('greatMMN_%s.mat',sides{i})),'greatMMN','Y','I','J','MMN_cmb');
    end
    
%     for i = 1:length(sides)
% %         load(fullfile(Wdir,sprintf('P3a_%s.mat',sides{i}))); maxChan = maxP3a.chan
%         load(fullfile(Wdir,sprintf('MMN_%s.mat',sides{i}))); maxChan = maxMMN.chan;
%         
%         cfg                 = [];
%         cfg.latency         = [0.200 0.300];
%         cfg.channel         = maxChan;        
%         P3a_cmb        = ft_selectdata(cfg, mmn);
%         [Y,I]       = findpeaks(P3a_cmb.avg);
%         [X,J]       = max(Y);
%         greatP3a.amp = X;
%         greatP3a.lat = (I(J)/fsample) - 1/fsample + P3a.time(1);
%         greatP3a.chan = maxChan;
% 
%         cfg_sing = [];
%         cfg_sing.channel = greatP3a.chan;
%         cfg_sing.ylim = ylims(1,:);
%         figure, ft_singleplotER(cfg_sing, std, dev, mmn);
%         hold on, plot(greatP3a.lat, greatP3a.amp, 'ro', 'linewidth', 1);
%         text(greatP3a.lat+.01, greatP3a.amp+.04e-12, [num2str(greatP3a.lat*1000) ' ms'])
%         title(sprintf('P3a cmb - %s: %s - GM',sides{i}, greatP3a.chan))      
%         legend({'std','bsl-same-oppCmb','mmnCmb'})
%         saveas(gcf,fullfile(Wdir,sprintf('greatP3a_%s_GM.pdf',sides{i})),'pdf');
%         
%         save(fullfile(Wdir,sprintf('greatP3a_%s.mat',sides{i})),'greatP3a','Y','I','J','P3a_cmb');
%     end
end

