function out = custom_mu_plotting(input)

out = [];
%% NB! CUSTOM plotting (see http://eelkespaak.nl/blog/customizing-common-m-eeg-plots-part-1-the-event-related-potential-field-erp-f/ for more details)

% STD-DEV-DEV-...-DEV plots for ALL deviants
% OR STD-STD-...-STD plots for ALL standards (in MuMufe)
% we use these plots to plot all the conditions/standards against each other

addpath /projects/MINDLAB2011_14-music-autism/scripts/MEG/fieldtrip-20161231/
ft_defaults;

% Savedir = '/Volumes/projects/MINDLAB2011_14-music-autism/scratch/MEG/plots/ft';
Savedir = '/projects/MINDLAB2011_14-music-autism/scratch/MEG/plots/ft';
Wdir = '/projects/MINDLAB2011_14-music-autism/scratch/MEG/grandmean/ft';
Rootdir = '/projects/MINDLAB2011_14-music-autism/scratch/MEG/';
addpath(genpath('~/matlab/boundedline'))
% datas = {'std', 'dev', 'mmn'};
femto = 1e-15;
conds = {'w-std', 'no-std'};
groups = {'CTRL', 'ASD'};
sides = {'left', 'right'};
% chans = {{'MEG2241','MEG2431'}, {'MEG0421','MEG0431'};
%     {'MEG0242+0243','MEG0232+0233'}, {'MEG1332+1333','MEG1612+1613'}};
chan_types = {'MAGS', 'GRADS'};
chan_type = input(3);
g = input(2);
p = input(1);

grouplist{1} = [1 3 5 6 7 8 11 12 15 41 44 48 49 50 52 55 62 65]; % All CTRL participants, identified by their exam nr (subj0059 is excluded from MuMUFE due to technical problems, see the database (wiki) for the full report of exclusions)
grouplist{2} = [19 21 24 27 30 37 38 39 43 45 46 53 54 56 57 60 61 63 64];                             % All ASD participants, identified by their exam nr (subj0058 is excluded due to medical history, see the database (wiki) for the full report of exclusions)

% NB! variance plotting not fully fixed yet
err = 0; % if 1 = variance on diff wave; if 2 = Cousineau/Morey variance on std and dev; 0 = no variance

% for plotting DEVs
devNames = {{
    'std', ...
    'pitch_m', ...
    'timbre_m', ...
    'loc_m', ...
    'intens_m', ...
    'slide_m', ...
    'duration_m', ...
    'rhythm_m'}
    
    {'std2_ns', ...
    'pitch_ns', ...
    'timbre_ns', ...
    'loc_ns', ...
    'intens_ns', ...
    'slide_ns', ...
    'duration_ns', ...
    'rhythm_ns'}};

% for plotting STDs
stdNames = {{
    'std', ...
    'std1', ...
    'std2', ...
    'std4', ...
    'std124', ...
    'std1_all', ...
    'std_key', ...
    'std_meter', ...
    'std_key_meter'}, ...

    {'std', ...
    'std1_ns', ...
    'std2_ns', ...
    'std4_ns' ...
    'std124_ns', ...
    'std1_all_ns', ...
    'std_key_ns', ...
    'std_meter_ns', ...
    'std_key_meter_ns'}};

ylims = [-80 80;  % mags
    -2000 4000];     % grads    

%% GRANDMEAN stds + std+devs

if ismember(1,input(1))
    for m = 1:length(sides)
        temp = load(fullfile(Wdir,sprintf('P50_%s_%s_%s.mat',groups{g}, chan_types{chan_type}, sides{m})));
        topP50.(sides{m}) = temp.topP50;
        maxP50.(sides{m}) = temp.maxP50;
    end
    for i = 1:length(conds) % looping over w-std and no-std
        if g == 2
            Sdir = fullfile(Savedir,'stds','ASD');
        elseif g == 1
            Sdir = fullfile(Savedir,'stds');  
        end
        % Sdir = Savedir;
        if ~exist(Sdir,'dir')
            mkdir(Sdir)
        end
        
        for j = 1:size(stdNames{i},2)    % based on number of rows in the relevant stdNames-cell we'll do X number of std traces
            load(fullfile(Wdir, sprintf('%s_GM_%s.mat', groups{g}, stdNames{i}{j})))
            cfg = [];
            GM = ft_combineplanar(cfg, GM);
            cfg = [];
            cfg.baseline = [-0.1 0];
            tlck_data.(stdNames{i}{j}) = ft_timelockbaseline(cfg, GM);
            clear GM
        end
        datatime = tlck_data.(stdNames{i}{1}).time;
        for m = 1:length(sides) % make a plot for the "peak" chans in each hemi
%             load(fullfile(Wdir,sprintf('P50_%s_%s_%s.mat',groups{g}, chan_types{chan_type}, sides{m})))
            % list of top 4 chans
            chaninds = match_str(tlck_data.(stdNames{i}{1}).label, topP50.(sides{m}).chan);
            %                 chaninds = match_str(tlck_data.(stdNames{i}{1}).label, chans{chan_type, m});
            % average over top 4 chans + converting to femtoTesla
            for j = 1:size(stdNames{i},2)
                erf.(stdNames{i}{j}).(sides{m}) = ...
                    mean(tlck_data.(stdNames{i}{j}).avg(chaninds, :), 1) ./ femto;
            end
            %             topP50.(sides{m}).amp_femto = topP50.(sides{m}).amp ./ femto;
            
            % see info regarding suitable colors here:
            % http://jfly.iam.u-tokyo.ac.jp/color/
            linecolors = [0 0 0 1;    % black
                230/255 159/255 0 1;   % orange
                86/255 180/255 233/255 1;      % sky blue
                0 158/255 115/255 1;    % bluish green
                213/255 94/255 0 1; % red (vermillion)
                240/255 228/255 66/255 1; % yellow
                0/255 114/255 178/255 1; % blue
                204/255 121/255 167/255 1; % reddish purple
                0 0 0 0.5; %grey
                0 0 0 0.05]; % black transparaent (aka. grey)
            
            figure(),
            hold on,
            
            %                     for k = 1:2     % only plot individual participants for the std and dev
            %                         plot(datatime, erf_struct.(datas{k}).(sides{m})', 'Color', linecolors(k, :))
            %                     end
            
            %             % draw the shaded area (outlining the GM-timewin for the MMN
            %             ww = fill([timewin_MMN(1)/1000-0.1 timewin_MMN(1)/1000-0.1 ...
            %                 timewin_MMN(end)/1000-0.1 timewin_MMN(end)/1000-0.1], ...
            %                 [ylims(1,1) ylims(1,2) ylims(1,2) ylims(1,1)], [.85 .85 .85]);
            %             ww.LineStyle = 'none';
            
            if err == 0
                for j = 1:size(stdNames{i},2)-3  % (excl. std-key etc.)   % size(stdNames{i},2)
                    nn(j) = plot(datatime, erf.(stdNames{i}{j}).(sides{m}), ...
                        'Color', linecolors(j,:), 'linewidth', 3); % plotting all stds
                end
            elseif err == 1
                % only plot individual participants for the diff wave
                plot(datatime, erf_struct.(sides{m})', 'Color', linecolors(3, :))
                % plot bounds for only mmn
                [ll, ~] = boundedline(datatime, erf.(sides{m}), bounds.(sides{m}), 'k--', 'alpha', 'transparency', 0.1);
                ll.LineWidth = 2;
                mm = plot(datatime, erf.(sides{m}), 'b', 'linewidth', 2);
                nn = plot(datatime, erf.(sides{m}), 'r', 'linewidth', 2);
                
            elseif err == 2
                % plot bounds for only std and dev
                [ll, ~] = boundedline(datatime, [erf.(sides{m}); erf.(sides{m})], ...
                    cat(3, bounds.(sides{m}), bounds.(sides{m})), 'alpha', 'transparency', 0.1);
                ll(1).LineWidth = 2;
                ll(2).LineWidth = 2;
                % plot the mmn in plain black and broken line
                mm = plot(datatime, erf.(sides{m}), 'k--', 'linewidth', 2);
            end
            
            %             if any(abs(diff(topP50.(sides{m}).amp_femto))<100)
            %                 topP50.(sides{m}).amp_femto(abs(diff(topP50.(sides{m}).amp_femto))<100) = topP50.(sides{m}).amp_femto(abs(diff(topP50.(sides{m}).amp_femto))<100) + 50;
            %             end
            %             oo = plot(topP50.(sides{m}).lat, topP50.(sides{m}).amp_femto, 'ko', 'linewidth', 1);
            
            ylim(ylims(chan_type,:))
            %                 xlim([-0.396 0.423])
            xlim([-0.086 0.423])  % adjusted by 13.6 ms
            
            xlabel('Time (s)')
            ylabel('RMS (fT/cm)')
            ax = gca();
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
            ax.TickDir = 'out';
            box off
            ax.XLabel.Position(2) = -ylims(chan_type,1)*0.25;
            ax.YLabel.Position = [-0.1 ylims(chan_type,2)*0.99 1];
            ax.YLabel.HorizontalAlignment = 'center';
            ax.Layer = 'top';
            
            title(sprintf('All standards  |  %s  |  %s  |  %s', conds{i}, sides{m}, chan_types{chan_type})); % custom title
            
            if err == 0
                legend(nn, strrep(stdNames{i}(1:length(nn)), '_', '-'), 'Location', 'NorthEast') % custom legend 'SouthOutside'; 'NorthEast'
                clear nn
            elseif err == 1
                legend([mm, nn, ll], 'std','dev','diff') % custom legend
                clear mm nn ll
            elseif err == 2
                legend([ll', mm], 'std','dev','diff') % custom legend
                clear ll mm
            end
            legend boxoff
            set(gca,'FontSize',12);
            %                     hline(0,'k');
            saveas(gcf, fullfile(Sdir, sprintf('CUSTOM_allstds_%s_GM_%s_%s_%s.pdf', groups{g}, chan_types{chan_type}, sides{m}, conds{i})),'pdf'); % save with custom title
            
            % individual plots for each standard
            for j = 1:size(stdNames{i},2)
                figure,
                mm = plot(datatime, erf.(stdNames{i}{j}).(sides{m}), ...
                    'Color', linecolors(j+(i-1),:), 'linewidth', 3); % plotting all stds
                ylim(ylims(chan_type,:))
                %                     xlim([-0.396 0.423])
                xlim([-0.086 0.423])  % adjusted by 13.6 ms
                
                xlabel('Time (s)')
                ylabel('RMS (fT/cm)')
                ax = gca();
                ax.XAxisLocation = 'origin';
                ax.YAxisLocation = 'origin';
                ax.TickDir = 'out';
                box off
                ax.XLabel.Position(2) = -500;
                ax.YLabel.Position = [-0.1 ylims(chan_type,2)*0.99 1];
                ax.YLabel.HorizontalAlignment = 'center';
                ax.Layer = 'top';
                
                title(sprintf('%s  |  %s  |  %s  |  %s', strrep(stdNames{i}{j}, '_', '-'), conds{i}, sides{m}, chan_types{chan_type})); % custom title
                legend(mm, strrep(stdNames{i}{j}, '_', '-'), 'Location', 'NorthEast') % custom legend 'SouthOutside'; 'NorthEast'
                legend boxoff
                set(gca,'FontSize',12);
                saveas(gcf, fullfile(Sdir, sprintf('CUSTOM_%s_GM_%s_%s_%s.pdf', groups{g}, chan_types{chan_type}, sides{m}, stdNames{i}{j})),'pdf'); % save with custom title
            end
        end
        clear tlck_data
        close all
    end
end

%% INDIVS stds + std+devs

ylims = [-80 80;  % mags
    -2000 6000];     % grads

if ismember(2,input(1))
    for m = 1:length(sides)
        temp = load(fullfile(Wdir,sprintf('P50_%s_%s_%s.mat',groups{g}, chan_types{chan_type}, sides{m})));
        topP50.(sides{m}) = temp.topP50;
        maxP50.(sides{m}) = temp.maxP50;
    end
    for s = 1:length(grouplist{g})
        for i = 1:length(conds) % looping over w-std and no-std
            if input(2) == 2
                Sdir = fullfile(Savedir,'stds','indivs','ASD');
            elseif input(2) == 1
                Sdir = fullfile(Savedir,'stds','indivs');
            end
            % Sdir = Savedir;
            if ~exist(Sdir,'dir')
                mkdir(Sdir)
            end
            Filedir = fullfile(Rootdir, groups{g}, 'ns/ft', sprintf('%04d', grouplist{g}(s)));
            
            for j = 1:size(stdNames{i},2)    % based on number of rows in the relevant stdNames-cell we'll do X number of std traces
                load(fullfile(Filedir, sprintf('%04d_%s_evoked.mat', grouplist{g}(s), stdNames{i}{j})))
                cfg = [];
                data_struct_tlck = ft_combineplanar(cfg, data_struct_tlck);
                cfg = [];
                cfg.baseline = [-0.1 0];
                tlck_data.(stdNames{i}{j}) = ft_timelockbaseline(cfg, data_struct_tlck);
                clear data_struct_tlck
            end
            datatime = tlck_data.(stdNames{i}{1}).time;
            for m = 1:length(sides) % make a plot for the "peak" chans in each hemi
%                 load(fullfile(Wdir,sprintf('P50_%s_%s_%s.mat',groups{g}, chan_types{chan_type}, sides{m})))
                % list of top 4 chans
                chaninds = match_str(tlck_data.(stdNames{i}{1}).label, topP50.(sides{m}).chan);
                %                 chaninds = match_str(tlck_data.(stdNames{i}{1}).label, chans{chan_type, m});
                % average over top 4 chans + converting to femtoTesla
                for j = 1:size(stdNames{i},2)
                    erf.(stdNames{i}{j}).(sides{m}) = ...
                        mean(tlck_data.(stdNames{i}{j}).avg(chaninds, :), 1) ./ femto;
                end
                %             topP50.(sides{m}).amp_femto = topP50.(sides{m}).amp ./ femto;
                
                % see info regarding suitable colors here:
                % http://jfly.iam.u-tokyo.ac.jp/color/
                linecolors = [0 0 0 1;    % black
                    230/255 159/255 0 1;   % orange
                    86/255 180/255 233/255 1;      % sky blue
                    0 158/255 115/255 1;    % bluish green
                    213/255 94/255 0 1; % red (vermillion)
                    240/255 228/255 66/255 1; % yellow
                    0/255 114/255 178/255 1; % blue
                    204/255 121/255 167/255 1; % reddish purple
                    0 0 0 0.5; %grey
                    0 0 0 0.05]; % black transparaent (aka. grey)
                
                figure(),
                hold on,
                
                %                     for k = 1:2     % only plot individual participants for the std and dev
                %                         plot(datatime, erf_struct.(datas{k}).(sides{m})', 'Color', linecolors(k, :))
                %                     end
                
                %             % draw the shaded area (outlining the GM-timewin for the MMN
                %             ww = fill([timewin_MMN(1)/1000-0.1 timewin_MMN(1)/1000-0.1 ...
                %                 timewin_MMN(end)/1000-0.1 timewin_MMN(end)/1000-0.1], ...
                %                 [ylims(1,1) ylims(1,2) ylims(1,2) ylims(1,1)], [.85 .85 .85]);
                %             ww.LineStyle = 'none';
                
                if err == 0
                    for j = 1:size(stdNames{i},2)-3  % (excl. std-key etc.)   % size(stdNames{i},2)
                        nn(j) = plot(datatime, erf.(stdNames{i}{j}).(sides{m}), ...
                            'Color', linecolors(j,:), 'linewidth', 3); % plotting all stds
                    end
                elseif err == 1
                    % only plot individual participants for the diff wave
                    plot(datatime, erf_struct.(sides{m})', 'Color', linecolors(3, :))
                    % plot bounds for only mmn
                    [ll, ~] = boundedline(datatime, erf.(sides{m}), bounds.(sides{m}), 'k--', 'alpha', 'transparency', 0.1);
                    ll.LineWidth = 2;
                    mm = plot(datatime, erf.(sides{m}), 'b', 'linewidth', 2);
                    nn = plot(datatime, erf.(sides{m}), 'r', 'linewidth', 2);
                    
                elseif err == 2
                    % plot bounds for only std and dev
                    [ll, ~] = boundedline(datatime, [erf.(sides{m}); erf.(sides{m})], ...
                        cat(3, bounds.(sides{m}), bounds.(sides{m})), 'alpha', 'transparency', 0.1);
                    ll(1).LineWidth = 2;
                    ll(2).LineWidth = 2;
                    % plot the mmn in plain black and broken line
                    mm = plot(datatime, erf.(sides{m}), 'k--', 'linewidth', 2);
                end
                
                %             if any(abs(diff(topP50.(sides{m}).amp_femto))<100)
                %                 topP50.(sides{m}).amp_femto(abs(diff(topP50.(sides{m}).amp_femto))<100) = topP50.(sides{m}).amp_femto(abs(diff(topP50.(sides{m}).amp_femto))<100) + 50;
                %             end
                %             oo = plot(topP50.(sides{m}).lat, topP50.(sides{m}).amp_femto, 'ko', 'linewidth', 1);
                
                ylim(ylims(chan_type,:))
                %                 xlim([-0.396 0.423])
                xlim([-0.086 0.423])  % adjusted by 13.6 ms
                
                xlabel('Time (s)')
                ylabel('RMS (fT/cm)')
                ax = gca();
                ax.XAxisLocation = 'origin';
                ax.YAxisLocation = 'origin';
                ax.TickDir = 'out';
                box off
                ax.XLabel.Position(2) = -ylims(chan_type,1)*0.25;
                ax.YLabel.Position = [-0.1 ylims(chan_type,2)*0.99 1];
                ax.YLabel.HorizontalAlignment = 'center';
                ax.Layer = 'top';
                
                title(sprintf('All standards (%04d)  |  %s  |  %s  |  %s', grouplist{g}(s), conds{i}, sides{m}, chan_types{chan_type})); % custom title
                
                if err == 0
                    legend(nn, strrep(stdNames{i}(1:length(nn)), '_', '-'), 'Location', 'NorthEast') % custom legend 'SouthOutside'; 'NorthEast'
                    clear nn
                elseif err == 1
                    legend([mm, nn, ll], 'std','dev','diff') % custom legend
                    clear mm nn ll
                elseif err == 2
                    legend([ll', mm], 'std','dev','diff') % custom legend
                    clear ll mm
                end
                legend boxoff
                set(gca,'FontSize',12);
                %                     hline(0,'k');
                saveas(gcf, fullfile(Sdir, sprintf('CUSTOM_allstds_%04d_%s_%s_%s.pdf', grouplist{g}(s), chan_types{chan_type}, sides{m}, conds{i})),'pdf'); % save with custom title
                
%                 % individual plots for each standard - not for indivs
%                 for j = 1:size(stdNames{i},2)
%                     figure,
%                     mm = plot(datatime, erf.(stdNames{i}{j}).(sides{m}), ...
%                         'Color', linecolors(j+(i-1),:), 'linewidth', 3); % plotting all stds
%                     ylim(ylims(chan_type,:))
%                     %                     xlim([-0.396 0.423])
%                     xlim([-0.086 0.423])  % adjusted by 13.6 ms
%                     
%                     xlabel('Time (s)')
%                     ylabel('RMS (fT/cm)')
%                     ax = gca();
%                     ax.XAxisLocation = 'origin';
%                     ax.YAxisLocation = 'origin';
%                     ax.TickDir = 'out';
%                     box off
%                     ax.XLabel.Position(2) = -500;
%                     ax.YLabel.Position = [-0.1 ylims(chan_type,2)*0.99 1];
%                     ax.YLabel.HorizontalAlignment = 'center';
%                     ax.Layer = 'top';
%                     
%                     title(sprintf('%s  |  %s  |  %s  |  %s', strrep(stdNames{i}{j}, '_', '-'), conds{i}, sides{m}, chan_types{chan_type})); % custom title
%                     legend(mm, strrep(stdNames{i}{j}, '_', '-'), 'Location', 'NorthEast') % custom legend 'SouthOutside'; 'NorthEast'
%                     legend boxoff
%                     set(gca,'FontSize',12);
%                     saveas(gcf, fullfile(Sdir, sprintf('CUSTOM_%s_GM_%s_%s_%s.pdf', groups{g}, chan_types{chan_type}, sides{m}, stdNames{i}{j})),'pdf'); % save with custom title
%                 end
            end
        end
        clear tlck_data
        close all
    end
end