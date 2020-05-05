%% NB! CUSTOM plotting (see http://eelkespaak.nl/blog/customizing-common-m-eeg-plots-part-1-the-event-related-potential-field-erp-f/ for more details)

addpath /projects/MINDLAB2011_14-music-autism/scripts/MEG/fieldtrip-20161231/
ft_defaults

% requires that boundedline has been downloaded from https://se.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m
% specify path to boundedline
addpath(genpath('~/matlab/boundedline')) 

% specify name of the fields containing the standard, deviant and
% difference ERPs
datas = {'std', 'dev', 'mmn'};

err = 1; % if 1 > variance on diff wave; if 2 > Cousineau/Morey variance on std and dev; if 3 > C/M correction + betw variance [NB! Not functional yet!]

% for plotting MMNs (incl. std and dev)
% this is a very specific customization for the project which this script 
% was originally intended - it is made up of 3 cells referring to 3 overall
% conditions which were of different "sizes", hence 2 cells with only 1
% condition (i.e. 1 std, 1 dev and 1 MMN) and 1 cell with 4 conditions
% (i.e. 4 std, 4 dev, 4 MMN)
stdNames = {
    {'bslCmb','stdbslCmb','MMN_bsl'}, ...
    {'visCmb','stdvisCmb','MMN_vis'}, ...
    {'bslCmb','stdbslCmb','MMN_bsl';
    'sameCmb','stdsameCmb','MMN_same';
    'oppCmb','stdoppCmb','MMN_opp';
    'visCmb','stdvisCmb','MMN_vis'}}; % cell 1 = aud, cell 2 = vis, and cell 3 with several rows = 4 MMNs in AV

% setting different ylimits for different conditions and settings
% NB! Needs to be changed - is currently adapted to MEG (should probably be
% something [-4 4] when plotting EEG data
if err == 1
    ylims = [-2150 2150     % aud cond
    -2150 2150          % vis cond
    -2150 2150];        % AV cond

elseif err == 2
    ylims = [-1000 5000     % aud cond
    -1000 5000         % vis cond
    -1000 5000];        % AV cond
end


for i = 1:length(conds) % looping over aud, vis, and AV
    % specifying the first part of where to save the plots
    Sdir = fullfile(Savedir,'std-dev-diff',conds{i});
    % creating the folder if doesn't already exist
    if ~exist(Sdir,'dir')
        mkdir(Sdir)
    end

    % specifying exactly where to save the plots
    Sdir_single = fullfile(Sdir,'single','groups');
    % creating the folder if doesn't already exist
    if ~exist(Sdir_single,'dir')
        mkdir(Sdir_single)
    end

    for j = 1:size(stdNames{i},1)    % based on number of rows in the relevant mmn_names-cell we'll do 1 or 4 plots 
        for g = 1 %:length(groups)  % possibility of specifying more than one group and looping over them
            % load grandmean-file
            load(fullfile(Wdir, sprintf('GM%s_%s_%s_tlck.mat', groups{g}, conds{i}, stdNames{i}{j,1})))
            % extracting the time vector
            % NB! this has to be updated with respect to the EEGLAB
            % dataformat
            datatime = data.(datas{1}).time;
            % loading data structure with individual ERPs for std, dev, MMN
            load(fullfile(Wdir, sprintf('GM_struct_%s_%s_tlck.mat', conds{i}, stdNames{i}{j,2})))
            load(fullfile(Wdir, sprintf('GM_struct_%s_%s_tlck.mat', conds{i}, stdNames{i}{j,1})))
            load(fullfile(Wdir, sprintf('GM_struct_%s_%s_tlck.mat', conds{i}, stdNames{i}{j,3})))
            % loading saved summary values (or "extracted values") defining
            % peak channels
            % NB! probably not relevant for your project, Quan
            load(fullfile(Pdir,'scratch/grandaverage/THIRD',sprintf('P50_%s.mat',sides{m})))
            % specifying indices of peak channels - could just be the index of 'Fz'
            % NB! Needs to be adapted to EEGLAB
            chaninds = match_str(data.(datas{1}).label, topP50.chan);
            % creating an average of the peak channels
            %%%%% NB! if only 1 peak channel specified (e.g. Fz), then this
            %%%%% just exctracts the data into an erf_struct-structure
            % NB! Should be adapted to EEGLAB
            for p = 1:length(GM_struct_std)
                erf_struct.(datas{1})(p,:) = mean(GM_struct_std{p}.avg(chaninds, :), 1) ./ femto;
                erf_struct.(datas{2})(p,:) = mean(GM_struct_dev{p}.avg(chaninds, :), 1) ./ femto;
                erf_struct.(datas{3})(p,:) = mean(GM_struct_mmn{p}.avg(chaninds, :), 1) ./ femto;
            end
            
            %%% COUSINEAU/MOREY CORRRECTION - NB! For plotting purposes only! %%%
            erf_CM.data(:,:,1) = erf_struct.(datas{1}); % std
            erf_CM.data(:,:,2) = erf_struct.(datas{2}); % dev
            erf_CM.new_data = zeros(size(erf_CM.data)); % zeros for new std + dev
            gm = mean(mean(erf_CM.data,3));         % mean across participants (dim=1) and conditions (dim=3)
            for k = size(erf_CM.data,1):-1:1
                pm(k,:) = mean(erf_CM.data(k,:,:),3);   % mean for each participant across conditions
                for c = 1:size(erf_CM.data,3)
                    erf_CM.new_data(k,:,c) = erf_CM.data(k,:,c) - pm(k,:) + gm;      % adjusting by the difference between pm and gm
                end
            end
            for c = 1:size(erf_CM.data,3)
                erf_CM.(datas{c}).mean = mean(erf_CM.new_data(:,:,c));
                erf_CM.(datas{c}).var = var(erf_CM.new_data(:,:,c));
                % multiplying the sample variance by ncond/(ncond-1) - see Morey (2008), eq. 2.
                erf_CM.(datas{c}).var_corr = erf_CM.(datas{c}).var.*(size(erf_CM.new_data,3)/(size(erf_CM.new_data,3)-1));
                erf_CM.(datas{c}).sd = sqrt(erf_CM.(datas{c}).var_corr);
                erf_CM.(datas{c}).sem = erf_CM.(datas{c}).sd/size(erf_CM.new_data,1);
            end
            
            bounds.(datas{1}) = std(erf_struct.(datas{1}), [], 1);
            bounds.(datas{2}) = std(erf_struct.(datas{2}), [], 1);
            bounds.(datas{3}) = std(erf_struct.(datas{3}), [], 1);
            
            % clearing the loaded data structures for the individual data
            clear GM_struct_std GM_struct_dev GM_struct_mmn 

            load(fullfile(Pdir,'scratch/grandaverage/THIRD/P50.mat'))
            % list of top 4 chans
            chaninds = match_str(data.(datas{1}).label, topP50.chan);
            % see info regarding suitable colors here:
            % http://jfly.iam.u-tokyo.ac.jp/color/
            linecolors = [230/255 159/255 0 1;   % orange
                86/255 180/255 233/255 1;      % sky blue
                0 158/255 115/255 1;    % bluish green
                213/255 94/255 0 1; % red (vermillion)
                0 0 0 0.05]; % black transparaent (aka. grey)
            %                     linecolors = [1 0 0 0.05;   % red transparent
            %                         0 0 1 0.05;      % blue transparent
            %                         0 0 0 0.05;];    % black transparent (aka. grey)
            
            % average over top 4 chans
            for n = 1:length(datas)
                erf.(datas{n}) = mean(data.(datas{n}).avg(chaninds, :), 1);
            end
            topP50.amp_femto = topP50.amp;
            
            figure(),
            hold on,
            
            % draw the shaded area (outlining the GM-timewin for the MMN
            rectangle('Position', [timewin_MMN(1)/1000-0.1 ylims(i,1) length(timewin_MMN)/1000 ylims(i,2)-ylims(i,1)], 'EdgeColor', [.85 .85 .85], 'FaceColor', [.85 .85 .85])
            
            if err == 1
                % only plot individual participants for the diff wave
                plot(datatime, erf_struct.(datas{3})', 'Color', linecolors(3, :))
                % plot bounds for only mmn
                [ll, ~] = boundedline(datatime, erf.(datas{3}), bounds.(datas{3}), 'k--', 'alpha', 'transparency', 0.1);
                ll.LineWidth = 2;
                mm = plot(datatime, erf.(datas{1}), 'b', 'linewidth', 2);
                nn = plot(datatime, erf.(datas{2}), 'r', 'linewidth', 2);
                set(gca, 'Color', 'None')
            elseif err == 2
                % plot bounds for only std and dev
                [ll, ~] = boundedline(datatime, [erf.(datas{1}); erf.(datas{2})], ...
                    cat(3, bounds.(datas{1}), bounds.(datas{2})), 'alpha', 'transparency', 0.1);
                %                             cat(3, erf_CM.(datas{1}).sd, erf_CM.(datas{2}).sd), 'alpha', 'transparency', 0.1);
                ll(1).LineWidth = 2;
                ll(2).LineWidth = 2;
                % plot the mmn in plain black and broken line
                mm = plot(datatime, erf.(datas{3}), 'k--', 'linewidth', 2);
                set(gca, 'Color', 'None')
                
%             %%%% NB! not fully functional - too many lines,
%             %%%% legend not working, etc.
%             elseif err == 3
%                 % plot bounds for only std and dev - but for both
%                 % between and withing
%                 [ll, ~] = boundedline(datatime, repmat([erf.(datas{1}); erf.(datas{2})],2,1), ...
%                     cat(3, erf_CM.(datas{1}).sd, erf_CM.(datas{2}).sd, ...
%                     bounds.(datas{1}), bounds.(datas{2})), 'alpha', 'transparency', 0.1);
%                 ll(1).LineWidth = 2;
%                 ll(2).LineWidth = 2;
%                 % plot the mmn in plain black and broken line
%                 mm = plot(datatime, erf.(datas{3}), 'k--', 'linewidth', 2);
%                 set(gca, 'Color', 'None')
            end
                        
            ylim(ylims(i,:))
            xlim([-0.1 0.416])
            
            xlabel('Time (s)')
            ylabel('RMS (fT/cm)')
            ax = gca();
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
            ax.TickDir = 'out';
            box off
            ax.XLabel.Position(2) = -1000;
            ax.Layer = 'top';
            
            title(sprintf('%s - %s (single): ERFs %s', groups{g}, conds{i}, stdNames{i}{j,1})); % custom title
            
            if err == 1
                legend([mm, nn, ll], 'std','dev','diff') % custom legend
            elseif err == 2
                legend([ll', mm], 'std','dev','diff') % custom legend
            end
            set(gca,'FontSize',12);
            %                     hline(0,'k');
            saveas(gcf, fullfile(Sdir_single, sprintf('CUSTOM_diff_betw_GM%s_%s_%s_ERFs_single.pdf', groups{g}, conds{i}, stdNames{i}{j,1})),'pdf'); % save with custom title
            
            %                 close all  % closes all open figures
        end
    end
end
% close all  % closes all open figures