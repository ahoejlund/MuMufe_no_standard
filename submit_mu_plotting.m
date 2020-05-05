%% CUSTOM PLOTTING
% submitting MEG-music-autism for custom plotting on the cluster


clusterconfig('long_running', 0, 'slot', 2); % the job shouldn't take more than 1 hr

groups = {'CTRL', 'ASD'};
chans = {'MAGS', 'GRADS'};

plots = [2]; % 1=grandmean; 2=indivs
group = [1 2];  %{'CTRL', 'ASD'};
chan = [2];  %{'MAGS', 'GRADS'};

input = cell(length(plots)*length(group)*length(chan),1);
for h = 1:length(plots)
    for i = 1:length(group)
        for j = 1:length(chan)
%             input{i+(h-1)*length(plots)+j+(i-1)*length(group)} = [plots(h) group(i) chan(j)];
%             input{j+(i-1)*length(group)} = [plots(h) group(i) chan(j)];        
            input{j+(i-1)*length(chan)} = [plots(h) group(i) chan(j)];
        end
    end
end

cell2mat(input(:))

jobid = job2cluster(@custom_mu_plotting, input);