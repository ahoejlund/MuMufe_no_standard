%% CUSTOM PLOTTING
% submitting MEG-music-autism for custom plotting on the cluster


clusterconfig('long_running', 0, 'slot', 4); % the job shouldn't take more than 1 hr

groups = {'CTRL', 'ASD'};
chans = {'MAGS', 'GRADS'};

plots = [1];
group = [1 2];
chan = [1 2];

input = cell(length(plots)*length(group)*length(chan),1);
for h = 1:length(plots)
    for i = 1:length(group)
        for j = 1:length(chan)
%             input{i+(h-1)*length(plots)+j+(i-1)*length(group)} = [plots(h) group(i) chan(j)];        % where i indicates the participant's sequential number and group{h}(i) indicates their actual participant id number
            input{j+(i-1)*length(group)} = [plots(h) group(i) chan(j)];        % where i indicates the participant's sequential number and group{h}(i) indicates their actual participant id number
        end
    end
end

cell2mat(input(:))

jobid = job2cluster(@custom_mu_plotting, input);