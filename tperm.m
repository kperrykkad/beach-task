%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION = tperm - Performs sign flip maximum pseudo t on data
% Sent to Kelsey Perrykkad on 18/05/21 by Jonathan Robinson who created it
% for use with the Beach Task output from goodparticipantsNull.m

% INPUT
% - data = parts x time x cond
% - perm[default=1000]= number of permutation to perform
% OUTPUT
% - stat = data structure contains:
%   - stat.GAdata - grand average of data across parts (time x cond)
%   - stat.DIFFdata - diff between conds (parts x time)
%   - stat.real_t  - array of t-values for real data (size of time)
%   - stat.null_t - array of null t-values (size of perm) based on maximum psuedo-t
%   - stat.thresh - t-value for significance
%     at .05, .01, .001
%   - stat.sig - logical values for significance in each timepoint 
%     at .05, .01, .001
%   - stat.ESTp = estimated p value for significance in each timepoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stat = tperm(data,perm)
    
    stat.GAdata = squeeze(nanmean(data,1));
    % default permutations
    if nargin < 2
        perm = 1000;
    end

    % make difference waveforms
    data_diff = squeeze((data(:,:,1)-data(:,:,2)));
    stat.DIFFdata= data_diff;
    % real t value
    [~,~,~,stats] = ttest(data_diff,0);
    stat.real_t = stats.tstat;

    % make vars
    partN = size(data,1);
    null = nan(1,perm);
    for thisPerm = 1:perm

       % random permute but not all or non
       flipN = 0;
       while flipN == 0 || flipN == partN 
           flips = rand(1,partN);
           flips(flips > 0.5)= 1;
           flips(flips <= 0.5)= -1;
           flipN = sum(flips == 1);
       end
        
       %generate sign flips
       data_null = data_diff.*flips';

       % run ttest
       [~,~,~,stats] = ttest(data_null,0);
       null(1,thisPerm) = max(stats.tstat);
       
       disp(['Permutations complete = ',num2str(thisPerm,'%03.f')])
    end
    
    % sort values
    stat.null_t = sort(null);
    
    % add thresholds
    stat.thresh(1,:) = [0.05,stat.null_t(1,floor(prctile(1:perm,95)))];
    stat.thresh(2,:) = [0.01,stat.null_t(1,floor(prctile(1:perm,99)))];
    stat.thresh(3,:) = [0.001,stat.null_t(1,floor(prctile(1:perm,99.9)))];

    % add significance flags
    stat.sig(1,:) =  stat.real_t > stat.thresh(1,2);
    stat.sig(2,:) =  stat.real_t > stat.thresh(2,2);
    stat.sig(3,:) =  stat.real_t > stat.thresh(3,2);
    
    % add estimated p-value of timepoints
    for t = 1:size(stat.real_t,2)
        [~,col] = find(stat.null_t > stat.real_t(1,t));
        if col > 0 
            stat.ESTp(t) = min(col);
        else
            stat.ESTp(t) = perm-1;
        end
        stat.ESTp(t) = (perm-stat.ESTp(t))/perm;
    end
end