function [leadin_erpe, leadout_erpe, prehist, posthist,reclassifyHyp] = eventRelatedPredictionError (trigger, PE, water, hypSelected, params, badtrials, leadinlength, leadoutlength, splitlogic, relatedevent)

% This function will create an ERP style analysis for one participant using triggers set in function (hyp switch or env switch - can split these by particular kinds when calling). These can
% later be made into a grand average. It will make two separate ERPEs for
% lead in and lead out and will shift the event to the nearest selected hypothesis (deletes nans around the trigger). Based on
% eventRelatedPredictionError.m from VarVol1. NOT QUITE AS UNIVERSALLY
% FLEXIBLE AS INTENDED IN THE END
%Edited by Kelsey Perrykkad Feb 2021, Created by Kelsey Perrykkad

%% Set up variables
epoch = [leadinlength,leadoutlength]; %frames before and after event
nTotalTrials = params.nTotalTrials;
nParticipants = params.nParticipants;
if ~isempty(splitlogic)
    conditions = 2; %number of splits in split logic - which splits triggers into two categories based on 1 and 0 (could adjust things in the future to have more flags)
else
    conditions = 1;
end
if ~isempty(relatedevent) %create empty structure for participant wise histograms of hypothesis switches before and after environment switches
    prehist = zeros(nParticipants, leadinlength+1);
    posthist = zeros(nParticipants, leadoutlength+1);
else
    prehist = nan;
    posthist = nan;
end
for c=1:conditions
    leadin_erpe{c} = nan(nParticipants,epoch(1)+1);
    leadout_erpe{c}  = nan(nParticipants,epoch(2)+1);
end
reclassifyHyp = cell(nParticipants,nTotalTrials);

%% Go through and create matrices to average
for p=1:nParticipants
    for c=1:conditions
        leadin_tomean = nan(1,epoch(1)+1);
        leadout_tomean =  nan(1,epoch(2)+1);
        for t = 1:nTotalTrials
            if badtrials(p,t) ==0
                for s = 1:length(trigger{p,t})
                    if isempty(splitlogic) || splitlogic{p,t}(s)==c-1 %only includes relevant triggers based on switch logic - if empty, then everything goes into c cell (1) anyway
                        %define lead in lead out event times, removing nans
                        %surrounding event
                        inTrigger = trigger{p,t}(s);
                        outTrigger = trigger{p,t}(s);
                        if outTrigger<length(PE{p}{t}) %if its the last frame, nan the lead out
                            if isnan(hypSelected{p,t}(outTrigger+1))%change the lead out trigger location if nan after
                                try
                                    outTrigger= min(find(~isnan(hypSelected{p,t}(outTrigger+1:end))))+outTrigger; %now the lead out trigger
                                catch %if there is no following selected frames
                                    outTrigger=nan;
                                end
                            end
                        else
                            outTrigger=nan;
                        end
                        
                        if isnan(hypSelected{p,t}(inTrigger-1))%change the lead in trigger location if nan before
                            try
                                inTrigger = max(find(~isnan(hypSelected{p,t}(1:inTrigger-1)))); %now the lead in trigger
                            catch %if there is no preceding selected frames
                                inTrigger=nan;
                            end
                        end
                        
                        %check for possible reclassification
                        if ~isempty(relatedevent) %only if we're trying to remove hypSwitches
                            if ~isnan(inTrigger) & ~isnan(outTrigger)
                                if any(ismember(inTrigger,relatedevent{p,t})) %check inTriggers for matching location hyp switches
                                    if water{p,t}(inTrigger)~= water{p,t}(outTrigger) && hypSelected{p,t}(inTrigger)~= hypSelected{p,t}(outTrigger) %this is a misclassified hyp switch
                                        reclassifyHyp{p,t}(end+1)= trigger{p,t}(s); %mark it for reclassification
                                        inTrigger = nan;
                                        outTrigger = nan;%nan the triggers so it doesn't epoch
                                    end
                                end
                            end                                    
                        end
                        
                        %create epochs for averaging
                         if ~isnan(inTrigger) & ~isnan(outTrigger)
                            if (inTrigger > epoch(1)) && (outTrigger <= length(PE{p}{t})-epoch(2)) %if both exist and don't need to be buffed out by nans
                                leadin_tomean(size(leadin_tomean,1)+1,:) = PE{p}{t}(inTrigger-epoch(1):inTrigger); % add a row to the to be meaned matrix that contains the epoch for this switch
                                leadout_tomean(size(leadout_tomean,1)+1,:) = PE{p}{t}(outTrigger:outTrigger+epoch(2));% add a row to the to be meaned matrix that contains the epoch for this switch
                            elseif inTrigger <= epoch(1) & outTrigger > length(PE{p}{t})-epoch(2)
                                leadin_tomean(size(leadin_tomean,1)+1,:) = [nan(1,epoch(1)-inTrigger+1), PE{p}{t}(1:inTrigger)];
                                leadout_tomean(size(leadout_tomean,1)+1,:) =[PE{p}{t}(outTrigger:length(PE{p}{t})), nan(1,epoch(2)+outTrigger-length(PE{p}{t}))];
                            elseif inTrigger <= epoch(1) % if it's too early in the trial, pad with nans at the beginning
                                leadin_tomean(size(leadin_tomean,1)+1,:) = [nan(1,epoch(1)-inTrigger+1), PE{p}{t}(1:inTrigger)];
                                leadout_tomean(size(leadout_tomean,1)+1,:) = PE{p}{t}(outTrigger:outTrigger+epoch(2));
                            elseif outTrigger > length(PE{p}{t})-epoch(2) %if it's too late in the trial, pad with nans at the end
                                leadin_tomean(size(leadin_tomean,1)+1,:) = PE{p}{t}(inTrigger-epoch(1):inTrigger);
                                leadout_tomean(size(leadout_tomean,1)+1,:) =[PE{p}{t}(outTrigger:length(PE{p}{t})), nan(1,epoch(2)+outTrigger-length(PE{p}{t}))];
                            else
                                fprintf('erpe error, p%i, t%i,s%i',p,t,s);
                            end
                        elseif ~isnan(inTrigger)
                            %only define inTrigger epoch
                            if inTrigger > epoch(1)
                                leadin_tomean(size(leadin_tomean,1)+1,:) = PE{p}{t}(inTrigger-epoch(1):inTrigger); % add a row to the to be meaned matrix that contains the epoch for this switch
                            else  % if it's too early in the trial, pad with nans at the beginning
                                leadin_tomean(size(leadin_tomean,1)+1,:) = [nan(1,epoch(1)-inTrigger+1), PE{p}{t}(1:inTrigger)];
                            end
                        elseif ~isnan(outTrigger)
                            %only define outTrigger epoch
                            if outTrigger < length(PE{p}{t})-epoch(2)
                                leadout_tomean(size(leadout_tomean,1)+1,:) = PE{p}{t}(outTrigger:outTrigger+epoch(2));% add a row to the to be meaned matrix that contains the epoch for this switch
                            else %if its too late, pad with nans at the end
                                leadout_tomean(size(leadout_tomean,1)+1,:) =[PE{p}{t}(outTrigger:length(PE{p}{t})), nan(1,epoch(2)+outTrigger-length(PE{p}{t}))];
                            end
                        end
                        
                        %hist related events around this event - primarily for
                        %plotting hypothesis switches around environment
                        %switches
                        if ~isempty(relatedevent)
                            if ~isnan(inTrigger)
                                if any(relatedevent{p,t}<=inTrigger & relatedevent{p,t}>inTrigger-epoch(1))
                                    frames = relatedevent{p,t}(find(relatedevent{p,t}<=inTrigger & relatedevent{p,t}>inTrigger-epoch(1)));
                                    for h = 1:length(frames)
                                        f = epoch(1)+1-(inTrigger-frames(h));
                                        prehist(p,f) =  prehist(p,f)+1;
                                    end
                                end
                            end
                            if ~isnan(outTrigger)
                                if any(relatedevent{p,t}>=outTrigger & relatedevent{p,t}<outTrigger+epoch(2))
                                    frames = relatedevent{p,t}(find(relatedevent{p,t}>=outTrigger & relatedevent{p,t}<outTrigger+epoch(2)));
                                    for h = 1:length(frames)
                                        f = frames(h)-outTrigger+1;
                                        posthist(p,f) = posthist(p,f)+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        %% Average Matrices into per participant erpes
        leadin_erpe{c}(p,:) = mean(leadin_tomean,1,'omitnan');
        leadout_erpe{c}(p,:) = mean(leadout_tomean,1,'omitnan');
    end
end

