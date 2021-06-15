function [slopebyAA, slopebyDomEnv, slopebyEnv] = slopeAnalysis(PE, accuracy, judgedAgency, domWater, water, badtrials, keep)
%% Slope calculations including buffering trial length
%created by Kelsey Perrykkad on 21/01/21

%param extraction
nParticipants = length(PE);
nTotalTrials = length(PE{1});

%set up output
slopebyAA = nan(nParticipants, 4);
slopebyDomEnv = nan(nParticipants,2);
slopebyEnv = nan(nParticipants,2);

for p=1:nParticipants
    if ~isempty(PE{p})
        % pad out PE in each trial to 450 frames
        %min(cellfun('size',PE{p},2)); %if you decide to downsample instead
        %of upsample this gives the shortest trial length for each
        %participant
        for t=1:nTotalTrials
            PEUp(t,:) = resample(PE{p}{t},450,length(PE{p}{t}),0); %upsampling PE with nearest neighbour antialliasing
            waterUp(t,:) = resample(water{p,t},450,length(water{p,t}),0); %upsampled logical about water location
        end
        % average over all trials in the different conditions to generate
        % slopes per condition of interest
        %% by Accuracy and Agency
        if length(find(accuracy(p,:)==1 & judgedAgency(p,:)==1& badtrials(p,:)==0))>1
            PEcorrectAgency(p,:) = mean(PEUp(find(accuracy(p,:)==1 & judgedAgency(p,:)==1 & badtrials(p,:)==0),:));
            fit = polyfit(1:length(PEcorrectAgency(p,find(~isnan(PEcorrectAgency(p,:))))), PEcorrectAgency(p,find(~isnan(PEcorrectAgency(p,:)))),1); %fit a polynomial of order 1 (a line) ignoring nan values, and put the gradient in the output (much simpler than gradientforstats.m)
            slopebyAA(p,1) = fit(1,1);
        elseif length(find(accuracy(p,:)==1 & judgedAgency(p,:)==1& badtrials(p,:)==0))==1
            PEcorrectAgency(p,:) = PEUp(find(accuracy(p,:)==1 & judgedAgency(p,:)==1& badtrials(p,:)==0),:);
            fit = polyfit(1:length(PEcorrectAgency(p,find(~isnan(PEcorrectAgency(p,:))))), PEcorrectAgency(p,find(~isnan(PEcorrectAgency(p,:)))),1); %fit a polynomial of order 1 (a line), and put the gradient in the output (much simpler than gradientforstats.m)
            slopebyAA(p,1) = fit(1,1);
        elseif length(find(accuracy(p,:)==1 & judgedAgency(p,:)==1& badtrials(p,:)==0))==0
        end
        
        
        if length(find(accuracy(p,:)==1& judgedAgency(p,:)==0& badtrials(p,:)==0))>1
            PEcorrectNoAgency(p,:) = mean(PEUp(find(accuracy(p,:)==1 & judgedAgency(p,:)==0& badtrials(p,:)==0),:));
            fit = polyfit(1:length(PEcorrectNoAgency(p,find(~isnan(PEcorrectNoAgency(p,:))))), PEcorrectNoAgency(p,find(~isnan(PEcorrectNoAgency(p,:)))),1);
            slopebyAA(p,2) = fit(1,1);
        elseif length(find(accuracy(p,:)==1& judgedAgency(p,:)==0& badtrials(p,:)==0))==1
            PEcorrectNoAgency(p,:) = PEUp(find(accuracy(p,:)==1& judgedAgency(p,:)==0& badtrials(p,:)==0),:);
            fit = polyfit(1:length(PEcorrectNoAgency(p,find(~isnan(PEcorrectNoAgency(p,:))))), PEcorrectNoAgency(p,find(~isnan(PEcorrectNoAgency(p,:)))),1);
            slopebyAA(p,2) = fit(1,1);
        elseif length(find(accuracy(p,:)==1& judgedAgency(p,:)==0& badtrials(p,:)==0))==0
        end
        
        
        if length(find(accuracy(p,:)==0 & judgedAgency(p,:)==1& badtrials(p,:)==0))>1
            PEincorrectAgency(p,:) = mean(PEUp(find(accuracy(p,:)==0 & judgedAgency(p,:)==1& badtrials(p,:)==0),:));
            fit = polyfit(1:length(PEincorrectAgency(p,find(~isnan(PEincorrectAgency(p,:))))), PEincorrectAgency(p,find(~isnan(PEincorrectAgency(p,:)))),1);
            slopebyAA(p,3) = fit(1,1);
        elseif length(find(accuracy(p,:)==0 & judgedAgency(p,:)==1& badtrials(p,:)==0))==1
            PEincorrectAgency(p,:) = PEUp(find(accuracy(p,:)==0 & judgedAgency(p,:)==1& badtrials(p,:)==0),:);
            fit = polyfit(1:length(PEincorrectAgency(p,find(~isnan(PEincorrectAgency(p,:))))), PEincorrectAgency(p,find(~isnan(PEincorrectAgency(p,:)))),1);
            slopebyAA(p,3) = fit(1,1);
        elseif length(find(accuracy(p,:)==0 & judgedAgency(p,:)==1& badtrials(p,:)==0))==0
        end
        
        
        if length(find(accuracy(p,:)==0& judgedAgency(p,:)==0& badtrials(p,:)==0))>1
            PEincorrectNoAgency(p,:) = mean(PEUp(find(accuracy(p,:)==0& judgedAgency(p,:)==0& badtrials(p,:)==0),:));
            fit = polyfit(1:length(PEincorrectNoAgency(p,find(~isnan(PEincorrectNoAgency(p,:))))), PEincorrectNoAgency(p,find(~isnan(PEincorrectNoAgency(p,:)))),1);
            slopebyAA(p,4) = fit(1,1);
        elseif length(find(accuracy(p,:)==0& judgedAgency(p,:)==0& badtrials(p,:)==0))==1
            PEincorrectNoAgency(p,:) = PEUp(find(accuracy(p,:)==0& judgedAgency(p,:)==0& badtrials(p,:)==0),:);
            fit = polyfit(1:length(PEincorrectNoAgency(p,find(~isnan(PEincorrectNoAgency(p,:))))), PEincorrectNoAgency(p,find(~isnan(PEincorrectNoAgency(p,:)))),1);
            slopebyAA(p,4) = fit(1,1);
        elseif length(find(accuracy(p,:)==0& judgedAgency(p,:)==0& badtrials(p,:)==0))==0
        end
        
        
        %% by Dom Env
        if length(find(domWater(p,:)==1& badtrials(p,:)==0))>1
            PEwaterDom(p,:) = mean(PEUp(find(domWater(p,:)==1& badtrials(p,:)==0),:));
            fit = polyfit(1:length(PEwaterDom(p,find(~isnan(PEwaterDom(p,:))))),PEwaterDom(p,find(~isnan(PEwaterDom(p,:)))),1);
            slopebyDomEnv(p,1) = fit(1,1);
        elseif length(find(domWater(p,:)==1& badtrials(p,:)==0))==1
            PEwaterDom(p,:) = PEUp(find(domWater(p,:)==1& badtrials(p,:)==0),:);
            fit = polyfit(1:length(PEwaterDom(p,find(~isnan(PEwaterDom(p,:))))),PEwaterDom(p,find(~isnan(PEwaterDom(p,:)))),1);
            slopebyDomEnv(p,1) = fit(1,1);
        elseif length(find(domWater(p,:)==1& badtrials(p,:)==0))==0
        end
        
        
        if length(find(domWater(p,:)==0& badtrials(p,:)==0))>1
            PEsandDom(p,:) = mean(PEUp(find(domWater(p,:)==0& badtrials(p,:)==0),:));
            fit = polyfit(1:length(PEsandDom(p,find(~isnan(PEsandDom(p,:))))),PEsandDom(p,find(~isnan(PEsandDom(p,:)))),1);
            slopebyDomEnv(p,2) = fit(1,1);
        elseif length(find(domWater(p,:)==0& badtrials(p,:)==0))==1
            PEsandDom(p,:) = PEUp(find(domWater(p,:)==0& badtrials(p,:)==0),:);
            fit = polyfit(1:length(PEsandDom(p,find(~isnan(PEsandDom(p,:))))),PEsandDom(p,find(~isnan(PEsandDom(p,:)))),1);
            slopebyDomEnv(p,2) = fit(1,1);
        elseif length(find(domWater(p,:)==0& badtrials(p,:)==0))==0
        end
        
        %% By Environment on each frame
        for t=1:nTotalTrials
            if badtrials(p,t) ==0
                PESandtoMean(t,1:450) = nan(1,450);
                PEWatertoMean(t,1:450) = nan(1,450);
                PEWatertoMean(t,find(waterUp(t,:)==1))=PEUp(t,find(waterUp(t,:)==1));
                PESandtoMean(t,find(waterUp(t,:)==0))=PEUp(t,find(waterUp(t,:)==0));
            end
        end
        PEwater(p,:) = nanmean(PEWatertoMean);
        fit = polyfit(1:length(PEwater(p,find(~isnan(PEwater(p,:))))),PEwater(p,find(~isnan(PEwater(p,:)))),1);
        slopebyEnv(p,1) = fit(1,1);
        
        PEsand(p,:) = nanmean(PESandtoMean);
        fit = polyfit(1:length(PEsand(p,find(~isnan(PEsand(p,:))))),PEsand(p,find(~isnan(PEsand(p,:)))),1);
        slopebyEnv(p,2) = fit(1,1);
        
    end
end

%% Optional Plots

%% Accuracy by Agency
figure
hold on
plot(1:450,nanmean(PEcorrectAgency(keep,:)),'Color',[8/255,72/255,135/255])
fit = polyfit(1:450,nanmean(PEcorrectAgency(keep,:)),1);
pts = polyval(fit,1:450);
plot(1:450,pts,'--','Color',[8/255,72/255,135/255], 'LineWidth', 3);
plot(1:450,nanmean(PEcorrectNoAgency(keep,:)),'Color',[144/255,156/255,194/255])
fit = polyfit(1:450,nanmean(PEcorrectNoAgency(keep,:)),1);
pts = polyval(fit,1:450);
plot(1:450,pts,'--','Color',[144/255,156/255,194/255], 'LineWidth', 3);
plot(1:450,nanmean(PEincorrectAgency(keep,:)),'Color',[245/255,138/255,7/255])
fit = polyfit(1:450,nanmean(PEincorrectAgency(keep,:)),1);
pts = polyval(fit,1:450);
plot(1:450,pts,'--','Color',[245/255,138/255,7/255], 'LineWidth', 3);
plot(1:450,nanmean(PEincorrectNoAgency(keep,:)),'Color',[249/255,171/255,85/255])
fit = polyfit(1:450,nanmean(PEincorrectNoAgency(keep,:)),1);
pts = polyval(fit,1:450);
plot(1:450,pts,'--','Color',[249/255,171/255,85/255], 'LineWidth', 3);
legend('Correct Agency','Correct Agency','Correct No Agency','Correct No Agency','Incorrect Agency','Incorrect Agency','Incorrect No Agency','Incorrect No Agency')


%% Dominant Env
figure
hold on
plot(1:450,nanmean(PEwaterDom(keep,:)),'Color',[108/255,166/255,193/255])
fit = polyfit(1:450,nanmean(PEwaterDom(keep,:)),1);
pts = polyval(fit,1:450);
plot(1:450,pts,'--','Color',[108/255,166/255,193/255], 'LineWidth', 3);
plot(1:450,nanmean(PEsandDom(keep,:)),'Color',[245/255,204/255,0/255])
fit = polyfit(1:450,nanmean(PEsandDom(keep,:)),1);
pts = polyval(fit,1:450);
plot(1:450,pts,'--','Color',[245/255,204/255,0/255], 'LineWidth', 3);
legend('Water Dominant Trials','Water Dominant Trials','Sand Dominant Trials', 'Sand Dominant Trials')

%% Env per frame
figure
hold on
plot(1:450,nanmean(PEwater(keep,:)),'Color',[23/255,37/255,90/255])
fit = polyfit(1:450,nanmean(PEwater(keep,:)),1);
pts = polyval(fit,1:450);
plot(1:450,pts,'--','Color',[23/255,37/255,90/255], 'LineWidth', 3);
plot(1:450,nanmean(PEsand(keep,:)),'Color',[189/255,153/255,117/255])
fit = polyfit(1:450,nanmean(PEsand(keep,:)),1);
pts = polyval(fit,1:450);
plot(1:450,pts,'--','Color',[189/255,153/255,117/255], 'LineWidth', 3);
legend('Water Frames','Water Frames','Sand Frames', 'Sand Frames')
