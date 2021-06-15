function [trialtable, participanttable] = TrialwiseResults(permresults)
%% Output for Statistics Program
%Created by Kelsey Perrykkad Dec 2020
%Loads in results file from BeachAnalysis.m and structures it trialwise and outputs
%it as a csv that can be loaded into a stats program.
%Only keeps non-eliminated participants
%Comparable to longform_trialwise.m

%Choose results file to convert through dialogue
folder = 'G:\My Drive\Experiments\JoA_squarestask\BeachTask\Data and Analysis\Analysis\Results\';
resultsfile = uigetfile(folder, 'Pick results file to convert');
suffix = erase(resultsfile,"results");
paramsfile = ['params' suffix];
load([folder resultsfile]); %will load a struct named 'results'
try
load([folder paramsfile]);
catch
    paramsfile = uigetfile(folder, 'Automatic params matching failed. Pick matching params file:');
    load([folder paramsfile]);
end

%% set up empty trialwise table
dataheaders = {'workerId','anonymousId','trial','nFrames','lags','nolagfreq','totalLagDuration','badtrial','noControl','accuracy','judgedAgency','chosenSquare','SoA','confidence','percFramesMoving','nHypSwitches','domWater','percWater', 'nWaves','modeWaveTime','sdWaveTime','percentChosenSelected','percentCorrectSelected', 'nTotalEnvSwitches','nTraverseEnvSwitches','nHypSwitchEnvSwitches','nQrEnvSwitches','avPE','speed','acceleration','jerk'};
sz = [length(results.keep)*params.nTotalTrials length(dataheaders)];
varTypes = {'string','double','double','double','double','double','double','logical','logical','logical','logical','double','double','double','double','double','logical','double','double','double','double','double','double','double','double','double','double','double','double','double','double'};
trialtable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',dataheaders);

%fill in trialwise variables
r = 1; %row counter
w = waitbar(0,'Starting trialwise output...');
for i = 1:length(results.keep)
    p = results.keep(i);
    waitbar(i/length(results.keep),w,sprintf('Trial-wise Progress: %d %%', floor(i/length(results.keep)*100)));
    for t = 1:params.nTotalTrials
        trialtable.workerId(r) = results.keepId(i);
        trialtable.anonymousId(r) = p;
        trialtable.trial(r) = t;
        trialtable.nFrames(r) = results.nFrames(p,t);
        trialtable.lags(r) = results.lagCount(p,t);
        trialtable.nolagfreq(r) = results.noLagFreq(p,t);
        trialtable.totalLagDuration(r) = results.totalLagDuration(p,t);
        trialtable.badtrial(r) = results.badtrials(p,t);
        trialtable.noControl(r) = results.noControl(p,t);
        trialtable.accuracy(r) = results.accuracy(p,t);
        trialtable.judgedAgency(r) = results.judgedAgency(p,t);
        trialtable.chosenSquare(r) = results.chosenSquare(p,t);
        trialtable.SoA(r) = results.senseOfAgency(p,t);
        trialtable.confidence(r) = results.confidence(p,t);
        trialtable.percFramesMoving(r) = results.percentFramesSpentMoving(p,t);
        trialtable.nHypSwitches(r) = results.nHypSwitches(p,t);
        trialtable.domWater(r) = results.domWater(p,t);
        trialtable.percWater(r) = results.percWater(p,t);
        trialtable.nWaves(r) = results.nWaves(p,t);
        trialtable.sdWaveTime(r) = results.sdWaveTime(p,t);
        trialtable.modeWaveTime(r) = results.modeWaveTime(p,t);
        trialtable.percentChosenSelected(r) = results.percentChosenSelected(p,t);
        trialtable.percentCorrectSelected(r) = results.percentCorrectSelected(p,t);
        trialtable.nTotalEnvSwitches(r) = results.nTotalEnvSwitches(p,t);
        trialtable.nTraverseEnvSwitches(r) = results.nTraverseEnvSwitches(p,t);
        trialtable.nHypSwitchEnvSwitches(r) = results.nHypSwitchEnvSwitches(p,t);
        trialtable.nQrEnvSwitches(r) = results.nQrEnvSwitches(p,t);
        trialtable.avPE(r) = results.avPE(p,t);
        trialtable.speed(r) = results.speed(p,t);
        trialtable.acceleration(r) = results.acceleration(p,t);
        trialtable.jerk(r) = results.jerk(p,t);
        
        r=r+1;
    end
end


%% set up empty participantwise table
dataheaders = {'workerId','anonymousId','nFrames','sdWaveTime','avAccuracy','dprime','criterion','avPercWater','percDomWater','dprimeWater','criterionWater','dprimeSand','criterionSand','avSoA','avConfidence', 'B1WaterDom','B3WaterDom','B1MeanPercentWater','B3MeanPercentWater'};
sz = [length(results.keep) length(dataheaders)];
varTypes = {'string','string','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double'};
participanttable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',dataheaders);

%participantwise variables
r = 1;
for i = 1:length(results.keep)
    p = results.keep(i);
    waitbar(i/length(results.keep),w,sprintf('Participant-wise Progress: %d %%', floor(i/length(results.keep)*100)));
    participanttable.workerId(r) = results.keepId(i);
    participanttable.anonymousId(r) = p;
    participanttable.nFrames(r) = mean(results.nFrames(p,find(~results.badtrials(p,:))));
    participanttable.sdWaveTime(r) = mean(results.sdWaveTime(p,find(~results.badtrials(p,:))));
    participanttable.avAccuracy(r) = results.avAccuracy(p);
    participanttable.dprime(r) = results.dprime(p);
    participanttable.criterion(r) = results.criterion(p);
    participanttable.avPercWater(r) = results.avPercWater(p);
    participanttable.percDomWater(r) = results.percDomWater(p);
    participanttable.dprimeWater(r) = results.dprimeWaterDom(p);
    participanttable.criterionWater(r) = results.criterionWaterDom(p);
    participanttable.dprimeSand(r) = results.dprimeSandDom(p);
    participanttable.criterionSand(r) = results.criterionSandDom(p);
    participanttable.avSoA(r) = results.avSoA(p);
    participanttable.avConfidence(r) = results.avConfidence(p);
    participanttable.B1WaterDom(r) = results.B1WaterDom(p);
    participanttable.B3WaterDom(r) = results.B3WaterDom(p);
    participanttable.B1MeanPercentWater(r) = results.B1MeanPercentWater(p);
    participanttable.B3MeanPercentWater(r) = results.B3MeanPercentWater(p);
    
    r=r+1;
end

%%  set up empty conditionwise tables
%set up acc/agency slopes
dataheaders = {'workerId','anonymousId','accuracy','agency','slope'};
sz = [length(results.keep)*4 length(dataheaders)];
varTypes = {'string','string','logical','logical','double'};
aaslopetable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',dataheaders);
%add variables
r = 1;
for i = 1:length(results.keep)
    p = results.keep(i);
    waitbar(i/length(results.keep),w,sprintf('AAslope Progress: %d %%', floor(i/length(results.keep)*100)));
    for c=1:4
        aaslopetable.workerId(r) = results.keepId(i);
        aaslopetable.anonymousId(r) = p;
        if c==1
            aaslopetable.accuracy(r)=1;
            aaslopetable.agency(r)=1;
            aaslopetable.slope(r) = results.slopebyAA(p,1);
            r=r+1;
        elseif c==2
            aaslopetable.accuracy(r)=1;
            aaslopetable.agency(r)=0;
            aaslopetable.slope(r) = results.slopebyAA(p,2);
            r=r+1;
        elseif c==3
            aaslopetable.accuracy(r)=0;
            aaslopetable.agency(r)=1;
            aaslopetable.slope(r) = results.slopebyAA(p,3);
            r=r+1;
        elseif c==4
            aaslopetable.accuracy(r)=0;
            aaslopetable.agency(r)=0;
            aaslopetable.slope(r) = results.slopebyAA(p,4);
            r=r+1;
        end
    end
end

%set up env slopes
dataheaders = {'workerId','anonymousId','waterframes','slope'};
sz = [length(results.keep)*2 length(dataheaders)];
varTypes = {'string','string','logical','double'};
envslopetable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',dataheaders);
%add variables
r = 1;
for i = 1:length(results.keep)
    p = results.keep(i);
    %waitbar(i/length(results.keep),w,sprintf('Env Slope Progress: %d %%', floor(i/length(results.keep)*100)));
    for c=1:2
        envslopetable.workerId(r) = results.keepId(i);
        envslopetable.anonymousId(r) = p;
        if c==1
            envslopetable.waterframes(r) = 1;
            envslopetable.slope(r) = results.slopebyEnv(p,1);
            r = r+1;
        elseif c==2
            envslopetable.waterframes(r) = 0;
            envslopetable.slope(r) = results.slopebyEnv(p,2);
            r = r+1;
        end
    end
end

%set up env erpe tables - old version pre-perm testing
dataheaders = {'workerId','anonymousId','timebin','directionToS1W2','avPE','dMouse','hypcount'};
sz = [length(results.keep)*12 length(dataheaders)];
varTypes = {'string','string','double','double','double','double','double'};
enverpetable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',dataheaders);
%add variables
r = 1;
if isempty(permresults) %old results
    for i = 1:length(results.keep)
        p = results.keep(i);
        waitbar(i/length(results.keep),w,sprintf('Env ERPE Progress: %d %%', floor(i/length(results.keep)*100)));
        for tb = 1:6 %5 frame bins
            for direction = 1:2
                enverpetable.workerId(r) = results.keepId(i);
                enverpetable.anonymousId(r) = p;
                enverpetable.timebin(r) = tb;
                enverpetable.directionToS1W2(r) = direction;
                if tb<4 %lead in
                    enverpetable.avPE(r)= nanmean(results.ESnoHS_leadin_erpe{direction}(p,tb*5-4:tb*5));
                else %lead out
                    enverpetable.avPE(r)= nanmean(results.ESnoHS_leadout_erpe{direction}(p,tb*5-19:tb*5-15));
                end
                if tb==1
                    enverpetable.hypcount(r)= nanmean(results.envHyp_prehist(p,tb*5-3:tb*5)); %eliminates that first 0 column
                elseif    tb<4 %lead in
                    enverpetable.hypcount(r)= nanmean(results.envHyp_prehist(p,tb*5-4:tb*5));
                elseif tb==6
                    enverpetable.hypcount(r)= nanmean(results.envHyp_posthist(p,tb*5-19:tb*5-16)); %eliminates that last 0 column
                else %lead out
                    enverpetable.hypcount(r)= nanmean(results.envHyp_posthist(p,tb*5-19:tb*5-15));
                end
                if tb<4 %lead in
                    enverpetable.dMouse(r)= nanmean(results.ESnoHS_leadin_erSpeed{direction}(p,tb*5-4:tb*5));
                    r=r+1;
                else %lead out
                    enverpetable.dMouse(r)= nanmean(results.ESnoHS_leadout_erSpeed{direction}(p,tb*5-19:tb*5-15));
                    r=r+1;
                end
            end
        end
    end
else
     for i = 1:length(results.keep)
        p = results.keep(i);
        waitbar(i/length(results.keep),w,sprintf('Env ERPE Progress: %d %%', floor(i/length(results.keep)*100)));
        for tb = 1:3 %bins based on permutation significance period (manually delineated - could be pretty easily automated in future)
            for direction = 1:2
                enverpetable.workerId(r) = results.keepId(i);
                enverpetable.anonymousId(r) = p;
                enverpetable.timebin(r) = tb;
                enverpetable.directionToS1W2(r) = direction;
                if tb==1 %pre
                    enverpetable.avPE(r)= nanmean(results.ESnoHS_leadin_erpe{direction}(p,1:9));
                elseif tb==2 %sig peak cluster (across lead in and lead out)
                    enverpetable.avPE(r)= nanmean([results.ESnoHS_leadin_erpe{direction}(p,10:15) results.ESnoHS_leadout_erpe{direction}(p,1:8)]);
                elseif tb ==3 %post
                    enverpetable.avPE(r)= nanmean(results.ESnoHS_leadout_erpe{direction}(p,9:15));
                end
                if tb==1 %pre
                    enverpetable.hypcount(r)= nanmean(results.envHyp_prehist(p,1:9)); 
                elseif    tb==2 %sig peak cluster
                    enverpetable.hypcount(r)= nanmean([results.envHyp_prehist(p,10:15) results.envHyp_posthist(p,1:8)]);
                elseif tb==3 %post
                    enverpetable.hypcount(r)= nanmean(results.envHyp_posthist(p,9:15)); 
                end
                if tb==1 %pre
                    enverpetable.dMouse(r)= nanmean(results.ESnoHS_leadin_erSpeed{direction}(p,1:9));
                    r=r+1;
                elseif tb == 2%sig peak cluster
                    enverpetable.dMouse(r)= nanmean([results.ESnoHS_leadin_erSpeed{direction}(p,10:15) results.ESnoHS_leadout_erSpeed{direction}(p,1:8)]);
                    r=r+1;
                elseif tb==3 %post
                    enverpetable.dMouse(r)= nanmean(results.ESnoHS_leadout_erSpeed{direction}(p,9:15));
                    r=r+1;
                end
            end
        end
    end
end

dataheaders = {'workerId','anonymousId','timebin','directionToS1W2','avPE','dMouse','hypcount'};
sz = [length(results.keep)*12 length(dataheaders)];
varTypes = {'string','string','double','double','double','double','double'};
enverpepeaktable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',dataheaders);
%add variables
r = 1;
if ~isempty(permresults) %new results
     for i = 1:length(results.keep)
        p = results.keep(i);
        waitbar(i/length(results.keep),w,sprintf('Env ERPE Progress: %d %%', floor(i/length(results.keep)*100)));
        for tb = 1:14 %bins based on permutation significance period (manually delineated - could be pretty easily automated in future)
            for direction = 1:2
                enverpepeaktable.workerId(r) = results.keepId(i);
                enverpepeaktable.anonymousId(r) = p;
                enverpepeaktable.timebin(r) = tb;
                enverpepeaktable.directionToS1W2(r) = direction;
                if tb<=6 %pre (10:15)
                    enverpepeaktable.avPE(r)= results.ESnoHS_leadin_erpe{direction}(p,9+tb);
                elseif tb>6 %post (16:23)
                    enverpepeaktable.avPE(r)= results.ESnoHS_leadout_erpe{direction}(p,1+(tb-7));
                end
                if tb<=6 %pre
                    enverpepeaktable.hypcount(r)= results.envHyp_prehist(p,9+tb); 
                elseif tb>6  %post
                    enverpepeaktable.hypcount(r)= results.envHyp_posthist(p,1+(tb-7)); 
                end
                if tb<=6 %pre
                    enverpepeaktable.dMouse(r)= results.ESnoHS_leadin_erSpeed{direction}(p,9+tb);
                    r=r+1;
                elseif tb>6 %post
                    enverpepeaktable.dMouse(r)= results.ESnoHS_leadout_erSpeed{direction}(p,1+(tb-7));
                    r=r+1;
                end
            end
        end
    end
end

%set up hyp switch erpe table
dataheaders = {'workerId','anonymousId','timebin','avPE'};
sz = [length(results.keep)*9 length(dataheaders)];
varTypes = {'string','string','double','double'};
hyperpetable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',dataheaders);
%add variables
r = 1;
for i = 1:length(results.keep)
    p = results.keep(i);
    waitbar(i/length(results.keep),w,sprintf('Hyp Switch ERPE Progress: %d %%', floor(i/length(results.keep)*100)));
    for tb = 1:9 %5 frame bins
        hyperpetable.workerId(r) = results.keepId(i);
        hyperpetable.anonymousId(r) = p;
        hyperpetable.timebin(r) = tb;
        if tb<7 %lead in
            hyperpetable.avPE(r)= nanmean(results.HS_leadin_erpe{1}(p,tb*5-4:tb*5));
            r=r+1;
        else %lead out
            hyperpetable.avPE(r)= nanmean(results.HS_leadout_erpe{1}(p,tb*5-34:tb*5-30));
            r=r+1;
        end
    end
end

close(w)

filename = [folder 'trialwise' erase(suffix,".mat") '.csv'];
writetable(trialtable,filename);
filename = [folder 'participantwise' erase(suffix,".mat") '.csv'];
writetable(participanttable,filename);
filename = [folder 'aaslopes' erase(suffix,".mat") '.csv'];
writetable(aaslopetable,filename);
filename = [folder 'envslopes' erase(suffix,".mat") '.csv'];
writetable(envslopetable,filename);
filename = [folder 'enverpe' erase(suffix,".mat") '.csv'];
writetable(enverpetable,filename);
filename = [folder 'enverpepeak' erase(suffix,".mat") '.csv'];
writetable(enverpepeaktable,filename);
filename = [folder 'hyperpe' erase(suffix,".mat") '.csv'];
writetable(hyperpetable,filename);