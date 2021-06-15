function [results,params] = BeachAnalysis(dataset,savedata)

%% Analysis of "Online Beach Task" %%
%Created by Kelsey Perrykkad Nov 2020
%Project coauthors: Jonathan Robinson and Jakob Hohwy

%% Description %%
%Moving the var/vol Squares task online to deal with the horrors of 2020.
%50:50 no control and agentive trials half (left/right) the screen has
%predictable variability (pos or neg 'waves' every .5s) and half has
%unpredictable variability (sand). All blocks are the 'same' statistically.
%Programmed to be presented at 30hz -
%note: expected variability in timing of the experiment because it was run
%online - all timing needs to be independently calculated per trial, and
%not 'assumed'
%Variability width and trial orders etc is set by WholeExperimentVariability.m and is
%predefined prior to running the experiment to reduce computational running
%costs. The only thing determined online is mouse movement components of
%each frame and the side the environments will be on.

%This program starts with data parsed from pavlovia using Jonno's
%parse_outMATread.py anaconda/python script, which processes a whole folder
%as downloaded from pavlovia.

%Notes for analysis: any frame based timestamp (eg. for object selection
%button presses) are in the full 900 frames of 15s at 60hz. So to match
%with the 30hz measures, do round(frame/2). Basically anything associated with responses is based on 900 frames all objects and movement variables are every second frame. Order = colour allocation.

%In general, data structure should be participants as rows, and trials as
%columns.

% Will want to eventually do most measures split by 1) Env 2) Accuracy 3)
% Judged Agency

%% Define parameters of the experiment
%This can be carried through to subfunctions, and has general experiment
%details - defined before the for loop through participants
params.datafolder = 'G:\My Drive\Experiments\JoA_squarestask\BeachTask\Data and Analysis\Analysis\Dataset';
params.nPractice = 16;
params.nBlocks = 3;
params.nTrialperBlock = 16;
params.nTotalTrials = params.nBlocks*params.nTrialperBlock; %not including practice
params.maxFrames = 450;
params.trialDuration=15; %in seconds
params.slothSelected=0.0873; %divider from experiment for different elements is called 'sloth'
params.slothDistracter=0.0436;
params.noSelectionSloth=10000;
params.nSquares=4;
params.waveFreq = .5; %sec per wave
params.programmedHz = 30; %how frequently frames should be presented in hz as programmed


%% Set up Results and run through each participant to extract relevant dependent variables into Results
params.filelist = dir(fullfile(params.datafolder,'*experiment.mat')); %creates a struct of file info for all files in the current directory
if dataset~='all'
    params.filelist=params.filelist(1:dataset);
end
params.nParticipants = length(params.filelist);
%set up at least one field in the struct - any others can be used to
%specify a desired order. If the var isn't in this line, it gets added at
%the end as it works through the script
results = struct('workerId', [], 'keep', [],'keepId',[], 'errorId', [], 'errorIndex',[], 'errorMessage', [], 'warns',[],'avFrames',[],'nFrames',[],'noControl',[],'judgedAgency',[],'chosenSquare', [],'accuracy',[],'avAccuracy',[],'confidence',[],'senseOfAgency',[],'percentSpentMoving',[],'speed',[],'acceleration',[],'jerk',[],'hypSwitches',[],'percWater',[],'domWater',[],'envSwitches',[],'PE', [], 'avPE',[],'slopebyAA',[],'slopebyEnv',[]); %if this is incomplete it's fine - should add undefined variables

for p=1:params.nParticipants
    % p = 7%for testing with one participant
    filename = [params.filelist(p).folder,'\',params.filelist(p).name];
    data = load(filename);
    
    %% Participant info
    results.os{p} = data.OS{1};
    results.browser{p} = data.browserDet{1};
    results.workerId{p} = data.workerId{1};
    try
        %% Trial Details
        for t = params.nPractice+1:params.nTotalTrials+params.nPractice
            if data.corrAns(t)==0
                results.noControl(p,t-params.nPractice) = 1;
            else
                results.noControl(p,t-params.nPractice) = 0;
            end
            if data.key_resp_whichkeys(t) == 0
                results.judgedAgency(p,t-params.nPractice) = 0;
            elseif data.key_resp_whichkeys(t) > 0
                results.judgedAgency(p,t-params.nPractice) = 1;
            end
        end
        
        %% Number of Warnings for Inactivity (not using selection buttons at all during a trial)
        %              if isfield(data,'warn')
        %                  results.nWarns(p) = nansum(data.warn);
        %                  results.warns(p,:) = data.warn;
        %              else %for pilots before warns was a thing - this still skips
        %              some for pilot -- recalculates warns - warns is not properly
        %              counted in existing data
        results.nWarns(p) = length(find(cellfun(@isempty,data.track_Model_key(params.nPractice+1:end))));
        findinactivity = find(cellfun(@isempty,data.track_Model_key(params.nPractice+1:end)));
        for t=1:params.nPractice
            data.warn(t)=nan;
        end
        for t = params.nPractice+1:params.nTotalTrials+params.nPractice
            if (find(t==findinactivity+params.nPractice)>0)==1
                data.warn(t)=1;
            else
                data.warn(t)=0;
            end
        end
        results.warns(p,:)=data.warn(params.nPractice+1:end);
        %              end
       
               
        %% Count screen minimisation from windims
        %set expected window dimensions
        win1 = data.winDims{1};
        win2 = data.winDims{2};
        counter =3;
        if win1 ~= win2 & counter<params.nTotalTrials+params.nPractice
            while win1 ~= win2 & counter<params.nTotalTrials+params.nPractice
                win1 = win2;
                win2 = data.winDims{counter};
                counter = counter+1;
            end
        end
        %count resizes
        results.screenResizeCount(p) = 0;
        for t = 1:params.nTotalTrials
            if data.winDims{t+params.nPractice} ~= win1
                results.screenResizeCount(p)= results.screenResizeCount(p)+1;
            end
        end
        %% Framerate
        for t = params.nPractice+1:params.nTotalTrials+params.nPractice
            results.nFrames(p,t-params.nPractice) = length(data.track_logrowT{t});
        end
        %optional output of figures
        %      figure
        %      plot(results.nFrames(p,:));
        
        %% Define Bad Trials
        for t = params.nPractice+1:params.nTotalTrials+params.nPractice %data trial numbers (17 is trial 1)
            [results.lagLocation{p,t-params.nPractice},results.noLagFreq(p,t-params.nPractice),results.totalLagDuration(p,t-params.nPractice)] = frameDiagnosis(data.track_logrowT{t}, params.maxFrames, params.trialDuration);%count lagged frames
            if length(data.track_logrowT{t})< params.maxFrames*.5 % <50% frame presentation
                results.badtrials(p,t-params.nPractice) = 1;
                %elseif %lagcount
            elseif results.noLagFreq(p,t-params.nPractice)>2*params.trialDuration/params.maxFrames %if the non lag times are more than twice as slow as they should be
                results.badtrials(p,t-params.nPractice) = 1;
            elseif results.warns(p,t-params.nPractice)==1
                results.badtrials(p,t-params.nPractice) = 1;
            elseif data.warnX(t)==1 && data.warnY(t)==1
                results.badtrials(p,t-params.nPractice) = 1;
            else
                results.badtrials(p,t-params.nPractice) = 0; %updated for percent spent moving just below
            end
            results.lagCount(p,t-params.nPractice) = length(results.lagLocation{p,t-params.nPractice});
        end
        
        %Update badtrials after trial spent moving duration 
        goodtrials = find(results.badtrials(p,:)==0); %this variable is not saved, is just used for averaging across vars at the end of this loop
        results.countBad(p) = sum(results.badtrials(p,:));
        
        %% Wave regularity per trial - in seconds (with and without lags?)
        %number of waves
        %av. sec per wave
        framesPerWave = params.programmedHz*params.waveFreq;
        for t = 1:params.nTotalTrials
            waveTimes = [];
            for w = 1:floor((length(data.track_logrowT{t+params.nPractice})-1)/framesPerWave) %for every complete wave in that trial
                if w > 1
                    waveTimes(end+1) = data.track_logrowT{t+params.nPractice}((framesPerWave*w)+1)-sum(waveTimes); %
                else
                    waveTimes(end+1) = data.track_logrowT{t+params.nPractice}(framesPerWave);
                end
            end
            results.modeWaveTime(p,t) = mode(waveTimes);
            results.sdWaveTime(p,t) = std(waveTimes);
            results.nWaves(p,t) = length(waveTimes);
        end
        
        %% Confidence Judgement
        for t = params.nPractice+1:params.nTotalTrials+params.nPractice
            results.confidence(p,t-params.nPractice) = data.key_resp_confkeys(t);
        end
        
        %% Sense of Agency Ratings
        for t = params.nPractice+1:params.nTotalTrials+params.nPractice
            results.senseOfAgency(p,t-params.nPractice) = data.key_resp_agentkeys(t);
        end
        results.avSoA(p) = nanmean(results.senseOfAgency(p,:));
        
        %% Chosen Square
        for t = params.nPractice+1:params.nTotalTrials+params.nPractice
            results.chosenSquare(p,t-params.nPractice) = data.key_resp_whichkeys(t);
        end
        results.histChosenSq(p,:)=hist(results.chosenSquare(p,:),5);
        
        %% Accuracy
        results.accuracy(p,1:params.nTotalTrials)=data.key_resp_whichcorr(params.nPractice+1:end);
        results.corrAns(p,1:params.nTotalTrials)=data.corrAns(params.nPractice+1:end);
        
        %% Replacement data variables
        fprintf('Participant %i: Recalculating unsaved variables for analysis \n\n',p);
        for t=params.nPractice+1:params.nTotalTrials+params.nPractice %skipping practice trials
            for f= 2:length(data.poly1x{t})
                
                %% Hypothesis Selection on each frame
                % results.hypSelected{p,t-params.nPractice}(f) = str2double(data.track_Model_key{t}(find(double(data.track_Model_Frame{t})/2<=f, 1, 'last' ))); %trying to do it based on button press timing
                %based on track_sloth
                if all(data.track_sloth{t}(f,:)==params.noSelectionSloth)%track sloth is all equal and therefore no square is selected
                    results.hypSelected{p,t-params.nPractice}(f) = nan;
                else %a square is selected
                    results.hypSelected{p,t-params.nPractice}(f) = find(max(data.track_sloth{t}(f,:))==data.track_sloth{t}(f,:));
                end
                
                %% Distance travelled on by Mouse and each square each frame
                if data.mouseX{t}(f)>0.88 || data.mouseX{t}(f)<-0.88 || data.mouseY{t}(f)>0.48 || data.mouseY{t}(f)<-0.48 %Nan distance measurements close to screen borders - wall
                    data.dSquares{t}(f,1:params.nSquares) = nan; %mouse too close to screen wall to measure intended distance, so don't count these frames
                    data.dMouse{t}(f) = nan;
                else
                    if data.poly1y{t}(f)<1 && data.poly1y{t}(f)>-1 && data.poly1x{t}(f)<1 && data.poly1x{t}(f)>-1
                        %set basic distances traveled
                        data.dMouse{t}(f) = sqrt(abs((data.mouseX{t}(f+1)-data.mouseX{t}(f)))^2+abs((data.mouseY{t}(f+1)-data.mouseY{t}(f)))^2); %mouse counts from a frame ahead?
                        data.dSquares{t}(f,1) = sqrt(abs((data.poly1x{t}(f)-data.poly1x{t}(f-1)))^2+abs((data.poly1y{t}(f)-data.poly1y{t}(f-1)))^2);
                        data.dSquares{t}(f,2) = sqrt(abs((data.poly2x{t}(f)-data.poly2x{t}(f-1)))^2+abs((data.poly2y{t}(f)-data.poly2y{t}(f-1)))^2);
                        data.dSquares{t}(f,3) = sqrt(abs((data.poly3x{t}(f)-data.poly3x{t}(f-1)))^2+abs((data.poly3y{t}(f)-data.poly3y{t}(f-1)))^2);
                        data.dSquares{t}(f,4) = sqrt(abs((data.poly4x{t}(f)-data.poly4x{t}(f-1)))^2+abs((data.poly4y{t}(f)-data.poly4y{t}(f-1)))^2);
                        %deal with overestimates because of screen wrapping
                        if any(data.track_wrap{t}(f,1:params.nSquares)~= 0) %if any of the squares wrapped
                            wrappedSq=find(data.track_wrap{t}(f,:)>0);
                            if any(find(results.hypSelected{p,t-params.nPractice}(f)==wrappedSq))%if any of the wrapped squares was the selected square
                                data.dSquares{t}(f,wrappedSq)= data.dMouse{t}(f)/params.slothDistracter;
                                data.dSquares{t}(f,wrappedSq(wrappedSq==results.hypSelected{p,t-params.nPractice}(f)))=data.dMouse{t}(f)/params.slothSelected; %replaces the ones that are selected
                            else %none of the wrapped squares are the selected square
                                data.dSquares{t}(f,wrappedSq)= data.dMouse{t}(f)/params.slothDistracter;
                            end
                        end
                    else %if data is outside normal limits - due to minimisation
                        data.dSquares{t}(f,1:params.nSquares) = nan; %mouse too close to screen wall to measure intended distance, so don't count these frames
                        data.dMouse{t}(f) = nan;
                    end
                end
            end
            results.dMouse{p}{t-params.nPractice}=data.dMouse{t};
        end
       
        fprintf('Participant %i: Variables added successfully. Continuing analysis... \n\n',p);
        
        %small script for visual checking of a trial for adjustments (1:4 is square distance, 5 is mouse distance, 6 is selected square)
        %     d =    data.dSquares{t}
        %     d(:,5) = data.dMouse{t}'
        %     d(:,6) = results.hypSelected{end,t-params.nPractice}
        
%         %% Proportion of Trial with chosen square or me square selected
%         -MOVED BELOW FOR GOODTRIAL UPDATING      
%         righttomeanS = [];
%         wrongtomeanS = [];
%         righttomeanMe = [];
%         wrongtomeanMe = [];
%         for t = 1:params.nTotalTrials
%             if results.chosenSquare(p,t)>0 %if they didn't say no control
%                 results.percentChosenSelected(p,t) = length(find(results.hypSelected{p,t}==results.chosenSquare(p,t)))/length(find(~isnan(results.hypSelected{p,t})))*100; %percent of frames where a square was selected and it was the finally chosen square
%                 results.percentCorrectSelected(p,t) = length(find(results.hypSelected{p,t}==data.corrAns(t+params.nPractice)))/length(find(~isnan(results.hypSelected{p,t})))*100; %percent of frames where a square was selected and it was the correct square
%                 if length(find(t==goodtrials))==1 %only add to average if it's a good trial
%                     if results.accuracy(p,t) ==1 %if they were right
%                         righttomeanS(end+1)=results.percentChosenSelected(p,t);
%                         righttomeanMe(end+1)=results.percentCorrectSelected(p,t);
%                     else %if they were wrong
%                         wrongtomeanS(end+1)=results.percentChosenSelected(p,t);
%                         wrongtomeanMe(end+1)=results.percentCorrectSelected(p,t);
%                     end
%                 end
%             else
%                 results.percentChosenSelected(p,t) = 0;
%             end
%         end
%         results.avCorrectPercentChosenSelected(p)= mean(righttomeanS);
%         results.avIncorrectPercentChosenSelected(p) = mean(wrongtomeanS);
%         results.avCorrectPercentCorrectSelected(p)= mean(righttomeanMe);
%         results.avIncorrectPercentCorrectSelected(p) = mean(wrongtomeanMe);
        
        %% Speed derivatives (AQ and Acceleration Replication)
        for t = 1:params.nTotalTrials
            deltad = min(data.dSquares{params.nPractice+t}'); %selected square is always the lowest distance, this lists that value for every frame
            for f = 2:length(data.dSquares{params.nPractice+t})
                d2(f) = deltad(f)-deltad(f-1);% acceleration
            end
            for f = 3:length(data.dSquares{params.nPractice+t})
                d3(f) = d2(f)-d2(f-1); %jerk
            end
            results.speed(p,t) = nanmean(deltad);
            results.acceleration(p,t) = nanmean(d2);
            results.jerk(p,t) = nanmean(d3);
        end
        
       %% Percent of Trial Spent Moving - needs to be run with speed derivs for distance measure(frames spent moving)
        for t = 1:params.nTotalTrials
            results.trialTime(p,t) = params.trialDuration-results.totalLagDuration(p,t); % THIS SHOULD BE REDONE IF USED - some lags were longer than the actual trial!
            results.framesMoving(p,t)= sum((~isnan(results.hypSelected{p,t})) & (data.dSquares{t+params.nPractice}(:,1)>0)');
            timemoving = 0;
            timemoving2 = nan;
            for f=2:length(data.track_logrowT{t+params.nPractice})
                if ~isnan(results.hypSelected{p,t}(f)) & data.dSquares{t+params.nPractice}(f,1)>0
                    timemoving = timemoving + (data.track_logrowT{t+params.nPractice}(f)-data.track_logrowT{t+params.nPractice}(f-1)); %add duration of that moving frame to time moving variable - though this will include lags while moving...
                end
                frametime(f-1) = data.track_logrowT{t+params.nPractice}(f)-data.track_logrowT{t+params.nPractice}(f-1);%alternative way of calcing time moving
            end
            timemoving2 = results.framesMoving(p,t)*results.noLagFreq(t);% or mean(frametime)
            results.avFrameDuration(p,t) = mean(frametime);
            results.percentSpentMoving(p,t)=timemoving/15*100; %change between timemoving and timemoving2 depending on how you want to calc it - can be vastly different... timemoving will inflate because lags are counted in num but not denom if (results.trialTime(p,t)) is used, timemoving2 will inflate because lags are part of the average frame duration (not if you use nolagfreq)
            results.percentFramesSpentMoving(p,t)= results.framesMoving(p,t)/results.nFrames(p,t)*100;
            if results.percentFramesSpentMoving(p,t)<1 %less than 1% of trial spent moving
                results.badtrials(p,t) = 1;
            end
        end
        
        %Update badtrials after trial spent moving duration 
        goodtrials = find(results.badtrials(p,:)==0); %this variable is not saved, is just used for averaging across vars at the end of this loop
        results.countBad(p) = sum(results.badtrials(p,:));
        
                %% Proportion of Trial with chosen square or me square selected
        righttomeanS = [];
        wrongtomeanS = [];
        righttomeanMe = [];
        wrongtomeanMe = [];
        for t = 1:params.nTotalTrials
            if results.chosenSquare(p,t)>0 %if they didn't say no control
                results.percentChosenSelected(p,t) = length(find(results.hypSelected{p,t}==results.chosenSquare(p,t)))/length(find(~isnan(results.hypSelected{p,t})))*100; %percent of frames where a square was selected and it was the finally chosen square
                results.percentCorrectSelected(p,t) = length(find(results.hypSelected{p,t}==data.corrAns(t+params.nPractice)))/length(find(~isnan(results.hypSelected{p,t})))*100; %percent of frames where a square was selected and it was the correct square
                if length(find(t==goodtrials))==1 %only add to average if it's a good trial
                    if results.accuracy(p,t) ==1 %if they were right
                        righttomeanS(end+1)=results.percentChosenSelected(p,t);
                        righttomeanMe(end+1)=results.percentCorrectSelected(p,t);
                    else %if they were wrong
                        wrongtomeanS(end+1)=results.percentChosenSelected(p,t);
                        wrongtomeanMe(end+1)=results.percentCorrectSelected(p,t);
                    end
                end
            else
                results.percentChosenSelected(p,t) = 0;
            end
        end
        results.avCorrectPercentChosenSelected(p)= mean(righttomeanS);
        results.avIncorrectPercentChosenSelected(p) = mean(wrongtomeanS);
        results.avCorrectPercentCorrectSelected(p)= mean(righttomeanMe);
        results.avIncorrectPercentCorrectSelected(p) = mean(wrongtomeanMe);
        
        %% N Selected Square Switches per Trial
        results.nHypSwitches(p,t) = 0;
        results.hypSwitches{p,t}=[];
        results.triggerHypSwitches{p,t}=[];
        for t = 1:params.nTotalTrials
            if ~isempty(data.track_Model_key{t+params.nPractice})
                for   b = 1:length(data.track_Model_key{t+params.nPractice}) %go through selection button markers
                    if ~isnan(str2double(data.track_Model_key{t+params.nPractice}(b))) %if it's not a button off recording
                        if ~isempty(results.hypSwitches{p,t}) %if it's not the first selection recorded
                            if str2double(data.track_Model_key{t+params.nPractice}(b))~= results.hypSwitches{p,t}(1,end) %if this selection is not the same as the last selection ie. its truly a 'switch'
                                results.nHypSwitches(p,t) = results.nHypSwitches(p,t)+1; %increase counter for the trial
                                %record switch number and frame
                                results.hypSwitches{p,t}(1,end+1)=str2double(data.track_Model_key{t+params.nPractice}(b)); %which selection
                                results.hypSwitches{p,t}(2,end)=round(data.track_Model_Frame{t+params.nPractice}(b)/2)+1; %what frame? in 30hz presentation, matching with d variables
                                results.triggerHypSwitches{p,t}(end+1)=round(data.track_Model_Frame{t+params.nPractice}(b)/2)+1; %frames for ERPEs
                            end
                        else
                            results.hypSwitches{p,t}(1,end+1)=str2double(data.track_Model_key{t+params.nPractice}(b)); %which selection
                            results.hypSwitches{p,t}(2,end)=nan; %since it's the first press, it only really counts for comparison or counting which ones get pressed at all, so we nan the "time of switch"
                        end
                    end
                end
            end
        end
        
        %% Percent of Trial Spent in Water
        %Side = 1 then any x value>0 (right) is water, Side = 2 x<0 (left) is water
        
        for t = 1:params.nTotalTrials
            for f = 1:length(data.poly1x{t+params.nPractice})
                %pull out the x value for the selected object on each frame
                if results.hypSelected{p,t}(f)==1
                    selectedX=data.poly1x{t+params.nPractice}(f);
                elseif results.hypSelected{p,t}(f)==2
                    selectedX=data.poly2x{t+params.nPractice}(f);
                elseif results.hypSelected{p,t}(f)==3
                    selectedX=data.poly3x{t+params.nPractice}(f);
                elseif results.hypSelected{p,t}(f)==4
                    selectedX=data.poly4x{t+params.nPractice}(f);
                else
                    selectedX=nan; %if not selected, nan
                end
                
                if isnan(selectedX)
                    results.water{p,t}(f)= nan;
                else %if selectedX is not nan
                    if data.side(t+params.nPractice)==1 %right is water
                        results.water{p,t}(f)= selectedX >= 0; %logicald
                    else %left is water
                        results.water{p,t}(f)= selectedX <= 0;
                    end
                end
                results.xPosSelected{p,t}(f) = selectedX;
            end
            results.percWater(p,t) =  nansum(results.water{p,t})/sum(~isnan(results.water{p,t}))*100; %water as a proportion of all frames where an object was selcted
            %results.percWater(p,t) = nansum(results.water{p,t})/results.framesMoving(p,t)*100; %as a proportion of frames moving
        end
        %% Dominant Environment
        results.domWater(p,:) = results.percWater(p,:)>50; % mark trials where water was dominant with a 1
        results.avPercWater(p) = nanmean(results.percWater(p,goodtrials));
        results.countDomWater(p) = sum(results.domWater(p,goodtrials));
        results.countDomSand(p) = sum(~results.domWater(p,goodtrials));
        results.percDomWater(p) = results.countDomWater(p)/(results.countDomWater(p)+results.countDomSand(p));
        
        %% Environment Switch Count
        for t=1:params.nTotalTrials
            results.nTotalEnvSwitches(p,t)=0;
            results.envSwitches{p,t}=[];
            results.directionToWaterEnvSwitches{p,t}=[];
            results.nHypSwitchEnvSwitches(p,t)=0;
            results.hypSwitchEnvSwitch{p,t}=[];
            results.nQrEnvSwitches(p,t)=0;
            results.qrEnvSwitch{p,t}=[];
            results.nTraverseEnvSwitches(p,t)=0;
            results.traverseEnvSwitch{p,t}=[];
            curEnv = results.water{p,t}(min(find(~isnan(results.water{p,t})))); %start with first env so it doesn't count as a switch
            prevEnv = NaN;
            curHyp = NaN;
            prevHyp = NaN;
            for f = 1:length(data.poly1x{t+params.nPractice})
                if ~isnan(results.hypSelected{p,t}(f)) %only count numerical hyp - nans are no selected sq
                    prevHyp = curHyp;
                    curHyp = results.hypSelected{p,t}(f);
                    prevEnv = curEnv;
                    curEnv = results.water{p,t}(f);
                end
                if f>1
                    if isnan(results.water{p,t}(f)) & ~isnan(results.water{p,t}(f-1)) & ~isnan(curEnv) & ~isnan(curHyp)
                        prevEnv = curEnv; %stops weird things where you get a switch right before no selection
                        prevHyp = curHyp;
                    end
                end
                if (prevEnv ~= curEnv) & all(~isnan([curEnv, curHyp, prevEnv, prevHyp]))%its a switch of some variety
                    results.envSwitches{p,t}(end+1) = f;
                    results.nTotalEnvSwitches(p,t) = results.nTotalEnvSwitches(p,t)+1;
                    if curEnv == 1 %if the switch is to water from sand
                        results.directionToWaterEnvSwitches{p,t}(end+1)=1;
                    else
                        results.directionToWaterEnvSwitches{p,t}(end+1)=0;
                    end
                    if any(results.triggerHypSwitches{p,t}==f)%||any(results.triggerHypSwitches{p,t}==f-1)||any(results.triggerHypSwitches{p,t}==f+1)%||curHyp ~= prevHyp %if the switch is within 1 frame either side of a hyp switch
                        results.hypSwitchEnvSwitch{p,t}(end+1) = f;
                        results.nHypSwitchEnvSwitches(p,t) = length(results.hypSwitchEnvSwitch{p,t});
                        if results.nHypSwitches(p,t)<results.nHypSwitchEnvSwitches(p,t)
                            fprintf('WARNING: number of hypothesis switch environment switches exceeds number of hypothesis switches in general, p=%i t=%i \n',p,t)
                        end
                    end
                end
            end
            %Split by QR vs movement
            if ~isempty(data.track_Env_Frame{t+params.nPractice})
                for    s=1:length(results.envSwitches{p,t})
                    if ~ismember(results.envSwitches{p,t}(s),results.hypSwitchEnvSwitch{p,t})%~any(results.hypSwitchEnvSwitch{p,t}==results.envSwitches{p,t}(s)) %if not a selection switch
                        try
                            if ismember(results.envSwitches{p,t}(s),round(data.track_Env_Frame{t+params.nPractice}/2))%QR switch
                                results.qrEnvSwitch{p,t}(end+1) = results.envSwitches{p,t}(s);
                                results.nQrEnvSwitches(p,t) = results.nQrEnvSwitches(p,t)+1;
                            elseif (results.envSwitches{p,t}(s)==max(find(~isnan(results.water{p,t}(1:round(data.track_Env_Frame{t+params.nPractice}(s)/2))))))|(results.envSwitches{p,t}(s)==min(find(~isnan(results.water{p,t}(round(data.track_Env_Frame{t+params.nPractice}(s)/2):end))))+round(data.track_Env_Frame{t+params.nPractice}(s)/2)-1) %capture non-exact qr switches during nan periods - if the switch time equals the first or last nan of that
                                results.qrEnvSwitch{p,t}(end+1) = results.envSwitches{p,t}(s);
                                results.nQrEnvSwitches(p,t) = results.nQrEnvSwitches(p,t)+1;
                            end
                        catch %this will catch when no min index of real values after nans or no max before nans
                        end
                    end
                end
            end
            %movement switch (not hyp switch, not qr switch)
            results.traverseEnvSwitch{p,t} = setdiff(setdiff(results.envSwitches{p,t},results.hypSwitchEnvSwitch{p,t}),results.qrEnvSwitch{p,t}); %all leftover env switches that are not hyp switches are traverse switches
            results.nTraverseEnvSwitches(p,t) = length(results.traverseEnvSwitch{p,t});
            
            results.nonHypEnvSwitch{p,t}= union(results.qrEnvSwitch{p,t},results.traverseEnvSwitch{p,t}); %env switches that are not related to hyp switches for triggering later
            
        end
        
        %env switch by direction
        for t=1:48
results.toWaterEnvSwitchPercent(p,t)=sum(results.directionToWaterEnvSwitches{p,t})/length(results.directionToWaterEnvSwitches{p,t});
end
        %% Average 'Prediction Error'
        [results.PE{p}] = predictionerror(p, data, params, results);
        for t = 1:params.nTotalTrials
            results.avPE(p,t)=nanmean(results.PE{p}{t});
        end
        
        %% Averages removing bad trials
        results.avFrames(p) = nanmean(results.nFrames(p,goodtrials));
        results.modeFrames(p) = mode(results.nFrames(p,goodtrials));
        results.avTrialTime(p) =nanmean(results.trialTime(p,goodtrials));
        results.avLagCount(p) = nanmean(results.lagCount(p,goodtrials));
        results.ALLavLagCount(p) = nanmean(results.lagCount(p,:));
        results.avConfidence(p) = nanmean(results.confidence(p,goodtrials));
        results.avAccuracy(p)=sum(results.accuracy(p,goodtrials))/(params.nTotalTrials-results.countBad(p));
        results.avSoA(p) = nanmean(results.senseOfAgency(p,goodtrials));
        results.avEnvSwitches(p) = nanmean(results.nTotalEnvSwitches(p,goodtrials));
        results.avHypSwitches(p) = nanmean(results.nHypSwitches(p,goodtrials));
        results.avPercentSpentMoving(p)= nanmean(results.percentSpentMoving(p,goodtrials));
        results.avPercentFramesSpentMoving(p)= nanmean(results.percentFramesSpentMoving(p,goodtrials));
        
        %% SDT removing bad trials
        results.SDTcounts(p,1:4)=zeros(1,4); %order: True pos, false pos, true neg, false neg
        results.SDTcountsW(p,1:4)=zeros(1,4); %restricted to Water dom trials
        results.SDTcountsS(p,1:4)=zeros(1,4); %restricted to Sand dom trials
        nSigW = 0; %number of control trials dominant water
        nNoiseW = 0;
        nSigS = 0; %number of control trials dominant sand
        nNoiseS = 0;
        for t = 1:params.nTotalTrials
            if results.badtrials(p,t)==0 %if not a bad trial
                if results.accuracy(p,t)== 1
                    if results.judgedAgency(p,t) == 1 %true pos
                        results.SDTcounts(p,1)=results.SDTcounts(p,1)+1;
                        if results.domWater(p,t)==1
                            results.SDTcountsW(p,1)=results.SDTcountsW(p,1)+1;
                        else
                            results.SDTcountsS(p,1)=results.SDTcountsS(p,1)+1;
                        end
                    elseif results.noControl(p,t)==1 %true neg
                        results.SDTcounts(p,3)=results.SDTcounts(p,3)+1;
                        if results.domWater(p,t)==1
                            results.SDTcountsW(p,3)=results.SDTcountsW(p,3)+1;
                        else
                            results.SDTcountsS(p,3)=results.SDTcountsS(p,3)+1;
                        end
                    else
                        sprintf('SDT ERROR: Cannot place correct trial, check analysis code, participant: %i, trial: %i \n',p,t);
                    end
                elseif results.judgedAgency(p,t) == 1 %false alarms (in both control (picked wrong square) and no control trials)
                    results.SDTcounts(p,2)=results.SDTcounts(p,2)+1;
                    if results.domWater(p,t)==1
                        results.SDTcountsW(p,2)=results.SDTcountsW(p,2)+1;
                    else
                        results.SDTcountsS(p,2)=results.SDTcountsS(p,2)+1;
                    end
                elseif results.judgedAgency(p,t) == 0 %false negative
                    results.SDTcounts(p,4)=results.SDTcounts(p,4)+1;
                    if results.domWater(p,t)==1
                        results.SDTcountsW(p,4)=results.SDTcountsW(p,4)+1;
                    else
                        results.SDTcountsS(p,4)=results.SDTcountsS(p,4)+1;
                    end
                else
                    sprintf('SDT ERROR: Cannot place incorrect trial, check analysis code, participant: %i, trial: %i\n',p,t);
                end
                if results.noControl(p,t)==1
                    if results.domWater(p,t)==1
                        nNoiseW = nNoiseW+1;
                    else
                        nNoiseS = nNoiseS+1;
                    end
                else
                    if results.domWater(p,t)==1
                        nSigW = nSigW+1;
                    else
                        nSigS = nSigS+1;
                    end
                end
            end
        end
        
        %use function to calc dprime and criterion across all trials
        nSignal = length(find(~results.noControl(p,goodtrials)));
        nNoise = length(find(results.noControl(p,goodtrials)));
        pHit = results.SDTcounts(p,1)/(results.SDTcounts(p,1)+results.SDTcounts(p,4)); %hit/hit+miss
        pFA = results.SDTcounts(p,2)/(results.SDTcounts(p,2)+results.SDTcounts(p,3)); %fa/fa+correct rej
        [results.dprime(p),results.criterion(p)]= dprime(pHit,pFA,nSignal,nNoise); %(pHit,pFA,nTarget,nDistract)
        % Note that number of noise trials ~= number of fa+corr rej -
        % because false alarms are also incorrect agentive responses in a
        % control trial
        
        %calc dprime and C for water dom trials
        pHit = results.SDTcountsW(p,1)/(results.SDTcountsW(p,1)+results.SDTcountsW(p,4)); %hit/hit+miss
        pFA = results.SDTcountsW(p,2)/(results.SDTcountsW(p,2)+results.SDTcountsW(p,3)); %fa/fa+correct rej
        [results.dprimeWaterDom(p),results.criterionWaterDom(p)]= dprime(pHit,pFA,nSigW,nNoiseW);
        
        %calc dprime and C for sand dom trials
        pHit = results.SDTcountsS(p,1)/(results.SDTcountsS(p,1)+results.SDTcountsS(p,4)); %hit/hit+miss
        pFA = results.SDTcountsS(p,2)/(results.SDTcountsS(p,2)+results.SDTcountsS(p,3)); %fa/fa+correct rej
        [results.dprimeSandDom(p),results.criterionSandDom(p)]= dprime(pHit,pFA,nSigS,nNoiseS);
        
        %% Comparing Dominance B1 vs B3
        results.B1MeanPercentWater(p) = nanmean(results.percWater(p,~results.badtrials(p,1:16))');
        results.B1WaterDom(p) = results.B1MeanPercentWater(p)>50;
        results.B1DomAmount(p) = abs(results.B1MeanPercentWater(p)-50);
        results.B3MeanPercentWater(p) = nanmean(results.percWater(p,~results.badtrials(p,32:48))');
        results.B3WaterDom(p) = results.B3MeanPercentWater(p)>50;
        results.B3DomAmount(p) = abs(results.B3MeanPercentWater(p)-50);
        results.B1B3DomDiff(p) = results.B3DomAmount(p)-results.B1DomAmount(p);
        results.B1B3MeanPercWaterDiff(p) = results.B3MeanPercentWater(p)- results.B1MeanPercentWater(p);
        
        %% List Indices of Participants to keep
        if results.avAccuracy(p)>=.25
            if results.countBad(p)<=.5*params.nTotalTrials
                if max(results.histChosenSq(p,:))<=40
                    results.keep(end+1)=p;
                end
            end
        end
        
    catch error
        results.errorId{end+1} = results.workerId(p);
        results.errorIndex(end+1) = p;
        results.errorMessage{end+1}= error;
    end
end

%% Analyses that don't require looping through participants but can be done straight from results file
try
    %% Prediction Error Slope and ERPEs
    fprintf('--Performing Prediction Error Slope Analysis-- \n\n');
    %% Prediction Error Slope
    [results.slopebyAA, results.slopebyDomEnv, results.slopebyEnv] = slopeAnalysis(results.PE, results.accuracy, results.judgedAgency, results.domWater, results.water, results.badtrials, results.keep);
    
    %optional plots 
%     figure
%     boxplot(results.slopebyAA(results.keep,:));
%     xticklabels({'Correct Agency','Correct No Agency','Incorrect Agency','Incorrect No Agency'});
%     ylabel('gradient of line of best fit');
%     figure
%     boxplot(results.slopebyDomEnv(results.keep,:));
%     xticklabels({'Water Dominant','Sand Dominant'});
%     ylabel('gradient of line of best fit');
%     boxplot(results.slopebyEnv(results.keep,:));
%     xticklabels({'Water Frames','Sand Frames'});
%     ylabel('gradient of line of best fit');
    
    %% Environment Switch ERPE
    fprintf('--Performing Environment ERPE Analysis-- \n\n');
    [results.ES_leadin_erpe, results.ES_leadout_erpe,~,~,~] = eventRelatedPredictionError (results.envSwitches, results.PE, results.water, results.hypSelected, params, results.badtrials,14,14,results.directionToWaterEnvSwitches,[]);
    %subset split logic for only nonhypswitcherpe
    for p = 1:params.nParticipants
        for t = 1:params.nTotalTrials
            directionToWater{p,t} = results.directionToWaterEnvSwitches{p,t}(find(ismember(results.envSwitches{p,t},results.nonHypEnvSwitch{p,t}')));
        end
    end
    [results.ESnoHS_leadin_erpe, results.ESnoHS_leadout_erpe,results.envHyp_prehist,results.envHyp_posthist,results.reclassHypSwitchEnvSwitch] = eventRelatedPredictionError (results.nonHypEnvSwitch, results.PE, results.water, results.hypSelected, params, results.badtrials,14,14,directionToWater,results.triggerHypSwitches);
   %reclassify env switches if necessary
   if isfield(results,'reclassHypSwitchEnvSwitch')
       for p = 1:params.nParticipants
           for t = 1:params.nTotalTrials
               if ~isempty(results.reclassHypSwitchEnvSwitch{p,t})
                   results.hypSwitchEnvSwitch{p,t} = union(results.reclassHypSwitchEnvSwitch{p,t}, results.hypSwitchEnvSwitch{p,t});
                   results.nHypSwitchEnvSwitches(p,t) = length(results.hypSwitchEnvSwitch{p,t});
                   if results.nHypSwitches(p,t)<results.nHypSwitchEnvSwitches(p,t)
                       fprintf('WARNING: number of hypothesis switch environment switches exceeds number of hypothesis switches in general, p=%i t=%i \n',p,t)
                   end
                   results.traverseEnvSwitch{p,t} = setdiff(results.traverseEnvSwitch{p,t},results.reclassHypSwitchEnvSwitch{p,t});
                   results.nTraverseEnvSwitches(p,t) = length(results.traverseEnvSwitch{p,t});
                   results.qrEnvSwitches{p,t} = setdiff(results.qrEnvSwitch{p,t},results.reclassHypSwitchEnvSwitch{p,t});
                   results.nQrEnvSwitches(p,t) = length(results.qrEnvSwitch{p,t});
                   results.nonHypEnvSwitch{p,t}= union(results.qrEnvSwitch{p,t},results.traverseEnvSwitch{p,t});
               end
           end
       end
   end
   %resubset split logic for only nonhypswitcherpe
    for p = 1:params.nParticipants
        for t = 1:params.nTotalTrials
            directionToWater2{p,t} = results.directionToWaterEnvSwitches{p,t}(find(ismember(results.envSwitches{p,t},results.nonHypEnvSwitch{p,t}')));
        end
    end
    %event related speed ESnoHS
    [results.ESnoHS_leadin_erSpeed, results.ESnoHS_leadout_erSpeed,~,~] = eventRelatedPredictionError (results.nonHypEnvSwitch, results.dMouse, results.water, results.hypSelected, params, results.badtrials,14,14,directionToWater2,[]);
    
    %optional plot
    %plot
%     figure;
%     subplot(2,2,1);
%     plot(nanmean(results.ESnoHS_leadin_erpe{2}(results.keep,:)));
%     hold on;
%     plot(nanmean(results.ESnoHS_leadin_erpe{1}(results.keep,:)));
%    % plot_areaerrorbar(nanmean(results.ESnoHS_leadin_erpe(results.keep,:)));
%     title('Non-Hypothesis Environment Switch ERPEs')
%      ylim([0 .13]); 
%      xlim([1 15]);
%     subplot(2,2,2);
%     plot(nanmean(results.ESnoHS_leadout_erpe{2}(results.keep,:)));
%     hold on;
%     plot(nanmean(results.ESnoHS_leadout_erpe{1}(results.keep,:)));
%    % plot_areaerrorbar(nanmean(results.ESnoHS_leadout_erpe(results.keep,:)));
%      ylim([0 .13]);
%      xlim([1 15]);
%     legend('Into Water Switches','Into Sand Switches') %plotted water first (2=1 in logic) so it was the blue default
%     subplot(2,2,3);
%     plot(mean(results.envHyp_prehist(results.keep,:)));
%          ylim([0 6]);
%      xlim([1 15]);
%    title('Mean Count Hyp Switch Per Frame in Epoch')
%     subplot(2,2,4);
%     plot(mean(results.envHyp_posthist(results.keep,:)));
%         ylim([0 6]);
%      xlim([1 15]);

    %% Null Distribution from Faux Hyp Switches
    %% this takes too long with the whole dataset, so use goodparticipantsNull.m to do this bit for only good participants after the other results are generated
%     %Trigger Generation
%     [results.nullTriggers, results.nullBoundaryID, results.fewerFauxCount] = nullCrossingGeneration (results.nonHypEnvSwitch, results.xPosSelected, results.triggerHypSwitches)
%     
%     %null event related prediction error
%     [results.NULL_ESnoHS_leadin_erpe, results.NULL_ESnoHS_leadout_erpe,~,~,~] = eventRelatedPredictionError (results.nullTriggers, results.PE, results.water, results.hypSelected, params, results.badtrials,14,14,directionToWater,results.triggerHypSwitches);

    %% Hypothesis Switch ERPE
     fprintf('--Performing Hypothesis ERPE Analysis-- \n\n');
    [results.HS_leadin_erpe, results.HS_leadout_erpe,~,~] = eventRelatedPredictionError (results.triggerHypSwitches, results.PE, results.water, results.hypSelected, params, results.badtrials,29,14,[],[]);
    %optional plot
    figure;
    subplot(1,2,1);
    plot(nanmean(results.HS_leadin_erpe{1}(results.keep,:)));
   % plot_areaerrorbar(nanmean(results.HS_leadin_erpe(results.keep,:)));
     ylim([0.02 0.055]); 
    xlim([1 30]);
    title('Hypothesis Switch ERPEs')
    hold on;
    subplot(1,2,2);
    plot(nanmean(results.HS_leadout_erpe{1}(results.keep,:)));
   % plot_areaerrorbar(nanmean(results.HS_leadout_erpe(results.keep,:)));
     ylim([0.02 .055]);
    xlim([1 15]);
    
catch error
    results.errorId{end+1}= results.workerId(p);
    results.errorMessage{end+1}= error;
end

%% Remove Errored participants from Keep list
results.keep = setdiff(results.keep,results.errorIndex);
results.keepId = results.workerId(results.keep);

%% Save!
if savedata == 1
    disp('---SAVING---');
    folder='\Results\';
    if dataset=='all'
        resultname = sprintf('results%s.mat',datestr(now,'ddmmyy'));
        paramsname = sprintf('params%s.mat',datestr(now,'ddmmyy'));
    else
        resultname = sprintf('results_short%i_%s.mat',dataset,datestr(now,'ddmmyy'));
        paramsname = sprintf('params_short%i_%s.mat',dataset,datestr(now,'ddmmyy'));
    end
    save([pwd folder resultname],'results');
    save([pwd folder paramsname],'params');
end

disp('---ANALYSIS COMPLETE ---');

