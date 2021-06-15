%% Create Null Distribution based on Faux Boundary Crossing Data

%Created by Kelsey Perrykkad on 06/05/21 for use with BeachAnalysis.m
%project in collaboration with Jakob Hohwy and Jonno Robinson

%primary purpose is to create a trigger matrix for random non-Environment
%Switch location crosses which will have the same structure as
%results.nonHypEnvSwitches



function [nullTriggers, nullBoundaryID, fewerFauxCount] = nullCrossingGeneration (nonHypEnvSwitch, xPosSelected, triggerHypSwitches)

%% xpos true boundary definitions
leftedge = -0.88;
rightedge = 0.88;
center = 0;
trueBound = [leftedge rightedge center];

%% set minimum distance from boundary
%minDistanceFromBoundary recommendation = .5*total epoch frame length*mean speed
%of kept participants in dist/frame [as long as this does not discount >50%
%available screen space]
%calc for n=89 on 06/05/21
%mean speed of kept participants = 0.0274
%epoch length lead in = 15
%based on above, this would put minDist at .4110
%in one half of the screen, total available distance = 0.88
%which would only leave 0.058 x pos range for faux boundary
%new plan. keep boundary selection at least 5 frames average speed distance on average
%away from true boundaries = 0.137, and double check each selected
%is >5 frames away from real boundary crossing (one timebin in previous ERPE exps)

bufferDistance = 0.137; %spatial buffer distance from true boundaries
availableRange = [leftedge+bufferDistance center-bufferDistance center+bufferDistance rightedge-bufferDistance];
rng('shuffle'); %shuffle random number generator


%% Fill out nullTriggers
nullTriggers = cell(size(nonHypEnvSwitch));
nullBoundaryID = cell(size(nonHypEnvSwitch));
fewerFauxCount = zeros(size(nonHypEnvSwitch));

for p = 1:size(nonHypEnvSwitch,1)
    fprintf('Participant %i Faux Trigger Generation Started \n', p);
    for t = 1:size(nonHypEnvSwitch,2)
        for s = 1:length(nonHypEnvSwitch{p,t}) %match the number of switches in original dataset for each trial
            passCheck = 0;
            tic
            while toc<2 & passCheck ==0  %have it time out after 2s of trying (in cases like p=1 t=5 where the only possible crossing is the real one)
                %% generate and record boundary ID
                screenHalf = randi(2);
                randLocPerc = rand(1);
                if screenHalf == 1;
                    availDistance = availableRange(2)-availableRange(1);
                    genBound = availableRange(1)+(randLocPerc*availDistance);
                else
                    availDistance = availableRange(4)-availableRange(3);
                    genBound = availableRange(3)+(randLocPerc*availDistance);
                end
                
                %% count nullBoundary Crossings based on script from BeachAnalysis for env crossings
                curEnv = xPosSelected{p,t}(min(find(~isnan(xPosSelected{p,t}))))>genBound; %start with first env so it doesn't count as a switch (right of boundary is 1, left is 0)
                prevEnv = NaN;
                candidateSwitch = [];
                for f = 1:length(xPosSelected{p,t})
                    if ~isnan(xPosSelected{p,t}(f)) %only count numerical hyp - nans are no selected sq
                        prevEnv = curEnv; %shuffle it all down one position
                        curEnv = xPosSelected{p,t}(f)>genBound;
                    end
                    if f>1
                        if isnan(xPosSelected{p,t}(f)) & ~isnan(xPosSelected{p,t}(f-1)) & ~isnan(curEnv)
                            prevEnv = curEnv; %stops weird things where you get a switch right before no selection
                        end
                    end
                    if (prevEnv ~= curEnv) & all(~isnan([curEnv, prevEnv]))%its a switch of some variety
                        candidateSwitch(end+1) = f;
                    end
                end
                if length(candidateSwitch)>1 %catch if none available and repeat process
                    whichSwitch = randi(length(candidateSwitch));
                    if ~(any(candidateSwitch(whichSwitch)==triggerHypSwitches{p,t})) && ~(any(candidateSwitch(whichSwitch)==nonHypEnvSwitch{p,t})) && ~(any(abs(nonHypEnvSwitch{p,t}-candidateSwitch(whichSwitch))<5)) && ~any(candidateSwitch(whichSwitch)==nullTriggers{p,t}) %acceptable triggers must not match existing faux trigger locations, hyp switch locations or real boundary crossings and be 5f away from real boundary crossings
                        passCheck=1;
                    end
                end
            end
            %% Plot (turn off!)
            %             plot(xPosSelected{p,t})
            %             hold on
            %             yline(genBound)
            %             plot(candidateSwitch(whichSwitch),genBound,'r*')
            
            %% select from available nullBoundary Crossings and record trigger time
            if passCheck ==1 % only record the switch if it passes criteria
                nullBoundaryID{p,t}(end+1)= genBound; %records generated and used bound location
                nullTriggers{p,t}(end+1) = candidateSwitch(whichSwitch); %records trigger location
            end
        end
       % fprintf('Trial %i Faux Trigger Generation Completed \n', t);
        if ~(length(nullTriggers{p,t}) == length(nonHypEnvSwitch{p,t}))
            fprintf('WARNING: Trial %i has mismatched numbers of Triggers. Real: %d Faux: %d \n', t, length(nonHypEnvSwitch{p,t}),length(nullTriggers{p,t}));
            fewerFauxCount(p,t) = length(nonHypEnvSwitch{p,t})-length(nullTriggers{p,t});
        end
    end
    fprintf('Participant %i of %i Faux Trigger Generation Completed \n', p, size(nonHypEnvSwitch,1));
end
