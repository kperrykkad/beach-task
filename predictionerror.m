function [PE] = predictionerror(p, data, params, results)
%% Calculate Prediction Error from XY and reconstructed distance vars - no Angle info saved

% This function was created by Kelsey Perrykkad on 20/11/2020 for the
% beach task. It will output a cell for each trial containing a vector of
% the prediction error on each frame.

% This script was entirely rewritten since VarVol as the online
% environment meant that angle and distance data on each frame was not
% captured in the dataset automatically. To avoid complex quadrant based
% trigonometric reconstructions of angles, this function only uses
% previous and current x and y for the selected object, the previous and
% current mouse positions and the reconstructed distance variable to
% measure the prediction error.

% Derivation notes for this algorithm are in Kelsey's onenote page called
% "Prediction Error" in the section "Experiment 4 - Beach Task" on 20 Nov
% 2020.

moveThresh = 0; %how much does the mouse have to move in x and y to count as having moved?

for t = 1:params.nTotalTrials %this is 'results' trial counter - add nPractice to get data trial number
    PE{t}(1:results.nFrames(p,t))=0;
    for f = 2:results.nFrames(p,t)
        selectedObj = results.hypSelected{p,t}(f);
        dMouse = data.dMouse{t+params.nPractice}(f); %how far did the mouse move?
        mx1 = data.mouseX{t+params.nPractice}(f); %previous
        mx2 = data.mouseX{t+params.nPractice}(f+1); %current (mouse is one ahead of poly)
        my1 = data.mouseY{t+params.nPractice}(f);
        my2 = data.mouseY{t+params.nPractice}(f+1);
        %current x and y FOR THE SELECTED OBJECT and previous x and y FOR THE SELECTED OBJECT (don't want hypothesis switches to impact on P.E.)
        if ~isnan(selectedObj) && dMouse>moveThresh %if either no object is selected or the mouse doesn't move, then prediction error is 0 (there is no prediction error if you don't move)
            if selectedObj == 1
                x2 = data.poly1x{t+params.nPractice}(f); %current
                y2 = data.poly1y{t+params.nPractice}(f);
                x1 = data.poly1x{t+params.nPractice}(f-1); %previous
                y1 = data.poly1y{t+params.nPractice}(f-1);
            elseif selectedObj == 2
                x2 = data.poly2x{t+params.nPractice}(f); %current
                y2 = data.poly2y{t+params.nPractice}(f);
                x1 = data.poly2x{t+params.nPractice}(f-1); %previous
                y1 = data.poly2y{t+params.nPractice}(f-1);
            elseif selectedObj == 3
                x2 = data.poly3x{t+params.nPractice}(f); %current
                y2 = data.poly3y{t+params.nPractice}(f);
                x1 = data.poly3x{t+params.nPractice}(f-1); %previous
                y1 = data.poly3y{t+params.nPractice}(f-1);
            elseif selectedObj == 4
                x2 = data.poly4x{t+params.nPractice}(f); %current
                y2 = data.poly4y{t+params.nPractice}(f);
                x1 = data.poly4x{t+params.nPractice}(f-1); %previous
                y1 = data.poly4y{t+params.nPractice}(f-1);
            end
            dSelectSquare = data.dSquares{t+params.nPractice}(f,selectedObj);
            multiplier = (dMouse+(dSelectSquare-dMouse))/dMouse; %how much should mouse movement be multiplied for expected value estimation
            expectedX = x1+((mx2-mx1)*multiplier);
            expectedY = y1+((my2-my1)*multiplier);
            PE{t}(f) = sqrt((x2-expectedX)^2+(y2-expectedY)^2);
            
            %if expected X or expected Y ends up off
            %screen and the object would have wrapped, nan PE -
            %expectations for wrapping would presumably be imprecise
            if abs(x2-x1)>=1.5*.75 %if the square wrapped in x (1.5 is distance across movable area)
                PE{t}(f)= nan;
            end
            if abs(y2-y1)>=0.7*.75 %if the square wrapped in y (.7 is distance across movable area)
                PE{t}(f)= nan;
            end
        elseif isnan(selectedObj)
            PE{t}(f)= nan;
        end
        %nan overestimates because of q and r - essentially non-estimatable
        %PE, not just 0
        if any((data.track_Env_Frame{t+params.nPractice}./2)+1==f)
            PE{t}(f)= nan;
        end
    end
end

