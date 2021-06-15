function [lagLocation, noLagFreq, totalLagDuration] = frameDiagnosis(trialdata, maxFrames, trialDuration)
%Identify hemorrhage lags (dropping MANY frames at once), output location
%to use in future analysis, take the average slope of frame time for each
%trial to use for bad trial definition (along with number of lags)

%trial data is the logrow times for that trial (data.track_logrowT)
%trial duration is in seconds
%max frames is logrow count expected based on programming (every other frame for more consistency)

%% Frame to frame time
expectDuration = trialDuration/maxFrames;
lagLocation = [];
lagDuration = [];
unacceptableDrops = 5; %how many dropped frames is too many dropped frames/ should count as a lag

for f = 2:length(trialdata)
    duration(f-1) = trialdata(f)-trialdata(f-1);
    if duration(f-1)>unacceptableDrops*expectDuration %if more than five frames were dropped in a row
        lagLocation(end+1)= f-1;
        lagDuration(end+1)=duration(f-1);
    end
end

goodDurations = find(duration<unacceptableDrops*expectDuration);
noLagFreq = mean(duration(goodDurations));
noLagVar = std(duration(goodDurations));
totalLagDuration = sum(lagDuration);