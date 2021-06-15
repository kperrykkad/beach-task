%% Some graphs from pilot for determining thresholds for participant elimination

%% Accuracy by Average Frame Count
%threshold at 255f average minimum
%threshold at .25 accuracy 
%chuck is a vector of indexes for a particular condition for removal
%keep is the successful participants

keep = find(results.countBad<=24);
accuracychuck = find(results.avAccuracy<.25);
warnschuck = find(results.nWarns>5);

figure;
hold on
scatter(results.avFrames,results.avAccuracy);
scatter(results.avFrames(keep),results.avAccuracy(keep),'filled','MarkerFaceColor',[0 .7 .7]);
scatter(results.avFrames(accuracychuck),results.avAccuracy(accuracychuck),'filled','MarkerFaceColor',[.7 0 .1]);
scatter(results.avFrames(accuracychuck),results.avAccuracy(accuracychuck),'filled','MarkerFaceColor',[.8 0 .05]);