%% creating the null events FOR ONLY GOOD PARTICIPANTS etc using just a results file

%load results file
 load('G:\My Drive\Experiments\JoA_squarestask\BeachTask\Data and Analysis\Analysis\Results\results050521_withNull.mat')
 load('G:\My Drive\Experiments\JoA_squarestask\BeachTask\Data and Analysis\Analysis\Results\params050521.mat')
 load('G:\My Drive\Experiments\JoA_squarestask\BeachTask\Data and Analysis\Analysis\Results\keeplist.mat')
 
 %generate a set of null triggers for only good participants
 %[nullTriggers, nullBoundaryID, fewerFauxCount] = nullCrossingGeneration (results.nonHypEnvSwitch(results.keep,:), results.xPosSelected(results.keep,:), results.triggerHypSwitches(results.keep,:));
        %load('G:\My Drive\Experiments\JoA_squarestask\BeachTask\Data and Analysis\Analysis\Results\null2.mat')
 %redefine keep participants to match the survey kept participants. expect
 %n = 84
        
%change number of participants in parameters to only good participants
 results.keep = keeplist; %remove the extra participants from the survey analysis
 params.nParticipants=length(results.keep); %this line was mistakenly above 15 until 1/06/21
 
 %do faux erpe
[results.FAUX_ESnoHS_leadin_erpe, results.FAUX_ESnoHS_leadout_erpe, ~,~,~] = eventRelatedPredictionError (results.nullTriggers(results.keep(:),:), results.PE(results.keep), results.water(results.keep,:), results.hypSelected(results.keep,:), params, results.badtrials(results.keep,:),14,14,[],[]);
 %do faux event related speed
 [results.FAUX_ESnoHS_leadin_erSPEED, results.FAUX_ESnoHS_leadout_erSPEED, ~,~,~] = eventRelatedPredictionError (results.nullTriggers(results.keep(:),:), results.dMouse(results.keep), results.water(results.keep,:), results.hypSelected(results.keep,:), params, results.badtrials(results.keep,:),14,14,[],[]);
%non direction split real erpe
[results.ESnoHS_noDirection_leadin_erpe, results.ESnoHS_noDirection_leadout_erpe,~,~,~] = eventRelatedPredictionError (results.nonHypEnvSwitch(results.keep(:),:), results.PE(results.keep), results.water(results.keep,:), results.hypSelected(results.keep,:), params, results.badtrials(results.keep,:),14,14,[],[]);
 %non direction split real speed 
 [results.ESnoHS_noDirection_leadin_erSpeed, results.ESnoHS_noDirection_leadout_erSpeed,~,~,~] = eventRelatedPredictionError (results.nonHypEnvSwitch(results.keep(:),:), results.dMouse(results.keep), results.water(results.keep,:), results.hypSelected(results.keep,:), params, results.badtrials(results.keep,:),14,14,[],[]);

 
 %create set for permutation testing
 realfauxPerm = nan(length(results.keep),30,2);
 for p = 1:length(results.keep)
     for c = 1:2
         if c == 1 %real
             for tpre = 1:15
                 realfauxPerm(p,tpre,c) = results.ESnoHS_noDirection_leadin_erpe{1}(p,tpre);
             end
             for tpost = 1:15
                 realfauxPerm(p,tpost+15,c) = results.ESnoHS_noDirection_leadout_erpe{1}(p,tpost);
             end
         elseif c == 2 %fake
             for tpre = 1:15
                 realfauxPerm(p,tpre,c) = results.FAUX_ESnoHS_leadin_erpe{1}(p,tpre);
             end
             for tpost = 1:15
                 realfauxPerm(p,tpost+15,c) = results.FAUX_ESnoHS_leadout_erpe{1}(p,tpost); %this line said leadin until 1/06/21
             end
         end
     end
 end
 
 %perform permutation testing
 permstat = tperm(realfauxPerm,10000);

 %save
 %save('results050521_withNull_withGoodErpe_n84.mat','results') %WARNING: careful
 %doing this - results.keep needs to match results.keepId, which is changed
 %above -- did this manually, need to fix if this is needed in results file
 %going forward
 save('results050521_permStat_n84.mat','permstat')
 
 %% plots
 
 %histogram boundary position null
 nullBoundaryID = results.nullBoundaryID(results.keep(:),:); %only the good participants
 hist([nullBoundaryID{:}],55);
 
 % t-values
plot(permstat.real_t)
hold on
t1=ones(1,30)*permstat.thresh(1,2)
t2=ones(1,30)*permstat.thresh(2,2)
t3=ones(1,30)*permstat.thresh(3,2)
plot(t1)
plot(t2)
plot(t3)
 
 %faux vs real erpe
figure
subplot(1,2,1);
plot(nanmean(results.FAUX_ESnoHS_leadin_erpe{1})); hold on;
plot(nanmean(results.ESnoHS_noDirection_leadin_erpe{1}))
ylim([0 .2]);
subplot(1,2,2);
plot(nanmean(results.FAUX_ESnoHS_leadout_erpe{1})); hold on;
plot(nanmean(results.ESnoHS_noDirection_leadout_erpe{1}))
ylim([0 .2]);
legend('fake','real')

%with sem shading
figure
subplot(1,2,1);
semshade(results.FAUX_ESnoHS_leadin_erpe{1},.25,'r'); hold on;
semshade(results.ESnoHS_noDirection_leadin_erpe{1},.25,'g');
ylim([0 .2]);
subplot(1,2,2);
semshade(results.FAUX_ESnoHS_leadout_erpe{1},.25,'r'); hold on;
semshade(results.ESnoHS_noDirection_leadout_erpe{1},.25,'g');
ylim([0 .2]);

%with sem shading and sig periods marked for paper
figure
semshade([results.FAUX_ESnoHS_leadin_erpe{1} results.FAUX_ESnoHS_leadout_erpe{1}],.25,'r'); hold on;
semshade([results.ESnoHS_noDirection_leadin_erpe{1} results.ESnoHS_noDirection_leadout_erpe{1}],.25,'b');
xlim([1 30])
plot(10:23,0.025*ones(size(10:18)),'LineWidth', 4)
ylim([0 .14])
legend('SEM','Control: faux boundaries','SEM','True: boundary crossings', 'p<0.05')
xline(15:16)

%speed vs PE faux erpe
figure
subplot(1,2,1);
plot(nanmean(results.FAUX_ESnoHS_leadin_erSPEED{1})/params.slothSelected); hold on;
plot(nanmean(results.FAUX_ESnoHS_leadin_erpe{1}))
ylim([0 .2]);
title('Faux Boundaries Only')
subplot(1,2,2);
plot(nanmean(results.FAUX_ESnoHS_leadout_erSPEED{1})/params.slothSelected); hold on;
plot(nanmean(results.FAUX_ESnoHS_leadout_erpe{1}))
ylim([0 .2]);
legend('speed','prediction error')

%speed vs PE real erpe
figure
subplot(1,2,1);
plot(nanmean(results.ESnoHS_noDirection_leadin_erSpeed{1})/params.slothSelected); hold on;
plot(nanmean(results.ESnoHS_noDirection_leadin_erpe{1}))
%ylim([0 .2]);
title('Real Boundaries Only')
subplot(1,2,2);
plot(nanmean(results.ESnoHS_noDirection_leadout_erSpeed{1})/params.slothSelected); hold on;
plot(nanmean(results.ESnoHS_noDirection_leadout_erpe{1}))
%ylim([0 .2]);
legend('speed','prediction error')