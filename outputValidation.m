% output results for comparison of hypothesis selection methods across
% beach and varvol 1 experiments

%created by Kelsey Perrykkad on 12/02/21

%% Load data
beach = load('G:\My Drive\Experiments\JoA_squarestask\BeachTask\Data and Analysis\Analysis\Results\results110221.mat');
varvol1 = load('G:\My Drive\Experiments\JoA_squarestask\VarVolFactorial\Results\FINAL Results_20200309T180926\Results_20200909T163050\Results_20200909T163050.mat');
load('G:\My Drive\Experiments\JoA_squarestask\BeachTask\Data and Analysis\Analysis\Results\beachidlist84.mat');
id = table2array(beachidlist1);

% %% set up table
% dataheaders = {'id','study','CorrPercOnSquare','IncorrPercOnCorr', 'IncorrPercOnChosen'};
% sz = [length(id)+length(varvol1.participantList) length(dataheaders)];
% varTypes = {'double','double','double','double','double'};
% hypValidation = table('Size',sz,'VariableTypes',varTypes,'VariableNames',dataheaders);
%
% %% Build Table
% r=1;
% for k = 1:length(id)
%     p = beach.results.keep(find(id(k)==beach.results.keepId));
%     hypValidation.id(r) = p;
%     hypValidation.study(r)=2;
%     hypValidation.CorrPercOnSquare(r)= beach.results.avCorrectPercentChosenSelected(p);
%     hypValidation.IncorrPercOnCorr(r)= beach.results.avIncorrectPercentCorrectSelected(p);
%     hypValidation.IncorrPercOnChosen(r)=beach.results.avIncorrectPercentChosenSelected(p);
%     r=r+1;
% end
%
% for p = 1:length(varvol1.participantList)
%     hypValidation.id(r) = p;
%     hypValidation.study(r)=1;
%     hypValidation.CorrPercOnSquare(r)= varvol1.behaviour.avlookingtimeRacc(p,1)/15*100;
%     hypValidation.IncorrPercOnCorr(r)= varvol1.behaviour.avlookingtimeCacc(p,2)/15*100;
%     hypValidation.IncorrPercOnChosen(r)= varvol1.behaviour.avlookingtimeRacc(p,2)/15*100;
%     r=r+1;
% end
%
% %Save
% writetable(hypValidation,'G:\My Drive\Experiments\JoA_squarestask\VarVolFactorial\Results\HypValidation.csv');

%% set up table
dataheaders = {'id','study','Acc','Corr_Chosen', 'PercentDwell'};
sz = [(length(id)+length(varvol1.participantList))*4 length(dataheaders)];
varTypes = {'double','double','double','double','double'};
hypValidation = table('Size',sz,'VariableTypes',varTypes,'VariableNames',dataheaders);

%% Build Table
r=1;
keeplist=[];
for k = 1:length(id)
    p = beach.results.keep(find(id(k)==beach.results.keepId));
    keeplist(end+1)=p;
    for acc=1:2
        for square=1:2
            hypValidation.id(r) = p;
            hypValidation.study(r)=2;
            hypValidation.Acc(r)=acc-1;
            hypValidation.Corr_Chosen(r)=square;
            if acc == 2 %correct trials
                hypValidation.PercentDwell(r)= beach.results.avCorrectPercentChosenSelected(p);
                r=r+1;
            else %incorrect trials
                if square == 1
                    hypValidation.PercentDwell(r)= beach.results.avIncorrectPercentCorrectSelected(p);
                    r=r+1;
                else
                    hypValidation.PercentDwell(r)=beach.results.avIncorrectPercentChosenSelected(p);
                    r=r+1;
                end
            end
        end
    end
end

for p = 1:length(varvol1.participantList)
    for acc=1:2
        for square=1:2
            hypValidation.id(r) = p;
            hypValidation.study(r)=1;
            hypValidation.Acc(r)=acc-1;
            hypValidation.Corr_Chosen(r)=square;
            if acc == 2 %correct trials
                hypValidation.PercentDwell(r)= varvol1.behaviour.avlookingtimeRacc(p,1)/15*100;
                r=r+1;
            else %incorrect trials
                if square == 1
                    hypValidation.PercentDwell(r)=varvol1.behaviour.avlookingtimeCacc(p,2)/15*100;
                    r=r+1;
                else
                    hypValidation.PercentDwell(r)=varvol1.behaviour.avlookingtimeRacc(p,2)/15*100;
                    r=r+1;
                end
            end
        end
    end
end

%Save
%writetable(hypValidation);

%plot
figure
subplot(2,3,1)
boxplot(beach.results.avCorrectPercentChosenSelected(keeplist));
ylim([0 100])
title('BEACH - Correct Trials')
ylabel('% Trial Spent on Correct/Chosen Square')
subplot(2,3,2)
boxplot(beach.results.avIncorrectPercentCorrectSelected(keeplist));
ylim([0 100])
title('BEACH - Incorrect Trials')
ylabel('% Trial Spent on Correct Square')
subplot(2,3,3)
boxplot(beach.results.avIncorrectPercentChosenSelected(keeplist));
ylim([0 100])
title('BEACH - Incorrect Trials')
ylabel('% Trial Spent on Chosen Square')
subplot(2,3,4)
boxplot(varvol1.behaviour.avlookingtimeRacc(:,1)/15*100);
ylim([0 100])
title('Perrykkad et al. 2021 - Correct Trials')
ylabel('% Trial Spent on Correct/Chosen Square')
subplot(2,3,5)
boxplot(varvol1.behaviour.avlookingtimeCacc(:,2)/15*100);
ylim([0 100])
title('Perrykkad et al. 2021 - Incorrect Trials')
ylabel('% Trial Spent on Correct Square')
subplot(2,3,6)
boxplot(varvol1.behaviour.avlookingtimeRacc(:,2)/15*100);
ylim([0 100])
title('Perrykkad et al. 2021 - Incorrect Trials')
ylabel('% Trial Spent on Chosen Square')