# Concatonate datasets for Beach Task as downloaded from qualtrics and Cloudresearch and Prolific and cross check keep id list from matlab output
# Written by Kelsey Perrykkad Dec 2020, adapted from script written June 2020 for use with "Self-Processing and Personality Traits" Qualtrics study

# load libraries
library(tidyverse)
#library(jmv)
#library(gamlj)
#library(ggplot2)
#library(ggthemr)
library(qualtRics)

#file locations
q1loc <- "G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Data Storage/FinalSet/Qualtrics/Causal Inference and Cognition - Beach Task - PART 1_November 30, 2020_14.53.csv"
q2loc <- "G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Data Storage/FinalSet/Qualtrics/Causal Inference and Cognition - Beach Task - PART 2_November 30, 2020_14.54.csv"
paymentinfoloc1 <- "G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Data Storage/MTurkDataset/Was it you A Cognitive Task [Was it you A Cognitive Task(~ 35 minutes)] (259360).csv"
paymentinfoloc2 <- "G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Data Storage/ProlificDataset/prolific_export_5facb810868a7d2edf03efd9.csv"


#find and load in datafiles
q1 <- read_survey(q1loc)
q2 <- read_survey(q2loc)
pay1 <- read.csv(paymentinfoloc1)
pay2 <- read.csv(paymentinfoloc2)

#rename for merging
pay1 <- rename(pay1, workerId=AmazonIdentifier)
pay2 <- rename(pay2, workerId=participant_id)

#join payment datasets
pay1 <- filter(pay1, ApprovalStatus=="Approved")
pay2 <- filter(pay2, status=="APPROVED")
paidId <- as.factor(c(as.character(pay1$workerId),as.character(pay2$workerId)))


#replace missing worker Ids in Q2
workerIdq2.missing.replaceattempted <- which(is.na(q2$workerId))
for (row.index in which(is.na(q2$workerId))){
  if(any(which(q1$IPAddress==q2$IPAddress[row.index]))){
    q2$workerId[row.index] <- q1$workerId[which(q1$IPAddress==q2$IPAddress[row.index])]
  }
}
workerIdq2.missing.notfound <- which(is.na(q2$workerId))

#join across all three, only keeping those with all three datasets
bothq <- inner_join(q1,q2,by="workerId", suffix=c(".q1",".q2"))

#data cleaning
#bothq$Q2.5 <- NULL #remove self reported TurkID
bothq$na_count <- apply(bothq, 1, function(x) sum(is.na(x))-36) #create a variable for how many missing answers there are (negatives may indicate postive optional answers to eg. needing glasses, 32 expected nas)
bothq <- filter(bothq, RecordedDate.q1>="2020-11-10") #remove pilots collected before 10 Nov
bothq <- filter(bothq, na_count<(.10*(ncol(bothq)-36))) #remove participants who skipped more than 10% of the questions

initial_duplicates <- which(duplicated(bothq$workerId))

#remove duplicates
bothq <- filter(bothq, Finished.q1==TRUE, Finished.q2==TRUE) #removes duplicates that happen from restarting one of the qualtrics surveys
to.delete.dup <- vector()
final_duplicates <- which(duplicated(bothq$workerId)) #update duplicates variable
precleancount <- nrow(bothq)

#make sure everyone who we paid we have data for
paidnotdata <- setdiff(paidId, bothq$workerId)

#remove participants for various reasons
#no head injury, no eye issues of note. 
bothq <- filter(bothq, workerId!="PUT ID HERE") #repeated qualtrics after reporting disqualifier -- ID in txt file. ignored by git for anonymity

#LOAD MATLAB FOR ID CROSS CHECK
matlab <-  read.csv('G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/participantwise050521_withNull_withGoodErpe.csv')
trialmatlab <- read.csv('G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/trialwise050521_withNull_withGoodErpe.csv')
slope1 <- read.csv('G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/envslopes050521_withNull_withGoodErpe.csv')
slope2 <- read.csv('G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/aaslopes050521_withNull_withGoodErpe.csv')
erpe1 <- read.csv('G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/enverpe050521_withNull_withGoodErpe.csv')
erpe2 <- read.csv('G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/hyperpe050521_withNull_withGoodErpe.csv')
erpep <- read.csv('G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/enverpepeak050521_withNull_withGoodErpe.csv')

#set up simpleq
simpleq <- data.frame(bothq$workerId,bothq$SC9,bothq$SC0,bothq$SC3)
simpleq <- rename(simpleq, workerId=bothq.workerId, SATQ=bothq.SC9, AQ=bothq.SC0, SCC=bothq.SC3)
bothq <- rename(bothq, SATQ=SC9, AQ=SC0, SCC=SC3)

#join datasets
participantdata <- inner_join(bothq,matlab,by="workerId", suffix=c(".qual",".matlab"))
simpleparticipantdata <- inner_join(simpleq,matlab,by="workerId", suffix=c(".qual",".matlab"))
trialdata <- inner_join(simpleparticipantdata, trialmatlab, by="workerId",suffix=c(".qual",".matlab"))
slope1data <-  inner_join(simpleparticipantdata, slope1, by="workerId",suffix=c(".qual",".matlab"))
slope2data <-  inner_join(simpleparticipantdata, slope2, by="workerId",suffix=c(".qual",".matlab"))
erpe1data <-  inner_join(simpleparticipantdata, erpe1, by="workerId",suffix=c(".qual",".matlab"))
erpe2data <-  inner_join(simpleparticipantdata, erpe2, by="workerId",suffix=c(".qual",".matlab"))
erpepdata <-  inner_join(simpleparticipantdata, erpep, by="workerId",suffix=c(".qual",".matlab"))

#check payment of final dataset
notpaid <- setdiff(participantdata$workerId, paidId) #returns items that are in the first set but not in the second

#deidentification
beachfinalidlist <- participantdata$workerId 
participantdata$workerId <- NULL #remove computer collected TurkID
simpleparticipantdata$workerId <- NULL
trialdata$workerId <- NULL
slope1data$workerId <- NULL
slope2data$workerId <- NULL
erpe1data$workerId <- NULL
erpe2data$workerId <- NULL
erpepdata$workerId <- NULL

#save struct
write_csv(participantdata, "G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/participantdata.csv")
write_csv(simpleparticipantdata, "G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/simpleparticipantdata.csv")
write_csv(trialdata, "G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/trialdata.csv")
write_csv(slope1data,"G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/envslopesdata.csv")
write_csv(slope2data,"G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/aaslopesdata.csv")
write_csv(erpe1data,"G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/enverpedata.csv")
write_csv(erpe2data,"G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/hyperpedata.csv")
write_csv(erpepdata,"G:/My Drive/Experiments/JoA_squarestask/BeachTask/Data and Analysis/Analysis/Results/enverpepeakdata.csv")