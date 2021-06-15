# beach-task
Analysis scripts and associated files related to the Beach Task by Perrykkad, Robinson and Hohwy

## Pipeline for Analysis

1. Run BeachAnalysis.m -- this script is the batch for processing the raw data from Pavlovia task

2. Run TrialwiseResults.m -- this will output tables that can be used for statistical analysis (more than are used for manuscript)

3. Run SurveyData.R (update matlab file paths)  -- add in any manually removed participants in line 66, this will combine the matlab processed data from the Pavlovia task with the survey data

   3a. This step may remove other participants. Create a new results.keep.participants variable to use for other analyses in matlab that were not part of the batch script.

   3b. Some manual looping through 2 and 3 will be required to match permutation testing participant list with output from trialwiseresults.m for which the windows are hard coded. For permutation analysis, use goodparticipantsNull.m 

   3c. Note: selecting null triggers for full dataset takes a long time. It is commented out in BeachAnalysis.m and can be run outside the loop -- one of these was generated for the full dataset (including rejected participant data), and this is what was used for the manuscript. There is some alternative script to do this for just kept participants in goodparticipantsNull.m, which was used to develop other scripts, but this was not used in the manuscript. 

4. Open trialdata.csv, simpleparticipantdata.csv (means across participants) and enverpepeakdata.csv (time window based off of permutation analysis) in Jamovi for stats analysis

5. Demographic info.R for demographic info

