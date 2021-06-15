#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 13:30:46 2020

@author: Jonathan E Robinson
"""
import pandas as pd
import numpy
import scipy.io as sio
import os


###############################################################################
# Function checks whether csv has already been created or needs creating.
###INPUTS###
#f = filename.
#participant = the prefix to he output mat you want
###OUPTUTS###
#pracDict = a dictionary containing keys and values ready to be saved as .mat files
# P###_practice.mat = saved out .mat in script dir of pracDict
#expDict = a dictionary containing keys and values ready to be saved as .mat files
# P###_experiment.mat = saved out .mat in script dir of expDict
###############################################################################
def csvParse(f,participant):
    dfRead = pd.DataFrame() #empty data fram
    try:
        dfRead = pd.read_csv(f) #read a csv file
        pracDict = {} # empty dicts to be filled
        expDict = {} # empty dicts to be filled
        
        for i in dfRead: #for each 'column' in structure. 
            #Here these take the form of keys (column name) and values (values under column name) in a dict
            irep = i.replace('.','') #remove dots for python save
            if ('practice' in i) or ('p_' in i and 'resp_' not in i) or ('cp' in i): #isolate practice items for practice dictionary
                valNum = 0
                for index,val in enumerate(dfRead[i]):
                    if type(val) == str:
                        val = val.replace('null','100000') #give null an impossible value 
                        dfRead[i][index] = val 
                try:
                    eval(dfRead[i][0]) # set to first practice trial index
                    valNum = 1 # can this array be evaluate set to 1 if yes
                except:
                    pass
                if valNum == 1:
                    for ind,j in enumerate(dfRead[i]): #move through trials
                        if type(j) == str:
                            tmp = numpy.array(eval(j)) #convert item to matrix
                        else:
                            tmp = numpy.array(j) # dont make matrix is eval is not possible
                        dfRead[i][ind] = tmp #put in shape
                    pracDict[irep] = pd.Series.to_numpy(dfRead[i]) #make CSV column into call array
                else:
                    pracDict[irep] = pd.Series.to_numpy(dfRead[i]) #make CSV column into call array
            
            else: #trial items
                valNum = 0
                for index,val in enumerate(dfRead[i]):
                    if type(val) == str:
                        val = val.replace('null','100000') #give null an impossible value
                        dfRead[i][index] = val 
                try:
                    eval(dfRead[i][16]) #set to first real trial (17) index [16]
                    valNum = 1 # can this array be evaluate 
                except:
                    pass
                if valNum == 1:
                    for ind,j in enumerate(dfRead[i]): #move through trials
                        if type(j) == str:
                            tmp = numpy.array(eval(j)) #convert item to matrix
                        else:
                            tmp = numpy.array(j) # dont make matrix is eval is not possible
                        dfRead[i][ind] = tmp #put in shape
                    expDict[irep] = pd.Series.to_numpy(dfRead[i]) #make CSV column into call array
                else:
                    expDict[irep] = pd.Series.to_numpy(dfRead[i]) #make CSV column into call array
        # save function call   
        saveMat(pracDict, participant+'_practice.mat')
        saveMat(expDict, participant+'_experiment.mat')
    except:
        print('Warning file empty: '+f)
        pracDict = {} # empty dicts to be filled
        expDict = {} # empty dicts to be filled
    return(pracDict, expDict)  
    
    
###############################################################################    
# Function saves Mat.
###INPUTS###
#dataDict = containing variables (keys) from columns and matrices (value) from columns filled by the program
#name = filename
###OUPTUTS###
#saved mat file in script location
###############################################################################
def saveMat(dataDict,name):
    sio.savemat(name, dataDict)


###############################################################################    
# Function calls CSVs in folder.
###INPUTS###
#folder = directory to be read
###OUPTUTS###
#saved mat files in script location
###############################################################################
def folderRead(folder):
    count = 0 #give paricipant a number based on load in order
    files = os.listdir(folder)
    files.sort() #numerically order
    for file in files: #iterate through files in folder
        if file.endswith(".csv"): #if they are csvs
            count += 1 #add to participant count   
            scount = str(count) #string to be used for padding
            print(file)
            pracDict, expDict = csvParse(folder+file,'P'+scount.zfill(3)) #carry out function
    return(pracDict, expDict)
    
#function call
x,y = folderRead('G:\My Drive\Experiments\JoA_squarestask\BeachTask\Data and Analysis\Data Storage\FinalSet\\')
#folderRead('<UserDirectory>/<DatacsvDirectory>/')

