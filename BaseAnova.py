#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:16:34 2023

@author: Aaron Mooney
"""

import pandas as pd
import statistics as stats
import scipy

### Calculate Anova Instructions ###

### Use this script to calculate the F-value from a series of raw test data with multiple levels.
### The data must be contained in a simple table with the following column format:
### Table Heading - [Sample #] [Level_Var_1] [Level_Var_2] [Level_Var_3]...[Level_Var_n]

### Remove extra or unused columns and rows. Save the file as a single sheet .csv file.
### If you do not know how to find the path to the data, save the data and the script
### into the same folder. 

### Be sure to change the file name in the next line to the name of your file.

### *************************************************************************

# *** Add your file name here! *** 
# Example: BaseTestData_2Lvls.csv

inFile = 'BaseAnova_TestData_4Lvls.csv'

### **************************************************************************

# Reads the data set into a Pandas DataFrame
df = pd.read_csv(inFile)

# Your raw experimental data may be of type 'Object' internally. This function
# converts it to type 'float' for use in calculations. 
def convertToFloat(df):
    """Reads in a DataFrame containing columns with float values set as
    'object' data types (df.dtype). Loops through the list and converts
    each to type float for use in calculations. Returns None"""
    
    colList = df.columns
    listLen = len(colList)
    
    for item in range(listLen):
        try:
            df[colList[item]] = df[colList[item]].astype(float)
        except:
            continue
        
    return None

# Step 1: Calculates the sums for each column. 
def getColSums(df):
    """Reads in a cleaned DataFrame. Reads the values in each column and
    calculates the sums for each column.
    Returns a list containing the column sums."""
    
    colList = df.columns
    listLen = len(colList)
    colSum = 0
    colSumsList = []
    
    for item in range(listLen):
        
        if colList[item] != 'Sample':
            colSum = df[colList[item]].sum()
            colSumsList.append(colSum)
        
    return colSumsList

# Step 2: Calculates the sum of squares values for each column.
def getSumSqdVals(df):
    """Reads in a DataFrame containing columns of values. Calculates the
    sum of the squared values for each column and sum of column squared values.
    Returns a list of the squared column values and total squared value."""
            
    sumSqd = []
    colList = df.columns
    listLen = len(colList)
    ySqd = 0
    
    for i in range(1,listLen):
        axSqd = sum(df[colList[i]] ** 2)
        sumSqd.append(axSqd)
        ySqd = ySqd + axSqd
   
    sumSqd.append(ySqd)
    
    return sumSqd


# Step 3: Calculates the mean values for each column.
def getMeanVals(df):
    """Reads in a DataFrame containing columns of values. Calculates the
    mean for each column. Returns the mean values for each column."""
            
    meanValsList = []
    colList = df.columns
    listLen = len(colList)
    meanVal = 0
    
    for i in range(1,listLen):
        meanVal = round((stats.mean(df[colList[i]])), 2)
        meanValsList.append(meanVal)
           
    return meanValsList

# Step 4a: Calculates the variance values for each column. While not strictly
# necessary for calculating the Anova, these values are used to determine
# whether the assumption of homogeneity of variance has been satisfied.
def getVarVals(df):
    """Reads in a DataFrame containing columns of values. Calculates the
    variance for each column. Returns the variance values for each column."""
            
    varianceValsList = []
    colList = df.columns
    listLen = len(colList)
    varianceVal = 0
    
    for i in range(1,listLen):
        varianceVal = round((stats.variance(df[colList[i]])), 2)
        varianceValsList.append(varianceVal)
           
    return varianceValsList

# Step 4b: Determines the FMax value, an indicator of whether the homogeneity
# of variance assumption has been met.
def getFMaxVal(df):
    """Reads in the DataFrame, checks the variance, and then calculates the
    ratio used for the FMax calculation for homogeneity of variance.
    Returns the FMaxVal."""
    
    varianceVals = getVarVals(df)
    minVarVal = min(varianceVals)
    maxVarVal = max(varianceVals)
    
    FMaxVal = round((maxVarVal / minVarVal), 2)
    
    return FMaxVal

# Step 4c: Evaluates whether homogeneity of variance assumption is met.
def checkFMaxVal(df):
    """Reads in the DataFrame. Calculates the FMaxVal which determines 
    whether the samples pass the homogeneity of variance assumption.
    Returns FMaxList."""
    
    FMaxVal = getFMaxVal(df)
    FMaxTest = 0
    
    print('\nHomogeneity of Variance Test Results:')
    
    if FMaxVal < 3:
        print('FMax:', FMaxVal, '< 3.0.')
        print('Samples meet homogeneity of variance assumption.')
        FMaxTest = 1
        
    else:
        print('FMax:', FMaxVal, '> 3.0.')
        print('Samples do not meet homogeneity of variance assumption.')
        print('Test is not appropriate for this sample.')
        FMaxTest = 0
    
    FMaxList = [FMaxVal, FMaxTest]
    
    return FMaxList

# Step 5: Calculates the 'Basic Ratios', the values necessary for determining
# the SS_A, SS_SA, and SS_T values.
def getBasicRatios(df):
    """Reads in a DataFrame. Calculates the statistics necessary to 
    return the Basic Ratios used in Anovas. 
    Returns a list object containing the Basic Ratio values."""
    
    colSumList = getColSums(df)
    meanValsList = getMeanVals(df)
    sumSqdVals = getSumSqdVals(df)
    
    listLen = len(colSumList)
    sumASqd = 0
    colSums = 0
    
    for i in range(listLen):
        sumASqd = sumASqd + (colSumList[i] ** 2)
        colSums = colSums + (colSumList[i])
    
       
    numConds = (df.shape[1] - 1)    # (a) in Basic Ratios.
    sampSize = df.shape[0]          # (n) in Basic Ratios.
        
    brY = sumSqdVals[-1]            # [Y] in Basic Ratios.
    brA = round((sumASqd / sampSize), 2) # [A] in Basic Ratios.
    brT = round(((colSums ** 2) / (numConds * sampSize)), 2) # [T] in Basic Ratios.
    
    basicRatios = [brY, brA, brT]
    
    return basicRatios

# Step 6: Calculate the degrees of freedom for the MS calculations.
def getDegFree(df):
    """Reads in a DataFrame. Calculates the degrees of freedom needed for 
    use in Anova calculations. 
    Returns a list object containing values for df_A, df_SA, df_T."""
    
    sampSize = (df.shape[0])
    numConds = (df.shape[1] - 1)

    df_A = (numConds - 1)
    df_SA = (numConds * (sampSize - 1))
    df_T = ((numConds * sampSize) - 1)
    
    degFree = [df_A, df_SA, df_T]
    
    return degFree

# Step 7: Calculate the SS_A, SS_SA, and SS_T values.
def getSSVals(df):
    """Reads in a DataFrame object and calculates the Sum of Squares
    values associated with each of the basic ratio components.
    Returns a list of the SS values."""

    basicRatios = getBasicRatios(df)
    
    brY = basicRatios[0]
    brA = basicRatios[1]
    brT = basicRatios[2]
    
    SS_A = round((brA - brT), 2)
    SS_SA = round((brY - brA), 2)
    SS_T = round((brY - brT), 2)
    
    SS_Vals = [SS_A, SS_SA, SS_T]
    
    return SS_Vals

# Step 8: Calculate the MS_A, MS_SA, and MS_T values.
def getMSVals(df):
    """Reads in a DataFrame and calculates the Mean Squares values 
    for each component of the Anova.
    Returns the a list of the MS values."""
    
    SS_Vals = getSSVals(df)
    degFree = getDegFree(df)
    
    SS_A = SS_Vals[0]
    SS_SA = SS_Vals[1]
    SS_T = SS_Vals[2]
    
    df_A = degFree[0]
    df_SA = degFree[1]
    df_T = degFree[2]
    
    MS_A = round((SS_A / df_A), 2)
    MS_SA = round((SS_SA / df_SA), 2)
    MS_T = round((SS_T / df_T), 2)
    
    MS_Vals = [MS_A, MS_SA, MS_T]
    
    return MS_Vals

# Step 9: Calculate the F-statistic for the samples.
def getFCalcVal(df):
    """Reads in a DataFrame. Calculates the associated values necessary
    to calculate the F-statistic value.
    Returns a list containing the F-statistic."""
    
    MS_Vals = getMSVals(df)
    
    MS_A = MS_Vals[0]
    MS_SA = MS_Vals[1]
    
    FCalcVal = round((MS_A / MS_SA), 2)
    
    return FCalcVal

# Step 10: Look up the tabled value of F based on the alpha and type of test.
def getFTabledVal(df, alpha, tailTest):
    """Reads in a DataFrame, the significance level (alpha), and the 
    type of test to use (one or two-tailed). Calculates the F-Tabled value
    based on those parameters.
    Returns the F-Tabled value."""
    
    q = (1 - (alpha / tailTest))
    degFree = getDegFree(df)
    
    df_A = degFree[0]
    df_SA = degFree[1]
    
    rawFTabVal = scipy.stats.f.ppf(q, df_A, df_SA)
    FTabVal = round((rawFTabVal), 2)
    
    return FTabVal       

# Step 11: Compare the calculated vs. the tabled F-values to determine significance.    
def getResultType(df, alpha, tailTest):
    """Reads in the DataFrame, significance level (alpha), and tailTest values.
    Compares the values and prints a comparison string.
    Returns None."""
    
    FCalcVal = getFCalcVal(df)
    FTabVal = getFTabledVal(df, alpha, tailTest)
    
    if FCalcVal > FTabVal:
        print('\nResults:')
        print('F-calculated:', FCalcVal, '\nF-tabled:', FTabVal, '\n')
        #print(FCalcVal, '>', FTabVal, 'at', alpha)
        print('Significant at the', alpha, 'level.')
        
    else:
        print('\nResults:')
        print('F-calculated:', FCalcVal, '\nF-tabled:', FTabVal, '\n')
        #print(FCalcVal, '<', FTabVal, 'at', alpha)
        print('Not significant at the', alpha, 'level.')
                
    return None

# This procedure calls steps 1 - 11 and prints the result. 
# Note: The tests will not run if the homogeneity of variance assumption isn't met.
def testData(df, alpha, tailTest):
    """Reads in a DataFrame, the significance level (alpha), and the type
    of test (one or two-tailed). Checks homogeneity of variance assumption,
    prints results. Returns None."""    
    
    varianceVals = getVarVals(df)
    FMaxList = checkFMaxVal(df)
    
    if FMaxList[1] == 1:
        FCalcVal = getFCalcVal(df)
        FTabVal = getFTabledVal(df, alpha, tailTest)
        getResultType(df, alpha, tailTest)
    
    else:
        print('Assumptions for homogeneity of variance not met.')
        print('This test is not appropriate for these samples.')
        
    return None

# Converts column variables into float values for use in calculations.
convertToFloat(df)

# Test significance level.
alpha = 0.05

# Type of test (1-tailed or 2-tailed).
tailTest = 1

# Calculate the F-Value.
testData(df, alpha, tailTest)




