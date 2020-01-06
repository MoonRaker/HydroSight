
import numpy as np


def doDataQualityAnalysis(headData, boreDepth, surface_elevation, casing_length, constuction_date, 
                          checkMinSartDate, checkMaxEndDate, chechDuplicateDates, checkMinHead, 
                          checkMaxHead, RateofChangeThreshold, ConstHeadThreshold, outlierNumStDevs, 
                          outlierForwardBackward):


    #EXPORTDATATABLE Summary of this function goes here
    #   Detailed explanation goes here

    # Handle situation where outlierForwardBackward is not set.
    if nargin<14:
        outlierForwadBackward = True

    # Minimum number of non-errorous observations required for undertaking
    # ARAM outlier detection
    minObsforOutlierDetection = 12

    # Duration for constant head error checkd min ob 
    constHeadThreshold_minObs = 3
    
    # Assign plausible dates for water level obs    
    plausibleEndDate = now()    
    
    errCode = -9999.99
    
    # Check there is enough data to run the analysis
    if len(headData)<=1:
        return
       
    # Sort by date 
    headData = np.sort(headData,1)

    # Convert the table data to arrays
    if istable(headData)==True:
        headData = headData[:,[1:2]]

    # Filter for plausible dates
    filt_date = False(len(headData),1)
    if checkMinSartDate==True:
        filt_date = headData[headData[:,1]<constuction_date]
    if checkMaxEndDate==True:
        filt_date = headData[(headData[:,1]>now()) | (headData[:,1]==filt_date)
    isErrorObs = filt_date

    # Filter date duplicates
    filt_duplicates = False(len(headData),1)
    if checkDuplicateDates==True:
        timeStep = [np.diff(headData[:,1]), inf]
        filt_duplicates = np.abs(timeStep)<np.sqrt(eps)    
        filt_duplicates[isErrorObs] = False        
    isErrorObs = filt_date | filt_duplicates # <<< ?

    # Check head is above the bottom of the bore.
    filt_minHead = False(len(headData),1)
    if checkMinHead==True:  
        filt_minHead_tmp = headData[~isErrorObs,2]<surface_elevation-boreDepth
        filt_minHead[~isErrorObs] = filt_minHead_tmp
        clear filt_minHead_tmp    
    isErrorObs = filt_date | filt_duplicates | filt_minHead # <<< ?

    
    # Check the head is below the top of the casing (assumes the aquifer is
    # unconfined)
    filt_maxHead = False(len(headData),1)
    if checkMaxHead==True:
        filt_maxHead_tmp = headData[~isErrorObs,2]>surface_elevation+casing_length        
        filt_maxHead[~isErrorObs] = filt_maxHead_tmp
        clear filt_maxHead_tmp
    isErrorObs = filt_date | filt_duplicates | filt_minHead | filt_maxHead # <<< ?

    
    # Filter for rapd change in headData
    filt_rapid = False(len(headData),1)    
    d_headData_dt = np.diff(headData[~isErrorObs,2]) ./ np.diff(headData[~isErrorObs,1])
    filt_rapid[~isErrorObs] = [False np.abs(d_headData_dt)>=RateofChangeThreshold]
    isErrorObs = filt_date | filt_duplicates | filt_minHead | filt_maxHead | filt_rapid # <<< ?
            
    # Filter out bore with a constant head for > 'ConstHeadThreshold' days. First the
    # duration of 'flat' periods is assessed.
    filt_flatExtendedDuration = False(len(isErrorObs),1)
    if RateofChangeThreshold>0:
        headData_tmp = headData[~isErrorObs,:]
        delta_headData_fwd = [False headData_tmp[2:end, 2]-headData_tmp[1:end-1, 2]]
        delta_headData_rvs = [headData_tmp[1:end-1, 2]-headData_tmp[2:end, 2] False]
        filt_flat = delta_headData_fwd==0 | delta_headData_rvs==0 # <<< ?
        if any(filt_flat)==True:

            filt_flatExtendedDuration_tmp = False(np.sum(~isErrorObs),1)

            # Filt out prior identified errors
            headData_filt = headData[~isErrorObs,:]

            startDate = 0
            endDate = 0
            startheadData = np.nan
            for j in range(2, len(headData_filt)):
                try:
                    if ((filt_flat[j]) & (~filt_flat[j-1])) | ((j==2) & (filt_flat[j])):
                        if (j==2) & (filt_flat[j-1]):
                            startDate = headData_filt[j-1,1]-np.sqrt(eps)
                            startheadData = headData_filt[j-1, 2]
                        else:
                            startDate = headData_filt[j,1]
                            startheadData = headData_filt[j,2]
                        endDate = 0                            
                    elif (startDate>0) & ((~filt_flat[j]) | (j==len(headData_filt)) | (headData_filt[j,2]~=startheadData)):
                        endDate = headData_filt[j,1]
                        if (j==len(headData_filt)) & (headData_filt[j,2]==headData_filt[j-1,2]):
                            endDate = endDate+np.sqrt(eps)

                        # Assess if the zero period is >60 days long
                        # and are > constHeadThreshold_minObs
                        filt_tmp = [(headData_filt[:,1]>= startDate) & (headData_filt[:,1] < endDate)]
                        consHead_dates =  headData_filt(filt_tmp,1)
                        if (consHead_dates(end)-consHead_dates[1]>=ConstHeadThreshold) & (np.sum(filt_tmp)>=constHeadThreshold_minObs):      
                            filt_flatExtendedDuration_tmp(filt_tmp) = True

                        # Reset markers
                        startDate = 0
                        endDate = 0   
                        startheadData = np.nan                            
                except:
                   print 'err'
            filt_flatExtendedDuration = False(len(isErrorObs),1)
            filt_flatExtendedDuration(~isErrorObs) = filt_flatExtendedDuration_tmp

    # Aggregate error filters
    isErrorObs = filt_date | filt_duplicates | filt_minHead | filt_maxHead | filt_rapid | filt_flatExtendedDuration  

    # Detect remaining outliers using a calibrated ARMA(1) model.            
    if (np.sum(~isErrorObs)>minObsforOutlierDetection) & (outlierNumStDevs>0):
        try:
            # Analyse outliers in forward time.
            isOutlierObs_forward, noise_sigma, ARMA_params, exp_model = outlierDetection(headData, isErrorObs, outlierNumStDevs)
            isOutlierObs = isOutlierObs_forward 
            
            if outlierForwardBackward:
                # Analyse outliers in reverse time.
                headData_reverse = headData[::-1]
                isErrorObs_reverse = isErrorObs[::-1]
                headData_reverse[:,1] = headData[end,1] - headData_reverse[:,1]
                isOutlierObs_reverse = outlierDetection(headData_reverse, isErrorObs_reverse, outlierNumStDevs)
                isOutlierObs_reverse = isOutlierObs_reverse[::-1]

                # Define as outlier if detected forward and reverse in time.
                isOutlierObs = (isOutlierObs_forward) & (isOutlierObs_reverse)
        except:
            print '    WARNING: Outlier detection failed.'
            isOutlierObs = False(len(isErrorObs_reverse))
            noise_sigma = []
            ARMA_params = []           
            exp_model = []
    else:
        isOutlierObs = False(len(isErrorObs))
        noise_sigma = []
        ARMA_params = []                
        exp_model = []

    # Delete calibration data files
    from os import remove
    remove('*.dat')
    
    # Aggregate the logical data from the analyis into a table and combine with the observed data.
    headData = table(year(headData(:,1)), month(headData(:,1)),  day(headData(:,1)), 
                     hour(headData(:,1)), minute(headData(:,1)), headData(:,2), 
                     filt_date, filt_duplicates, filt_minHead, filt_maxHead, filt_rapid, 
                     filt_flatExtendedDuration, isOutlierObs, 
                     'VariableNames',['Year', 'Month', 'Day', 'Hour', 'Minute', 'Head', 'Date_Error', 
                     'Duplicate_Date_Error', 'Min_Head_Error','Max_Head_Error','Rate_of_Change_Error',
                     'Const_Hear_Error','Outlier_Obs'])

    return headData,noise_sigma, ARMA_params, exp_model