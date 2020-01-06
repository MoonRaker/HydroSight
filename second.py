
import numpy as np
from datetime import datetime

    
def second(date_asNum):
    
    #SECOND Summary of this function goes here

    if isalpha(date_asNum):
        print 'date_asNum must be a date vector, not character.'
        return
     
    # Get date vectors
    c = datetime.strptime(date_asNum)
        
    # Get year and reformat to the same shape as input data.
    s = np.reshape(c[:,6], np.shape(date_asNum))
    
    # from finance toolbox second.m
    s = round(1000.*s)./1000
    
    return s

