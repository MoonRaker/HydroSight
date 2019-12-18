
import numpy as np
from datetime import datetime


#HOUR Summary of this function goes here
def hour(date_asNum):

    if isalpha(date_asNum):
        print 'date_asNum must be a date vector, not character.'
        return
     
    # Get date vectors
    c = datetime.strptime(date_asNum)
        
    # Get year and reformat to the same shape as input data.
    h = np.reshape(c[:,4], np.shape(date_asNum))
    
    return h
