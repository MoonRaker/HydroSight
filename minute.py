
import numpy as np
from datetime import datetime

    
#MINUTE Summary of this function goes here
def minute(date_asNum):
    
    if isalpha(date_asNum):
        print 'date_asNum must be a date vector, not character.'
        return 
     
    # Get date vectors
    c = datetime.strptime(date_asNum)
        
    # Get year and reformat to the same shape as input data.
    m = np.reshape(c[:,5], np.shape(date_asNum))
    
    return m

