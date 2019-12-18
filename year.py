
import numpy as np
from datetime import datetime
    

#YEAR Summary of this function goes here
#   Detailed explanation goes here
def year(date_asNum):
    
    if isalpha(date_asNum):
        print 'date_asNum must be a date vector, not character.'
        return 
     
    # Get date vectors
    c = datetime.strptime(date_asNum)
        
    # Get year and reformat to the same shape as input data.
    y = reshape(c[:,1], np.shape(date_asNum))
    
    return y

