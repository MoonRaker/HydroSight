
import numpy as np


class responseFunction_Pearsons:


    # Pearson's type III impulse response transfer function class. 


    def __init__(self, bore_ID, forcingDataSiteID, siteCoordinates, options, params):
    
        # Call inherited model constructor.
        #obj = responseFunction_Pearsons(bore_ID, forcingDataSiteID, siteCoordinates, options)

        # Assign model parameters if input.
        if nargin>5:
            setParameters(obj, params)
            
            
    def theta(self, t):
    
        # Calculate impulse-response function.
        # Call the source model theta function and change the sign of the output.
        return -responseFunction_Pearsons(self, t).theta


    def intTheta_lowerTail(self, t):

        # Calculate integral of impulse-response function from 0 to 1.
        # Call the source model intTheta function and change the sign of
        # the output.
        return -intTheta_lowerTail@responseFunction_Pearsons(self, t)                              

    
    def intTheta_upperTail2Inf(self, t):

        # Calculate integral of impulse-response function from t to inf.
        # This is used to minimise the impact from a finit forcign data
        # set.
        # Call the source model intTheta function and change the sign of
        # the output.
        return -intTheta_upperTail2Inf@responseFunction_Pearsons(self, t)                              
           

    def delete(self):

        # delete class destructor
        #
        # Syntax:
        #   delete(self)
        #
        # Description:
        #   Loops through parameters and, if not an object, empties them. Else, calls
        #   the sub-object's destructor.
        #
        # Input:
        #   obj -  model object
        #
        # Output:  
        #   (none)
        #
        # Author: 
        #   Dr. Tim Peterson, The Department of Infrastructure Engineering, 
        #   The University of Melbourne.
        #
        # Date:
        #   24 Aug 2016
        ##            
        
        propNames = properties(self)
        for i in range(len(propNames)):
           if np.isempty(self.propNames[i])):
               continue               
           if isobject(self.propNames[i])):
               delete(self.propNames[i]))
           else:
               self.propNames[i]) = [] 
