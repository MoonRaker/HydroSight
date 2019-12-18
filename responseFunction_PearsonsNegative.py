
import numpy as np


# Pearson's type III impulse response transfer function class. 
class responseFunction_Pearsons():


    def __init__(self, bore_ID, forcingDataSiteID, siteCoordinates, options, params):
        # Call inherited model constructor.
        #obj = responseFunction_Pearsons(bore_ID, forcingDataSiteID, siteCoordinates, options)

        # Assign model parameters if input.
        if nargin>5:
            setParameters(obj, params)
            
            
    # Calculate impulse-response function.
    def theta(self, t):
        # Call the source model theta function and change the sign of the output.
        return -responseFunction_Pearsons(self, t).theta


    # Calculate integral of impulse-response function from 0 to 1.
    def intTheta_lowerTail(self, t):
        # Call the source model intTheta function and change the sign of
        # the output.
        return -intTheta_lowerTail@responseFunction_Pearsons(self, t)                              

    
    # Calculate integral of impulse-response function from t to inf.
    # This is used to minimise the impact from a finit forcign data
    # set.
    def intTheta_upperTail2Inf(self, t):
        # Call the source model intTheta function and change the sign of
        # the output.
        return -intTheta_upperTail2Inf@responseFunction_Pearsons(self, t)                              
           

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
     def delete(self):
         propNames = properties(self)
         for i in range(len(propNames)):
            if np.isempty(self.propNames[i])):
                continue               
            if isobject(self.propNames[i])):
                delete(self.propNames[i]))
            else:
                self.propNames[i]) = [] 
