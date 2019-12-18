
import numpy as np


# Pearson's type III impulse response transfer function class. 
class responseFunction_FerrisKnowlesJacobs():
    
    
    def __init__(self):
        obj = responseFunction_FerrisKnowles(bore_ID, forcingDataSiteID, siteCoordinates, options)
        #obj = responseFunction_JacobsCorrection() # How to choose which?
        #if nargin==5: # <<< ?
        return 
         
         
    # Set parameters
    def setParameters(self, params):
        if len(params)==3:
            responseFunction_FerrisKnowles(self, params[1:2, :])
            #setParameters@responseFunction_JacobsCorrection(self, params[3, :])
        elif len(params)==2:
            responseFunction_FerrisKnowles(self, params)
        elif len(params)==1:
            responseFunction_JacobsCorrection(self, params)
    
    
    # Get model parameters
    def getParameters(self):           
        params, param_names = responseFunction_FerrisKnowles(self)
        params[3,:], param_names[3,:] = responseFunction_JacobsCorrection.getParameters(self)            
        return params, param_names
    
    
    def getParameterValidity(self, params, param_names):                       
        isValidParameter = responseFunction_FerrisKnowles.getParameterValidity(self, params[1:2,:], param_names[1:2])
        isValidParameter[3,:] = responseFunction_JacobsCorrection.getParameterValidity(self, params[3,:], param_names[3])    
        return isValidParameter
    
    # Return fixed upper and lower bounds to the parameters.
    def getParameters_physicalLimit(self):
        params_upperLimit, params_lowerLimit = responseFunction_FerrisKnowles.getParameters_physicalLimit(self)
        params_upperLimit[3], params_lowerLimit[3] = responseFunction_JacobsCorrection.getParameters_physicalLimit(self)
        return params_upperLimit, params_lowerLimit
    
    
    # Return fixed upper and lower plausible parameter ranges. 
    # This is used to define reasonable range for the initial parameter sets
    # for the calibration. These parameter ranges are only used in the 
    # calibration if the user does not input parameter ranges.
    def getParameters_plausibleLimit(self):
        params_upperLimit, params_lowerLimit] = responseFunction_FerrisKnowles.getParameters_plausibleLimit(self)
        params_upperLimit[3], params_lowerLimit[3] = responseFunction_JacobsCorrection.getParameters_plausibleLimit(self)
        return params_upperLimit, params_lowerLimit
    
    
    # Transform the estimate of the response function * the forcing.
    # This undertakes the Jacob's correction for an unconfined aquifer.
    # If, in solution to the quadratic equation, a complex number is
    # produced, then the input h_star value is returned. Peterson Feb
    # 2013.
    def transform_h_star(self, h_star_est):       
       return responseFunction_JacobsCorrection.transform_h_star(self, h_star_est)
    
    
    # Extract the estimates of aquifer properties from the values of
    # alpha, beta and zeta.
    def getDerivedParameters(self):            
        params, param_names = responseFunction_FerrisKnowles.getDerivedParameters(self)
        Ksat = T/10.**self.zeta   
        params = [params[1,:] params[2,:] Ksat]
        param_names[1], param_names[2] = 'Ksat', 'Lateral conductivity'
        return params, param_names
    
    
    # delete class destructor
    #
    # Syntax:
    #   delete(obj)
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
           if isobject(self.propNames[i]):
               delete(self.propNames[i])
           else:
               self.propNames[i] = []
    