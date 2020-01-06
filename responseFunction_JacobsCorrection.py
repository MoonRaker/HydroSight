
import numpy as np


class responseFunction_JacobsCorrection:
    

    #RESPONSEFUNCTION_JACOBSCORRECTION Summary of this class goes here
    #   Detailed explanation goes here


    def __init__(self, params):

        if nargin==0:
            params = np.log10(1.)        
        # Set parameters for transfer function.
        setParameters(self, params)                 
                  
   
    def setParameters(self, params):

        # Set parameters
        self.zeta = params[1,:]
    
    
    def getParameters(self):

        # Get model parameters
        params[1,:] = self.zeta
        param_names = ['zeta']
        return params, param_names
    

    def getParameters_physicalLimit(self):

        # Return fixed upper and lower bounds to the parameters.
        params_upperLimit = np.inf
        params_lowerLimit = np.log10(eps())
        return params_upperLimit, params_lowerLimit

    
    def getParameters_plausibleLimit(self):

        # Return fixed upper and lower plausible parameter ranges. 
        # This is used to define reasonable range for the initial parameter sets
        # for the calibration. These parameter ranges are only used in the 
        # calibration if the user does not input parameter ranges.
        params_upperLimit = np.log10(100.)
        params_lowerLimit = np.log10(np.sqrt(eps()))-5.
    	return params_upperLimit, params_lowerLimit

            
    def getParameterValidity(self, params, param_names):                   
        
        zeta_filt = param_names[param_names=='zeta']
        zeta = params[zeta_filt,:]
        
        # Back transform to get Sat. Thickness.
        SatThickness = 10.**zeta

        # Get physical bounds.
        params_upperLimit, params_lowerLimit = getParameters_physicalLimit(self)

        # Check parameters are within bounds and T>0 and 0<S<1.
        if (SatThickness > 0) & (params >= params_lowerLimit(:, np.ones([1, np.shape(params, 2)])) & (params <= params_upperLimit(:, np.ones([1, np.shape(params, 2)])):            
            isValidParameter = True
        else:
            isValidParameter = False
	
    	return isValidParameter

    
    def transform_h_star(self, h_star_est):

        # Transform the estimate of the response function * the forcing.
        # This undertakes the Jacob's correction for an unconfined aquifer.
        # If, in solution to the quadratic equation, a complex number is
        # produced, then the input h_star value is returned. Peterson Feb
        # 2013.
        temp = 1.-2.*h_star_est[:,2]./10.**self.zeta  # <<< ./ ?
        filt = temp[temp>=0.]
        result = h_star_est[:,2]
        result[filt] = -10.**self.zeta .* (-1.+np.sqrt(1.-2.*h_star_est[filt,2]./10.*self.zeta)) # <<< .* and ./ ?
        return result


    def delete(self):

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
                  
        propNames = properties(self)
        for i in range(len(propNames)):
            if np.isempty(self.(propNames[i])):
                continue               
            if np.isobject(self.(propNames[i])): # <<< isselfect ?
                delete(self.propNames[i])
            else:               
                self.propNames[i] = [] 
