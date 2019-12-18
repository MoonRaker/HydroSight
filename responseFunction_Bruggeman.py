
import numpy as np


# Pearson's type III impulse response transfer function class. 
class responseFunction_Bruggeman(self):


    def __init__(self, bore_ID, forcingDataSiteID, siteCoordinates, options, params):           
        # Define default parameters 
        if nargin==4:
            params=[10., 10., 10.]
            
        # Set parameters for transfer function.
        setParameters(self, params)                 
        
        # No settings are required.
        self.settings = []
       
       
        # Set parameters
        def setParameters(self, params):
            self.alpha = params[1,:]
            self.beta  = params[2,:]
            self.gamma = params[3,:]            
        
        
        # Get model parameters
        def getParameters(self):
            params[1,:] = self.alpha
            params[2,:] = self.beta
            params[3,:] = self.gamma       
            param_names = ['alpha', 'beta', 'gamma']
            return params, param_names
        
        
        def getParameterValidity(self, params, param_names):
            # Get physical bounds.
            params_upperLimit, params_lowerLimit = getParameters_physicalLimit(self)
            # Check parameters are within bounds and T>0 and 0<S<1.
            isValidParameter = (params[params>=params_lowerLimit[:,np.ones([1, np.shape(params[2]))) & 
                               (params[params<=params_upperLimit[:,np.ones([1, np.shape(params)[2])) 
            return isValidParameter


        # Return fixed upper and lower bounds to the parameters.
        def getParameters_physicalLimit(self):
            params_upperLimit = np.inf*np.ones([3,1])
            params_lowerLimit = np.zeros([3,1])
            return params_upperLimit, params_lowerLimit
        
        
        # Return fixed upper and lower plausible parameter ranges. 
        # This is used to define reasonable range for the initial parameter sets
        # for the calibration. These parameter ranges are only used in the 
        # calibration if the user does not input parameter ranges.
        def getParameters_plausibleLimit(self):
            params_upperLimit = [100., 100., 100.]
            params_lowerLimit = [  0.,   0.,   0.]
            return params_upperLimit, params_lowerLimit
        
        
        # Calculate impulse-response function.
        def theta(self, t):
            result = -self.gamma ./ np.sqrt(np.pi*self.beta**2./self.alpha**2. .* t .** 3.) .* np.exp(-self.alpha**2. ./(self.beta .^ 2. .* t)-self.beta .** 2.*t)
            # Set theta at first time point to zero. NOTE: the first time
            # point is more accuratly estimated by intTheta_lowerTail().
            result[t==0, :] = 0
            return result 
        
        
        # Calculate integral of impulse-response function from t to inf.
        # This is used to minimise the impact from a finit forcign data
        # set.
        # TODO: IMPLEMENTED integral of theta
        def intTheta_upperTail2Inf(self, t):
            result = 0. 
            return result 


        # Calculate integral of impulse-response function from t to inf.
        # This is used to minimise the impact from a finit forcign data set.
        # TODO: IMPLEMENTED integral of theta
        def intTheta_lowerTail(self, t):
            result = 0. 
            return result 
        
        
        # Transform the estimate of the response function * the forcing.
        def transform_h_star(self, h_star_est):
           result = h_star_est[:, end]
           return result 
        
        
        # Return the derived variables.
        def getDerivedParameters(self):
            params = []
            param_names = cell(0,2)  
            return params, param_names


        def getDerivedDataTypes(self):
            derivedData_types = 'weighting'
            return derivedData_types 
             
             
        # Return the theta values for the GUI 
        def getDerivedData(self, derivedData_variable, t, axisHandle):
            params, param_names = getParameters(self)
            nparamSets = np.shape(params)[2]
            setParameters(self, params[:,1])
            derivedData_tmp = theta(self, t)            
            if nparamSets>1:
                derivedData = np.zeros(np.shape(derivedData_tmp)[1], nparamSets)
                derivedData[:,1] = derivedData_tmp            
                for i in range(2, nparamSets):
                    setParameters(self,params[:,i])
                    derivedData[:,i] = theta(self, t)
                setParameters(self, params)
                
                # Calculate percentiles
                derivedData_prctiles = np.percentile(derivedData,[5. 10. 25. 50. 75. 90. 95.])
                
                # Plot percentiles
                #XFill = [t' fliplr(t')]
                #YFill = [derivedData_prctiles(:,1)', fliplr(derivedData_prctiles(:,7)')]                   
                #fill(XFill, YFill,[0.8 0.8 0.8],'Parent',axisHandle)
                #hold(axisHandle,'on')                    
                #YFill = [derivedData_prctiles(:,2)', fliplr(derivedData_prctiles(:,6)')]                   
                #fill(XFill, YFill,[0.6 0.6 0.6],'Parent',axisHandle)                    
                #hold(axisHandle,'on')
                #YFill = [derivedData_prctiles(:,3)', fliplr(derivedData_prctiles(:,5)')]                   
                #fill(XFill, YFill,[0.4 0.4 0.4],'Parent',axisHandle)                    
                #hold(axisHandle,'on')
                #clear XFill YFill     

                # Plot median
                plt.plot(axisHandle, t, derivedData_prctiles[:,4], 'b-')
                #hold(axisHandle,'off')                
                
                ind = find(np.abs(derivedData_prctiles[:,4])>np.max(np.abs(derivedData_prctiles[:,4]))*0.01, 1, 'last')
                if isempty(ind):
                    ind = len(t)
                xlim(axisHandle, [1, t[ind]])
                
                # Add legend
                plt.legend(axisHandle, '5-95th#ile', '10-90th#ile', '25-75th#ile', 'median', 'Location', 'northeastoutside')   
                
                # Add data column names
                derivedData = [t, derivedData]
                derivedData_names = cell(nparamSets+1, 1)
                derivedData_names[1, 1] = 'Time lag (days)'
                #derivedData_names[2:end, 1] = strcat(repmat({'Weight-Parm. Set'},1,nparamSets )',num2str([1:nparamSets ]'))                
            else:
                plt.plot(axisHandle, t, derivedData_tmp, 'b-')                                   
                ind = find(np.abs(derivedData_tmp)>np.max(np.abs(derivedData_tmp))*0.05, 1, 'last')
                if isempty(ind):
                    ind = len(t)
                xlim(axisHandle, [1, t[ind]])
                
                derivedData_names = ['Time lag (days)', 'Weight']                
                derivedData = [t, derivedData_tmp]

            plt.xlabel(axisHandle, 'Time lag (days)')
            plt.ylabel(axisHandle, 'Weight')            
            #box(axisHandle,'on')
            return derivedData, derivedData_names
        
        
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
            if isempty(self.propNames[i]):
               continue
            if isobject(self.propNames[i]):
                delete(self.propNames[i])
            else:               
                self.propNames[i] = [] 
