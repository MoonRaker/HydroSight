
import numpy as np


# Pearson's type III impulse response transfer def class. 
class responsedef_Hantush(self):


        def __init__(self, bore_ID, forcingDataSiteID, siteCoordinates, options, params):
            
            # Define default parameters 
            if nargin==4:
                self.gamma = 0.1
                
            # Set parameters for transfer def.
            #setParameters(self, params)                 
       
       
        # Set parameters
        def setParameters(self, params):
            if np.shape(params)[1]==3:
                setParameters @responsedef_FerrisKnowles(self, params(1:2,:)) # <<< How to handle inheritance?       
                self.gamma = params[3,:]            
            elif np.shape(params)[1]==2:
                setParameters@responsedef_FerrisKnowles(self, params(1:2,:)) # <<< Ditto ?
        
        
        # Get model parameters
        def getParameters(self):
            params, param_names = responsedef_FerrisKnowles.getParameters(self) # <<< ?
            params[3,:] = self.gamma       
            param_names[3,1] = 'gamma'
            return params, param_names
        
        
        def getParameterValidity(self, params, param_names):
            # Initialise output.
            isValidParameter = True(np.shape(params))
            
            alpha_filt = param_names[param_names=='alpha']
            beta_filt  = param_names[param_names=='beta']
            gamma_filt = param_names[param_names=='gamma']
            
            alpha = params[alpha_filt,:]
            beta  = params[beta_filt,:]
            gamma = params[gamma_filt,:]
            
            # Calculate T, leakage and S.
            T = 1./(4.*np.pi .* 10.**alpha)
            S = 4.*10.**beta .* T    
            Leakage = 1./(S .* 10.**gamma)            
                        
            # Get physical bounds.
            params_upperLimit, params_lowerLimit = getParameters_physicalLimit(self)

            # Check gamma is within bounds.
            isValidParameter = ((repmat(S >=0 & S <1 & T>= 0,size(params,1),1)) & 
                                (params >= params_lowerLimit(1:size(params,1),ones(1,size(params,2)))) & 
                                (params <= params_upperLimit(1:size(params,1),ones(1,size(params,2)))))       
            
            return isValidParameter 


        # Return fixed upper and lower bounds to the parameters.
        def getParameters_physicalLimit(self):
            params_upperLimit, params_lowerLimit = getParameters_physicalLimit@responsedef_FerrisKnowles(self)
            params_upperLimit[3,1] =  np.inf
            params_lowerLimit[3,1] = -np.inf
            return params_upperLimit, params_lowerLimit
        
        
        # Return fixed upper and lower plausible parameter ranges. 
        # This is used to define reasonable range for the initial parameter sets
        # for the calibration. These parameter ranges are only used in the 
        # calibration if the user does not input parameter ranges.
        def getParameters_plausibleLimit(self):
            params_upperLimit, params_lowerLimit = getParameters_plausibleLimit@responsedef_FerrisKnowles(self)
            params_upperLimit[3,1] =  10.
            params_lowerLimit[3,1] = -10.                        
            return params_upperLimit, params_lowerLimit
        
        
        # Calculate impulse-response def.
        def theta(self, t):           
            # Loop though each production bore and, if image wells exist,
            # account for them in the drawdown.
            result = np.zeros(np.shape(t)[1], np.shape(self.settings.pumpingBores)[1])
            for i in range(np.shape(self.settings.pumpingBores)[1]):
                # Calc. distance to obs well.
                pumpDistancesSqr = ((self.settings.obsBore.Easting -self.settings.pumpingBores[i,1].Easting) .** 2. + 
                                    (self.settings.obsBore.Northing-self.settings.pumpingBores[i,1].Northing) .** 2.)
                
                if isfield(self.settings.pumpingBores[i,1], 'imageBoreID'):

                    # Calculate the distance to each image bore.
                    imageDistancesSqr = ((self.settings.obsBore.Easting -self.settings.pumpingBores[i,1].imageBoreEasting) .** 2. +
                                         (self.settings.obsBore.Northing-self.settings.pumpingBores[i,1].imageBoreNorthing) .** 2.)
                    
                    imageWellMultiplier = np.zeros([np.shape(self.settings.pumpingBores[i,1].imageBoreType)[1], 1])
                    
                    # create filter for recharge image wells
                    filt = cellfun(@(x)strcmp(x,'Recharge'),self.settings.pumpingBores[i,1].imageBoreType) # <<< ?
                    imageWellMultiplier[filt] = 1
                    
                    # create filter for no flow image wells
                    filt = cellfun(@(x)strcmp(x,'No flow'),self.settings.pumpingBores[i,1].imageBoreType) # <<< ?
                    imageWellMultiplier[filt] = -1
    
                    # Calculate the drawdown from the production well plus the influence from each image well.
                    result[:,i] = (bsxfun(@plus, -10.**self.alpha ./ t .* np.exp(-10.**self.beta*(pumpDistancesSqr ./ t)-10.**self.gamma .* t), 
                                   np.sum(bsxfun(@times, imageWellMultiplier, 
                                   bsxfun(@times, 10.**self.alpha ./ t, np.exp(bsxfun(@plus, -10.**self.beta* 
                                   bsxfun(@rdivide, imageDistancesSqr, t), -10.**self.gamma.*t)))),2)))
                else:
                    result[:,i] = -10.**self.alpha ./ t .* np.exp(-10.**self.beta*(pumpDistancesSqr ./ t)-10.**self.gamma .* t)
            
            # Set theta at first time point to zero. NOTE: the first time
            # point is more accuratly estimated by intTheta_lowerTail().
            result[t==0, :] = 0.
            return result 
            
            
        # Calculate integral of impulse-response def from t to inf.
        # This is used to minimise the impact from a finite forcing data
        # set.
        # TODO: IMPLEMENTED integral of theta
        def intTheta_upperTail2Inf(self, t):
            result = np.zeros([np.shape(self.settings.pumpingBores)[1], 1]) 
            return result


        # Numerical integration of impulse-response def from 0 to 1.
        # This is undertaken to ensure the first time step is accuratly
        # estimated. This was found to be important for highly transmissive aquifers.
        def intTheta_lowerTail(self, t):             
            # Calculate Theta at 1 minute time steps.
            delta_t = 1./(60.*24.)
            t_0to1 = [eps:delta_t:t]
            theta_0to1 = theta(self, t_0to1)           
            # Undertake Simpson's 3/8 composite integratraion for 1 minute time steps.            
            result = (3.*delta_t/8 .* [theta_0to1[1,:]+np.sum(3.*theta_0to1[[2:3:end-3],:]+3.*
                                       theta_0to1[[3:3:end-2],:]+2.*theta_0to1[[4:3:end-1],:],1]+
                                       theta_0to1[end,:]])
            return result
            
        
        # Extract the estimates of aquifer properties from the values of
        # alpha, beta and zeta.
        def get_AquiferProperties(self):
            T = 1./(4.*np.pi .* 10.**self.alpha)
            S = 4.*10.**self.beta .* T    
            Leakage = 1./(S .* 10.**self.gamma)            
            return T, S, Leakage
       
        def getDerivedParameters(self):
            T, S, Leakage = get_AquiferProperties(self)
            params = [T, S, Leakage]
            param_names = ['T: Transmissivity (Head units^2/day)', 'S: Storativity', 'G : Leakage Param.']
            return params, param_names

        def getDerivedDataTypes(self):
            derivedData_types = cell(size(self.settings.pumpingBores,1),1)
            for i=1:size(self.settings.pumpingBores,1)
                derivedData_types[i,1] = [self.settings.pumpingBores[i,1].BoreID,' weighting']
            return derivedData_types 
        
        
        # Return the theta values for the GUI 
        def getDerivedData(self, derivedData_variable, t, axisHandle):           
            # Find which bore to extract data for.
            ind = []
            for i in range(np.shape(self.settings.pumpingBores)[1]):
                if strcmp([self.settings.pumpingBores[i,1].BoreID,' weighting'], derivedData_variable): # <<< ?
                    ind = i
                    break
            if isempty(ind):
                print 'The following derived variable could not be found: ', derivedData_variable
            
            params, param_names = getParameters(self)
            nparamSets = np.shape(params)[2]
            setParameters(self, params[:,1])
            derivedData_tmp = theta[self, t]
            if nparamSets>1:
                derivedData = np.zeros(np.shape(derivedData_tmp)[1], nparamSets)                
                derivedData[:,1] = derivedData_tmp[:,ind]
                for i in range(2, nparamSets):
                    setParameters(self, params[:,i])
                    derivedData_tmp = theta[self, t]
                    derivedData[:,i] = derivedData_tmp[:,ind]
                setParameters(self, params)
                
                # Calculate percentiles
                derivedData_prctiles = np.percentile(derivedData,[5., 10., 25., 50., 75., 90., 95.])
                
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
                
                ind = find(np.abs(derivedData_prctiles[:,4])>np.max(np.abs(derivedData_prctiles[:,4]))*0.05, 1, 'last')
                if isempty(ind):
                    ind = len(t)   
                plt.xlim(axisHandle, [1, t[ind]])
                                                
                # Add legend
                plt.legend(axisHandle, '5-95th#ile', '10-90th#ile', '25-75th#ile', 'median', 'Location', 'northeastoutside')   
                
                # Add data column names
                derivedData_names = cell(nparamSets+1,1)
                derivedData_names[1,1] = 'Time lag (days)'
                derivedData_names[2:end, 1] = strcat(repmat({'Weight-Parm. Set '},1,nparamSets )',num2str([1:nparamSets ]'))                
                
                derivedData = [t, derivedData]
            else:
                plt.plot(axisHandle, t, derivedData_tmp[:,ind], 'b-')      
                t_ind = find(np.abs(derivedData_tmp[:,ind])>np.max(np.abs(derivedData_tmp[:,ind]))*0.05, 1, 'last')
                if isempty(ind):
                    ind = len(t)
                plt.xlim(axisHandle, [1, t[ind]])
                
                derivedData_names = ['Time lag (days)', 'Weight']
                derivedData = [t,derivedData_tmp[:,ind]]
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
