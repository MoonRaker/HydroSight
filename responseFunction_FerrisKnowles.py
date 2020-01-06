
import numpy as np


class responsedef_FerrisKnowles:


    # Pearson's type III impulse response transfer def class. 


    def __init__(self, bore_ID, forcingDataSiteID, siteCoordinates, options, params):
        
        # Get the obs bore easting and northing.
        filt = cellfun(@(x)strcmp(x, bore_ID), siteCoordinates[:, 1]) # <<< ?
        self.settings.obsBore.BoreID   = bore_ID
        self.settings.obsBore.Easting  = siteCoordinates[filt, 2]
        self.settings.obsBore.Northing = siteCoordinates[filt, 3]
                                      
        # Get the number of pumping bores and loop through each to get
        # their easting and northing.               
        if iscell(forcingDataSiteID):
            nForcingSites = len(forcingDataSiteID)
        else:
            nForcingSites = 1
            forcingDataSiteID = [forcingDataSiteID]
        for j in range(nForcingSites):           
            filt = cellfun(@(x)strcmp(x, forcingDataSiteID[j]), siteCoordinates[:, 1]) # <<< ?
            self.settings.pumpingBores[j,1].BoreID   = siteCoordinates[filt, 1]
            self.settings.pumpingBores[j,1].Easting  = siteCoordinates[filt, 2]
            self.settings.pumpingBores[j,1].Northing = siteCoordinates[filt, 3]
                    
        # Set the image well options.
        nOptions = np.shape(options)[1]
        
        if nOptions>0:
            
            # Get the list of available options
            defaultOptions = responsedef_FerrisKnowles.modelOptions(bore_ID, forcingDataSiteID, siteCoordinates)                

            # Check that the options is a cell selfect of Nx3
            if (~isempty(options)) & (len(options)!=len(defaultOptions)):
                print 'The input options is inconsistent with that for responsedef_FerrisKnowles.'
                                            
            # Extract the available types of image wells.
            availableImageTypes = defaultOptions[1].colFormats[end]
            
            # Check the image well type is valid.
            for i in range(np.shape(options[1])[1]):
                filt = cellfun(@(x)strcmp(x, options[1][i, 3]), availableImageTypes) # <<< ?
                if ~any(filt):
                    print 'The image well types specified within the third column of the input options cell array can only be "Recharge" or "No flow".'
            
            # Check the first column contains only production bore IDs and
            # the second column does not contain production bore IDs (or obs bore ID). 
            for i in range(np.shape(options[1])[1]):
                # Check the left column ID is a forcingDataSiteID
                isSiteIDError = True
                for j in range(nForcingSites):
                    if forcingDataSiteID[j]==options[1][i, 1]:
                        isSiteIDError = False                        
                    elif bore_ID==options[1][i,1]:
                        isSiteIDError = True
                        break
                if isSiteIDError:
                    print 'The left column of the input data for image wells must contain only production bore IDs and cannot contain the obs. bore ID.'

                # Check the right column ID is a forcingDataSiteID
                isSiteIDError = False
                for j in range(nForcingSites):
                    if (forcingDataSiteID[j]==options[1][i,2]) | (bore_ID==options[1][i,2]):                            
                        isSiteIDError = True
                        break
                if isSiteIDError: 
                    print 'The right column of the input data for image wells cannot contain production bore IDs or the observation bore ID.'
            
            # Cycle through each production bore and get the image well
            # site IDs for the production bore.
            for j in range(np.shape(self.settings.pumpingBores)[1]):
                if ~isempty(options[1]):
                    filt = cellfun(@(x)strcmp(x, self.settings.pumpingBores[j,1].BoreID), options[1][:, 1])
                    if any(filt):
                        self.settings.pumpingBores[j,1].imageBoreID   = options[1][filt, 2]
                        self.settings.pumpingBores[j,1].imageBoreType = options[1][filt, 3]
            
            # Now cycle though each production bore and get the
            # coordinates for each image well.
            for i in range(np.shape(self.settings.pumpingBores)[1]):                                                            
                if isfield(self.settings.pumpingBores[i,1], 'imageBoreID'):
                    # Cycle though each image bore for current production bore and 
                    # find the coordinates.                        
                    nImageBores = np.shape(self.settings.pumpingBores[i,1].imageBoreID)[1])
                    for j in range(nImageBores):
                        # Get the image bore easting and northing.
                        filt = cellfun(@(x)(strcmp(x, self.settings.pumpingBores[i,1].imageBoreID(j,1))), siteCoordinates[:,1])
                        self.settings.pumpingBores[i,1].imageBoreEasting[j,1]  = siteCoordinates[filt, 2]
                        self.settings.pumpingBores[i,1].imageBoreNorthing[j,1] = siteCoordinates[filt, 3] 
        
        # Set the search radius option
        if ~isempty(options):
            if options[2][1,2]=='true':
                dynPropMetaData[1] = addprop(self, 'searchRadiusFrac') # <<< ?
                self.searchRadiusFrac = 0.5     

                # Find the maximum distance to pumps 
                for i in range(np.shape(self.settings.pumpingBores)[1]):
                    # Calc. distance to obs well.
                    pumpDistances[i] = (np.sqrt((self.settings.obsBore.Easting -self.settings.pumpingBores[i,1].Easting) .** 2.+
                                                 self.settings.obsBore.Northing-self.settings.pumpingBores[i,1].Northing) .** 2.))                
                self.settings.pumpingBoresMaxDistance = np.max(pumpDistances)
            if options[2][2,2]=='true':
                dynPropMetaData[2] = addprop(self, 'searchRadiusIsotropicRatio') # <<< ?
                self.searchRadiusIsotropicRatio = 0.5                
        else:
            self.settings.pumpingBoresMaxDistance = np.inf                
        
        # Define default parameters 
        if nargin==4:
            params = [np.log10(1e-4), np.log10(0.01)]
            
            if isprop(self, 'searchRadiusFrac'):
                params = [params, 0.5]
            if isprop(self, 'searchRadiusIsotropicRatio'):
                params = [params, 0.5]
           
        # Set parameters for transfer def.
        setParameters(self, params)                 
        
        
        def setParameters(self, params):

            # Set parameters
            self.alpha = params[1,:]
            self.beta  = params[2,:]        
            if isprop(self, 'searchRadiusFrac'):
                self.searchRadiusFrac = params[3]
            if isprop(self, 'searchRadiusIsotropicRatio')
                self.searchRadiusIsotropicRatio = params[4]
        
        
        def getParameters(self):

            # Get model parameters
            params[1,:] = self.alpha
            params[2,:] = self.beta    
            param_names = ['alpha', 'beta']
            if isprop(self, 'searchRadiusFrac')
                params[3,:] = self.searchRadiusFrac 
                param_names[3] = 'searchRadiusFrac'
            if isprop(self, 'searchRadiusIsotropicRatio')
                params[4,:] = self.searchRadiusIsotropicRatio 
                param_names[4] = 'searchRadiusIsotropicRatio'
            return params, param_names
        
        
        def getParameterValidity(self, params, param_names):

            alpha_filt = param_names=='alpha'
            beta_filt  = param_names=='beta'
            
            alpha = params[alpha_filt,:]
            beta  = params[beta_filt, :]
            
            # Calculate hydraulic transmissivity and S.
            T = 1./(4.*np.pi .* 10.**alpha)
            S = 4.*10.**beta .* T    
            
            # Get physical bounds.
            params_upperLimit, params_lowerLimit = getParameters_physicalLimit(self)

            # Check parameters are within bounds and T>0 and 0<S<1.
            isValidParameter = repmat((S>=1e-6) & (S<1) & (T>0), np.shape(params)[1], 1) & 
                                      (params>=params_lowerLimit[1:np.shape(params)[1]),ones([1, np.shape(params)[2])]))) & 
                                      (params<=params_upperLimit[1:np.shape(params)[1]),ones([1, np.shape(params)[2])])))
           
            return isValidParameter 
            
        
        def getParameters_physicalLimit(self):

            # Return fixed upper and lower bounds to the parameters.
            params_upperLimit = inf(2,1)
            params_lowerLimit = [log10(sqrt(eps)) log10(eps)]            
            
            if isprop(self,'searchRadiusFrac')
                params_upperLimit = [params_upperLimit 1]
                params_lowerLimit = [params_lowerLimit 0]
            if isprop(self,'searchRadiusIsotropicRatio')
                params_upperLimit = [params_upperLimit 1]
                params_lowerLimit = [params_lowerLimit 0]
 
            return params_upperLimit, params_lowerLimit
        
        
        def getParameters_plausibleLimit(self):

            # Return fixed upper and lower plausible parameter ranges. 
            # This is used to define reasonable range for the initial parameter sets
            # for the calibration. These parameter ranges are only used in the 
            # calibration if the user does not input parameter ranges.
            params_upperLimit = [ 1., -1.]
            params_lowerLimit = [-5., -7.]            

            if isprop(self, 'searchRadiusFrac'):
                params_upperLimit = [params_upperLimit, 1.]
                params_lowerLimit = [params_lowerLimit, 0.]
            if isprop(self, 'searchRadiusIsotropicRatio'):
                params_upperLimit = [params_upperLimit, 1.]
                params_lowerLimit = [params_lowerLimit, 0.]
                
            return params_upperLimit, params_lowerLimit
        
        
        def theta(self, t):  

            # Calculate impulse-response def for each pumping bore.
            # Loop though each production bore and, if image wells exist,
            # account for them in the drawdown.
            result = np.zeros([np.shape(t)[1]), np.shape(self.settings.pumpingBores)[1])])
            for i in range(np.shape(self.settings.pumpingBores)[1]):
                # Calc. distance to obs well.
                pumpDistancesSqr = (((self.settings.obsBore.Easting -self.settings.pumpingBores[i,1].Easting) .** 2.) +
                                    ((self.settings.obsBore.Northing-self.settings.pumpingBores[i,1].Northing) .** 2.))
                
                if isprop(self, 'searchRadiusFrac'):
                   if np.sqrt(pumpDistancesSqr)>(self.searchRadiusFrac*self.settings.pumpingBoresMaxDistance):
                       result[:,i] = 0
                       continue
                
                if isfield(self.settings.pumpingBores[i,1], 'imageBoreID'):

                    # Calculate the distance to each image bore.
                    imageDistancesSqr = (((self.settings.obsBore.Easting -self.settings.pumpingBores[i,1].imageBoreEasting) .**2.) +
                                         ((self.settings.obsBore.Northing-self.settings.pumpingBores[i,1].imageBoreNorthing) .** 2.))

                    imageWellMultiplier=zeros([np.shape(self.settings.pumpingBores[i,1].imageBoreType)[1], 1])
                    
                    # create filter for recharge image wells
                    filt =  cellfun(@(x)strcmp(x, 'Recharge'), self.settings.pumpingBores[i,1].imageBoreType)
                    imageWellMultiplier(filt) = 1
                    
                    # create filter for no flow image wells
                    filt =  cellfun(@(x)strcmp(x, 'No flow'), self.settings.pumpingBores[i,1].imageBoreType)
                    imageWellMultiplier(filt) = -1
    
                    # Calculate the drawdown from the production well plus
                    # the influence from each image well.
                    result[:,i] = 10.**self.alpha ./ t .* bsxfun(@plus, -np.exp(-10.**self.beta*(pumpDistancesSqr ./ t)), 
                                                   np.sum(bsxfun(@times, imageWellMultiplier, np.exp(-10.**^self.beta*
                                                          bsxfun(@rdivide, imageDistancesSqr, t))), 2))
                else:
                    result[:,i] = -10.**self.alpha ./ t .* np.exp(-10.**self.beta*(pumpDistancesSqr ./ t))
            
            # Set theta at first time point to zero. NOTE: the first time
            # point is more accuratly estimated by intTheta_lowerTail().
            result[t==0,:] = 0
            
            return result 
        
        
        def intTheta_upperTail2Inf(self, t):

            # Calculate integral of impulse-response def from t to inf.
            # This is used to minimise the impact from a finit forcign data
            # set.
            # TODO: IMPLEMENTED integral of theta
            result = np.zeros([np.shape(self.settings.pumpingBores)[1], 1]) 
            return result


        def intTheta_lowerTail(self, t):

            # Numerical integration of impulse-response def from 0 to 1.
            # This is undertaken to ensure the first time step is accuratly
            # estimated. This was found to be important for highly transmissive aquifers.
            # Loop though each production bore and, if image wells exist,
            # account for them in the drawdown.
            result = np.zeros([np.shape(t)[1], np.shape(self.settings.pumpingBores)[1]])
            for i in range(np.shape(self.settings.pumpingBores)[1]):
                # Calc. distance to obs well.
                pumpDistancesSqr = (((self.settings.obsBore.Easting -self.settings.pumpingBores[i,1].Easting) .** 2.) +
                                    ((self.settings.obsBore.Northing-self.settings.pumpingBores[i,1].Northing) .** 2.)
                
                if isfield(self.settings.pumpingBores[i,1], 'imageBoreID'):

                    # Calculate the distance to each image bore.
                    imageDistancesSqr = (((self.settings.obsBore.Easting -self.settings.pumpingBores[i,1].imageBoreEasting) .**2.) +
                                         ((self.settings.obsBore.Northing-self.settings.pumpingBores[i,1].imageBoreNorthing) .** 2.)
                    
                    imageWellMultiplier = np.zeros([np.shape(self.settings.pumpingBores[i,1].imageBoreType)[1], 1])
                    
                    # create filter for recharge image wells
                    filt =  cellfun(@(x)strcmp(x, 'Recharge'),self.settings.pumpingBores[i,1].imageBoreType)
                    imageWellMultiplier(filt) = 1
                    
                    # create filter for no flow image wells
                    filt =  cellfun(@(x)strcmp(x, 'No flow'),self.settings.pumpingBores[i,1].imageBoreType)
                    imageWellMultiplier(filt) = -1
    
                    # Calculate the drawdown from the production well plus
                    # the influence from each image well.
                    from scipy.special import expi
                    result[:,i] = (-10.**self.alpha*expi(10.**self.beta*(pumpDistancesSqr ./ t)) +
                                   +np.sum(imageWellMultiplier .* 10.**self.alpha .* expi(10.**self.beta*((imageDistancesSqr) ./ t))))
                else:
                    result[:,i] = -10.**self.alpha .* expi(10.**self.beta*(pumpDistancesSqr ./ t))
                
                # TEMP: CHECK integral using trapz
                # NOTE: As n approaches zero, theta(0) approaches inf. Trapz
                # integration of such a def produces a poor numerical estimate.
                #t_0to1 = 10.^([-10:0.0001:0])'
                #theta_0to1 = theta(self, t_0to1)
                #result_trapz = trapz(t_0to1, theta_0to1)                
                
        return result 
        
        
        def transform_h_star(self, h_star_est):

           # Transform the estimate of the response def * the forcing.
           result = h_star_est[:,2]
           return result
        
        
        def getDerivedParameters(self):

            # Extract the estimates of aquifer properties from the values of
            # alpha, beta and gamma.
            T = 1 ./ (4.*np.pi .* 10.**self.alpha)
            S = 4 .* 10.**self.beta .* T            
            params = [T, S]
            param_names = ['T : Transmissivity (Head units^2/day)', 'S : Storativity']
            return params, param_names            
        
        
        def getDerivedDataTypes(self):

            derivedData_types = cell(np.shape(self.settings.pumpingBores)[1], 1)
            for i in range(np.shape(self.settings.pumpingBores)[1]):
                derivedData_types[i,1] = [self.settings.pumpingBores[i,1].BoreID, ' weighting']
            return derivedData_types 
        
        
        def getDerivedData(self,derivedData_variable,t,axisHandle):
            
            # Return the theta values for the GUI 
            import matplotlib as mpl
            import matplotlib.pyplot as plt
            
            # Find which bore to extract data for.
            ind = []
            for i in range(np.shape(self.settings.pumpingBores)[1]):
                if [self.settings.pumpingBores[i,1].BoreID,' weighting']==derivedData_variable:
                    ind = i
                    break
            if isempty(ind):
                print 'The following derived variable could not be found: ', derivedData_variable
            
            params, param_names = getParameters(self)
            nparamSets = np.shape(params)[2]
            setParameters(self, params[:,1])
            derivedData_tmp = theta(self, t)
            if nparamSets>1:
                derivedData = np.zeros([np.shape(derivedData_tmp)[1], nparamSets])                
                derivedData[:,1] = derivedData_tmp[:, ind]
                for i in range(2, nparamSets):
                    setParameters(self, params[:,i])
                    derivedData_tmp = theta(self, t)
                    derivedData[:,i] = derivedData_tmp[:,ind]
                setParameters(self,params)
                
                # Calculate percentiles
                derivedData_prctiles = np.percentile(derivedData,[5., 10., 25., 50., 75., 90., 95.], 2)
                
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
                plt.set_xlim(axisHandle, [1., t[ind]])
                                                
                # Add legend
                plt.legend(axisHandle, '5-95th#ile', '10-90th#ile', '25-75th#ile', 'median', 'Location', 'northeastoutside')   
                
                # Add data column names
                derivedData_names = cell(nparamSets+1, 1)
                derivedData_names[1,1] = 'Time lag (days)'
                derivedData_names(2:end,1) = strcat(repmat(['Weight-Parm. Set '],1,nparamSets )',num2str([1:nparamSets ]'))                
                
                derivedData = [t, derivedData]
            else:
                plt.plot(axisHandle, t, derivedData_tmp[:,ind], 'b-')      
                t_ind = find(np.abs(derivedData_tmp[:,ind])>np.max(np.abs(derivedData_tmp[:,ind]))*0.05, 1, 'last')
                if isempty(t_ind):
                    t_ind  = len(t)
                plt.set_xlim(axisHandle, [1., t[t_ind]])
                
                derivedData_names = ['Time lag (days)', 'Weight']
                derivedData = [t, derivedData_tmp[:,ind]]
            plt.set_xlabel(axisHandle, 'Time lag (days)')
            plt.set_ylabel(axisHandle, 'Weight')
            #box(axisHandle,'on')
            
            return derivedData, derivedData_names
        
        
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
        #   self -  model selfect
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
            if isempty(self.propNames[i]):
               continue
            if isobject(self.propNames[i]):
                delete(self.propNames[i])
            else:               
                self.propNames[i] = [] 
