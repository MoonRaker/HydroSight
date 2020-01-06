
import numpy as np


class responsedef_Pearsons:


    # Pearson's type III impulse response transfer def class. 


    def __init__(self, bore_ID, forcingDataSiteID, siteCoordinates, options, params):
                        
        # Define default parameters 
        if nargin==4:
            params = [np.log10(1.), np.log10(0.01), np.log10(1.5)]
        
        # Set parameters for transfer def.
        setParameters(self, params)     
        
        # Initialise the private property, t_limit, 
        # defining the lower limit to an exponetial 
        # repsonse def.
        self.settings.t_limit = np.nan
        self.settings.weight_at_limit = np.nan
        
        if (~isempty(options)) & (iscell(options)):  
            self.settings.params_lowerPhysicalLimit = cell2mat(options[:,2])
            self.settings.params_upperPhysicalLimit = cell2mat(options[:,3])
       
       
        def setParameters(self, params)

            # Set parameters
            self.A = params[1,:]
            self.b = params[2,:]
            self.n = params[3,:]            
        
        
        def getParameters(self):

            # Get model parameters
            params[1,:] = self.A
            params[2,:] = self.b
            params[3,:] = self.n       
            param_names = ['A', 'b', 'n']
            return params, param_names
        
        
        def getParameterValidity(self, params, param_names):

            isValidParameter = True(np.shape(params))

            A_filt = [param_names=='A']
            b_filt = [param_names=='b']
            n_filt = [param_names=='n']
            
            A = params[A_filt, :]
            b = params[b_filt, :]
            n = params[n_filt, :]
                                    
            # Get physical bounds.
            params_upperLimit, params_lowerLimit = getParameters_physicalLimit(self)

    	    # Check parameters are within bounds.
            isValidParameter = ((params>=params_lowerLimit[:, np.ones([1, np.shape(params)[2]]]) & 
                                (params<=params_upperLimit[:, np.ones([1, np.shape(params)[2]]]))
                 
            # Check the b and n parameters will not produce np.nan or Inf
            # values when theta() is integrated to -inf (see intTheta()).
            isValidParameter(n_filt, gamma(10.**n)==np.inf) = False
            isValidParameter(b_filt | n_filt, (10.**b) .** (10.**n)<=0) = False            

            return isValidParameter

        
        def getParameters_physicalLimit(self):

            # Return fixed upper and lower bounds to the parameters.
            # NOTE: The upper limit for 'b' is set to that at which
            # exp(-10^b * t) <= sqrt(eps) where t = 1 day.
            #params_upperLimit = [inf log10(-log(sqrt(eps))) inf]
            #params_upperLimit = [log10(1/1000) log10(-log(sqrt(eps))) inf]
            #params_lowerLimit = [log10(sqrt(eps)) log10(sqrt(eps)) log10(sqrt(eps))]
            #params_upperLimit = [log10(1/1000/1e-6) log10(-log(sqrt(eps))) inf]
            #params_lowerLimit = [log10(1/1000) log10(sqrt(eps)) log10(sqrt(eps))]    
            if isfield(self.settings, 'params_lowerPhysicalLimit'):
                params_lowerLimit = self.settings.params_lowerPhysicalLimit
            else:
                params_lowerLimit = [np.log10(np.sqrt(eps)), np.log10(np.sqrt(eps)), np.log10(np.sqrt(eps))]                     
            if isfield(self.settings, 'params_upperPhysicalLimit'):
                params_upperLimit = self.settings.params_upperPhysicalLimit
            else:
                params_upperLimit = [np.inf, np.log10(-np.log(np.sqrt(eps))), np.inf]
            return params_upperLimit, params_lowerLimit
        
        
        def getParameters_plausibleLimit(self):
        
            # Return fixed upper and lower plausible parameter ranges. 
            # This is used to define reasonable range for the initial parameter sets
            # for the calibration. These parameter ranges are only used in the 
            # calibration if the user does not input parameter ranges.
            #params_upperLimit = [np.log10(1./1000./1e-4), np.log10(0.1), np.log10(10.)]
            params_upperLimit = [np.log10(1.), np.log10(0.1), np.log10(10.)]
            #params_lowerLimit = [np.log10(np.sqrt(eps))+2., np.log10(np.sqrt(eps))+2., np.log10(np.sqrt(eps))+2.]
            params_lowerLimit = [-5., -5., -2.]
            #params_lowerLimit = [np.log10(1./1000./0.1), -5., -2.]
            return params_upperLimit, params_lowerLimit
        
        
        def theta(self, t):

            # Calculate impulse-response def.
            # Back-transform parameter 'n' from natural log domain to
            # non-transformed domain. This transformation of n was
            # undertaken because ocassionally the optimal value of n was
            # very large eg 200. By transforming it as below, such large
            # values do not occur in the transformed domain and this allows
            # for a greater parameter range to be more easily accessed by
            # the gradient based calibration.
            n_backTrans = 10.**self.n
            
            # Back transform 'b'
            b_backTrans = 10.**self.b
            
            # Back transform 'A'
            A_backTrans = 10.**self.A
            
            # If n>1, then theta will have a value of 't' at which the
            # gradient is zero. By putting this value of 't' into theta
            # the value at the peak can be determined. This substitution
            # is done inside the theta calculation to avoid problems of
            # theta=inf at time time point.
            #
            # Also, the original version of the theta def (as given within von Asmuth 2005) 
            # can produce Inf values when n is large and b is also large. If
            # these are filter out then major discontinuities can be
            # produced in the response surface. If this occurs, then an
            # algabraically rearranged version is used that minimises the
            # emergence of inf.            
            # 
            # Lastly, the theta def below are a modified version of
            # von Asmuth 2005. The modification was undertaken to remove
            # the considerable covariance between 'b' and 'A' parameters
            if n_backTrans>1:
                
                t_peak = (n_backTrans-1.)/b_backTrans                              
                         
                # Calculate Pearsons value.
                result = A_backTrans ./ (t_peak .** (n_backTrans-1.)*np.exp(-b_backTrans*t_peak)) .* t .** (n_backTrans-1.) .* np.exp(-b_backTrans .* t)                
                                
                # New algebraic version to which minimise Inf and np.nan. The
                # equation is identical to that above but rearranged to
                # minimise Inf values. This was found to be essential for reliable calibration when n is large.           
                # Importantly, the scaling by the peak value of theta is
                # undertaken below inside the bracketted term set to a power of
                # (n_backTrans-1). This was undertaken to reduce the liklihood that
                # any time point have a value of inf.
                if any((np.isnan(result)) | (np.isinf(result))):            
                    result = A_backTrans .* (t .* b_backTrans ./ (n_backTrans-1.)*np.exp(1.) .* np.exp(-b_backTrans .* t ./ (n_backTrans-1.))) .** (n_backTrans-1.)                      
            else:
                # If the lower limit to an exponential response def
                # has not been defined, the set it to 100 years prior to the start
                # of the forcing data record.
                if np.isnan(self.settings.t_limit):
                    self.settings.t_limit = np.max(t)+365.*100.                                    
                self.settings.weight_at_limit = self.settings.t_limit .** (n_backTrans-1.) .* np.exp(-b_backTrans .* self.settings.t_limit)                               

                # Calculate the integral from t=0 to weight_at_limit. This
                # is used to normalise the weights by the integral. This is
                # undertaken to minimise the covariance of parameters b and
                # n with A.
                #integral0toInf = 1./(1-self.settings.weight_at_limit) .* (gamma(n_backTrans)./(b_backTrans^n_backTrans) .* (gammainc(0 ,n_backTrans ,'upper') - gammainc(b_backTrans .* self.settings.t_limit ,n_backTrans ,'upper')) - self.settings.weight_at_limit.* self.settings.t_limit)
                
                # Calculate the weighting def with the weight at the
                # limit subtracted and rescaled to between zero and one and
                # then multiplied by A.
                #result = A_backTrans./(integral0toInf.*(1-self.settings.weight_at_limit)) .* ( t.^(n_backTrans-1) .* exp( -b_backTrans .* t ) - self.settings.weight_at_limit)
                result = A_backTrans ./ (1.-self.settings.weight_at_limit) .* (t .** (n_backTrans-1.) .* np.exp(-b_backTrans .* t)-self.settings.weight_at_limit)
            
            # Set theta at first time point to zero. NOTE: the first time
            # point is more accuratly estimated by intTheta_lowerTail().
            result[t==0,:] = 0
            
            return result

        
        def theta_normalised(self, t):

            # Get non-normalised theta result
            result = theta(self, t)
            
            # Back transform 'A'
            A_backTrans = 10.**self.A
            
            # Normalised result by dividing theta by A_backTrans (ir peak
            # value)
            result  = result ./ A_backTrans
            
            return result, A_backTrans

            
        def intTheta_upperTail2Inf(self, t):

            # Calculate integral of impulse-response def from t to inf.
            # This is used to minimise the impact from a finit forcign data
            # set.
            # Back-transform parameter 'n' from natural log domain to
            # non-transformed domain. This transformation of n was
            # undertaken because ocassionally the optimal value of n was
            # very large eg 200. By transforming it as below, such large
            # values do not occur in the transformed domain and this allows
            # for a greater parameter range to be more easily accessed by
            # the gradient based calibration.
            n_backTrans = 10.**self.n
            
            # Back transform 'b'
            b_backTrans = 10.**self.b
            
            # Back transform 'A'
            A_backTrans = 10.**self.A
            
            # If n>1, then theta will have a value of 't' at which the
            # gradient is zero. By putting this value of 't' into theta
            # then value at the peak can be determined. 
            if n_backTrans>1:
                t_peak = (n_backTrans-1.)/b_backTrans
                theta_peak = t_peak**(n_backTrans-1.)*np.exp(-b_backTrans*t_peak)
                
                if any(np.isinf(theta_peak)):
                      theta_peak = ((n_backTrans-1.)/b_backTrans*np.exp(-1.))**(n_backTrans-1.)

                result = (A_backTrans .* gamma(n_backTrans) ./ (b_backTrans**n_backTrans .* theta_peak) .* 
                         (gammainc(b_backTrans .* t, n_backTrans, 'upper')-gammainc(b_backTrans .* np.inf, n_backTrans, 'upper')))
                
            else:
                # The lower limit to an exponential response def
                # should have been defined prior in a call to theta().
                if (np.isnan(self.settings.t_limit)) | (isnp.nan(self.settings.weight_at_limit)):
                    print 'When "exp(self.n) - 1<1", the method "theta()" must be called prior to "intTheta()" so that the lower time limit can be set.'
                delta_t_to_limit = self.settings.t_limit-t
                result = (A_backTrans ./ (1.-self.settings.weight_at_limit) .* (gamma(n_backTrans) ./ (b_backTrans**n_backTrans) .* 
                         (gammainc(b_backTrans .* t, n_backTrans, 'upper')-gammainc(b_backTrans .* self.settings.t_limit, n_backTrans, 'upper'))-
                         self.settings.weight_at_limit.*delta_t_to_limit)
            
            # Trials indicated that when tor (ie t herein) is very large,
            # result can equal np.nan.            
            if any((np.isnan(result)) | (np.isinf(result))):
                result[:] = np.nan                
           
            return result

        
        def intTheta_lowerTail(self, t):

            # Numerical integration of impulse-response def from 0 to 1.
            # This is undertaken to ensure the first time step is accurately
            # estimated.
            # Back-transform parameter 'n' from natural log domain to
            # non-transformed domain. This transformation of n was
            # undertaken because ocassionally the optimal value of n was
            # very large eg 200. By transforming it as below, such large
            # values do not occur in the transformed domain and this allows
            # for a greater parameter range to be more easily accessed by
            # the gradient based calibration.
            n_backTrans = 10.**self.n
            
            # Back transform 'b'
            b_backTrans = 10.**self.b
            
            # Back transform 'A'
            A_backTrans = 10.**self.A           
            
            # If n>1, then theta will have a value of 't' at which the
            # gradient is zero. By putting this value of 't' into theta
            # then value at the peak can be determined. 
            if n_backTrans>1:
                t_peak = (n_backTrans-1.)/b_backTrans
                theta_peak =  t_peak**(n_backTrans-1.)*np.exp(-b_backTrans*t_peak)
                
                # Recalculate theta_peak in a way to minimise rounding
                # errors.
                if any(np.isinf(theta_peak)):
                    theta_peak = ((n_backTrans-1.)/b_backTrans*np.exp(-1.))**(n_backTrans-1.)

                result = (A_backTrans .* gamma(n_backTrans) ./ (b_backTrans**n_backTrans .* theta_peak) .* 
                         (gammainc(0., n_backTrans, 'upper')-gammainc(b_backTrans .* t, n_backTrans, 'upper')))

                # If theta_peak still equals inf, then the integral from
                # 0-1 can be approximated as 0. To achieve this,                 
                if any(np.isnan(result)):
                    filt = np.isinf(theta_peak)
                    result[filt] = 0.
                
            else:
                
                # The lower limit to an exponential response def
                # should have been defined prior in a call to theta().
                if (np.isnan(self.settings.t_limit)) | (np.isnan(self.settings.weight_at_limit)):
                    print 'When "exp(self.n) - 1<1", the method "theta()" must be called prior to "intTheta()" so that the lower time limit can be set.'
                                
                # NOTE: As t -> 0, t^(n_backTrans-1) -> Inf. This can cause
                # the integral of the lower tail (it t from 0 to 1 day) to
                # be exceeding large when n_backTrans<0 (ie decaying
                # exponential like def).  This can cause the
                # contribution from climate to be implausible and (at least
                # for the Clydebank built in expamle) can cause the pumping
                # componp.nant to produce S<=10-6. To address this weakness,
                # the following excludes the first 1 hour from the
                # integration of the lower tail and the rescales it to the 
                # duration of the input t.
                t_to_omit =  1./24.
                result_to_omit = (A_backTrans ./ (1.-self.settings.weight_at_limit) .* (gamma(n_backTrans) ./ (b_backTrans**n_backTrans) .* 
                                 (gammainc(0., n_backTrans, 'upper')-gammainc(b_backTrans .* t_to_omit, n_backTrans, 'upper'))- 
                                 self.settings.weight_at_limit.*t_to_omit))
                result = (A_backTrans ./ (1.-self.settings.weight_at_limit) .* (gamma(n_backTrans) ./ (b_backTrans**n_backTrans) .* 
                         (gammainc(0., n_backTrans, 'upper')-gammainc(b_backTrans .* t, n_backTrans, 'upper'))-self.settings.weight_at_limit.*t))
                result = (result-result_to_omit) ./ (t-t_to_omit)
                result = np.zeros(np.shape(t))
            
            # TEMP: CHECK integral using trapz
            # NOTE: As n approaches zero, theta(0) approaches inf. Trapz
            # integration of such a def produces a poor numerical estimate.
#              t_0to1 = 10.^([-100:0.0001:0])'
#              theta_0to1 = theta(self, t_0to1)
#              result_trapz = trapz(t_0to1, theta_0to1)
#             if abs(result_trapz - result) > abs(0.05*result)
#                 display(['Pearsons tail integration error. Analytical est:', num2str(result),' Trapz:', num2str(result_trapz)])
#             elif( isnp.nan(result) || isinf(abs(result)))
#                     display('Pearsons tail integration error (inf or np.nan)')

        return result 
        
        
        # Transform the estimate of the response def * the forcing.
        def transform_h_star(self, h_star_est):
           result = h_star_est[:,end]
           return result 

        
        # Return the derived lag time (ie peak of def)
        def getDerivedParameters(self):
        
            # Back-transform parameter 'n' from natural log domain to
            # non-transformed domain. This transformation of n was
            # undertaken because ocassionally the optimal value of n was
            # very large eg 200. By transforming it as below, such large
            # values do not occur in the transformed domain and this allows
            # for a greater parameter range to be more easily accessed by
            # the gradient based calibration.
            n_backTrans = 10.**self.n  
            
            # Back transform 'b'
            b_backTrans = 10.**self.b
                        
            # If n>1, then theta will have a value of 't' at which the
            # gradient is zero. By putting this value of 't' into theta
            # then value at the peak can be determined. 
            t_peak = np.zeros(np.shape(n_backTrans))
            theta_peak = np.zeros(np.shape(n_backTrans))
            filt = n_backTrans>1.
            t_peak[filt] = (n_backTrans[filt]-1.) ./ b_backTrans[filt]
            theta_peak[filt] = t_peak[filt] .** (n_backTrans[filt]-1.) .* np.exp(-b_backTrans[filt] .* t_peak[filt])
                        
            params = [t_peak, theta_peak]
            param_names = ['Lag : Lag time from input to head (days)', 'Peak : Peak weighting of input to head']
            
            return params, param_names
            

        def getDerivedDataTypes(self):
            derivedData_types = 'weighting'
            return derivedData_types
            
        
        # Return the theta values for the GUI
        def getDerivedData(self,derivedData_variable,t,axisHandle):
           
            import matplotlib as mpl
            import matplotlib.pyplot as plt
           
            params, param_names = getParameters(self)
            nparamSets = np.shape(params)[2]
            setParameters(self, params[:,1])
            derivedData_tmp = theta(self, t)            
            if nparamSets>1:
                derivedData = np.zeros([np.shape(derivedData_tmp)[1], nparamSets])
                derivedData[:,1] = derivedData_tmp            
                for i in range(2, nparamSets):
                    setParameters(self, params[:,i])
                    derivedData[:,i] = theta(self, t)
                setParameters(self, params)
                
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
                if np.isempty(ind):
                    ind = len(t)
                plt.set_xlim(axisHandle, [1., t[ind]])
                
                # Add legend
                plt.legend(axisHandle, '5-95th#ile', '10-90th#ile', '25-75th#ile', 'median', 'Location', 'northeastoutside')   
                
                # Add data column names
                derivedData = [t, derivedData]
                derivedData_names = cell(nparamSets+1, 1)
                derivedData_names[1,1] = 'Time lag (days)'
                derivedData_names[2:end, 1] = strcat(repmat(['Weight-Parm. Set'], 1, nparamSets ), str([1:nparamSets ]))                
            else:
                plt.plot(axisHandle, t, derivedData_tmp, 'b-')                                   
                ind = find(np.abs(derivedData_tmp)>np.max(np.abs(derivedData_tmp))*0.05, 1, 'last')
                if np.isempty(ind):
                    ind = len(t)
                elif ind==1:
                    ind = np.ceil(len(t)*0.05)
                plt.set_xlim(axisHandle, [1., t[ind]])
                
                derivedData_names = ['Time lag (days)', 'Weight']                
                derivedData =[t, derivedData_tmp]

            plt.set_xlabel(axisHandle,'Time lag (days)')
            plt.set_ylabel(axisHandle,'Weight')            
            #box(axisHandle,'on')
            
            return derivedData, derivedData_names
            
        
    # delete class destructor
    #
    # Syntax:
    #   delete(self)
    #
    # Description:
    #   Loops through parameters and, if not an selfect, empties them. Else, calls
    #   the sub-selfect's destructor.
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
