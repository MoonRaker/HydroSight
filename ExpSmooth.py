
import numpy as np


#ExpSmooth Summary of this class goes here
#   Detailed explanation goes here
class ExpSmooth():
    
    
    def __init__(self, bore_ID, obsHead, forcingData_data,  forcingData_colnames, siteCoordinates, varargin):
            # Set observed head data
            self.inputData.head = obsHead
            
            # Initialise parameters
            self.parameters.alpha = -1            # Exponential smoothing parameter            
            self.parameters.beta  = -3            # Noise model parameter            
            self.parameters.gamma = -1            # Exponential trend parameter (optional)
                        
            # Set Parameter names
            self.variables.param_names = ['Auto-regressive', 'alpha''Moving-average', 'beta''Auto-regressive', 'gamma']
            
        
        def solve(self, time_points, tor_min, tor_max):
            # Check that the model has first been calibrated.
            #if ~isfield(self.variables, 'meanHead_calib') || ~isfield(self.parameters,'Const')
            if ~isfield(self.variables, 'meanHead_calib'):
                print 'The model does not appear to have first been calibrated. Please calibrate the model before running a simulation.'
            
            # Check the time points are all unique
            if length(unique(time_points))!=len(time_points):
                print 'The time points for simulation must be unique.'
            
            # Create logical matrix indicating if the time point is an
            # observation. If true, then the data point is used to update
            # the exponential smoothing.  Else, a forecast is made using
            # the exonential smoothing using the smooth terms from the 
            # previous observation. 
            # To calculate this vector, the following steps are undertaken:
            # 1. Unique time points are derived from the simulation time
            # points and the observed time points within the calibration
            # period.
            # 2. Find the time points within the unique list that are
            # observations.
            # 3. Create a logical vector with the time points from 2 as true
            # 4. Assign vector from 3 to the selfect for access within the
            # objective def.
            time_points_all = time_points self.variables.calibration_time_points
            time_points_all = np.unique(np.sort(time_points_all))
            dummy, ind = np.intersect(time_points_all, self.variables.calibration_time_points)
            self.variables.isObsTimePoints = False(np.shape(time_points_all))
            self.variables.isObsTimePoints[ind] = True                                    
            
            # Create vector of the time steps for only the time points with
            # observed heads.
            self.variables.delta_t = np.diff(time_points_all(self.variables.isObsTimePoints)) ./ 365.
            self.variables.meanDelta_t = np.mean(self.variables.delta_t)
            
            # Convert logical to double for MEX input
            self.variables.isObsTimePoints = np.double(self.variables.isObsTimePoints)
            
            # Set percentile for noise 
            Pnoise = 0.95            
            
            # Calc deterministic component of the head at 'time_points_all'.
            params = getParameters(self)
            self.variables.doingCalibration = False
            dummy, headtmp, self.variables.h_forecast = objectivedef(self, params, time_points_all)            
            
            # Filter 'head' to only those time points input to the def.
            dummy, ind = np.intersect(time_points_all, time_points)            
            headtmp = [time_points, headtmp[ind,:]]            
                        
            if np.shape(params)[2]>1:
                head  = np.zeros([np.shape(headtmp)[1], np.shape(headtmp)[2], np.shape(params)[2]])
                noise = np.zeros([np.shape(headtmp)[1], 3, np.shape(params)[2]])
                head[:,:,1] = headtmp
                for ii in range(np.shape(params)[2]):
                    
                    # Calc deterministic component of the head at 'time_points_all'.
                    params = getParameters(self)
                    self.variables.doingCalibration = False
                    dummy, headtmp, self.variables.h_forecast = objectivedef(self, params[:,ii], time_points_all)            

                    # Filter 'head' to only those time points input to the
                    # def.
                    dummy, ind = np.intersect(time_points_all, time_points)            
                    head[:,:,ii] = [time_points, head[ind,:]]                                
                
                    # Create noise component output.
                    if isfield(self.variables, 'sigma_n'):
                        noise[:,:,ii] = [head[:,1,ii],  np.ones([np.shape(head)[1], 2]) .* np.norminv(Pnoise, 0, 1) .* self.variables.sigma_n[ii]]
                    else:
                        noise[:,:,ii] = [head[:,1,ii], np.zeros([np.shape(head)[1], 2])]
            else:
                head = headtmp                
                
                # Create noise component output.
                if isfield(self.variables, 'sigma_n'):
                    noise[:,:] = [head[:,1], np.norminv(Pnoise, 0, 1) .* self.variables.sigma_n(np.ones([np.shape(head)[1], 2)])]
                else:
                    noise[:,:] = [head[:,1], np.zeros([np.shape(head)[1], 2])]
                        
            # Assign column names
            colnames = ['time', 'h_star']
            
            return head, colnames, noise
        
        
        def calibration_initialise(self, t_start, t_end):
            
            # Extract time steps
            t_filt = find((self.inputData.head[:,1]>=t_start) & (self.inputData.head[:,1]<=t_end))   
            self.variables.calibration_time_points = self.inputData.head[t_filt, 1]
            time_points = self.variables.calibration_time_points
            
            # Store time difference
            self.variables.delta_t = np.diff(time_points, axis=1) ./ 365.
            self.variables.meanDelta_t = np.mean(self.variables.delta_t)
            
            # Set a flag to indicate that calibration is complete.
            self.variables.doingCalibration = True
            
            # Calculate the mean head during the calibration period
            self.variables.meanHead_calib = np.mean(self.inputData.head[t_filt, 2])

            # Create logical matrix indicating if the time point is an
            # observation. If true, then the data point is used to update
            # the exponential smoothing.  Else, a forecast is made using
            # the exonential smoothing using the smooth terms from the 
            # previous observation. For the calibration, this vector is
            # true.
            self.variables.isObsTimePoints = True(np.shape(self.variables.calibration_time_points))
            
            # Convert logical to double for MEX input
            self.variables.isObsTimePoints = np.double(self.variables.isObsTimePoints)
            
            # Estimate the initial value and slope by fitting a smoothed spline.
            #spline_model = smooth(time_points, self.inputData.head(t_filt,2), 'rloess')
            spline_time_points = [0, 1, (time_points[2:end]-time_points[1])] ./ 365.
            spline_vals = csaps((time_points-time_points[1]) ./ 365, self.inputData.head[t_filt, 2]-self.variables.meanHead_calib, 0.1, spline_time_points)
            self.variables.initialHead = spline_vals[1]
            self.variables.initialTrend = (spline_vals[2]-spline_vals[1]) ./ (spline_time_points[2]-spline_time_points[1])
            
            # Find upper limit for alpha above which numerical errors arise.
            beta_upperLimit = 1000.
            delta_time = np.diff(self.inputData.head[:,1],1) ./ 365.
            while ((np.abs(np.sum(np.exp(-2.*beta_upperLimit .* delta_time)))<eps) | 
                   (np.exp(np.mean(np.log(1.-np.exp(-2.*beta_upperLimit .* delta_time))))< eps)):
                beta_upperLimit -= 0.1
                if beta_upperLimit<=eps:
                    break
            if beta_upperLimit<=eps:
                self.variables.beta_upperLimit = 3.
            else:
                # Transform alpha log10 space.
                self.variables.beta_upperLimit = np.log10(beta_upperLimit)
            
            # Find lower limit for alpha. This is based on the assumption
            # that if the best model (ie explains 95# of variance) then 
            # for such a model the estimate of the
            # moving average noise should be << than the observed
            # head standard dev..                        
            if ~np.isinf(self.variables.beta_upperLimit):
                beta_lowerLimit = -10.
                obsHead_std = np.std(self.inputData.head[t_filt, 2])   
                while np.sqrt(np.mean(0.05 .* obsHead_std**2. ./ (1.-np.exp(-2. .* 10.**beta_lowerLimit .* self.variables.delta_t))))>obsHead_std:
                    if beta_lowerLimit>=(self.variables.beta_upperLimit-2.):
                        break
                    beta_lowerLimit += 0.1
                self.variables.beta_lowerLimit = beta_lowerLimit
            else:
                self.variables.beta_lowerLimit = -5.           
            
            # Assign initial params to outputs
            params_initial = getParameters(self)
            
            # Clear estimate of constant
            if isfield(self.parameters, 'Const'):
                self.parameters.Const = []
            
            return params_initial, time_points
        
        
        def calibration_finalise(self, params, useLikelihood):
            
            # Re-calc objective def and deterministic component of the head and innocations.
            # Importantly, the drainage elevation (ie the constant term for
            # the regression) is calculated within 'objectivedef' and
            # assigned to the object. When the calibrated model is solved
            # for a different time period then this
            # drainage value will be used by 'objectivedef'.
            self.variables.objFn, self.variables.h_star, self.variables.h_forecast = objectivedef(params, self.variables.calibration_time_points, self, useLikelihood)                        
 
            # Re-calc residuals and assign to object                        
            t_filt = ((self.inputData.head[:,1]>=self.variables.calibration_time_points[1]) & 
                      (self.inputData.head[:,1]<=self.variables.calibration_time_points[end])                           
            self.variables.resid = self.inputData.head[t_filt, 2]-self.variables.h_forecast
            
            # Calculate mean of noise. This should be zero +/- eps
            # because the drainage value is approximated assuming n-bar = 0.
            self.variables.n_bar = np.real(np.mean(self.variables.resid))
            
            # Calculate innovations
            innov = self.variables.resid[2:end]-self.variables.resid[1:end-1] .* np.exp(-10.**self.parameters.beta .* self.variables.delta_t)            
            
            # Calculate noise standard deviation.
            self.variables.sigma_n = np.sqrt(np.mean(innov .** 2. ./ (1.-np.exp(-2. .* 10.**self.parameters.beta .* self.variables.delta_t))))
                        
            # Set a flag to indicate that calibration is complete.
            self.variables.doingCalibration = False
            
            return
            
        
        def objectivedef(params, time_points, self, getLikelihood):

            # Check the required object variables are set.
            if (~isfield(self.variables,'isObsTimePoints')) | (~isfield(self.variables,'meanHead_calib')):
                print 'The model does not appear to have been initialised for calibration.'
               
            # Check the input time points and the logical vector denoting
            # observation time points are of equal length
            if len(time_points)!=len(self.variables.isObsTimePoints):
                print 'The input time points vector must be the same length as the logical vector of observation time points.'
                        
            # Set model parameters         
            setParameters(self, params, ['alpha','beta','gamma'])            
            
            # Get model parameters and transform them from log10 space.
            alpha = 10.**self.parameters.alpha            
            beta  = 10.**self.parameters.beta
            gamma = 10.**self.parameters.gamma
            
            # Get time points
            t = self.inputData.head[:,1]

            # Create filter for time points and apply to t
            t_filt = find((t>=time_points[1]) & (t<=time_points[end]))                   
            
            # Setup time-varying weigting terms from input parameters
            # and their values at t0 using ONLY the time points with
            # observed head values!
            q = np.mean(self.variables.delta_t)
            alpha_i = 1.-(1.-alpha)**q
            gamma_i = 1.-(1.-gamma)**q
            
            # Initialise the output and subract the mean from the observed
            # head.
            h_ar = np.zeros([len(time_points), 1])
            h_obs = self.inputData.head[:,2]-self.variables.meanHead_calib            
            
            # Assign linear regression estimate the initial slope and intercept.
            #h_trend(1) = self.variables.initialTrend_calib
            #h_trend = self.variables.initialTrend_calib
            h_trend = self.variables.initialTrend
            h_ar[1] = self.variables.initialHead
            h_forecast = h_ar
            indPrevObsTimePoint = 1
            indPrevObs = 1
            
            # Undertake double exponential smoothing.
            # Note: It is based on Cipra T. and Hanzák T. (2008). Exponential 
            # smoothing for irregular time series. Kybernetika,  44(3), 385-399.
            # Importantly, this can be undertaken for time-points that have observed 
            # heads and those that are to be interpolated (or extrapolated). 
            # The for loop cycles though each time point to be simulated.
            # If the time point is an observation then the alpha, gamma and
            # h_trend are updates and an index is stored pointng to the last true obs point. 
            # If the time point is not an observation then a forecast is
            # undertaken using the most recent values of alpha, gamma and
            # h_trend
            for i in range(2, len(time_points)):
               
               delta_t = (time_points[i]-time_points[indPrevObsTimePoint])/365.
                
               if self.variables.isObsTimePoints(i)
                   # Update smoothing model using the observation at the current time point.
                   if indPrevObs==1:
                       gamma_weight = (1.-gamma)**delta_t
                   else:                                              
                       gamma_weight = delta_t_prev/delta_t*(1.-gamma)**delta_t            
                   alpha_weight = (1.-alpha) .** delta_t
                   alpha_i = alpha_i ./ (alpha_weight+alpha_i) 
                   gamma_i = gamma_i ./ (gamma_weight+gamma_i)
                   h_ar[i] = (1.-alpha_i)*(h_ar[indPrevObsTimePoint]+delta_t*h_trend)+alpha_i*h_obs[indPrevObs+1])
                   h_forecast[i] = h_ar[indPrevObsTimePoint]+delta_t*h_trend
                   h_trend = (1.-gamma_i)*h_trend+gamma_i*(h_ar[i]-h_ar[indPrevObsTimePoint]) ./ delta_t
                                                         
                   indPrevObsTimePoint = i
                   indPrevObs += 1
                   delta_t_prev = delta_t
               else:
                   # Make a forecast to the current non-observation time point.
                   h_forecast[i] = h_ar[indPrevObsTimePoint]+delta_t*h_trend
                   h_ar[i] = h_forecast[i]
            
            # Add the mean head onto the smoothed estimate and forecast
            h_ar += self.variables.meanHead_calib
            h_forecast += self.variables.meanHead_calib
             
            # Assign output for non-corrected head
            h_star = h_ar
            
            # Calculate the moving average componant.                        
            if ~self.variables.doingCalibration:
                objFn = []
                return

            # Create natrix ob observed and forecast heads
            h_mat = [self.inputData.head[t_filt, 2], h_forecast]

            # Calculate residuals (excluding input outliers).
            resid = np.diff(h_mat, 1, 2)

            # Calculate the innovations
            innov = resid[2:end, :]-resid[1:end-1, :] .* np.exp(bsxfun(@times, -beta, self.variables.delta_t))

            # Calculate the weighted least squares objective def
            selfFn = (np.sum(bsxfun(@rdivide, np.exp(np.mean(np.log(1.-np.exp(bsxfun(@times, -2.*beta ,self.variables.delta_t))))),
                      (1.-np.exp(bsxfun(@times, -2.*beta, self.variables.delta_t)))) .* innov .** 2.)   
                    
            # Calculate log likelihood    
            if getLikelihood:
                N = np.shape(resid)[1]
                objFn = -0.5*N*(np.log(2.*np.pi)+np.log(objFn ./ N)+1.) 
            
            return selfFn, h_star, h_forecast
            
        
        def setParameters(self, params, param_names):
            self.parameters.param_names[1] = params[1]
            self.parameters.param_names[2] = params[2]
            self.parameters.param_names[3] = params[3]
        
        
        def getParameters(self):
            params[1,:] = self.parameters.alpha
            params[2,:] = self.parameters.beta
            params[3,:] = self.parameters.gamma
            param_names = ['alpha', 'beta', 'gamma']
            return params, param_names
        
        
        def getParameters_physicalLimit(self):
            if isfield(self.variables, 'beta_upperLimit'):
                params_upperLimit = [0., self.variables.beta_upperLimit, 0.]
            else:
                params_upperLimit = [0 5  0]
            if isfield(self.variables, 'beta_lowerLimit'):
                params_lowerLimit = [-np.inf, self.variables.beta_lowerLimit, -np.inf]
            else:
                params_lowerLimit = [-np.inf, -5., -np.inf]
            return params_upperLimit, params_lowerLimit 
        
        
        def getParameters_plausibleLimit(self):
            params_upperLimit = [ 0.,  2.,  0.]
            params_lowerLimit = [-3., -2., -3.]
            return params_upperLimit, params_lowerLimit
        
        
        def getParameterValidity(self, params, time_points):
            # Get physical limits and test if parames are within the range
            params_upperLimit, params_lowerLimit = getParameters_physicalLimit(self)
            isValidParameter = (params>=repmat(params_lowerLimit, 1, np.shape(params)[2])) & (params<=repmat(params_upperLimit, 1, np.shape(params)[2]))
                        
            # Check the alphanoise parameter is large enough not to cause numerical
            # problems in the calcuation of the objective def.            
            filt_noiseErr = ((np.exp(np.mean(np.log(1.-np.exp(bsxfun(@times, -2 .* 10.**params[2,:], self.variables.delta_t))), 1))<=eps) |
                             (np.abs(np.sum(np.exp(bsxfun(@times, -2 .* 10.**params[2,:], self.variables.delta_t)), 1))<eps)                
            isValidParameter[2, filt_noiseErr]= False                
            return isValidParameter 
       
       
        # Get the forcing data from the model
        def getForcingData(self):
            forcingData = []
            forcingData_colnames = [] 
            return forcingData, forcingData_colnames
        
        
        # Set the forcing data from the model
        #def setForcingData(self, forcingData, forcingData_colnames):
        #    # do nothing. The model does not use forcing data.
            

        # Get the obs head data from the model
        def getObservedHead(self):
            head = self.inputData.head
            return head 
        
        
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
            if isempty(self.propNames[i]):
               continue
            if isobject(self.propNames[i]):
                delete(self.propNames[i])
            else:               
                self.propNames[i] = [] 
