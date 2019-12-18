

class(responseFunction_abstract):
    #RESPONSEFUNCTION_ABSTRACT Summary of this class goes here
    #Detailed explanation goes here
                  
    def __init__(self):
        setParameters(self, params)
            
        self.params, self.param_names = getParameters(self)
            
        # Check if each parameters is valid. This is primarily used to 
        # check parameters are within the physical bounds and sensible (eg 0<specific yield<1)
        self.isValidParameter = getParameterValidity(self, params, param_names)

        # Calculate impulse-response function.
        self.result = theta(self, t)
        
        # Calculate integral of impulse-response function from t to inf.
        # This is used to minimise the impact from a finit forcign data set.
        self.result = intTheta_upperTail2Inf(self, t)           

        # Calculate integral of impulse-response function from t to inf.
        # This is used to minimise the impact from a finit forcign data set.
        self.result = intTheta_lowerTail(self, t)           
                
        # Return pre-set physical limits to the function parameters.
        self.params_upperLimit, self.params_lowerLimit = getParameters_physicalLimit(self, param_name)
        
        # Transform the result of the response function multiplied by the
        # forcing. This method was included so that the groundwater pumping 
        # transfer function of Shapoori, Peterson, Western and Costelleo
        # 2013 could be corrected from a confined aquifer to an unconfined
        # aquifer. The method is called from model_IRF.get_h_star
        # Peterson Feb 2013
        self.result = transform_h_star(self, h_star_est)
