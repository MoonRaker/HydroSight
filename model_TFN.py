
import numpy as np


class model_TFN:


    # Class definition for Transfer def Noise (TFN) model for use with HydroSight
    #
    # Description: 
    #   The class defines a transfer def Noise (TFN) model for
    #   simulating time series of groundwater head. The model should be
    #   defined, calibrated and solved using HydroSight() or 
    #   the graphical user interface.
    #
    #   This class uses an object-oriented structure to provide a highly flexible 
    #   model structure and the efficient inclusion of new structural
    #   components. Currently, linear and nonlinear TFN model can be built.
    #   The linear models are based upon von Asmuth et al. 2002 and von Asmuth et al. 
    #   2008 and the nonlinear models use a nonlinear transform of the input
    #   climate data to account for runoff and nonlinear free drainage (see
    #   Peterson & Western, 2014) 
    #
    #   Additional key features of the model include:
    #
    #       - long term historic daily climate forcing is used to estimate the
    #       groundwater head at each time point via use of continuous transfer
    #       def. 
    #
    #       - the non-linear response of groundwater head to climate forcing
    #       can be accounted for by transforming the input forcing data.
    #       Currently, a simple 1-D vertically integrated soil model is
    #       available (see climateTransform_soilMoistureModels), allowing the
    #       construction of a soil model with between 1 and 7 parameters.
    #
    #       - Pumping drawdown can be simulated using the Ferris-Knowles 
    #       response def for confined aquifers, the Ferris-
    #       Knowles-Jacobs def for an unconfined aquifer
    #       or the Hantush def for a leaky aquifer. Each of these
    #       defs allow multiple production bores to be accounted for and
    #       each can have recharge or no-flow boundary conditions.
    #
    #       - The influence of streamflow can be approximated using the Bruggeman 
    #       reponse def. Importantly, this def is a prototype and
    #       should be used with caution. It is likely that the covariance
    #       between rainfall and streamflow complicates the estimation of the
    #       impacts of streamflow.
    #
    #       - a model can be fit to irregularly sampled groundwater head
    #       observations via use of an exponential noise def.
    #
    #       - new data transfer defs can easily be defined simply
    #       by the creation of new response def class definitions.
    #
    #       - calibration of the model is not to the observed head, but to the
    #       innovations between time steps. That is, the residual between the
    #       observed and modelled head is derived and the innovation is
    #       calculated as the prior residual minus the later residual
    #       multiplied by the exponental noise def.
    #
    #       - the contribution of a given driver to the head is calculated as
    #       the integral of the historic forcing times the weighting def.
    #       This is undertaken by get_h_star() and the mex .c compiled def
    #       doIRFconvolution(). If the forcing is instantaneous (for example a
    #       flux a soil moisture model) then Simpon's integration is
    #       undertaken. However, if the forcing is a daily integral (such as
    #       precipitation or daily pumping volumes) then daily trapazoidal
    #       integration of the weighting def is undertaken and then 
    #       multiplied by the daily flux.
    #
    #   Below are links to the details of the public methods of this class. See
    #   the 'model_TFN' constructor for details of how to build a model.
    #
    # See also:
    #   HydroSight: time_series_model_calibration_and_construction
    #   model_TFN: model_construction
    #   calibration_finalise: initialisation_of_model_prior_to_calibration
    #   calibration_initialise: initialisation_of_model_prior_to_calibration
    #   get_h_star: main_method_for_calculating_the_head_contributions.
    #   getParameters: returns_a_vector_of_parameter_values_and_names
    #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration
    #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
    #   solve: solve_the_model_at_user_input_sime_points
    #
    # Dependencies:
    #   responsedef_Pearsons.m
    #   responsedef_Hantush.m
    #   responsedef_Bruggeman.m
    #   responsedef_FerrisKnowles.m
    #   responsedef_Hantush.m
    #   responsedef_FerrisKnowlesJacobs.m
    #   derivedweighting_PearsonsNegativeRescaled.m
    #   derivedweighting_PearsonsPositiveRescaled.m
    #   climateTransform_soilMoistureModels.m
    #
    # References:
    #   von Asmuth J. R., Bierkens M. F. P., Mass K., 2002, Transfer
    #   dunction-noise modeling in continuous time using predefined impulse
    #   response defs.
    #
    #   von Asmuth J. R., Mass K., Bakker M., Peterson J., 2008, Modeling time
    #   series of ground water head fluctuations subject to multiple stresses.
    #   Groundwater, 46(1), pp30-40.
    #
    #   Peterson and Western (2014), Nonlinear time-series modeling of unconfined
    #   groundwater head, Water Resources Research, DOI: 10.1002/2013WR014800
    #
    # Author: 
    #   Dr. Tim Peterson, The Department of Infrastructure Engineering, 
    #   The University of Melbourne.
    #
    # Date:
    #   26 Sept 2014


    def __init__(self, bore_ID, obsHead, forcingData_data,  forcingData_colnames, siteCoordinates, varargin):
    
        # model_TFN constructs for linear and nonlinear Transfer def Noise model
        #
        # Syntax:
        #   model = model_TFN(bore_ID, obsHead, forcingData_data, ...
        #           forcingData_colnames, siteCoordinates, modelOptions)
        #
        # Description:
        #   Builds the model_TFN object, declares initial parameters and sets the
        #   observation and forcing data within the object. This method should be 
        #   called from HydroSight.
        #
        #   A wide range of models can be built within this method. See the inputs
        #   section below for details of the model components that can be built and 
        #   Peterson and Western (2014).
        #
        # Inputs:
        #   obsHead_input - matrix of observed groundwater elevation data. The
        #   matrix must comprise of two columns: date (as a double data type) 
        #   and water level elevation. The data must be of a daily time-step or 
        #   larger and can be non-uniform.
        #
        #   forcingData_data - M x N numeric matrix of daily foring data. The matrix 
        #   must comprise of the following columns: year, month, day and
        #   observation data. The observation data must include precipitation but
        #   may include other data as defined within the input model options. 
        #   The record must be continuous (no skipped days) and cannot contain blank 
        #   observations.
        #
        #   forcingData_colnames - 1 x N cell array of the column names within the
        #   input numerical array "forcingData_data". The first three columns must
        #   be: year month, day. The latter columns must define the forcing data
        #   column names. Imprtantly, each forcing data column name must be listed
        #   within the input "siteCoordinates". 
        #
        #   siteCoordinates - Mx3 cell matrix containing the following columns :
        #   site ID, projected easting, projected northing. Importantly,
        #   the input 'bore_ID' must be listed and all columns within
        #   the forcing data (excluding the year, month,day). When each response 
        #   def is created, the coordinates for the observation bore and the 
        #   corresponding forcing are are passed to the response def. This is
        #   required for simulating the impacts of pumping. Additionally,
        #   coordinates not listed in the input 'forcingData' can be input.
        #   This provides a means for image wells to be input (via the input 
        #   modelOptions).
        #
        #   modelOptions - cell matrix defining the model components, ie how the
        #   model should be constructed. The cell matrix must be three columns,
        #   where the first, second and third columns define the model component
        #   type to be set, the component property to be set, and the value to be
        #   set to the component property. Each model component can have the following
        #   property inputs ((i.e. 2nd column):
        #
        #       - 'weightingdef': Property value is to define the component type.
        #       Available component types include:
        #           - 'responsedef_Pearsons' (for climate forcing)
        #           - 'responsedef_Hantush' (for groundwater pumping)
        #           - 'responsedef_FerrisKnowles' (for groundwater pumping)
        #           - 'responsedef_Hantush' (for groundwater pumping)
        #           - 'responsedef_FerrisKnowlesJacobs' (for groundwater pumping)
        #           - 'responsedef_Bruggeman' (for streamflow)
        #           - 'derivedweighting_PearsonsNegativeRescaled' (see below)
        #           - 'derivedweighting_PearsonsPositiveRescaled' (see below)
        #            
        #       - 'forcingdata': Property value is the column number, column name 
        #       within forcingData (as input to HydroSight) or a cell 
        #       array containing the options for a forcing transformation
        #       object. The forcing transformation input must be a Nx2 cell array  
        #       with the following property (left column) and value (right column) 
        #       settings:
        #           - 'transformdef' property and the transformation class name 
        #           e.g. 'climateTransform_soilMoistureModels'
        #
        #           - 'forcingdata' property and a Nx2 cell array declaring the 
        #           input forcing variable required by the transformation def
        #           (left column) and the name of the input forcing data or column
        #           number. eg   ['precip',3'et',4].
        #
        #           -  'outputdata' property and the output flux to be output from
        #           the forcing transformation def. 
        #
        #           - 'options' propert and an input entirely dependent upon the
        #           forcing transformation def. For
        #           'climateTransform_soilMoistureModels' this is used to define
        #           the form of the soil moisture model i.e. which ODE parameters
        #           to fix and which to calibrate. See the documenttation for each
        #           transformation def for details.
        #
        #       - 'inputcomponent': Property value is the name of another model
        #       component (i.e.  the first column value). This option allows,
        #       for example, the parameterised Pearson's weighting def for 
        #       a precipitation component to also be used for, say, an ET
        #       component. When building such a model, the weighting def
        #       should be a derived def such as
        #       'derivedweighting_PearsonsNegativeRescaled'
        #
        #       - 'options': cell array options specific for a weighting def.
        #       The structure of the weighting def options are entirely
        #       dependent upon the weighting def. Currently, only the
        #       groundwater pumping weighting defs allow the input of options.
        #       These pumping options allow the simulation of multiple recharge or
        #       no-flow image wells per pumping bore. These options must be a Nx3 cell
        #       array. The first column must be a string for the site ID of a
        #       production bore. The second column must be the site of an image
        #       well. Coordinate for both the production bore ID and the image
        #       well must be listed within the input coordinated cell array.
        #       The third column gives the type of image well. The availabes
        #       type are: "recharge" for say a river and "no flow" for say a
        #       aquifer no-flow boundary.
        #
        # Outputs:
        #   model - model_TFN class object 
        #
        # Example:
        #   see HydroSight: time_series_model_calibration_and_construction
        #
        # See also:
        #   HydroSight: time_series_model_calibration_and_construction
        #   calibration_finalise: initialisation_of_model_prior_to_calibration
        #   calibration_initialise: initialisation_of_model_prior_to_calibration
        #   get_h_star: main_method_for_calculating_the_head_contributions.
        #   getParameters: returns_a_vector_of_parameter_values_and_names
        #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration
        #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
        #   solve: solve_the_model_at_user_input_sime_points
        #
        # References:
        #   Peterson and Western (2014), Nonlinear time-series modeling of unconfined
        #   groundwater head, Water Resources Research, DOI: 10.1002/2013WR014800
        #
        # Author: 
        #   Dr. Tim Peterson, The Department of Infrastructure
        #   Engineering, The University of Melbourne.
        #
        # Date:
        #   26 Sept 2014
                
        # Expand cell structure within varargin. It is required to
        # contain the model options

        varargin = varargin[1]
        
        # Check there are an even number model type options, ie for
        # option option type there is an option value.
        if len(varargin)<=0:
           print 'Default model of only precipitation forcing has been adopted because no model components were input'
        elif np.shape(varargin, 2)!=3:
           print 'Invalid model options. The inputs must of three columns in the following order: model_component, property, property_value'
        
        # Check the forcing data does not contain nans or infs
        if any(any(np.isnan(forcingData_data) | np.isinf(forcingData_data))):
            ind = find(any(np.isnan(forcingData_data) | np.isinf(forcingData_data),2), 1, 'first')
            print 'Forcing data contains empty, np.nan or inf values (line ', str(ind), ')'
        
        # Check that forcing data exists before the first head
        # observation.
        if np.floor(np.min(forcingData_data[:,1]))>=np.floor(np.min(obsHead[:,1]))
            print 'Forcing data must be input prior to the first water level observation. Ideally, the forcing data should start some years prior to the head.'
                        
        # Set bore ID, observed head and forcing data and site coordinates
        self.inputData.bore_ID = bore_ID            
        self.inputData.head = obsHead            
        self.inputData.siteCoordinates = siteCoordinates
        setForcingData(self, forcingData_data, forcingData_colnames)
        
        # Cycle though the model components and their option and
        # declare instances of objects were appropriate.
        valid_properties = ['weightingdef', 'forcingdata', 'options', 'inputcomponent']            
        valid_transformProperties = ['transformdef', 'forcingdata', 'outputdata','options', 'inputcomponent']
        
        # Check the model options.
        for ii in range(np.shape(varargin)[1]):                
        
            modelComponent = varargin[ii, 1]
            propertyType  = varargin[ii, 2]                
            propertyValue = varargin[ii, 3]            

            # Check the component properties is valid.
            ind = find(strncmpi(propertyType, valid_properties, len(propertyType)))
            if len(ind)!=1:
                print 'Invalid component property type: ', str(propertyType)

            # Check the component properties are valid.
            #----------------------------------------------------------
            if propertyType==valid_properties[1]:
                # If the propertyType is weightingdef, check it is a valid object.
                if np.isempty(meta.class.fromName(propertyValue)):
                    # Report error if the input is not a class definition.
                    print 'The following weighting def name does not appear to be a class object:', propertyValue

            elif strcmp(propertyType, valid_properties[2]):
                # If the propertyType is forcingdata, check if the
                # input is cell, character (ie forcing site name) or an
                # integer(s) (ie forcing data column number). If the
                # input is a cell, then check if it is a list of
                # forcing site names or a def name and inputs to
                # transform forcing data (eg using a soil moisture
                # model).
                if np.iscell(propertyValue):
                    # If 'propertyValue' is a vector then it should be a list of forcing sites.
                    if np.isvector(propertyValue):
                        # Check each forcing site name is valid.
                        for j in range(len(propertyValue)):
                            if (np.isnumeric(propertyValue[j])) & ((propertyValue[j]<0.) | (propertyValue[j]+2.>np.shape(forcingData_data)[2]):
                                print 'Invalid forcing column number for component:', modelComponent,'. It must be an intereger >0 and <= the number of forcing data columns.'
                            elif np.ischar(propertyValue[j]):
                                # Check the forcing column name is valid.
                                filt = cellfun(@(x)(strcmp(x,propertyValue[j])),forcingData_colnames)                                
                                if np.isempty(find[filt]):
                                    print 'Invalid forcing column name for component:', modelComponent,'. Each forcing name must be listed within the input forcing data column names.'

                    else:
                        # Check that the input is of the expected
                        # format for creating a class to derive the
                        # forcing data. The Expected format is an Nx2
                        # matrix (where N<=5) and the first colum
                        # should have rows for 'tranformationdef'
                        # and 'outputvariable'.
                        if (np.shape(propertyValue)[2]!=2) | (np.shape(propertyValue)[1]<2) | (np.shape(propertyValue)[1]>len(valid_transformProperties)):
                            print 'Invalid inputs for transforming the forcing data for component:', modelComponent,'. The input cell array must have >=2 and <=5 rows and 2 columns.'
                        
                        # Check it has the required minimum inputs to the first columns.
                        if (np.sum(propertyValue[:,1]==valid_transformProperties[1])!=1) | (np.sum(propertyValue[:,1]==valid_transformProperties[3])!=1):
                            print 'Invalid inputs for transforming the forcing data for component:', modelComponent,'. The first column of the input cell array must have rows for "transformdef" and "outputdata".' 
                        
                        # Check that the first column only contains
                        # inputs that are acceptable.
                        for j in range(np.shape(propertyValue)[1]):                                                                
                            if np.sum(valid_transformProperties==propertyValue[j,1])!=1):
                                print 'Invalid inputs for transforming the forcing data for component:', modelComponent,'. Read the help for "model_TFN" to see the accepted values.'
                        
                        # Check that the transformation def name
                        # is a valid class definition name.
                        filt = propertyValue[:,1]==valid_transformProperties[1]
                        if np.isempty(meta.class.fromName(propertyValue[filt,2])):
                            print 'The following forcing transform def name is not accessible:', propertyValue[filt,2]
                        
                        # Check that the forcing transformation
                        # def name is consistent with the abstract
                        # 'forcingTransform_abstract'.
                        try:
                            if (~any(strcmp(findAbstractName(propertyValue[filt,2]), 'forcingTransform_abstract'))) & (~any(strcmp(findAbstractName(propertyValue[filt,2]),'derivedForcingTransform_abstract'))):
                                print 'The following forcing transform def class definition is not derived from the "forcingTransform_abstract.m" abstract:',propertyValue[filt,2]])
                        except:
                            print '... Warning: Checking that the required abstract for the forcing transform class definition was used failed. This is may be because the version of matlab is pre-2014a.'
                        
                        # Get a list of required forcing inputs.
                        requiredForcingInputs, isOptionalInput = feval(strcat(propertyValue[filt,2], '.inputForcingData_required'), bore_ID, forcingData_data, forcingData_colnames, siteCoordinates)
                                                    
                        # Check that the input forcing cell array has the correct dimensions.
                        filt = propertyValue[:,1]==valid_transformProperties[2]
                        if any(filt):
                            if (((np.shape(propertyValue[filt,2])[2]!=2) | 
                                ((np.shape(propertyValue[filt,2])[1]!=len(requiredForcingInputs(~isOptionalInput))) & 
                                ((np.shape(propertyValue[filt,2])[1]!=len(requiredForcingInputs)) & 
                                ((np.shape(propertyValue[filt,2])[1]>1)):
                                print 'Invalid forcing data for the forcing transform def name for component:', modelComponent,'. It must be a cell array of two columns and ', str(len(requiredForcingInputs)), ' rows (one row for each required input forcing).'
                        
                            # Check that the required inputs are specified.
                            for j in range(len(requiredFocingInputs)):
                                if isOptionalInput[j]
                                    continue
                                if ~any(propertyValue[filt,2]==requiredForcingInputs[j]):
                                    print 'Invalid inputs for transforming the forcing data for component:', modelComponent,'. A forcing data site name must be specified for the following transform def input:', requiredForcingInputs[j]
                        
                        # Get a list of valid output forcing variable names.
                        filt = propertyValue[:,1]==valid_transformProperties[1]
                        optionalFocingOutputs = feval([propertyValue[filt,2],'.outputForcingdata_options'], bore_ID, forcingData_data,  forcingData_colnames, siteCoordinates)
                                                                                
                        # Check that the output variable a char or cell vector.
                        filt = propertyValue[:,1]==valid_transformProperties[3])                              
                        if np.isnumeric(propertyValue[filt,2]) | (np.iscell(propertyValue[filt,2]) & ~np.isvector(propertyValue[filt,2]))
                            print 'Invalid output forcing variable for the forcing transform def name for component:', modelComponent, ' . It must be a string or a cell vector of strings of valid output variable names.'
                        
                        # Check each output variable name is valid.
                        if np.ischar(propertyValue[filt,2]):
                            if ~any(strcmp(optionalFocingOutputs, propertyValue[filt,2])):
                                print 'Invalid output forcing variable for the forcing transform def name for component:', modelComponent,' . The output variable name must be equal to one of the listed output variables for the transformaion def.'

                        else:
                            for j=1:len(propertyValue[filt,2])
                                if ~any(strcmp(optionalFocingOutputs, propertyValue[filt,2][j])):
                                    print 'Invalid output forcing variable for the forcing transform def name for component:' , modelComponent, ' . Each output variable name must be equal to one of the listed output variables for the transformaion def.'
                        
                        # Get a list of np.unique model components that have transformed forcing data
                        # and the transformation model does not require
                        # the input of another component.
                        k = 0
                        for j in range(np.shape(varargin)[1]):
                            if (varargin[j,2]==valid_properties[2]) & (np.iscell(varargin[j,3])) & (~np.isvector(varargin[j,3])):
                        
                                # Check if the transformation model
                                # input as the required inputs and does
                                # not require the input of another
                                # component.                                    
                                if ((any(strcmp(varargin[j,3][:,1], valid_transformProperties[1]))) &
                                    (any(strcmp(varargin[j,3][:,1], valid_transformProperties[2]))) & 
                                    (any(strcmp(varargin[j,3][:,1], valid_transformProperties[3]))) & 
                                   (~any(strcmp(varargin[j,3][:,1], valid_transformProperties[5])))):
                                    k +=1
                                    # Record the component name.
                                    modelComponents_transformed[k,1] = varargin[j,1]
                                    
                                    # Add the transformation model name.
                                    filt = varargin[j,3][:,1][varargin[j,3][:,1]==valid_transformProperties[1]]
                                    modelComponents_transformed[k,2] = varargin[j,3][filt,2]

                        # Check that each complete transforamtion model
                        # is np.unique.
                        if np.shape(np.unique(modelComponents_transformed[:,2])[1])!=np.shape(modelComponents_transformed[:,2])[1]:
                            print 'Non-unique complete transformation models. Each complete transformation model must be input for only one model component. A complete model is that having at least the following inputs: transformdef, forcingdata, outputdata'
                        
                        # Get a np.unique list of model components with
                        # complete transforamtion models and a list of
                        # np.unique tranformation defs
                        modelComponents_transformed_uniqueComponents = np.unique(modelComponents_transformed[:,1])
                        modelComponents_transformed_uniqueTransFunc = np.unique(modelComponents_transformed[:,2])
                        
                        # Check the inputcomponent is a valid
                        # component.
                        filt = propertyValue[:,1]==valid_transformProperties[5]
                        if (any[filt]) & (np.ischar(propertyValue[filt,2])):
                            if ~any(strcmp(modelComponents_transformed_uniqueTransFunc[:,1], propertyValue[filt,2])):
                                print 'Invalid source def for the forcing transform def name for component:', modelComponent, ' . The input def must be a tranformation def name that is listed within the model.'

                        elif (any[filt]) & (np.iscell(propertyValue[filt,2])):
                            for j in range(len(propertyValue[filt,2])):
                                if ~any(strcmp(modelComponents_transformed_uniqueTransFunc[:,1], propertyValue[filt,2][j])):
                                    print 'Invalid source def for the forcing transform def name for component:', modelComponent, ' . The input def must be a tranformation def name that is listed within the model.'

                        elif any(filt):
                            print 'Invalid source def for the forcing transform def name for component:', modelComponent, ' . The input component must be a a tranformation def name (string or cell vector of strings) that is listed within the model input options and it must also have a def transformation.'
                    
                elif (~np.isnumeric(propertyValue)) & ((np.isnumeric(propertyValue)) & (any(propertyValue<0))) & (~ischar(propertyValue)):
                    print 'Invalid property value for property type ', valid_properties[2], '. It must be a forcing data site name or column number integer >0 (scalar or vector)'
                                        
            elif propertyType==valid_properties[4]:    # If the propertyType is inputcomponent, check the input is a string and a valid component name.
            
                if np.ischar(propertyValue):
                    # Check that it is valid component name.
                    filt = cellfun(@(x)(strcmp(x,propertyValue)), varargin[:,1])
                    if np.isempty[filt]:
                        print 'Invalid input component name for component:', modelComponent, '. The named input component must be a component within the model.'

                else:
                    print 'Invalid input component name for component:', modelComponent, '. It must be a string defining a component type.'

        # Record the source data column for each component or the
        # transformation model and output variables and input transformation component.
        filt = find(strcmpi(varargin[:,2],'forcingdata'))
        transformdefIndexes_noInputComponents = []
        transformdefIndexes_inputComponents = []
        for ii in filt:
           modelComponent = varargin[ii,1]
           if np.isnumeric(varargin[ii,3]):
              self.inputData.componentData.(modelComponent).dataColumn = varargin[ii,3]-2
           elif np.ischar(varargin[ii,3]):
              filt_forcingCols = cellfun(@(x,y)(strcmp(x,varargin[ii,3])), forcingData_colnames)
              self.inputData.componentData.(modelComponent).dataColumn = find(filt_forcingCols)
           elif np.iscell(varargin[ii,3]):
              # If it is a cell vector, then get the column number for
              # each element. Else, if it is a Nx2 cell array with the
              # required format, the component uses transformed forcing
              # data so store the transformation object name and the
              # required output variable names.
              if (np.isvector(varargin[ii,3])) & (~np.iscell(varargin[ii,3][1])):
                for j in range(np.shape(varargin[ii,3])[1]):
                    filt_forcingCols = cellfun(@(x,y)(strcmp(x,varargin[ii,3][j])), forcingData_colnames)
                    self.inputData.componentData.(modelComponent).dataColumn(j,1) = find(filt_forcingCols)                                    
              elif ((np.shape(varargin[ii,3],2)==2) &
                    (any(strcmp(varargin[ii,3][:,1], valid_transformProperties[1]))) &
                    (any(strcmp(varargin[ii,3][:,1], valid_transformProperties[3]))) |
                    (np.shape(varargin[ii,3][1])[[2]==2) & 
                    (any(strcmp(varargin[ii,3][1][:,1], valid_transformProperties[1]))) &
                    (any(strcmp(varargin[ii,3][1][:,1], valid_transformProperties[3])))):
          
                
                # Record the transformation model name and required
                # output forcing data from the transformation model for
                # input to the model component.
                if np.shape(varargin[ii,3])[2]==2:
                    filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[1])
                    self.inputData.componentData.(modelComponent).forcing_object = varargin[ii,3][filt,2]
                    filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[3])
                    self.inputData.componentData.(modelComponent).outputVariable = varargin[ii,3][filt,2]
                else:
                    for j in range(np.shape(varargin[ii,3])[1]):
                        filt = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[1])
                        self.inputData.componentData.(modelComponent).forcing_object = varargin[ii,3][j][filt,2]
                        filt = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[3])
                        self.inputData.componentData.(modelComponent).outputVariable[j,1] = varargin[ii,3][j][filt,2]                        
                # Record the input forcing data columns for the 
                # transforamtion model if provided. Else, seach the other
                # transforamtion components for the forcing data
                # columns.
                if np.shape(varargin[ii,3])[2]==2:
                    filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[2])                    
                    if any(filt):
                        self.inputData.componentData.(modelComponent).inputForcing = varargin[ii,3][filt,2]
                    else:

                        # Get the name of the source transforamtion
                        # model.
                        filt_transf = strcmp(varargin[ii,3][:,1], valid_transformProperties[1])
                        sourceTransformModel = varargin[ii,3](filt_transf,2)

                        # Create a filter for the forcingData property.
                        filt_transf = find(strcmp(varargin[:,2], valid_properties[2]))

                        # Loop through all forcingData inputs to find the
                        # transformationModel that is the source to the
                        # current transformation model.
                        for j in filt_transf:
                           if ((np.shape(varargin[j,3])[2]==2) & 
                               (np.iscell(varargin[j,3])) &   
                               (any(strcmp(varargin[j,3][:,1], valid_transformProperties[1]))) &
                               (any(strcmp(varargin[j,3][:,1], valid_transformProperties[2]))) &
                               (any(strcmp(varargin[j,3][:,1], valid_transformProperties[3]))) &
                               (any(strcmp(varargin[j,3][:,2], sourceTransformModel))):
                               # Get the row containing the input forcing data
                               filt = strcmp(varargin[j,3][:,1], valid_transformProperties[2])

                               # Assign source forcing names to the current
                               # model component.
                               self.inputData.componentData.(modelComponent).inputForcing = varargin[j,3][filt,2]                               

                else:
                    for j in range(np.shape(varargin[ii,3])[1]):  
                        filt = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[2])                    
                        if any(filt):
                            inputForcing = varargin[ii,3][j][filt,2]
                            filt = cellfun(@(x) ~np.isempty(x), inputForcing[:,2])
                            inputForcing = inputForcing[filt,:]                                
                            self.inputData.componentData.(modelComponent).inputForcing = inputForcing 
                        #else:
                        #    print 'Inheritated transformation models with multiple outputs is not operational.')
                        #    # Get the name of the source transforamtion
                        #    # model.
                        #    filt_transf = strcmp(varargin[ii,3][:,1], valid_transformProperties[1])
                        #    sourceTransformModel = varargin[ii,3](filt_transf,2)
 			#
                        #    # Create a filter for the forcingData property.
                        #    filt_transf = find(strcmp(varargin[:,2], valid_properties[2]))
 			#
                        #    # Loop through all forcingData inputs to find the
                        #    # transformationModel that is the source to the
                        #    # current transformation model.
                        #    for j=filt_transf'
                        #       if np.shape(varargin[j,3],2) == 2 ...
                        #       & np.iscell(varargin[j,3]) ...        
                        #       & any(strcmp(varargin[j,3][:,1], valid_transformProperties[1])) ...
                        #       & any(strcmp(varargin[j,3][:,1], valid_transformProperties[2])) ...
                        #       & any(strcmp(varargin[j,3][:,1], valid_transformProperties[3])) ...
                        #       & any(strcmp(varargin[j,3][:,2], sourceTransformModel))
                        #           # Get the row containing the input forcing data
                        #           filt = strcmp(varargin[j,3][:,1], valid_transformProperties[2])
 			#
                        #           # Assign source forcing names to the current
                        #           # model component.
                        #           self.inputData.componentData.(modelComponent).inputForcing = varargin[j,3][filt,2]                               

                # Check that forcing data was identified for the current component.
                if ~np.isfield(self.inputData.componentData.(modelComponent), 'inputForcing'):
                    print 'Invalid forcing data for component:', modelComponent, '. No forcing data names or column numbers were input or found by searching for complete transformation models of the same name as input for this component.'
                

                # Subtract '2' from the column numbers (if
                # not a string). This is done because when
                # the user inputs the column number the
                # first three columns are year,month,day.
                # However, when passed to model_TFN(), the
                # first three columns are replaced by a
                # single column of the date number.
                for j in range(np.shape(self.inputData.componentData.(modelComponent).inputForcing)[1]):
                    if np.isnumeric(self.inputData.componentData.(modelComponent).inputForcing[j,2]):
                        self.inputData.componentData.(modelComponent).inputForcing[j,2] = self.inputData.componentData.(modelComponent).inputForcing[j,2]-2
                
                # Record indexes to each component with a
                # transformation model and the type of model.
                # Additionally, if the transformation model uses inputs
                # from a transformation model, then also record this data.
                if np.shape(varargin[ii,3])[2])==2:
                    filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[5])
                    if any(filt):

                        # Only add an index for the creation the
                        # transformation object if input and output forcing
                        # is defined. If only output is defined, then the 
                        # transforamtion is assumed to be simply calling
                        # a transformation object created for another
                        # component.
                        filt_inputForcing = strcmp(varargin[ii,3][:,1], valid_transformProperties[2])
                        filt_outputForcing = strcmp(varargin[ii,3][:,1], valid_transformProperties[3])
                        if (any(filt_inputForcing)) & (any(filt_outputForcing)):
                            self.inputData.componentData.(modelComponent).isForcingModel2BeRun = True
                            self.inputData.componentData.(modelComponent).inputForcingComponent = varargin[ii,3][filt,2]
                            transformdefIndexes_inputComponents = transformdefIndexes_inputComponents, ii
                        else:
                            print 'Invalid forcing transformation for component:', modelComponent, '. When a transformation model uses another transformation model as an input then the forcing data and the output variable must be specified.'

                    else:
                        # Only add an index for the creation of the
                        # transformation object if input and output forcing
                        # is defined. If only output is defined, then the 
                        # transforamtion is assumed to be simply calling
                        # a transformation object created for another
                        # component.
                        filt_inputForcing = strcmp(varargin[ii,3][:,1], valid_transformProperties[2])
                        filt_outputForcing = strcmp(varargin[ii,3][:,1], valid_transformProperties[3])
                        if (any(filt_inputForcing)) & (any(filt_outputForcing)):
                            self.inputData.componentData.(modelComponent).isForcingModel2BeRun = True
                            transformdefIndexes_noInputComponents = transformdefIndexes_noInputComponents, ii                      
                        elif (~any(filt_inputForcing)) & (any(filt_outputForcing)):
                            self.inputData.componentData.(modelComponent).isForcingModel2BeRun = False

                else:
                    for j in range(np.shape(varargin[ii,3])[1]):  
                        filt = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[5])
                        if any(filt):

                            # Only add an index for the creation the
                            # transformation object if input and output forcing
                            # is defined. If only output is defined, then the 
                            # transforamtion is assumed to be simply calling
                            # a transformation object created for another
                            # component.
                            filt_inputForcing = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[2])
                            filt_outputForcing = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[3])
                            if (any(filt_inputForcing)) & (any(filt_outputForcing)):
                                self.inputData.componentData.(modelComponent).isForcingModel2BeRun = True
                                self.inputData.componentData.(modelComponent).inputForcingComponent = varargin[ii,3][j][filt,2]
                                transformdefIndexes_inputComponents = [transformdefIndexes_inputComponents ii]
                            else:
                                print 'Invalid forcing transformation for component:', modelComponent, '. When a transformation model uses another transformation model as an input then the forcing data and the output variable must be specified.'

                        else:

                            # Only add an index for the creation of the
                            # transformation object if input and output forcing
                            # is defined. If only output is defined, then the 
                            # transforamtion is assumed to be simply calling
                            # a transformation object created for another
                            # component.
                            filt_inputForcing = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[2])
                            filt_outputForcing = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[3])
                            if (any(filt_inputForcing)) & (any(filt_outputForcing)):
                                self.inputData.componentData.(modelComponent).isForcingModel2BeRun = True
                                transformdefIndexes_noInputComponents = transformdefIndexes_noInputComponents, ii                           
                            #elif ~any(filt_inputForcing) & any(filt_outputForcing)
                            #    self.inputData.componentData.(modelComponent).isForcingModel2BeRun = False
                
              else:
                print 'An unexpected error occured for component :', modelComponent, '. The input for "forcingdata" should be either a forcing data site name(s) or number(s) or a cell array for an input transformation model.'

           else:
              print 'The data colum for each model component must be either the column number within the forcing data matrix or a string of the forcing site name.'

        # Check that each component has a forcing column or a forcing 
        # transformation defined.
        modelComponent = np.unique(varargin[:,1])
        for i in range(len(modelComponent)):
            if ~np.isfield(self.inputData.componentData,modelComponent[i]):         
                print 'The input model options must specify a forcing data column or forcing transformation model for each model components. The following component has no such input:', modelComponent[i]

        # Get a list of rows in varargin defining the weighting
        # def names.
        weightingdefIndexes = find(strcmpi(varargin[:,2], valid_properties[1]))
                               
        # Build the objects for forcing transformations that DO NOT
        # require an input weighting def object.            
        for ii in transformdefIndexes_noInputComponents:
            
            modelComponent = varargin[ii,1]
            propertyValue  = varargin[ii,3]

            # Get the following inputs to build the model object:
            # transformation def name, input data column names,
            # the variables to which the input data is to be assigned
            # too, and additional transformation model options.
            if np.shape(varargin[ii,3])[2]==2:
                filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[1])
                transformobject_name = varargin[ii,3][filt,2]
                filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[2])
                transformobject_inputs = varargin[ii,3][filt,2]
                filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[4])
                if any(filt):
                    transformobject_options = varargin[ii,3][filt,2]
                else:
                    transformobject_options =[]

            else:
                for j in range(np.shape(varargin[ii,3])[1]):  
                    filt = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[1])
                    transformobject_name = varargin[ii,3][j][filt,2]
                    filt = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[2])
                    transformobject_inputs = varargin[ii,3][j][filt,2]
                    filt = strcmp(varargin[ii,3][j][:,1], valid_transformProperties[4])
                    if any(filt):
                        transformobject_options = varargin[ii,3][j][filt,2]
                        break

            # Get the names (or column numbers) of the input forcing data. 
            forcingData_inputs = self.inputData.componentData.(modelComponent).inputForcing[:,2]                            
            
            # Remove empty rows
            filt = cellfun(@(x) ~np.isempty(x) & ~strcmp(x, '(none)'), forcingData_inputs)
            forcingData_inputs = forcingData_inputs[filt]
            filt = cellfun(@(x) ~np.isempty(x) & ~strcmp(x, '(none)'), transformobject_inputs[:,2])
            transformobject_inputs = transformobject_inputs[filt,:]
            
            # Find the required columns in forcing data so that only the
            # required data is input.
            colNum = 1
            rowFilt = True(np.shape(transformobject_inputs)[1], 1)
            for j in range(np.shape(transformobject_inputs)[1]):
                # Skip if optional forcing input
                if transformobject_inputs[j,2]=='(none)':
                    rowFilt[j] = False
                    continue                        

                if np.isnumeric(transformobject_inputs[j,2]):
                    colNum = [colNum, transformobject_inputs[j,2]-2]                        
                elif np.ischar(transformobject_inputs[j,2]):
                    filt = find(strcmpi(forcingData_colnames, transformobject_inputs[j,2]))
                    if sum[filt]==0:
                        print 'An unexpected error occured for component :', modelComponent, '. Within the input cell array for "forcingdata", the second column contains a forcing site name that is not listed within the input forcing data column names.'
                    colNum = colNum, filt
                    
                else:
                    print 'An unexpected error occured for component :', modelComponent, '. Within the input cell array for "forcingdata", the second column contains a forcing data column number or site name that is listed within the input forcing data column names.'

                transformobject_inputs[j,2] = j+1
            
            # Filter transformed object input rows to remove those that
            # are input as options ie '(none)'.
            transformobject_inputs = transformobject_inputs[rowFilt,:]
            
            try:
                #colNum_extras = True(len(forcingData_colnames),1)
                #colNum_extras(colNum) = False
                #colNum_extras = find(colNum_extras)
                #forcingData_colnamesTmp = forcingData_colnames([colNumcolNum_extras])
                #forcingData_dataTmp = forcingData_data(:,[colNumcolNum_extras])
                self.parameters.(transformobject_name) = feval(transformobject_name, bore_ID, forcingData_data(:,colNum), forcingData_colnames(colNum), siteCoordinates, transformobject_inputs, transformobject_options)
                 
                # Check each output variable is valid. If not valid, then remove it. 
                # If valid then add coordinates for the output
                # variable. If there is no valid output variable then
                # throw an error.                    
                variableName = self.inputData.componentData.(modelComponent).outputVariable
                if np.ischar(variableName):
                    variableNameAsCell[1] = variableName
                    variableName = variableNameAsCell
                    clear variableNameAsCell

                isValidVariableName = False[1,len(variableName)]
                t = self.inputData.forcingData[:,1]
                setTransformedForcing(self.parameters.(transformobject_name), t, True)
                for j in range(len(variableName)):
                    try:                                                         
                        forcingData = getTransformedForcing(self.parameters.(transformobject_name), variableName[j])
                        isValidVariableName[j] = True
                    except:
                        isValidVariableName[j] = False

                if ~any(isValidVariableName):
                    print 'No transformed forcing data for component "',modelComponent,'" could be derived' 

                variableName = variableName(isValidVariableName)
                self.inputData.componentData.(modelComponent).outputVariable = variableName
                coordinatesTmp = getCoordinates(self.parameters.(transformobject_name), variableName)
                siteCoordinates(np.shape(siteCoordinates,1)+1: np.shape(siteCoordinates)[1] + np.shape(coordinatesTmp)[1], 1:3) = coordinatesTmp
                
            except:
                print 'ERROR: Invalid model component class object for forcing transform: ', transformobject_name
                #rethrow(exception)

                    
        # Build the objects for forcing transformations that DO
        # require an input weighting def object.
        # NOTE: if the component is derived from another component
        # object, then like the construction of the weighting def
        # object, the forcing for the source component is input in the
        # construction of the forcing transformation. If the source
        # component forcing is also a transformation, then the
        # transformtion object is passed. If not, then the input
        # forcing data is passed.            
        for ii in transformdefIndexes_inputComponents:
            
            modelComponent = varargin[ii,1]
            propertyValue  = varargin[ii,3]

            # Get the following inputs to build the model object:
            # transformation def name, input data column names,
            # the variables to which the input data is to be assigned
            # too, and additional transformation model options.
            filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[1])
            transformobject_name = varargin[ii,3][filt,2]
            filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[2])
            transformobject_inputs = varargin[ii,3][filt,2]
            filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[4])
            if any(filt):
                transformobject_options = varargin[ii,3][filt,2]
            else:
                transformobject_options = []

            filt = strcmp(varargin[ii,3][:,1], valid_transformProperties[5])
            transformobject_sourceName = varargin[ii,3][filt,2]
            
            ## Find the name of the input transformation object for the
            ## specified input component name.
            #transformobject_inputComponentName = self.inputData.componentData.(transformobject_inputComponentName).forcing_object
            
            # Find the required columns in forcing data so that only the
            # required data is input.
            colNum = 1
            for j in range(np.shape(forcingData_inputs)[1]):
                if np.isnumeric(forcingData_inputs[j,1]):
                    colNum = colNum, forcingData_inputs[j,1]             
                
                elif np.ischar(forcingData_inputs[j,1]):
                    filt = find(strcmpi(forcingData_colnames, forcingData_inputs[j,1]))
                    if np.sum[filt]==0:
                        print 'An unexpected error occured for component :', modelComponent, '. Within the input cell array for "forcingdata", the second column contains a forcing site name that is not listed within the input forcing data column names.'
                    colNum = [colNum filt]                                
                
                else:
                    print 'An unexpected error occured for component :', modelComponent,'. Within the input cell array for "forcingdata", the second column contains a forcing data column number or site name that is listed within the input forcing data column names.' 

            try:
                #colNum_extras = True(len(forcingData_colnames),1)
                #colNum_extras(colNum) = False
                #colNum_extras = find(colNum_extras)
                #forcingData_colnamesTmp = forcingData_colnames([colNumcolNum_extras])
                #forcingData_dataTmp = forcingData_data(:,[colNumcolNum_extras])                                                                                    
                self.parameters.(transformobject_name) = feval(transformobject_name, bore_ID, forcingData_data(:,colNum), forcingData_colnames(colNum), siteCoordinates, transformobject_inputs, self.parameters.(transformobject_sourceName), transformobject_options)
                
                # Add coordinates for the output variable
                variableName = self.inputData.componentData.(modelComponent).outputVariable
                coordinatesTmp = getCoordinates(self.parameters.(transformobject_name), variableName)
                siteCoordinates(np.shape(siteCoordinates)[1])+1: np.shape(siteCoordinates)[1] + np.shape(coordinatesTmp)[1],1:3) = coordinatesTmp
                                    
            except:
                print 'ERROR: Invalid model component class object for forcing transform: ', transformobject_name
                #rethrow(exception)

        # Build the component weighting object for each component that DO NOT require an input weighting def object
        for ii in weightingdefIndexes:
        
            modelComponent = varargin[ii,1]
            propertyValue  = varargin[ii,3]

            # Check if the current model object requires input of
            # another weighting def object.
            inputWeightingdefIndex = find(strcmpi(varargin[:,1],modelComponent) & strcmpi(varargin[:,2],valid_properties[4]))
            
            # Build weighting def for those NOT requiring the
            # input of other weighting def objects.
            if np.isempty(inputWeightingdefIndex):
                try:
                    # Get the column number for forcing data
                    if np.isfield(self.inputData.componentData.(modelComponent), 'dataColumn'):
                        colNum = self.inputData.componentData.(modelComponent).dataColumn
                        forcingData_inputs = forcingData_colnames(colNum)
                    elif np.isfield(self.inputData.componentData.(modelComponent),'inputForcing'):

                        # Get the names (or column numbers) of the input forcing data. 
                        forcingData_inputs = self.inputData.componentData.(modelComponent).outputVariable
                        
                        # Remove empty rows
                        if np.iscell(forcingData_inputs):
                            filt = cellfun(@(x) ~np.isempty(x), forcingData_inputs)
                            forcingData_inputs = forcingData_inputs[filt]
                        #for j=1:np.shape(forcingData_inputs,1)
                        #    # Skip if optional forcing input
                        #    if strcmp(forcingData_inputs[j,1],'(none)')
                        #        continue
                        #    end
                        #    if np.isnumeric(forcingData_inputs[j,1])
                        #        colNum = [colNum forcingData_inputs[j,1]]             
                        #    elif np.ischar(forcingData_inputs[j,1])
                        #        filt = find(strcmpi(forcingData_colnames, forcingData_inputs[j,1]))
                        #        if sum[filt]==0
                        #            print 'An unexpected error occured for component :', modelComponent,'. Within the input cell array for "forcingdata", the second column contains a forcing site name that is not listed within the input forcing data column names.' ])
                        #        end
                        #        colNum = [colNum filt]                                
                        #    else:
                        #        print 'An unexpected error occured for component :', modelComponent,'. Within the input cell array for "forcingdata", the second column contains a forcing data column number or site name that is listed within the input forcing data column names.' ])

                    # Check that the response 
                    # def name is consistent with the abstract
                    # 'responsedef_abstract'.
                    try:
                        if ~strcmp(findAbstractName(propertyValue),'responsedef_abstract'):
                            print 'The following response def def class definition is not derived from the "responsedef_abstract.m" anstract:', propertyValue
                    except:
                        print '... Warning: Checking that the required abstract for the response def transform class definition was used failed. This is may be because the version of matlab is pre-2014a.'

                    # Filter for options.
                    filt = find(strcmpi(varargin[:,1],modelComponent) & strcmpi(varargin[:,2],valid_properties[3]))                       
                    
                    # Convert input forcing data to cell if string
                    if np.ischar(forcingData_inputs):
                        forcingData_inputs_tmp[1] = forcingData_inputs
                        forcingData_inputs = forcingData_inputs_tmp
                        clear forcingData_inputs_tmp
                    
                    # Call the object 
                    if any(filt):
                        self.parameters.(modelComponent) = feval(propertyValue, bore_ID, forcingData_inputs, siteCoordinates, varargin[filt,3])                             
                    else:
                        self.parameters.(modelComponent) = feval(propertyValue, bore_ID, forcingData_inputs, siteCoordinates, []) 
                except:
                     print 'ERROR: Invalid weighting def class object name: ', char(propertyValue), '. The weighting def object could not be created.'
                     #rethrow(exception)
        
        # Build the component weighting object for each component that DO require an input weighting def object
        for ii in weightingdefIndexes:
        
            modelComponent = varargin[ii,1]
            propertyValue = varargin[ii,3]

            # Check if the current model object requires input of
            # another weighting def object.
            inputWeightingdefIndex = find(strcmpi(varargin[:,1],modelComponent) & strcmpi(varargin[:,2],valid_properties[4]))
            
            # Build weighting def for those that DO require the
            # input of other weighting def objects.
            if ~np.isempty(inputWeightingdefIndex):
                
                # Get the name of the component to be input to the
                # weighting def.
                inputWeightingdefName = varargin[inputWeightingdefIndex,3]                    
                
                # Add input component name to input data fields
                self.inputData.componentData.(modelComponent).inputWeightingComponent = inputWeightingdefName
                
                # Check that the input weighting def has been
                # created.
                if ~np.isfield(self.parameters,inputWeightingdefName):
                    print 'The input component name for component "', modelComponent, '" must be a component that itself is not derived from another component.'
                
                try:
                    # Get the column number for forcing data
                    if np.isfield(self.inputData.componentData.(modelComponent),'dataColumn'):
                        colNum = self.inputData.componentData.(modelComponent).dataColumn
                        forcingData_inputs = forcingData_colnames(colNum)
                    elif np.isfield(self.inputData.componentData.(modelComponent),'inputForcing'):

                        # Get the names (or column numbers) of the input forcing data. 
                        forcingData_inputs = self.inputData.componentData.(modelComponent).outputVariable
                        
                        # Remove empty rows
                        if np.iscell(forcingData_inputs):
                            filt = cellfun(@(x) ~np.isempty(x), forcingData_inputs)
                            forcingData_inputs = forcingData_inputs[filt]
                    
                    ## Get the column number for forcing data
                    #colNum = 1
                    #if np.isfield(self.inputData.componentData.(modelComponent),'dataColumn')
                    #    colNum = self.inputData.componentData.(modelComponent).dataColumn
                    #elif np.isfield(self.inputData.componentData.(modelComponent),'inputForcing')
 		    #
                    #    # Get the names (or column numbers) of the input forcing data. 
                    #    forcingData_inputs = self.inputData.componentData.(modelComponent).inputForcing[:,2]
                    #    
                    #    for j=1:np.shape(forcingData_inputs,1)
                    #        if np.isnumeric(forcingData_inputs[j,1])
                    #            colNum = [colNum forcingData_inputs[j,1]]             
                    #        elif np.ischar(forcingData_inputs[j,1])
                    #            filt = find(strcmpi(forcingData_colnames, forcingData_inputs[j,1]))
                    #            if sum[filt]==0
                    #                print 'An unexpected error occured for component :', modelComponent,'. Within the input cell array for "forcingdata", the second column contains a forcing site name that is not listed within the input forcing data column names.' ])
                    #            colNum = [colNum filt]                                
                    #        else:
                    #            print 'An unexpected error occured for component :', modelComponent,'. Within the input cell array for "forcingdata", the second column contains a forcing data column number or site name that is listed within the input forcing data column names.' ])

                    # Convert input forcing data to cell if string
                    if np.ischar(forcingData_inputs):
                        forcingData_inputs_tmp[1] = forcingData_inputs
                        forcingData_inputs = forcingData_inputs_tmp
                        clear forcingData_inputs_tmp

                    # Filter for options.
                    filt = find(strcmpi(varargin[:,1],modelComponent) & strcmpi(varargin[:,2],valid_properties[3]))
                    
                    # Call the object. 
                    # NOTE: these objects have an additional input
                    # (compared to those created above) for the input
                    # of a previously build model weighting def
                    # object.
                    self.parameters.(modelComponent) = feval(propertyValue, bore_ID,forcingData_inputs, siteCoordinates, self.parameters.(inputWeightingdefName), varargin(filt,3)) 
                except:
                     print 'ERROR: Invalid weighting def class object name: ', char(propertyValue), '. The weighting def object could not be created.'
                     #rethrow(exception)
        
        # Add noise component
        self.parameters.noise.type  = 'transferdef'
        self.parameters.noise.alpha = np.log10(0.1)    
        
        # Set the parameter names variable.   
        [junk, self.variables.param_names] = getParameters(self)            

        # Set variable declarign that calibration is not being
        # undertaken.
        self.variables.doingCalibration = False
        
        # Store only input forcing data required by the model. The
        # forcing data columns are removed to reduce RAM requirements.
        
        # Find required columns.
        modelComponentAll = fieldnames(self.inputData.componentData)
        forcingData_requiredColNames = []
        for i in range(len(modelComponentAll)):
            if np.isfield(self.inputData.componentData.(modelComponentAll[i]), 'dataColumn')
                forcingData_requiredColNames = [forcingData_requiredColNames[:], forcingData_colnames[self.inputData.componentData.(modelComponentAll[i]).dataColumn]]
            else:
                forcingData_requiredColNames = [forcingData_requiredColNames[:], self.inputData.componentData.(modelComponentAll[i]).inputForcing[:,2]]
        
        # Create filter for required colnames
        filt = cellfun(@(x) any(strcmp(x, forcingData_requiredColNames)) , forcingData_colnames)
        filt[1] = True
        
        forcingData_data = forcingData_data[:,filt]
        forcingData_colnames = forcingData_colnames[filt]
        setForcingData(self, forcingData_data, forcingData_colnames)


        # Get the observed head
        def getObservedHead(self):
            head = self.inputData.head
            return head
        
        
        # Get the forcing data from the model
        def getForcingData(self):
            forcingData = self.inputData.forcingData
            forcingData_colnames = self.inputData.forcingData_colnames            
            return forcingData, forcingData_colnames
        
        
        def getDerivedForcingData(self, t):

            # Initialise outputs
            forcingData = []
            forcingData_colnames = []                        
            
            # Get the derived forcing data in the sub-model objects
            if ~np.isempty(self.parameters):
                modelnames = fieldnames[self.parameters]
                for i in range(len(modelnames)):
                    if np.isobject(self.parameters.(modelnames[i])):
                        forcingData_colnames_tmp = []
                        forcingData_tmp = []
                        if ((any(strcmp(methods(self.parameters.(modelnames[i])), 'setTransformedForcing'))) & 
                            (any(strcmp(methods(self.parameters.(modelnames[i])), 'getTransformedForcing')))):
                            # Get the list of all possible forcing data outputs.                                
                            # In doing so, handle situations where older
                            # HydroSight models did not have the bore ID or
                            # site coordinates stored within the object.
                            if ~np.isfield(self.inputData, 'bore_ID'):
                                bore_ID = ''
                            else:
                                bore_ID = self.inputData.bore_ID

                            if ~np.isfield(self.inputData, 'siteCoordinates')
                                siteCoordinates = table()
                            else:
                                siteCoordinates = self.inputData.siteCoordinates

                            variable_names = feval([modelnames[i], '.outputForcingdata_options'], bore_ID , self.inputData.forcingData, self.inputData.forcingData_colnames, siteCoordinates)

                            # Set the forcing
                            setTransformedForcing(self.parameters.(modelnames[i]), t, False)
                            
                            # Get the forcing for each variable
                            for j in range(len(variable_names)):
                                try:
                                    forcingData = [forcingData, getTransformedForcing(self.parameters.(modelnames[i]),variable_names[j])]
                                    forcingData_colnames = [forcingData_colnames[:], variable_names[j]]
                                except:
                                    continue

            return forcingData, forcingData_colnames
            
        
        ## Set the forcing data from the model
        def setForcingData(self, forcingData, forcingData_colnames):
            self.inputData.forcingData = forcingData
            self.inputData.forcingData_colnames = forcingData_colnames
            
            # Update forcing data in the sub-model objects
            if ~np.isempty(self.parameters):
                modelnames = fieldnames[self.parameters]
                for i in range(len(modelnames)):
                    if np.isobject(self.parameters.(modelnames[i])):
                        try:
                            setForcingData(self.parameters.(modelnames[i]), forcingData, forcingData_colnames)
                        except:
                            pass

        
        def getStochForcingData(self):
            stochForcingData = []
            # Get derived forcing data in the sub-model objects
            if ~np.isempty(self.parameters):
                modelnames = fieldnames[self.parameters]
                for i in range(len(modelnames)):
                    if np.isobject(self.parameters.(modelnames[i])):
                        if any(strcmp('getStochForcingData', methods(self.parameters.(modelnames[i])))):
                            forcingData_tmp = getStochForcingData(self.parameters.(modelnames[i]))
                            stochForcingData.(modelnames[i]) = forcingData_tmp
            return stochForcingData

           
        def updateStochForcingData(self, stochForcingData, refineStochForcingMethod):
            # Set derived forcing data in the sub-model objects
            #isValidDerivedForcing =false
            componentNames = fieldnames[self.parameters]
            finishedStochForcing = False
            for i inn range(len(componentNames)):
                if np.isobject(self.parameters.(componentNames[i])):
                    if any(strcmp('updateStochForcingData',methods(self.parameters.(componentNames[i])))):
                        if nargin==1:
                            updateStochForcingData(self.parameters.(componentNames[i]))                                
                        else:
                            forcingDataComponent = fieldnames(stochForcingData)
                            filt = strcmp(forcingDataComponent,componentNames[i])
                            forcingDataComponent = forcingDataComponent[filt]
                            
                            if nargin==2:
                                updateStochForcingData(self.parameters.(componentNames[i]),stochForcingData.(forcingDataComponent))
                            elif nargin==3
                                finishedStochForcing = updateStochForcingData(self.parameters.(componentNames[i]),stochForcingData.(forcingDataComponent), refineStochForcingMethod)
            return finishedStochForcing
            
        
        def updateStochForcingParameters(self, stochForcingData, params):
            # Set the input parameters
            setParameters(self, params, self.variables.param_names)       
            
            # For those parameter compant objects with a method
            # 'updateStochForcingParameters', update the component
            # paramerters.
            componentNames = fieldnames[self.parameters]
            didUpdate = False
            for i in range(len(componentNames)):
                if np.isobject(self.parameters.(componentNames[i]))
                    if any(strcmp('updateStochForcingParameters',methods(self.parameters.(componentNames[i])))):
                        didUpdate = True
                        forcingDataComponent = fieldnames(stochForcingData)
                        filt = strcmp(forcingDataComponent,componentNames[i])
                        forcingDataComponent = forcingDataComponent[filt] 
                        updateStochForcingParameters(self.parameters.(componentNames[i]),stochForcingData.(forcingDataComponent))
            
            # If updated, then Get the new parameters.
            if didUpdate:
                params = getParameters(self)
            return params

        
        def setStochForcingState(self, doingCalibration, t_start_calib, t_end_calib):
            # Finalise the calibration of the stochastic forcing data.
            if ~np.isempty(self.parameters):
                modelnames = fieldnames[self.parameters]
                for i in range(len(modelnames)):
                    if np.isobject(self.parameters.(modelnames[i])):
                        if any(strcmp('setStochForcingState',methods(self.parameters.(modelnames[i])))):
                            setStochForcingState(self.parameters.(modelnames[i]), doingCalibration, t_start_calib, t_end_calib)

        
        def solve(self, time_points):
         
            ## Solve the model for the input time points
            # solve solves the model for the input time points.
            #
            # Syntax:
            #   [head, colnames, noise] = solve(self, time_points)
            #
            # Description:
            #   Solves the model using the model parameters and, depending upon
            #   the method's inputs, limits the groundwater head to the
            #   contribution from various periods of climate forcing and plots the
            #   results. The latter is achieved by inputting min and max values for
            #   tor but to date is not incorporated into the HydroSight()
            #   callign methods.
            #
            # Input:
            #   self -  model object
            #
            #   time_points - column vector of the time points to be simulated.
            #
            # Outputs:
            #   head - MxN matrix of simulated head with the following columns: date/time,
            #   head, head due to model component i.
            #
            #   colnames - Nx1 column names for matrix 'head'.
            #
            #   noise - Mx3 matrix of the estimated upper and lower magnitude of the
            #   time-series noise componat at M time steps.
            #
            # Example:
            #   see HydroSight: time_series_model_calibration_and_construction
            #
            # See also:
            #   HydroSight: time_series_model_calibration_and_construction
            #   model_TFN: model_construction
            #   calibration_finalise: initialisation_of_model_prior_to_calibration
            #   calibration_initialise: initialisation_of_model_prior_to_calibration
            #   get_h_star: main_method_for_calculating_the_head_contributions.
            #   getParameters: returns_a_vector_of_parameter_values_and_names
            #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration
            #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
            #
            # Author: 
            #   Dr. Tim Peterson, The Department of Infrastructure
            #   Engineering, The University of Melbourne.
            #
            # Date:
            #   26 Sept 2014

            # Check that the model has first been calibrated.
            if ~np.isfield(self.variables, 'd')
                print 'The model does not appear to have first been calibrated. Please calibrate the model before running a simulation.')
                        
            # Clear self.variables of temporary variables.
            if np.isfield(self.variables, 'theta_est_indexes_min'), self.variables = rmfield(self.variables, 'theta_est_indexes_min')
            if np.isfield(self.variables, 'theta_est_indexes_max'), self.variables = rmfield(self.variables, 'theta_est_indexes_max')
            #if np.isfield(self.variables, 'SMS_frac'), self.variables = rmfield(self.variables, 'SMS_frac')
            #if np.isfield(self.variables, 'recharge'), self.variables = rmfield(self.variables, 'recharge')
            #if np.isfield(self.variables, 'SMSC'), self.variables = rmfield(self.variables, 'SMSC')
            #if np.isfield(self.variables, 'SMSC_1'), self.variables = rmfield(self.variables, 'SMSC_1')
            #if np.isfield(self.variables, 'SMSC_2'), self.variables = rmfield(self.variables, 'SMSC_2')
            #if np.isfield(self.variables, 'f'), self.variables = rmfield(self.variables, 'f')
            #if np.isfield(self.variables, 'b'), self.variables = rmfield(self.variables, 'b')                 
            #if np.isfield(self.variables, 'ksat'), self.variables = rmfield(self.variables, 'ksat')
            #if np.isfield(self.variables, 'ksat_1'), self.variables = rmfield(self.variables, 'ksat_1')
            #if np.isfield(self.variables, 'ksat_2'), self.variables = rmfield(self.variables, 'ksat_2')
            
            # Set a flag to indicate that calibration is NOT being undertaken.
            # self.variables.doingCalibration = getInnovations
            
            # Setup matrix of indexes for tor at each time_points
            filt = self.inputData.forcingData[:,1]<=np.ceil(time_points[end])
            tor  = np.flipud([0:time_points[end]-self.inputData.forcingData[filt,1]+1])
            ntor = np.shape(tor)[1]                                                     
            clear tor
            
            self.variables.theta_est_indexes_min = np.zeros([1, len(time_points)])
            self.variables.theta_est_indexes_max = np.zeros([1, len(time_points)])
                        
            for ii in range(len(time_points)):
                ntheta = np.sum(self.inputData.forcingData[:,1]<=time_points[ii])
                self.variables.theta_est_indexes_min[ii] = ntor-ntheta
                self.variables.theta_est_indexes_max[ii] = np.max([1, ntor])
            
            # Free memory within mex def (just in case there'se been a
            # prior calibration that crashed prior to clearing the MEX
            # sttaic variables)
            # Free memory within mex def
            try:
                junk = doIRFconvolutionPhi([], [], [], [], False, 0)            
            except:
                pass               
            
            # Get the parameter sets (for use in resetting if >sets)
            params, param_names = getParameters(self)
            
            # If the number of parameter sets is >1 then temporarily apply 
            # only the first parameter set. This is done only to reduce the
            # RAM requirements for the broadcasting of the self variable in
            # the following parfor.
            if np.shape(params)[2]>1:
                setParameters(self, params[:,1], param_names)
            
            # Set percentile for noise 
            Pnoise = 0.95
            
            # Solve the modle using each parameter set.
            self.variables.delta_time = np.diff(time_points)
            headtmp = cell(1,np.shape(params,2))
            noisetmp = cell(1,np.shape(params,2))
            components = fieldnames[self.inputData.componentData]
            ncomponents = np.shape(components)[1]                             
            for ii in range(np.shape(params)[2]):
                # Get the calibration estimate of the mean forcing for the
                # current parameter set. This is a bit of a work around to
                # handle the issue of each parameter set having a np.unique
                # mean forcing (if a forcing transform is undertaken). The
                # workaround was required when DREAM was addded.
                for j in range(ncomponents):
                    calibData[ii,1].mean_forcing.(components[j]) = self.variables.(components[j]).forcingMean[:,ii]
                              
                # Add drainage elevation to the varargin variable sent to
                # objectivedef.                
                calibData[ii,1].drainage_elevation = self.variables.d[ii]
                
                # Add noise std dev
                if np.isfield(self.variables,'sigma_n'):
                    calibData[ii,1].sigma_n = self.variables.sigma_n(ii)
                else:
                    calibData[ii,1].sigma_n = 0
            
            # Solve model and add drainage constants
            dummy, headtmp[1], colnames = objectivedef(params[:,1], time_points, self, calibData[1])                        
            if np.shape(params)[2]>1:
                for jj in range(2:np.shape(params)[2]):
                    dummy, headtmp[jj] = objectivedef(params[:,jj], time_points, self, calibData[jj])
            
            # Add drainage constants and calculate total error bounds.
            for ii=1:np.shape(params,2)
                headtmp[ii][:,2] = headtmp[ii][:,2] + calibData(ii).drainage_elevation
            
                noisetmp[ii] = [headtmp[ii][:,1], np.ones(np.shape(headtmp[ii],1),2) .* norminv(Pnoise,0,1) .* calibData(ii).sigma_n]

            head = np.zeros(np.shape(headtmp[1],1),np.shape(headtmp[1],2), np.shape(params,2))
            noise = np.zeros(np.shape(headtmp[1],1),3, np.shape(params,2))            
            for ii in range(np.shape(params)[2]):
                head[:,:,ii]  = headtmp[ii]
                noise[:,:,ii] = noisetmp[ii]

            #colnames = colnames[1]
            clear headtmp, noisetmp
                             
            # Set the parameters if >1 parameter sets
            if np.shape(params)[2]>1:
                setParameters(self,params, param_names)
            
            # Clear matrix of indexes for tor at each time_points
            self.variables = rmfield(self.variables, 'theta_est_indexes_min')
            self.variables = rmfield(self.variables, 'theta_est_indexes_max')
            
            return head, colnames, noise
        
        
        def calibration_initialise(self, t_start, t_end):

            ## Initialise the model prior to calibration.
            # calibration_initialise initialises the model prior to calibration.
            #
            # Syntax:
            #   [params_initial, time_points] = calibration_initialise(self, t_start, t_end)
            #
            # Description:
            #   Sets up model variables required for calibration. Most imporantly, it
            #   calculates the mean observed head, h_bar, and row indexes, 
            #   theta_est_indexes,  for efficient calculation of the response
            #   defs. The calibration also requires the method to return a column 
            #   vector of the initial parameters and a column vector of time points for
            #   which observation data exists.
            #
            # Input:
            #   self -  model object
            #
            #   t_start - scaler start time, eg datenum(1995,1,1)
            #
            #   t_end - scaler end time, eg datenum(2005,1,1)
            #
            # Outputs:
            #   params_initial - column vector of the initial parameters.
            #
            #   time_points - column vector of time points forwhich observation data
            #   exists.
            #
            # Example:
            #   see HydroSight: time_series_model_calibration_and_construction
            #
            # See also:
            #   HydroSight: time_series_model_calibration_and_construction
            #   model_TFN: model_construction
            #   calibration_finalise: initialisation_of_model_prior_to_calibration
            #   get_h_star: main_method_for_calculating_the_head_contributions.
            #   getParameters: returns_a_vector_of_parameter_values_and_names
            #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration
            #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
            #   solve: solve_the_model_at_user_input_sime_points
            #
            # Author: 
            #   Dr. Tim Peterson, The Department of Infrastructure
            #   Engineering, The University of Melbourne.
            #
            # Date:
            #   26 Sept 2014

            # clear old variables
            self.variables = []

            # Set a flag to indicate that calibration is being undertaken
            # and the calibraion start and end dates            
            self.variables.doingCalibration = True
            self.variables.t_start = t_start
            self.variables.t_end = t_end
                        
            # Free memory within mex def (just in case there'se been a
            # prior calibration that crashed prior to clearing the MEX
            # sttaic variables)
            # Free memory within mex def
            try:
                junk = doIRFconvolutionPhi([], [], [], [], False, 0)            
            except:
                pass               
                        
            # Get parameter names and initial values
            params_initial, self.variables.param_names = getParameters(self)
            
            # Extract time steps
            t_filt = find(self.inputData.head[:,1]>=t_start) & (self.inputData.head[:,1]<=t_end)   
            time_points = self.inputData.head[t_filt, 1]
            self.variables.time_points = time_points
            
            # Set the stochastic forcing to be in a 'calibration' state.
            setStochForcingState(self, self.variables.doingCalibration, t_start, t_end)            
            
            # Set initial time for trend term
            self.variables.trend_startTime = self.variables.time_points[1]
            
            # Calc time step np.shapes. 
            self.variables.delta_time = np.diff(self.variables.time_points)
            
            # Calcaulte mean head.
            self.variables.h_bar = np.mean(self.inputData.head[t_filt, 2])                            
            
            # Setup matrix of indexes for tor at each time_points
            filt = self.inputData.forcingData[: ,1]<=np.ceil(time_points[end])
            tor  = np.flipud([0:time_points[end]-self.inputData.forcingData[filt,1]+1])
            ntor = np.shape(tor)[1]                                                     
            clear tor
            
            self.variables.theta_est_indexes_min = np.zeros([1, len(time_points)])
            self.variables.theta_est_indexes_max = np.zeros([1, len(time_points)])
                        
            for ii in range(len(time_points)):
                ntheta = np.sum(self.inputData.forcingData[:,1]<=time_points[ii])                                
                self.variables.theta_est_indexes_min(ii) = ntor-ntheta
                self.variables.theta_est_indexes_max(ii) = np.max([1, ntor])

            self.variables.nobjectivedef_calls = 0
            
            return params_initial, time_points
        
        
        def calibration_finalise(self, params, useLikelihood):           

            ## Finalise the model following calibration.
            # calibration_finalise finalises the model following calibration.
            #
            # Syntax:
            #   calibration_finalise(self, params)   
            #
            # Description:
            #   Finalises the model following calibration and assigns the final 
            #   parameters and additional variables to the object for later simulation. 
            #   Of the variables calculated, the most essential for the model 
            #   is the scalar drainage, self.variables.d. Other variables that are also 
            #   important include: 
            #       - a vector of innovations,  self.variables.innov, for detection 
            #         of serial correlation in the model errors 
            #       - the noise standard deviation, self.variables.sigma_n.
            #
            # Input:
            #   self -  model object
            #
            #   params - column vector of the optima parameters derived from
            #   calibration.
            #
            # Outputs:
            #   (none, the results are output to self.variables)
            #
            # Example:
            #   see HydroSight: time_series_model_calibration_and_construction
            #
            # See also:
            #   HydroSight: time_series_model_calibration_and_construction
            #   model_TFN: model_construction
            #   calibration_initialise: initialisation_of_model_prior_to_calibration
            #   get_h_star: main_method_for_calculating_the_head_contributions.
            #   getParameters: returns_a_vector_of_parameter_values_and_names
            #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration
            #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
            #   solve: solve_the_model_at_user_input_sime_points
            #
            # Author: 
            #   Dr. Tim Peterson, The Department of Infrastructure
            #   Engineering, The University of Melbourne.
            #
            # Date:
            #   26 Sept 2014            

            # Set the stochastic forcing to NOT be in a 'calibration' state.
            setStochForcingState(self, False, self.variables.t_start, self.variables.t_end)

            for j in range(np.shape(params)[2]:
                # Re-calc objective def and deterministic component of the head and innocations.
                # Importantly, the drainage elevation (ie the constant term for
                # the regression) is calculated within 'objectivedef' and
                # assigned to the object. When the calibrated model is solved
                # for a different time period (or climate data) then this
                # drainage value will be used by 'objectivedef'.
                try:
                    self.variables.selfFn[:,j], h_star, dummy , self.variables.d[j] = objectivedef(params[:,j], self.variables.time_points, self,useLikelihood)                        
                except:
                    if np.shape(params)[2]>1:
                        continue
                    else:
                        #print ME.message)
                        pass #?

                # Calculate the mean forcing rates. These mean rates are used
                # to calculate the contribution from the tail of the theta
                # weighting beyond the last observation. That is, the theta
                # def is integrated from the last time point of the
                # forcing to negative infinity. This integral is then
                # multiplied by the mean forcing rate. To ensure future
                # simulations use an identical mean forcing rate, the means are
                # calculated here and used in all future simulations.
                components = fieldnames[self.inputData.componentData]
                ncomponents = np.shape(components)[1]
                for i in range(ncomponents):
                    self.variables.(components[i]).forcingMean[:,j] = np.mean(self.variables.(components[i]).forcingData)

                t_filt = find(self.inputData.head[:,1]>=self.variables.time_points(1) & self.inputData.head[:,1]<=self.variables.time_points[end])               
                self.variables.resid[:,j] = self.inputData.head[t_filt,2]-h_star[:,2]+self.variables.d[j]

                # Calculate mean of noise. This should be zero +/- eps
                # because the drainage value is approximated assuming n-bar = 0.
                self.variables.n_bar[j] = np.real(np.mean(self.variables.resid[:,j]))

                # Calculate innovations
                innov = self.variables.resid[2:end, j]-self.variables.resid[1:end-1, j] .* np.exp(-10.**self.parameters.noise.alpha .* self.variables.delta_time)

                # Calculate noise standard deviation.
                self.variables.sigma_n[j] = np.sqrt(np.mean(innov .** 2 ./ (1.-np.exp(-2 .* 10.**self.parameters.noise.alpha .* self.variables.delta_time))))
            
            # Set a flag to indicate that calibration is complete.
            self.variables.doingCalibration = False
            
            # Set model parameters (if params are multiple sets)
            if np.shape(params)[2]>1:
                setParameters(self, params, self.variables.param_names)            
            
            # Free memory within mex def
            try:
                junk = doIRFconvolutionPhi([], [], [], [], False, 0)            
            except:
                pass


        def objectivedef(params, time_points, self, varargin):

            ## Calculate objective def vector. 
            # objectivedef calculates the objective def vector. 
            #
            # Syntax:
            #   [selfFn, h_star, colnames, drainage_elevation] = objectivedef(params,time_points, self)
            #
            # Description:
            #   Solves the model for the input parameters and calculates the objective
            #   def vector. Importantly, the objective def vector is not
            #   simply the difference between the observed and modelled heads. Because
            #   the model uses a noise model, the residual between the observed 
            #   and modelled head is first derived and then the innovation
            #   is calculated as the prior residual minus the later residual multiplied
            #   by the exponental noise def. Finally, the objective def
            #   weights this vector according to the time-step between observations.
            #
            #   Imporantly, the numerator of the weighting equation from von Asmuth et al 2002
            #   rounds to zero when the number of samples is very large. This occurs
            #   because it is effecively a geometric mean and its product term for n
            #   (where n is the number of observation for calibration minus 1) rounds
            #   to zero as a result of machine precision. This was overcome by adoption
            #   of a restructuring of the geometric meazn in term of exp and log terms. 
            #
            # Inputs:
            #   params - column vector of the optima parameters derived from
            #   calibration.
            #
            #   time_points - column vector of the time points to be simulated.  
            #
            #   self -  model object
            #
            # Outputs:
            #   selfFn - scalar objective def value.
            #
            #   h_star - matrix of the contribution from various model components and
            #   their summed influence. The matrix columns are in the order of:
            #   date/time, summed contribution to the head, contribution from
            #   component i.
            #
            #   colnames - column names for matrix 'head'.
            #
            #   drainage_elevation - drainage elevation constant.
            #
            # Example:
            #   see HydroSight: time_series_model_calibration_and_construction
            #
            # See also:
            #   HydroSight: time_series_model_calibration_and_construction
            #   model_TFN: model_construction
            #   calibration_finalise: initialisation_of_model_prior_to_calibration
            #   calibration_initialise: initialisation_of_model_prior_to_calibration
            #   get_h_star: main_method_for_calculating_the_head_contributions.
            #   getParameters: returns_a_vector_of_parameter_values_and_names
            #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
            #   solve: solve_the_model_at_user_input_sime_points
            #
            # References:
            #   von Asmuth J. R., Bierkens M. F. P., Mass K., 2002, Transfer
            #   dunction-noise modeling in continuous time using predefined impulse
            #   response defs.
            #
            # Author: 
            #   Dr. Tim Peterson, The Department of Infrastructure
            #   Engineering, The University of Melbourne.
            #
            # Date:
            #   26 Sept 2014    
            
            # Set model parameters
            setParameters(self, params, self.variables.param_names)
            
            # If varargin is a structural variable then the model
            # simulation is to use predefined values of the drainage
            # elevation and the mean forcing. Note, these inputs are only
            # to be provided if not doing simulation.
            getLikelihood = False
            drainage_elevation = []
            mean_forcing = []
            if ~np.isempty(varargin):
                if np.isstruct(varargin[1]):
                    drainage_elevation = varargin[1].drainage_elevation
                    mean_forcing = varargin[1].mean_forcing
                elif np.islogical(varargin[1]):
                    getLikelihood = varargin[1]
                else:
                    print 'The input varargin must be either a logical or a structural variable.'
            
            # Calc deterministic component of the head.
            if np.isempty(mean_forcing):
                h_star, colnames = get_h_star(self, time_points)                 
            else:
                h_star, colnames = get_h_star(self, time_points,mean_forcing)                 
            
            # Return of there are np.nan or inf value
            if (any(np.isnan(h_star[:,2])) | (np.isinf(h_star[:,2]))):
                selfFn = np.inf
                return
                            
            # If the results from this method call are not to be used for
            # summarising calibration results, then exit here. This is
            # required because the innovations can only be calculated at
            # time points for which there are observations. 
            if ~self.variables.doingCalibration:
                selfFn = []
                return
            
            # Calculate residual between observed and modelled.
            # Importantly, an approximation of the drainage level,d, is
            # required here to calculate the residuals. To achieve this
            # requires an assumption that the mean noise, n_bar, equals
            # zero. If this is not the case, then d_bar calculated below
            # may differ from d calculated within 'calibration_finalise'.
            t_filt = find(self.inputData.head[:,1]>=time_points[1] & self.inputData.head[:,1]<=time_points[end])          
            if np.isempty(drainage_elevation):
                drainage_elevation = self.variables.h_bar-np.mean(h_star[:,2])      

            resid = self.inputData.head[t_filt,2]-(h_star[:,2]+drainage_elevation)                 
            
            # Calculate innovations using residuals from the deterministic components.            
            innov = resid[2:end]-resid[1:end-1] .* np.exp(-10.**self.parameters.noise.alpha .* self.variables.delta_time)
            
            # Calculate objective def
            selfFn = (np.sum(np.exp(np.mean(np.log(1.-np.exp(-2.*10.**self.parameters.noise.alpha .* self.variables.delta_time)))) ./
                      (1.-np.exp(-2.*10.**self.parameters.noise.alpha .* self.variables.delta_time)) .* innov .** 2.))
  
            # Calculate log likelihood    
            if getLikelihood:
                N = np.shape(resid)[1]
                selfFn = -0.5*N*(np.log(2.*np.pi)+np.log(selfFn ./ N)+1.) 

            # Increment count of def calls
            self.variables.nobjectivedef_calls = self.variables.nobjectivedef_calls+np.shape(params)[2]
            
            return selfFn, h_star, colnames, drainage_elevation


        def setParameters(self, params, param_names):
           
            ## Set the model parameters to the model object from a vector.
            # setParameters set the model parameters to the model object from a vector.
            #
            # Syntax:
            #   setParameters(self, params, param_names)
            #
            # Description:
            #   Assigns a vector of parameter values to the model object using the
            #   input component and parameter names. This method is predominately used
            #   by the model calibration.
            #
            # Input:
            #   self -  model object
            #
            #   params - column vector of the parameter values.
            #
            #   param_names - two column n-row (for n parameters) cell matrix of
            #   compnp.nant name (column 1) and parameter name (column 2).
            #
            # Outputs:
            #   params_initial - column vector of the initial parameters.
            #
            #   time_points - column vector of time points forwhich observation data
            #   exists.
            #
            # Example:
            #   see HydroSight: time_series_model_calibration_and_construction
            #
            # See also:
            #   HydroSight: time_series_model_calibration_and_construction
            #   model_TFN: model_construction
            #   calibration_finalise: initialisation_of_model_prior_to_calibration
            #   calibration_initialise: initialisation_of_model_prior_to_calibration
            #   get_h_star: main_method_for_calculating_the_head_contributions.
            #   getParameters: returns_a_vector_of_parameter_values_and_names
            #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration  
            #   solve: solve_the_model_at_user_input_sime_points
            #
            # Author: 
            #   Dr. Tim Peterson, The Department of Infrastructure
            #   Engineering, The University of Melbourne.
            #
            # Date:
            #   26 Sept 2014            
            
            # Get np.unique component names and count the number of parameters
            # for this component.
            
            ncomponents = 1
            components[ncomponents ,1] = param_names[1,1]
            components[ncomponents ,2] = 1
            if np.shape(param_names)[1]>1:
                for ii in range(2:np.shape(param_names,1)
                    if ~strcmp(param_names[ii,1], param_names[ii-1,1])
                        ncomponents +=1
                        components[ncomponents ,1] = param_names[ii,1]
                        components[ncomponents ,2] = 1
                    elif ~strcmp(param_names[ii,2], 'type'):
                        components[ncomponents ,2] = components[ncomponents,2]+1
                        
            params_index = 0
            for ii in range(np.shape(components)[1]):
                currentField = char(components[ii,1]) 
                
                # Get parameter names for this component
                if np.isobject(self.parameters.(currentField)):
                    param_values, component_params = getParameters(self.parameters.(currentField))
                    clear param_values
                else:
                     component_params = fieldnames[self.parameters.(currentField)]

                # Scan though and remove those that are called 'type'
                param_index = range(np.shape(component_params)[1]) 
                for j in range(np.shape(component_params)[1]):
                    if (strcmp(component_params[j], 'type')) | (strcmp(component_params[j], 'variables')):
                        param_index[j] = 0

                component_params = component_params[param_index>0]
                
                # Check all parameters for this component are to be set
                # with new values.
                if components[ii, 2]!=np.shape(component_params)[1]:
                    print 'The number of parameters to be set for the following model', 'component must equal the total number of parameters for the component: ,',currentField
                
                # Get input parameter values for this component.
                component_param_vals = np.zeros([np.shape(component_params)[1],np.shape(params)[2]])
                for j in range(np.shape(component_params)[1]):
                    params_index +=1
                    component_param_vals[j,:] = params[params_index,:]
                
                # Get model parameters for each component.
                # If the component is an object, then call the objects
                # getParameters method.
                if np.isobject(self.parameters.(currentField)):                
                    # Call object method to set parameter values.
                    setParameters(self.parameters.(currentField), component_param_vals)                   
                else:
                    # Non-object component.
                    for j in range(np.shape(component_params)[1]):               
                        self.parameters.(currentField).(char(component_params[j])) = component_param_vals[j,:]

        
        def getParameters(self):

            ## Returns the model parameters from the model object.
            # getParameters returns the model parameters from the model object.
            #
            # Syntax:
            #   [params, param_names] = getParameters(self)
            #
            # Description:
            #   Cycles through all model components and parameters and returns a vector
            #   of parameter values and a cell matrix of their espective component
            #   and parameter names.
            #
            # Input:
            #   self -  model object.
            #
            # Outputs:
            #   params - column vector of the parameter values.
            #
            #   param_names - two column n-row (for n parameters) cell matrix of
            #   compnp.nant name (column 1) and parameter name (column 2).
            #
            # Example:
            #   see HydroSight: time_series_model_calibration_and_construction
            #
            # See also:
            #   HydroSight: time_series_model_calibration_and_construction
            #   model_TFN: model_construction
            #   calibration_finalise: initialisation_of_model_prior_to_calibration
            #   calibration_initialise: initialisation_of_model_prior_to_calibration
            #   get_h_star: main_method_for_calculating_the_head_contributions.
            #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration  
            #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
            #   solve: solve_the_model_at_user_input_sime_points
            #
            # Author: 
            #   Dr. Tim Peterson, The Department of Infrastructure
            #   Engineering, The University of Melbourne.
            #
            # Date:
            #   26 Sept 2014            
            
            param_names = []
            params = []
            nparams = 0
            # Get model components
            components = fieldnames[self.parameters]
            
            for ii in range(np.shape(components)[1]):
                currentField = char(components[ii]) 
                # Get model parameters for each component.
                # If the component is an object, then call the objects
                # getParameters method.
                if np.isobject(self.parameters.(currentField)):
                    # Call object method.
                    params_temp, param_names_temp = getParameters(self.parameters.(currentField))
                    
                    for j in range(np.shape(param_names_temp)[1]):
                        nparams +=1
                        param_names[nparams,1] = currentField
                        param_names[nparams,2] = param_names_temp[j]
                        params[nparams,:] = params_temp[j,:]                        
                    
                else:
                    # Non-object component.
                    component_params = fieldnames[self.parameters.(currentField)]
                    for j in range(np.shape(component_params)[1]):
                       if ~strcmp(component_params[j], 'type'):
                          nparams +=1
                          param_names[nparams,1] = currentField
                          param_names[nparams,2] = char(component_params[j])
                          params[nparams,:] = self.parameters.(currentField).(char(component_params[j]))
            
            return params, param_names
        

        def getDerivedParameters(self):
        
            ## Returns the model parameters from the model object.
            # getDerivedParameters returns the derived parameters from the model object.
            #
            # Syntax:
            #   [params, param_names] = getDerivedParameters(self)
            #
            # Description:
            #   Cycles through all model components and parameters and returns a vector
            #   of derived parameter values and a cell matrix of their respective component
            #   and parameter names. The derived parameters are calibrated parameters
            #   or constants in the components but variables derived from the
            #   calibrated parameters. For example, the drawdown response def
            #   (e.g. responsedef_FerrisKnowles) can calculate the T and S from
            #   the calibrated parameters.
            #
            # Input:
            #   self -  model object.
            #
            # Outputs:
            #   params - column vector of the parameter values.
            #
            #   param_names - two column n-row (for n parameters) cell matrix of
            #   compnp.nant name (column 1) and parameter name (column 2).
            #
            # Example:
            #   see HydroSight: time_series_model_calibration_and_construction
            #
            # See also:
            #   HydroSight: time_series_model_calibration_and_construction
            #   model_TFN: model_construction
            #   calibration_finalise: initialisation_of_model_prior_to_calibration
            #   calibration_initialise: initialisation_of_model_prior_to_calibration
            #   get_h_star: main_method_for_calculating_the_head_contributions.
            #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration  
            #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
            #   solve: solve_the_model_at_user_input_sime_points
            #
            # Author: 
            #   Dr. Tim Peterson, The Department of Infrastructure
            #   Engineering, The University of Melbourne.
            #
            # Date:
            #   26 Sept 2014            
            
            param_names = []
            params = []
            nparams = 0
            # Get model components
            components = fieldnames[self.parameters]
            
            for ii in range(np.shape(components)[1]):
                currentField = char(components[ii]) 
                # Get model parameters for each component.
                # If the component is an object, then call the objects
                # getParameters method.
                if np.isobject(self.parameters.(currentField)):
                    # Call object method.
                    params_temp, param_names_temp = getDerivedParameters(self.parameters.(currentField))
                    
                    for j in range(np.shape(param_names_temp)[1]):
                        nparams +=1
                        param_names[nparams,1] = currentField
                        param_names[nparams,2] = param_names_temp[j]
                        params[nparams,:] = params_temp[j,:]           
            
            return params, param_names

             
        def derivedData_types = getDerivedDataTypes(self):

            # Initialise outputs
            derivedData_types = cell(0,2)
            
            # Get model components
            components = fieldnames[self.parameters]
            
            derivedData_types_p1 = []
            derivedData_types_p2 = []
            for ii in range(np.shape(components)[1]):
                modelcomponent = char(components[ii]) 
                # Get derived data for each component.
                if ((np.isobject(self.parameters.(modelcomponent))) & 
                    (any(strcmp(methods(self.parameters.(modelcomponent)), 'getDerivedData'))) & 
                    (any(strcmp(methods(self.parameters.(modelcomponent)), 'getDerivedDataTypes')))):
            
                    dataTypes_tmp = getDerivedDataTypes(self.parameters.(modelcomponent))
                    if ~np.iscell(dataTypes_tmp):
                        dataTypes = cell(1,1)
                        dataTypes[1] = dataTypes_tmp
                    else:
                        dataTypes = dataTypes_tmp
                    k = np.shape(derivedData_types)[1]
                    for j in range(len(dataTypes)):
                        derivedData_types[k+j, 1] = modelcomponent
                        derivedData_types[k+j, 2] = dataTypes[j]

        
        def getDerivedData(self, modelcomponent, derivedData_variable, t, plot_axes):

            # Initialise outputs
            derivedData = []
            derivedData_names = []
            
            if ((np.isobject(self.parameters.(modelcomponent)) & 
                (any(strcmp(methods(self.parameters.(modelcomponent)),'getDerivedData')))):
                derivedData, derivedData_names = getDerivedData(self.parameters.(modelcomponent), derivedData_variable, t, plot_axes)  

            return derivedData, derivedData_names

        
        def acceptStochForcingSolution(self, selfFuncVal, selfFuncVal_prior, stochForcingData):

            # Initialse some outputs
            #stochForcingData_new = stochForcingData
            stochForcingData_new = []
            
            # Get model components
            components = fieldnames[self.parameters]
            accepted = True(np.shape(components)[1],1)
            for ii in range(np.shape(components)[1]):
                currentField = char(components[ii]) 
                # Get model parameters for each component.
                # If the component is an object, then call the objects
                # getParameters method.
                if ((np.isobject(self.parameters.(currentField))) &
                    (any(strcmp(methods(self.parameters.(currentField)), 'acceptStochForcingSolution')))):

                    # Check if the current component is listed within the
                    # input stochastic forcing data.
                    forcingDataComponent = fieldnames(stochForcingData)
                    filt = strcmp(forcingDataComponent,components[ii])
                    
                    # Test if the solution should be accepted. Pass the
                    # stochastic forcing data is required.
                    if np.isempty(filt):
                        accepted[ii] = acceptStochForcingSolution(self.parameters.(currentField), selfFuncVal, selfFuncVal_prior)
                    else:
                        forcingDataComponent = forcingDataComponent[filt]
                        stochForcingData_new.(forcingDataComponent), accepted(ii) = acceptStochForcingSolution(self.parameters.(currentField), selfFuncVal, selfFuncVal_prior, stochForcingData.(forcingDataComponent))

            if any(~accepted):
                accepted = False
            else:
                accepted = True
            
            return stochForcingData_new, accepted
            
        
        def getParameters_physicalLimit(self):

            # getParameters_physicalLimit returns the physical limits to each model parameter.
            #
            # Syntax:
            #   [params_upperLimit, params_lowerLimit] = getParameters_physicalLimit(self)
            #
            # Description:
            #   Cycles through all model components and parameters and returns a vector
            #   of the physical upper and lower parameter bounds as defined by the
            #   weighting defs.
            #
            # Input:
            #   self -  model object.
            #
            # Outputs:
            #   params_upperLimit - column vector of the upper parameter bounds.
            #
            #   params_lowerLimit - column vector of the lower parameter bounds
            #
            # Example:
            #   see HydroSight: time_series_model_calibration_and_construction
            #
            # See also:
            #   HydroSight: time_series_model_calibration_and_construction
            #   model_TFN: model_construction
            #   calibration_finalise: initialisation_of_model_prior_to_calibration
            #   calibration_initialise: initialisation_of_model_prior_to_calibration
            #   get_h_star: main_method_for_calculating_the_head_contributions.
            #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration  
            #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
            #   solve: solve_the_model_at_user_input_sime_points
            #
            # Author: 
            #   Dr. Tim Peterson, The Department of Infrastructure
            #   Engineering, The University of Melbourne.
            #
            # Date:
            #   26 Sept 2014   

            params_lowerLimit = []
            params_upperLimit = []
            components = fieldnames[self.parameters]            
            for iiin range(np.shape(components)[1]):
                currentField = char(components[ii]) 
                # Get model parameters for each component.
                # If the component is an object, then call the objects
                # getParameters_physicalLimit method for each parameter.
                if np.isobject(self.parameters.(currentField)):
                    # Call object method.
                    params_upperLimit_temp, params_lowerLimit_temp = getParameters_physicalLimit(self.parameters.(currentField))                
                    params_upperLimit = [params_upperLimit, params_upperLimit_temp]
                    params_lowerLimit = [params_lowerLimit, params_lowerLimit_temp]
                                        
                else:
                    params_plausibleUpperLimit, params_plausibleLowerLimit = getParameters_plausibleLimit(self)
                    
                    # This parameter is assumed to be the noise parameter 'alpha'.  
                    ind = len(params_upperLimit)+1
                    params_upperLimit = [params_upperLimit, params_plausibleUpperLimit(ind)]                                                        
                    params_lowerLimit = [params_lowerLimit, params_plausibleLowerLimit(ind)]        
                    
           return params_upperLimit, params_lowerLimit

        
    def getParameters_plausibleLimit(self):

        # getParameters_plausibleLimit returns the plausible limits to each model parameter.
        #
        # Syntax:
        #   [params_upperLimit, params_lowerLimit] = getParameters_plausibleLimit(self)
        #
        # Description:
        #   Cycles though all model components and parameters and returns a vector
        #   of the plausible upper and lower parameter range as defined by the
        #   weighting defs.
        #
        # Input:
        #   self -  model object.
        #
        # Outputs:
        #   params_upperLimit - column vector of the upper parameter plausible bounds.
        #
        #   params_lowerLimit - column vector of the lower parameter plausible bounds
        #
        # Example:
        #   see HydroSight: time_series_model_calibration_and_construction
        #
        # See also:
        #   HydroSight: time_series_model_calibration_and_construction
        #   model_TFN: model_construction
        #   calibration_finalise: initialisation_of_model_prior_to_calibration
        #   calibration_initialise: initialisation_of_model_prior_to_calibration
        #   get_h_star: main_method_for_calculating_the_head_contributions.
        #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration  
        #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
        #   solve: solve_the_model_at_user_input_sime_points
        #
        # Author: 
        #   Dr. Tim Peterson, The Department of Infrastructure
        #   Engineering, The University of Melbourne.
        #
        # Date:
        #   26 Sept 2014   
                        
        params_lowerLimit = []
        params_upperLimit = []
        components = fieldnames[self.parameters]            
        for ii in range(np.shape(components)[1]):
            currentField = char(components[ii]) 
            # Get model parameters for each component.
            # If the component is an object, then call the objects
            # getParameters_physicalLimit method for each parameter.
            if np.isobject(self.parameters.(currentField)):
                # Call object method.
                params_upperLimit_temp, params_lowerLimit_temp = getParameters_plausibleLimit(self.parameters.(currentField))
                params_upperLimit = [params_upperLimit, params_upperLimit_temp]
                params_lowerLimit = [params_lowerLimit, params_lowerLimit_temp]
            else:
                component_params = fieldnames[self.parameters.(currentField)]
                for j in range(np.shape(component_params)[1]):
                   if ~strcmp(component_params[j], 'type'):
                       if (strcmp(currentField, 'et')) & (strcmp(component_params[j], 'k')):
                            # This parameter is assumed to be the ET
                            # parameter scaling when the ET
                            # uses the precipitation transformation
                            # def.
                            params_upperLimit = [params_upperLimit, 1.]
                            params_lowerLimit = [params_lowerLimit, 0.]
                       elif ((strcmp(currentField, 'landchange')) & (strcmp(component_params[j], 'precip_scalar')) |
                             (strcmp(currentField, 'landchange')) & (strcmp(component_params[j], 'et_scalar'))):  
                           # This parameter is the scaling parameter
                           # for either the ET or precip transformation
                           # defs.
                           params_upperLimit = [params_upperLimit,  1.]
                           params_lowerLimit = [params_lowerLimit, -1.]
                       else:
                            # This parameter is assumed to be the noise parameter 'alpha'.  
                            alpha_upperLimit = 100 
                            while ((np.abs(np.sum(np.exp(-2.*alpha_upperLimit .* self.variables.delta_time)))<eps) |
                                   (np.exp(np.mean(np.log(1.-np.exp(-2.*alpha_upperLimit .* self.variables.delta_time))))<eps)):
                                alpha_upperLimit -=0.01
                                if alpha_upperLimit<=eps:
                                    break

                            if alpha_upperLimit<=eps:
                                alpha_upperLimit = np.inf
                            else:
                                # Transform alpha log10 space.
                                alpha_upperLimit = np.log10(alpha_upperLimit)
                            
                            params_upperLimit = [params_upperLimit, alpha_upperLimit]
                            params_lowerLimit = [params_lowerLimit, np.log10(np.sqrt(eps))+4.]
                            
        return params_upperLimit, params_lowerLimit

        
    def getParameterValidity(self, params, time_points):
    
        # isValidParameter returns a logical vector for the validity or each parameter.
        #
        # Syntax:
        #   isValidParameter = getParameterValidity(self, params, time_points)
        #
        # Description:
        #   Cycles though all model components and parameters and returns a logical 
        #   vector denoting if each parameter is valid as defined by each weighting
        #   def.
        #
        # Input:
        #   self -  model object.
        #
        #   params - vector of model parameters
        #
        #   time_points - vector of simulation time points
        #
        # Outputs:
        #   isValidParameter - column vector of the parameter validity.
        #
        # Example:
        #   see HydroSight: time_series_model_calibration_and_construction
        #
        # See also:
        #   HydroSight: time_series_model_calibration_and_construction
        #   model_TFN: model_construction
        #   calibration_finalise: initialisation_of_model_prior_to_calibration
        #   calibration_initialise: initialisation_of_model_prior_to_calibration
        #   get_h_star: main_method_for_calculating_the_head_contributions.
        #   objectivedef: returns_a_vector_of_innovation_errors_for_calibration  
        #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
        #   solve: solve_the_model_at_user_input_sime_points
        #
        # Author: 
        #   Dr. Tim Peterson, The Department of Infrastructure
        #   Engineering, The University of Melbourne.
        #
        # Date:
        #   26 Sept 2014         
        
        isValidParameter = True(np.shape(params))
        components = fieldnames[self.parameters]            
        
        junk, param_names = getParameters(self)
        
        for ii in range(np.shape(components)[1]):
            currentField = char(components[ii]) 
            
            # Find parameters for the component
            filt = strcmp(components[ii], param_names[:,1])

            # Get model parameters for each component.
            # If the component is an object, then call the objects
            # getParameterValidity method for each parameter.
            if np.isobject(self.parameters.(currentField)):
                
                # Call object method and pass it the parameter vector and
                # parameter names.
                # NOTE: the parameters are not set within the component
                # object. This was done to avoid a call to
                # setParameters().
                isValidParameter[filt,:] = getParameterValidity(self.parameters.(currentField), params[filt,:], param_names[filt,2])
             
            # Check the alphanoise parameter is large enough not to cause numerical
            # problems in the calcuation of the objective def.
            elif strcmp('noise', currentField):
                alpha = params[filt,:] 
                filt_noiseErr = ((np.exp(np.mean(np.log(1.-np.exp(bsxfun(@times,-2.*10.**alpha, self.variables.delta_time))) ,1))<=eps) |
                                 (np.abs(np.sum(np.exp(bsxfun(@times,-2.*10.*alpha, self.variables.delta_time)),1))<eps)):
                isValidParameter[filt,filt_noiseErr] = False                
                
            else:
                isValidParameter[filt,:] = True
            
            # Break if any parameter sets are invalid!
            if (np.shape(params)[2]==1) & (any(any(~isValidParameter))):
                return

        return isValidParameter

    
    def plot_transferdefs(self, t_max):
    
        # plot_transferdefs plot the weighting defs
        #
        # Syntax:
        #   plot_transferdefs(self, t_max)
        #
        # Description:
        #   Creates a plot of each weighting def. 
        #
        # Input:
        #   self -  model object
        #
        #   t_max - scaler number of the maximum duration to plot (in days)
        #
        # Output:  
        #   (none)
        #
        # See also:
        #   HydroSight: time_series_model_calibration_and_construction
        #   model_TFN: model_construction
        #
        # Dependencies
        #   model_TFN.m
        #
        # Author: 
        #   Dr. Tim Peterson, The Department of Infrastructure Engineering, 
        #   The University of Melbourne.
        #
        # Date:
        #   7 May 2012       
        
        import matplotlib.pyplot as plt
        
        # Get components of the model.
        components = fieldnames[self.parameters]

        # Derive filter for components with linear contribution to
        # head, eg remove index to soil moisture component.
        component_indexes = True[np.shape(components)]
        component_indexes(strcmp(components, 'soilmoisture')) = False
        component_indexes[end] = False
        component_indexes = find(component_indexes) 

        # Calculate components with transfer def.
        nTranfterdefs = 0
        tor = range(t_max)
        for ii in range(component_indexes):

            # Get model parameters for each component.
            # If the component is an object, then call the objects
            # getParameters method.
            if np.isobject(self.parameters.(char(components[ii]))):                               
                
                # Calcule theta for each time point of forcing data.
                try:
                    tmp = theta(self.parameters.(char(components[ii])), tor)
                    nTranfterdefs =+1
                    if np.shape(tmp)[2]>1:
                        theta_est[:, nTranfterdefs] = np.sum(tmp,2)
                    else:
                        theta_est[:, nTranfterdefs] = tmp
                except:
                    pass
        
        plt.figure()
        for ii in range(nTranfterdefs):
            plt.add_subplot(nTranfterdefs, 1, ii)
            plt.plot(tor, theta_est[:,ii])
            plt.xlabel('Duration into the past (days)')
            plt.ylabel('Transfer def weight')
            plt.title(['Transfer def for: ', char(components[ii])])
            #box on
    
    
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
        #   self -  model object
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
            if np.isobject(self.(propNames[i])):
                delete(self.(propNames[i]))
            else:
                self.(propNames[i]) = [] 

       
    def get_h_star(self, time_points, varargin):
        
        ## Main method calculating the contribution to the head from each model component.
        # get_h_star private method calculating the contribution to the head from each model component.
        #
        # Syntax:
        #   [h_star, colnames] = get_h_star(self, time_points)
        #
        # Description:
        #   This method performs the main calculates of the model. The method 
        #   is private and is called by either solve() or objectivedef() public
        #   methods. The method is highly flexible, allowing any number of model
        #   components and forcing transformations. The numerical integration of
        #   the convolution def is also using a highly efficient MEX .c
        #   'doIRFconvolution.c', which undertakes trapazoidal integration for
        #   integrated forcing (eg daily precipitation) or Simpson's 3/8 composite
        #   inegration for instantaneous fluxes.
        
        #   Depending upon the user set model components, this method 
        #   undertakes the following seqential steps:
        #
        #   1) Transform the forcing data. This is undertaken using for the nonlinear
        #   TFN models of Peterson and Westenr (2014). 
        #
        #   2) If model component 'i' is an impulse response def (IRF), the 
        #   component method 'theta' is called with all daily step time points less
        #   than or equal to the latest date head for simulation. Theta and the
        #   forcing data are then passed to doIRFconvolution() for the integration at
        #   each water level observation time point.
        #   
        #   3) If the component is not an IRF, then a user specified
        #   IRF matrix from step 2 (which is often precipitation or ET) is scaled 
        #   by an appropriate parameter and daily forcing data and then integrated
        #   using doIRFconvolution().
        #
        #   4) Finally, the contribution from all components are summed to produce
        #   h* (as per von Asmuth et al 2002).
        #
        # Inputs:
        #   self -  model object
        #
        #   time_points - column vector of the time points to be simulated.  
        #
        # Outputs:
        #
        #   h_star - matrix of the contribution from various model components and
        #   their summed influence. The matrix columns are in the order of:
        #   date/time, summed contribution to the head, contribution frolm
        #   component i.
        #
        #   colnames - column names for matrix 'head'.
        #
        # Example:
        #   see HydroSight: time_series_model_calibration_and_construction
        #
        # See also:
        #   HydroSight: time_series_model_calibration_and_construction
        #   model_TFN: model_construction
        #   calibration_finalise: initialisation_of_model_prior_to_calibration
        #   calibration_initialise: initialisation_of_model_prior_to_calibration
        #   getParameters: returns_a_vector_of_parameter_values_and_names
        #   setParameters: sets_model_parameters_from_input_vector_of_parameter_values
        #   solve: solve_the_model_at_user_input_sime_points
        #
        # References:
        #   von Asmuth J. R., Bierkens M. F. P., Mass K., 2002, Transfer
        #   dunction-noise modeling in continuous time using predefined impulse
        #   response defs.
        #
        # Author: 
        #   Dr. Tim Peterson, The Department of Infrastructure
        #   Engineering, The University of Melbourne.
        #
        # Date:
        #   26 Sept 2014    

        # Get components of the model.
        components = fieldnames[self.inputData.componentData]
        ncomponents = np.shape(components)[1]  
	
        # Calc theta_t for max time to initial time (NOTE: min tor = 0).                        
        filt = self.inputData.forcingData[:,1)<=np.ceil(time_points[end])
        tor = np.flipud([0:time_points[end]-self.inputData.forcingData[1,1]+1])
        tor_end = tor(self.variables.theta_est_indexes_min[1,:])
        t = self.inputData.forcingData[filt ,1]
	
        # Calculate the transformation models that DO NOT require the
        # input of a forcing model and the forcing model is denoted as the
        # complete model to be ran (ie not also used by another
        # component).
        for i in range(ncomponents):
            if ((np.isfield(self.inputData.componentData.(components[i]), 'forcing_object')) &
                (~np.isfield(self.inputData.componentData.(components[i]), 'inputForcingComponent')) &
                (np.isfield(self.inputData.componentData.(components[i]), 'isForcingModel2BeRun')) &
                (self.inputData.componentData.(components[i]).isForcingModel2BeRun)):
                setTransformedForcing(self.parameters.(self.inputData.componentData.(components[i]).forcing_object), t, True)                    
        
        # Calculate the transformation models that DO require the
        # input of a forcing model  and the forcing model is denoted as the
        # complete model to be ran (ie not also used by another
        # component).
        for i in range(ncomponents):
            if ((np.isfield(self.inputData.componentData.(components[i]), 'forcing_object')) &
                (np.isfield(self.inputData.componentData.(components[i]), 'inputForcingComponent')) &
                (np.isfield(self.inputData.componentData.(components[i]), 'isForcingModel2BeRun')) &
                (self.inputData.componentData.(components[i]).isForcingModel2BeRun)):
                setTransformedForcing(self.parameters.(self.inputData.componentData.(components[i]).forcing_object), t, True)
        
        # Assign the forcing data for each model component.
        # Also, if the model is being calibrated, then also
        # calculate the mean forcing. The later is essential to reduce
        # the calibration error when using a forcing series
        # significantly shorter than the point back in time (ie tor) at
        # which the transfer def is approx. zero.
        nOutputColumns = 0
        isForcingADailyIntegral = False[ncomponents,1]
        for i in range(ncomponents):
           
            if ((np.isfield(self.inputData.componentData.(components[i]), 'forcing_object')) &
                (np.isfield(self.inputData.componentData.(components[i]), 'outputVariable'))):                
	
                # Get te transformed forcing. 
                # NOTE: getTransformedForcing() must return a boolean
                # scaler deonting if the forcing is an instantaneous flux
                # or an integral over the day. For the former, Simpson's
                # 3/8 integration is used for the convolution while for
                # the latter trapzoidal integration of theta is undertaken and then 
                # the daily integration result multiplied by the integrated daily forcing.
                self.variables.(components[i]).forcingData, isForcingADailyIntegral[i] = getTransformedForcing(self.parameters.(self.inputData.componentData.(components[i]).forcing_object), self.inputData.componentData.(components[i]).outputVariable)           
                self.variables.(components[i]).forcingData_colnames = self.inputData.componentData.(components[i]).outputVariable
                                    
            else:
                
                # Non-transformed forcing is assumed to be a daily
                # integral. Hence, doIRFconvolution() undertakes daily
                # trapzoidal integration of theta and multiplies it by
                # the integrated daily forcing.
                isForcingADailyIntegral[i] = True
                self.variables.(components[i]).forcingData = self.inputData.forcingData(filt, self.inputData.componentData.(components[i]).dataColumn)
                self.variables.(components[i]).forcingData_colnames = self.inputData.forcingData_colnames(self.inputData.componentData.(components[i]).dataColumn)

            # Increase the number of output columns for method.
            nOutputColumns = nOutputColumns+np.shape(self.variables.(components[i]).forcingData,2)
        
        # Initialise ouput vector.
        h_star = np.zeros([np.shape(time_points)[1], nOutputColumns])       
        iOutputColumns = 0
        # Calculate each transfer def.
        for i in range(ncomponents):
                            
            # Calcule theta for each time point of forcing data.
            theta_est_temp = theta(self.parameters.(char(components[i])), tor)                
	
            # Get analytical estimates of lower and upper theta tails
            integralTheta_upperTail = intTheta_upperTail2Inf(self.parameters.(char(components[i])), tor_end)                           
            integralTheta_lowerTail = intTheta_lowerTail(self.parameters.(char(components[i])), 1)
	
            # Get the mean forcing.
            if (~np.isempty(varargin) & (np.isfield(varargin[1],components[i])):
                #forcingMean = self.variables.(components[i]).forcingMean                    
                forcingMean = varargin[1].(components[i])
            else:
                forcingMean = np.mean(self.variables.(components[i]).forcingData)
            
            # Integrate transfer def over tor.
            for j in range(np.shape(theta_est_temp)[2]):
                # Increment the output volumn index.
                iOutputColumns +=1
	
                # Try to call doIRFconvolution using Xeon Phi
                # Offload coprocessors. This will only work if the
                # computer has (1) the intel compiler >2013.1 and (2)
                # xeon phi cards. The code first tried to call the
                # mex def.  
                if ~np.isfield(self.variables, 'useXeonPhiCard'):
                    self.variables.useXeonPhiCard = True
                
                try:
                    if self.variables.useXeonPhiCard:
                        #print 'Offloading convolution algorithm to Xeon Phi coprocessor!')
                        h_star[:, iOutputColumns] = doIRFconvolutionPhi(theta_est_temp[:,j], self.variables.theta_est_indexes_min, self.variables.theta_est_indexes_max[1], self.variables.(components[i]).forcingData[:,j], isForcingADailyIntegral[i], integralTheta_lowerTail[j])+integralTheta_upperTail[j,:] .* forcingMean[j]
                    else:
                        h_star[:, iOutputColumns] = doIRFconvolution(theta_est_temp[:,j], self.variables.theta_est_indexes_min, self.variables.theta_est_indexes_max[1], self.variables.(components[i]).forcingData[:,j], isForcingADailyIntegral[i], integralTheta_lowerTail[j])+integralTheta_upperTail[j,:] .* forcingMean[j]
                        
                except:
                    #print 'Offloading convolution algorithm to Xeon Phi coprocessor failed - falling back to CPU!')
                    self.variables.useXeonPhiCard = False
                    h_star[:, iOutputColumns] = doIRFconvolution(theta_est_temp[:,j], self.variables.theta_est_indexes_min, self.variables.theta_est_indexes_max[1], self.variables.(components[i]).forcingData[:,j], isForcingADailyIntegral[i], integralTheta_lowerTail[j])+integralTheta_upperTail[j,:] .* forcingMean[j]

                # Transform the h_star estimate for the current
                # component. This feature was included so that h_star
                # estimate fro groundwater pumping could be corrected 
                # for an unconfined aquifer using Jacobs correction.
                # Peterson Feb 2013.
                h_star[:, iOutputColumns] = transform_h_star(self.parameters.(char(components[i])), [time_points, h_star[:, iOutputColumns]])
                    
                # Add output name to the cell array
                if np.ischar(self.variables.(components[i]).forcingData_colnames):
                    colnames[iOutputColumns] = components[i]
                else:
                    colnames[iOutputColumns] = [components[i], ' - ', self.variables.(components[i]).forcingData_colnames[j]]
        
        # Sum all components (excluding soil moisture) and add time
        # vector.
        if np.shape(h_star)[2]>1:
            h_star = [time_points, np.sum(h_star,2), h_star]
            colnames = ['time', 'Head', colnames[:]]
        else:
            h_star = [time_points, h_star]
            colnames = ['time', 'Head']               

        return h_star, colnames