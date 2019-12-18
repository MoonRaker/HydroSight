
import numpy as np


def findClassDefsUsingAbstractName( abstractName, model_file_name):

    # Set some file names to ignore. This is undertaken because
    # requiredFilesAndProducts() (using matlab 2014b) appears to give inconsistent
    # results and for to reduce GUI start up time.
    fnames_ignore = ['model_TFN',
                     'example_TFN_model',
                     'HydroSight',
                     'ipdm',
                     'cmaes',
                     'fminsearchbnd',
                     'variogram',
                     'variogramfit',
                     'forcingTransform_abstract',
                     'derivedForcingTransform_abstract',
                     'stochForcingTransform_abstract',
                     'responseFunction_abstract',
                     'derivedResponseFunction_abstract',
                     'doIRFconvolution',
                     'findAbstractName',
                     'model_abstract']
    fnames_ignore = np.unique([abstractName, fnames_ignore ])

    # Check 'abstractName' is an abstract class 
    classNames = []
    if abstractName=='handle':
        return
    
    metaclass = meta.class.fromName(abstractName)
    if any(properties(metaclass)=='Abstract'):
        try:
            abstractNames[1] = metaclass.Name
        except:
            abstractNames[1] = ''            
    else:
        # Loop though all methods. If all methods are abstract, then 
        # 'abstractName' is deemed as an abstract classdef.
        isAbstract = True
        for i in range(len(metaclass.Methods):
            if (~metaclass.Methods{i,1}.Abstract) & (metaclass.Methods[i,1].DefiningClass.Name==metaclass.Methods[i,1].Name)) & (metaclass.Methods[i,1].DefiningClass.Name=='handle'):
                isAbstract = False
                break
        if ~isAbstract:
            return
    
    # Find path to the specified model.
    if exists(model_file_name):
        modelPath = fileparts(which(model_file_name))        
    else:
        print 'The following model file does not exist: ', model_file_name
    
    # Get list of all .m files listed within the model folder
    from os import listdir
    allFoldersFiles = listdir(fullfile(modelPath, ['**',filesep,'*.m']))
    
    # Get just the file names
    all_m_Files = cell(len(allFoldersFiles))
    for i in range(len(allFoldersFiles):
        pathstr, name = fileparts(allFoldersFiles[i].name)
        all_m_Files[i] = name
    
    isFileFromAbstract = False(size(allFoldersFiles))
    
    # Determine which version of depfun to use. Post Matlab 2014b, depfun
    # was removed.
    useDepFun = year(version('-date'))<=2014
    
    # Loop through each file name and assess if the file is dependent on
    # the specified abstract
    for i in range(len(allFoldersFiles):

        if ~any(cellfun(@(x) ~isempty(x), strfind( fnames_ignore , all_m_Files{i}))):

           # Get list of dependent function                              
           if useDepFun:
               depfunlist = depfun(all_m_Files(i),'-quiet')               
           else:
               depfunlist = matlab.codetools.requiredFilesAndProducts(all_m_Files(i))
               depfunlist = depfunlist              

           # Find if there is dependence upon the required abstract name
           if any(cellfun(@(x) ~isempty(x), strfind( depfunlist, abstractName))):
               isFileFromAbstract(i) = True               
    
    # File the list of class names to those dependent upon the abstract
    classNames = all_m_Files(isFileFromAbstract)
    
    return classNames

