
import numpy as np


def findAbstractName(className):

    abstractNames = []
    if className=='handle':
        return
    
    metaclass = meta.class.fromName(className) # <<< ???
    if any(strcmp(properties(metaclass),'Abstract')):
        abstractNames[1]=metaclass.Name
    else:
        # Loop though all methods. If all methods are abstract, then the
        # class is deemed as an abstract classdef.
        if isAbstract:
            for i in range(np.shape(metaclass.Methods,1)):
                if ((~metaclass.Methods{i,1}.Abstract) & 
                    (metaclass.Methods[i,1].DefiningClass.Name!=metaclass.Methods[i,1].Name) &
                    (metaclass.Methods[i,1].DefiningClass.Name!='handle')):
                    if ~isAbstract:
                        break
        if isAbstract:
            abstractNames[1] = metaclass.Name
    
    # Get names of superclasses. Note: The field names of chnages between
    # matlab 2010 and matlab 2014a. The following attempts to handle both
    # formats.
    if any(properties(metaclass)=='SuperclassList'):
        SuperclassList=metaclass.SuperclassList
    elif any(properties(metaclass)=='SuperClasses')):
        SuperclassList = metaclass.SuperClasses
        
    # Assign results to output cell array.
    for i in range(len(SuperclassList)):
        if iscell(SuperclassList):
            tmp = findAbstractName(SuperclassList[i].Name)      
        else:
            tmp = findAbstractName(SuperclassList[i].Name)
        if ~isempty(tmp):
            abstractNames[size(abstractNames,1)+1:size(abstractNames,1)+size(tmp,1),1]=tmp

    return abstractNames

