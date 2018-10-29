'''SWMM5 tools for reading and modifying the input file.'''

# IMPORTS:
##########################################################################
##########################################################################
import pandas as pd
import numpy as np


# REFERENCE DICTIONARIES:
###########################################################################
###########################################################################

# Dictionary of section start and end strings:
section_dic = {'TITLE':'OPTIONS', 'OPTIONS':'EVAPORATION', 'EVAPORATION':'RAINGAGES', 
               'RAINGAGES':'SUBCATCHMENTS', 'SUBCATCHMENTS':'SUBAREAS', 'SUBAREAS':'INFILTRATION', 
               'INFILTRATION':'JUNCTIONS', 'JUNCTIONS':'OUTFALLS', 'OUTFALLS':'CONDUITS', 
               'CONDUITS':'XSECTIONS', 'XSECTIONS':'POLLUTANTS', 'POLLUTANTS':'LANDUSES',
               'LANDUSES':'COVERAGES', 'COVERAGES':'LOADINGS', 'LOADINGS':'BUILDUP', 
               'BUILDUP':'WASHOFF', 'WASHOFF':'TIMESERIES', 'TIMESERIES':'REPORT', 'REPORT':'TAGS', 
               'TAGS':'MAP', 'MAP':'COORDINATES', 'COORDINATES':'VERTICES', 'VERTICES':'Polygons', 
               'Polygons':'SYMBOLS'}


# CLASSES:
#########################################################################
#########################################################################
class Junction:
    '''Junction class for SWMM junctions (drainage system nodes where pipes join together)
    Junctions can receive external inflows, become pressurized, lose water, or experience ponding. 
    
    Parameters listed in SWMM User Manual:
    - invert elevation (elevation of bottom of junction above the reference level)
    - height to ground surface
    - ponded surface area when flooded (optional)
    - external inflow data (optional)
    
    Attributes:
    - Elevation
    - MaxDepth
    - InitDepth
    - SurDepth
    - Aponded
    '''
    
    def __init__(self, name, pars):
        '''Create a new junction with the attribute values provided.
        name = name string
        pars = [elev, maxdepth, initdepth, surdepth, Apond]'''
        self.Name = name
        self.Elevation = pars[0]
        self.MaxDepth = pars[1]
        self.InitDepth = pars[2]
        self.SurDepth = pars[3]
        self.Aponded = pars[4]
        
    def setElev(self,elev):
        'Change the elevation of the junction'
        self.Elevation = elev
    
    def setMaxDepth(self,maxdepth):
        'Change the max depth of the junction'
        self.MaxDepth = maxdepth
        
    def setInitDepth(self,initdepth):
        'Change the initial depth of the junction'
        self.InitDepth = initdepth
    
    def setSurDepth(self,surdepth):
        'Change the surcharge depth of the junction'
        self.SurDepth = surdepth
    
    def setAponded(self,Apond):
        'Change the ponded area of the junction'
        self.Aponded = Apond

########################################################################
class Coordinate:
    '''Coordinate class for SWMM coordinates ((x,y) values for locations of nodes (junctions or outfalls)) 
    
    Attributes:
    - Name
    - x coordinate
    - y coordinate
    '''
    
    def __init__(self, name, pars):
        '''Create a new coordinate with the attribute values provided.
        name = name string
        pars = [x,y]'''
        self.Name = name
        self.x = pars[0]
        self.y = pars[1]
        
    def setX(self,x):
        'Change the x coordinate of the node'
        self.x = x
    
    def setY(self,y):
        'Change the y coordinate of the node'
        self.y = y
        

#########################################################################        
class Conduit:
    '''Conduit class for SWMM conduit connections between nodes 
    
    Attributes:
    - Name           
    - FromNode        
    - ToNode          
    - Length     
    - Roughness  
    - InOffset   
    - OutOffset  
    - InitFlow   
    - MaxFlow 
    '''
    
    def __init__(self, name, pars):
        '''Create a new conduit with the attribute values provided.
        name = name string
        pars = [x,y]'''
        self.Name = name
        self.FromNode = pars[0]
        self.ToNode = pars[1]
        self.Length = pars[2]     
        self.Roughness = pars[3]  
        self.InOffset = pars[4]   
        self.OutOffset = pars[5]  
        self.InitFlow = pars[6]   
        self.MaxFlow = pars[7]
        
    def setFromNode(self,fnode):
        'Change the conduit origin node'
        self.FromNode = fnode
    
    def setToNode(self,tnode):
        'Change the conduit terminus node'
        self.ToNode = tnode
        
    def setLength(self,length):
        'Change the conduit length'
        self.Length = length
        
    def setRoughness(self,rough):
        'Change the conduit roughness'
        self.Roughness = rough
        
    def setInOffset(self,inoff):
        'Change the conduit in offset'
        self.InOffset = inoff
        
    def setOutOffset(self,outoff):
        'Change the conduit out offset'
        self.OutOffset = outoff
        
    def setInitFlow(self,initflow):
        'Change the initial flow to the conduit'
        self.InitFlow = initflow
        
    def setMaxFlow(self,maxflow):
        'Change the max flow through the conduit'
        self.MaxFlow = maxflow
    
        
########################################################################
def convert(string):
    
    '''Detects whether a string can be converted to a float. If yes, returns 
    the float. If no, returns the original string.'''
    
    try:
        float(string)
        return float(string)
    except ValueError:
        return string
    
    
##########################################################################
def batch_create(df,Class):
    
    '''Batch create objects from the information in the dictionary.
    
    Inputs:
    
    df: dataframe with item names and parameter values
    Class: class of object to create (ex: Junction, Coordinate)
    
    Example usage:
    Junction_objects = batch_create(Junctions,si.Junction)
    '''

    Objects = []                    #create empty list to store objects
    for i in range(len(df)):        #loop over objects
        name = df.iloc[i].Name      #get list of object names from df
        vals = df.iloc[i,1:]        #get list of object vals from df
        O = Class(name,vals)        #create object
        Objects.append(O)           #append to list
    return Objects    



##############################################################################
def add_objects(Objects, Class, newkeys, newdic):
    
    '''Add junction objects to the existing list.
    
    Inputs:
    
    Objects: list of existing objects (if none, pass an empty list)
    Class: class of object to create (ex: Junction, Coordinate)
    newkeys: strings for new item names
    newdic: dictionary of new parameter values (floats) ex:[Elevation,MaxDepth,InitDepth,SurDepth,Aponded]
    '''
    
    for key in newkeys:               #loop over dictionary keys
        O = Class(key,newdic[key])   #use selected class to create a new object from dictionary
        Objects.append(O)             #append the new object to the list
         
    return Objects                    #return list of objects


######################################################################
def plot_conduits(Coordinates,Conduits,coordinate_dic,label=True):
    
    '''Plot a map of the system nodes and conduits.
    
    Inputs:
    
    Coordinates: list of coordinate objects
    Conduits: list of conduit objects
    coordinate_dic: dictionary of coordinate parameters
    label: optional labels (defaults to on, set to False to turn off)
    '''
    
    from matplotlib import pyplot as plt    #import matplotlib 

    # Generate node x and y datasets:
    X = []                                  #create empty list for X coordinates
    Y = []                                  #create empty list for Y coordinates
    for i in range(len(Coordinates)):       #loop over coordinate objects
        X.append(Coordinates[i].x)          #pull & store x values
        Y.append(Coordinates[i].y)          #pull & store y values

    # Plot nodes:
    fig, ax = plt.subplots()                #create figure and axis objects
    ax.plot(X,Y,'.')                        #plot nodes as points
    if label == True:                       #if labeling option is selected
        offset = (max(X)-min(X))*.01        #calculate how much to offset labels (proportional to scale)
        for i in range(len(Coordinates)):   #loop over nodes & add node labels
            ax.annotate(Coordinates[i].Name, (X[i]+offset,Y[i]+offset)) 

    # Plot conduits:
    for i in range(len(Conduits)):                          #loop over each conduit
        fromx = coordinate_dic[Conduits[i].FromNode][0]     #pull coordinates for start point
        fromy = coordinate_dic[Conduits[i].FromNode][1]
        tox = coordinate_dic[Conduits[i].ToNode][0]         #pull coordinates for end point
        toy = coordinate_dic[Conduits[i].ToNode][1]
        X = [fromx,tox]                                     #group in a list
        Y = [fromy,toy]
        ax.plot(X,Y,'k--')                                  #plot connecting line
        if label == True:                                   #if label option is selected
            mpx = fromx + (tox-fromx)/2                     #calculate midpoint location
            mpy = fromy + (toy-fromy)/2
            ax.annotate(Conduits[i].Name, (mpx,mpy))        #add labels at midpoint