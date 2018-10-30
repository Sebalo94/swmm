'''SWMM5 tools for reading and modifying the input file.'''

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
    
        
# FUNCTIONS:
########################################################################
########################################################################
def read_infile(filename):
    
    '''Reads the SWMM input file and returns a list of line strings.'''
    
    lines = []                              #initialize empty list of line strings
    with open (filename, "r") as file:      #open file for reading
        for line in file:                   #loop over each line
            lines.append(line)              #add each line to the list of line strings 
    return lines                            #return a list of line strings


########################################################################
def write_infile(filename, lines):
    
    '''Writes a new SWMM input file from a list of strings.'''
    
    with open(filename, 'w') as f:  #create and open new text file for writing
        for line in lines:          #loop over lines to write
            write = f.write(line)   #write each line and then go to a new line


########################################################################
def convert(string):
    
    '''Detects whether a string can be converted to a float. If yes, returns 
    the float. If no, returns the original string.'''
    
    try:
        float(string)
        return float(string)
    except ValueError:
        return string
    
    
########################################################################
def read_section(lines, section):
    
    '''Extracts information from the specified section of the SWMM5 input file.
    
    Inputs:
    
    lines: list of line strings from input file (created by read_input function)
    section: string indicating the name of the section to read (ex: 'JUNCTIONS, 'CONDUITS', etc.)
    
    Returns:
    
    items: list of item strings of item info
    itemsdic: dictionary of floats of item info keyed to item name strings
    keys: item name key strings
    lnums: list of line numbers for each item
    ''' 
    
    start_str = section                     #pull section start string from input
    end_str = section_dic[section]          #pull section end string from dictionary
    
    for line_num,line in enumerate(lines):  #loop over each index,string pair for lines in list file
        if line.find(start_str) > -1:       #find start of desired section
            start = line_num + 3            #save line number for start of desired section
        if line.find(end_str) > -1:         #find start of next section
            end = line_num - 1              #save line number of end of desired section
    lnums = range(start,end+1)              #save list of line numbers in the section

    keys = []                               #create empty list for item info
    itemsdic = {}                           #create empty dictionary for item info
    for line in lines[start:end]:           #loop over lines in desired section
        splitline = line.split()            #split line into its components (break at spaces)
        floats = [convert(splitline[i+1]) for i in range(len(splitline)-1)] #get values & convert to floats
        keys.append(splitline[0])
        itemsdic[splitline[0]] = floats     #assign info to dictionary keyed by item name string
    
    return [itemsdic, keys]   #return item info (list & dic), key strings, & line numbers


##########################################################################
def batch_create(keys, dic, Class):
    
    '''Batch create junction objects from the information in the junction dictionary.
    
    Inputs:
    
    keys: strings for item names
    dic: dictionary of object parameter values (floats)
    Class: class of object to create (ex: Junction, Coordinate)
    '''
    
    Objects = []                    #create empty list to store objects
    for key in keys:                #loop over dictionary keys
        O = Class(key,dic[key])     #use selected class to create a new object from dictionary values
        Objects.append(O)           #append the new object to the list
         
    return Objects                  #return list of objects


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


###################################################################
def update_infile(filename,section,oldlines,Items):
    
    '''Write a new SWMM input file that includes added or updated items.
    
    Inputs:
    
    filename: name string for the new file to be written (must have extension .inp)
    section: string indicating the name of the section to update (ex: 'JUNCTIONS, 'CONDUITS', etc.)
    oldlines: list of line strings read in from old input file (using read_infile)
    Items: dictionary of new items to add with format Items = {'Name':[par1, par2, par3...]}
    '''
    
    #Find desired section lines:
    start_str = section                     #pull section start string from input
    end_str = section_dic[start_str]        #pull section end string from dictionary
    
    for line_num,line in enumerate(oldlines):  #loop over each index,string pair for lines in list file
        if line.find(start_str) > -1:       #find start of desired section
            start = line_num + 3            #save line number for start of desired section
        if line.find(end_str) > -1:         #find start of next section
            end = line_num - 1              #save line number of end of desired section
    lnums = range(start,end+1)              #save list of line numbers in the section
    
    #Create new line strings:
    newlines = []                           #create empty list of new lines to be added
    for I in Items:                         #loop over new items to be added
        #choose desired section & create new formatted line string for each item:
        if section == 'JUNCTIONS':
            newline = "{:<17}{:<11}{:<11}{:<11}{:<11}{:<11}\n".format(I.Name,I.Elevation,
                I.MaxDepth,I.InitDepth,I.SurDepth,I.Aponded)
        if section == 'CONDUITS':
            newline = "{:<17}{:<17}{:<17}{:<11}{:<11}{:<11}{:<11}{:<11}{:<11}\n".format(I.Name,
                I.FromNode,I.ToNode,I.Length, I.Roughness, I.InOffset, I.OutOffset, I.InitFlow, 
                I.MaxFlow)
        if section == 'COORDINATES':
            newline = "{:<17}{:<19}{:<19}\n".format(I.Name,I.x,I.y)
        if section not in ['JUNCTIONS','CONDUITS','COORDINATES']:
            print('Error: Section does not exist or is not yet enabled.')  #print error message
        
        newlines.append(newline)                            #append new line to list of new lines

    #Update lines in file:
    oldlines[lnums[0]:lnums[-1]] = newlines                #insert new lines to list of lines
    write_infile(filename, oldlines)                       #write new input file
    
    
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