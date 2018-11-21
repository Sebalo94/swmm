'''SWMM5 tools for reading and modifying the input file.
Note: input file must have only one line for each section header, and no empty columns.
Formats are based on the SWMM documentation.'''

# IMPORTS:
##########################################################################
##########################################################################
import pandas as pd
import numpy as np


# REFERENCE DICTIONARIES:
###########################################################################
###########################################################################

#List of section names:
section_names = ['junctions','outfalls','conduits','xsections','losses','inflows','timeseries','coordinates',
                 'polygons','symbols']

#Dictionary for section markers (what to look for in the .inp file):
key_marker_dic = {'junctions':'[JUNCTIONS]','outfalls':'[OUTFALLS]','conduits':'[CONDUITS]',
                'xsections':'[XSECTIONS]','losses':'LOSSES','inflows':'[INFLOWS]',
                'timeseries':'[TIMESERIES]','coordinates':'[COORDINATES]',
                'polygons':'[Polygons]','symbols':'[SYMBOLS]'}

#Dictionary for the markers indicating start of next section in .inp file:
start_end_dic = {'[TITLE]':'[OPTIONS]', '[OPTIONS]':'[EVAPORATION]', '[EVAPORATION]':'[RAINGAGES]', 
               '[RAINGAGES]':'[SUBCATCHMENTS]', '[SUBCATCHMENTS]':'[SUBAREAS]', '[SUBAREAS]':'[INFILTRATION]', 
               '[INFILTRATION]':'[JUNCTIONS]', '[JUNCTIONS]':'[OUTFALLS]', '[OUTFALLS]':'[CONDUITS]', 
               '[CONDUITS]':'[XSECTIONS]', '[XSECTIONS]':'[POLLUTANTS]', '[POLLUTANTS]':'[LANDUSES]',
               '[LANDUSES]':'[COVERAGES]', '[COVERAGES]':'[LOADINGS]', '[LOADINGS]':'[BUILDUP]', 
               '[BUILDUP]':'[WASHOFF]', '[WASHOFF]':'[TIMESERIES]', 
               '[TIMESERIES]':'[REPORT]', '[REPORT]':'[TAGS]', '[TAGS]':'[MAP]', '[MAP]':'[COORDINATES]',                    '[COORDINATES]':'[VERTICES]', '[VERTICES]':'[Polygons]', '[Polygons]':'[SYMBOLS]'}

#List of indices for section format dataframe:
section_indices = ['cnames','format_strings'] 

#List of format strings for sections:
section_formats = [[['Name','Elevation','MaxDepth','InitDepth','SurDepth','Aponded'], 
                     '{Name:<17}{Elevation:<11}{MaxDepth:<11}{InitDepth:<11}{SurDepth:<11}{Aponded:<11}'],
    [['Name','Elevation','Type','Gated'], 
      '{Name:<17}{Elevation:<11}{Type:<11}{Gated:<11}'],
    [['Name','FromNode','ToNode','Length','Roughness','InOffset','OutOffset','InitFlow','MaxFlow'], 
      '{Name:<17}{FromNode:<17}{ToNode:<17}{Length:<11}{Roughness:<11}{InOffset:<11}{OutOffset:<11}{InitFlow:<11}{MaxFlow:<11}'],
    [['Name','Elevation','MaxDepth','InitDepth','SurDepth','Aponded'], 
      '{:<17}{:<11}{:<11}{:<11}{:<11}{:<11}'],
    [['Name','Elevation','MaxDepth','InitDepth','SurDepth','Aponded'], 
      '{:<17}{:<11}{:<11}{:<11}{:<11}{:<11}'],   
    [['Name','Elevation','MaxDepth','InitDepth','SurDepth','Aponded'], 
      '{:<17}{:<11}{:<11}{:<11}{:<11}{:<11}'],
    [['Name','Elevation','MaxDepth','InitDepth','SurDepth','Aponded'], 
      '{:<17}{:<11}{:<11}{:<11}{:<11}{:<11}'],
    [['Name','Elevation','MaxDepth','InitDepth','SurDepth','Aponded'], 
      '{:<17}{:<11}{:<11}{:<11}{:<11}{:<11}'],
    [['Name','Elevation','MaxDepth','InitDepth','SurDepth','Aponded'], 
      '{:<17}{:<11}{:<11}{:<11}{:<11}{:<11}'],
    [['Name','Elevation','MaxDepth','InitDepth','SurDepth','Aponded'], 
      '{:<17}{:<11}{:<11}{:<11}{:<11}{:<11}']]

#Compiled dataframe of section names and formats:
section_info = pd.DataFrame(section_formats, columns=section_indices, index=section_names)

        
# FUNCTIONS:
########################################################################
########################################################################
def read_input(inputfile):
    
    ''' Read entire input file into a dataframe. 
    
    Inputs:
    inputfile: SWMM text input file name (ex: 'template.inp'). 
               If being used as a template, file must have section placeholder strings.
    
    Outputs:
    input_lines: dataframe of input file lines as strings''' 
    
    input_lines = pd.read_csv(inputfile,sep='\t',names=['cols'],skip_blank_lines=False)    #read text file to df
    return input_lines


########################################################################
def edit_section(input_lines, sections, data, new_filename='new_input.inp'):
    
    ''' Read in the desired sections from the template input file, and replace the placeholders
    with new information.
    
    Inputs:
    df: dataframe of line strings for the template input file
    section: list of strings indicating the section to edit 
             (must be one of: junctions, outfalls, conduits, xsections, losses, inflows, 
             timeseries, coordinates, polygons, symbols)    
    data: dic of lists of parameters for each item to add, for each section. 
          Example: data = {'junctions': [['J5',87.,4.,0.,0.,0.], ['J6',86.,4.,0.,0.,0.]],
                           'outfalls':  [['Out1', 85, 'FREE', 'NO','']]} 
    new_filename: optional string for name of new input file (should end in .inp)
    
    Outputs:
    Section_Data: dic of dataframes for the new items in each section. Each section's info can be
               accessed by section keyword, and each item in that section can then be accessed
               by name or by column
    new_input.inp: new csv SWMM input file with updated info
    '''
    
    Section_Data = {}                               #empty dic to store section items
    for i in range(len(sections)):                  #loop over sections
        #Create dataframe for new items:
        section = sections[i]                       #select current section
        section_data = data[section]                #select data for current section
        index = input_lines.index[(input_lines.values==section).any(1)][0]  #get section index
        new_indices = np.arange(index,index+1, 1/len(section_data)) #set new ind to insert between rows 
        colnames = section_info.cnames[section]                                 #get section column names
        Items = pd.DataFrame(section_data, columns=colnames, index=new_indices) #create new df for items 
    
        #Convert to strings:
        formatstr = section_info.format_strings[section]        #get format string for current section
        Istrings = []                                           #create empty list to store item strings
        for i in new_indices:                                   #loop over rows (indices) of new items to add
            I = Items.loc[i,:]                                  #get row for current item
            string = formatstr.format(**I)                      #convert row to a single string
            Istrings.append(string)                             #append string to list of line strings
        input_col = input_lines.columns[0]                      #get original dataframe column
        new_lines = pd.DataFrame(Istrings,columns=[input_col],index=new_indices)   #create new df with the strings to write

        #Update the df of lines from the original input file:
        input_lines = input_lines.drop(index)                               #remove section placeholder string
        input_lines = input_lines.append(new_lines, ignore_index=False)     #append new section lines
        input_lines = input_lines.sort_index().reset_index(drop=True)       #re-sort by index
        
        #Save new data for each section:
        Section_Data[section] = Items                                       #store current section's items
        
    #Write a new input file:    
    input_lines.to_csv(new_filename, sep='\t', header=False, index=False, na_rep=' ') #write the new text file
        
    return Section_Data


######################################################################
def read_sections(input_lines,sections):
    '''Read data from the input file for the specified sections, and return a 
    dataframe split into columns for each section. Relies on formatting dictionaries.
    
    Inputs:
    input_lines: dataframe of strings for each line of the input file.
    sections: list of strings indicating which sections to read. 
            (must be one of: junctions, outfalls, conduits, xsections, losses, inflows, 
             timeseries, coordinates, polygons, symbols)
             
    Outputs:
    Section_Data: dictionary of dataframes containing the data for each section,
                  keyed by section name'''   
    
    Section_Data = {}                               #create empty dic to store section data
    for i in range(len(sections)):                  #loop over sections
        section = sections[i]                       #select current section
        start_marker = key_marker_dic[section]      #get marker string for start of section
        start = input_lines.index[(input_lines.values==start_marker).any(1)][0]+4 #get index for start of section
        end_marker = start_end_dic[start_marker]    #get marker string for start of next section
        end = input_lines.index[(input_lines.values==end_marker).any(1)][0]-1     #get index for end of section
        
        data = input_lines[start:end]               #get data for current section
        cnames = section_info.cnames[section]       #get column names for current section
        data = pd.DataFrame(data.cols.str.split().tolist(),columns=cnames) #split line string into columns
        Section_Data[section] = data                #store current section dataframe in dictionary
    return Section_Data

######################################################################
def write_new_input(input_lines, sections, Section_Data, new_filename='new_input.inp'):
    '''Writes a new input file with updated sections based on the information in Section_Data.
    
    Inputs:
    input_lines: dataframe of line strings from old input file
    sections: list of string keys for sections to update
              (must be one of: junctions, outfalls, conduits, xsections, losses, inflows, 
              timeseries, coordinates, polygons, symbols)
    Section_Data: dictionary of dataframes with section data, keyed to section names
    new_filename: name string for new input file (should end in .inp)
    
    Outputs:
    new_input.inp: updated text SWMM input file'''

    for i in range(len(sections)):
        #Get section info:
        section = sections[i]                       #select current section
        data = Section_Data[section]                #select data for current section
        start_marker = key_marker_dic[section]      #get marker string for start of section
        start = input_lines.index[(input_lines.values==start_marker).any(1)][0]+4 #get index for section start 
        end_marker = start_end_dic[start_marker]    #get marker string for start of next section
        end = input_lines.index[(input_lines.values==end_marker).any(1)][0]-1     #get index for section end
        
        old_indices = list(input_lines.iloc[start:end].index.values) #get list of old indices
        new_indices = np.arange(start,start+1, 1/len(data)) #set new ind to insert between rows 
    
        #Convert to strings:
        formatstr = section_info.format_strings[section]      #get format string for current section
        Istrings = []                                         #create empty list to store item strings
        for i in data.index:                                  #loop over rows (indices) of new items to add
            I = data.loc[i,:]                                 #get row for current item
            string = formatstr.format(**I)                    #convert row to a single string
            Istrings.append(string)                           #append string to list of line strings
        input_col = input_lines.columns[0]                    #get original dataframe column
        new_lines = pd.DataFrame(Istrings,columns=[input_col],index=new_indices)   #create new df with the strings to write
        
        #Update the df of lines from the original input file:
        input_lines = input_lines.drop(old_indices)                         #remove section placeholder string
        input_lines = input_lines.append(new_lines, ignore_index=False)     #append new section lines
        input_lines = input_lines.sort_index().reset_index(drop=True)       #re-sort by index
    
    #Write a new input file:    
    input_lines.to_csv(new_filename, sep='\t', header=False, index=False, na_rep=' ') #write the new text file

    
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