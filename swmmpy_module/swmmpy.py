'''Python tools to interact with EPA SWMM.
Chloé Fandel, 2018.
cfandel@email.arizona.edu

Uses a template .inp file to enable editing data in different sections (junctions, conduits, etc.) of the input file. 
Data can be added from .csv files, or directly from pandas dataframes.
Can also run SWMM directly from Python.

See accompanying notebook for examples of how to use.'''


########################################################################################
import pandas as pd
import numpy as np
import subprocess as sp


#########################################################################################

def import_template(template_filename='template.inp'):
    
    '''Imports the template.inp file into a pandas dataframe.
    Inputs:
    template_filename: filename/path for the template.inp file
                       this must be a text file correctly formatted for SWMM, with placeholder strings in the sections to be edited
                       to create a sample template file, use the tutorial1.inp file from the SWMM documentation, and replace the data
                       in the [JUNCTIONS] section with the word 'junctions'
    Outputs:
    template:          pandas dataframe of the template file, with one column, and a line for each line of text.'''
    
    template = pd.read_csv(template_filename, header=None, skip_blank_lines=False)  #import template .inp file to dataframe
    
    return template


#########################################################################################

def import_data(data_filename):
    
    '''Imports data from a csv file into a pandas dataframe with columns.
    Inputs:
    data_filename:     filename/path for the data.csv file for a section
                       this must be a csv file of data formatted with the first row as a comment row, 
                       and the second row with the correct SWMM column headers for that section 
                       (see SWMM documentation for description of data in each section)
    Outputs:
    template:          pandas dataframe of the data file, with columns, and a line for each item'''
    
    data = pd.read_csv(data_filename, header=1)    #import data to dataframe that has separate columns
    
    return data


#########################################################################################

def insert_data(template, placeholder, data, show=False):
    
    '''Inserts data from a pandas dataframe into a template dataframe by replacing the section placeholder in the template.
    The template file and data file must be read in first. Best to make a copy of the template.
    
    Inputs:
    template:     pandas dataframe with one column, and one line per line of text in the template input file.
                  NOTE: each section that will be edited MUST have a section placeholder string at the desired insert location, 
                  and a section header string with dashes for the column spacing, like so: ;;---- ------ ------
    placeholder:  string corresponding to the placeholder string to be replaced
    data:         pandas dataframe with rows of data split into columns, to be inserted at the placeholder location
                  NOTE: the first column will get converted to an integer
    show:         True/False. If True, prints the updated section
    
    Outputs:
    new:          pandas dataframe with one column, and one line per line of text in the desired output file, 
                  with the new data inserted where the placeholder was in the original file
    '''
    
    #Get row index and location of placeholder:
    ind = template.index[template[0]==placeholder]              #df.index gets the index, df[0] looks in column 0 (returns an index object not just the name)
    loc = template.index.get_loc(ind[0])                        #get the integer position of the index specified (need to select first item (index name) in index object)

    #Convert data to formatted line strings:
    #Create format string based on template file:
    fmt_line = template[0].loc[loc-1]                           #get string to model formatting on (the row of dashes in the template file)
    cols = fmt_line.split()                                     #split string into chunks separated by whitespace
    l = [len(item) for item in cols]                            #get length of each column
    form = ''                                                   #create empty string to fill
    form = [form + '{d['+ str(i) + ']:<' + str(l[i]+1) + '}' for i in range(len(l))]  #concatenate formatting info and return a list of format strings (one per column)
    form = ''.join(form)                                        #join format strings into one for entire row

    #Insert values into format string
    data = data.round(6)                                #round data to 6 decimal places (to correct for artifacts caused by binary/float conversions)
    data_strings = pd.DataFrame(columns=[0])            #create empty dataframe to fill
    for ind in data.index:                              #loop over lines of data
        dl = data.loc[ind].tolist()                     #get line of data to be formatted
        if type(dl[0])!=type('str'):                    #if name isn't a string
            dl[0] = int(dl[0])                          #make sure name column is an integer not a float (WHY is this not automatic?)
        line = form.format(d=dl)                        #insert each item in list into format string
        data_strings.loc[ind] = line                    #insert line string into new dataframe

    #Replace placeholder with new data strings:
    #Split original df into two, one for everything before the placeholder, one for everything after:
    dfA = template[template.index < loc]                            #create df for everything above placeholder
    dfB = template[template.index > loc]                            #create df for everything below placeholder

    #Append the three dfs to each other (part above, part to insert, part below):
    new = dfA.append(data_strings, ignore_index=True)               #append additional part to top part
    new = new.append(dfB,          ignore_index=True)               #append bottom part to new df
    
    if show==True:
        print(new.iloc[loc-2:loc+len(data_strings)])               #print updated section if desired

    return new


#########################################################################################

def write_input(inputfile, placeholders, data_filenames, template_filename='template.inp'):
    
    '''Write a SWMM input file using a template.inp file with placeholders in each sections,
    and data.csv files for each section of data. The placeholders will be replaced with data.
    
    Inputs:
    inputfile:           string for name of input file to write (must end in .inp). Example: project.inp
    placeholders:        list of placeholder strings to be replaced with data. Ex: ['junctions', 'conduits']
    data_filenames:      list of file name strings for data.csv files to insert into the input file. Ex: ['junctions.csv', 'conduits.csv']
    template_filename:   string for name of template file (defaults to 'template.inp')
    
    Outputs:
    project.inp:         SWMM5 input text file
    '''
    
    template = import_template(template_filename)                               #import template file from txt
    for i in range(len(placeholders)):                                          #loop over list of placeholders
        data = import_data(data_filenames[i])                                   #import data to insert from csv
        template = insert_data(template,  placeholders[i], data)                #replace placeholder string with data
    template.to_csv(inputfile, header=False, index=False, quoting=3)            #write dataframe to .inp text file with specified name
        

#########################################################################################

def run(inputfile, reportfile, outputfile, exe='swmm5.exe'):
    
    '''Run SWMM using the specified input file, and create the specified output files.
    Inputs:
    inputfile:  name of .inp file to use with SWMM (must be formatted according to SWMM documentation)
    reportfile: name of .rpt file that SWMM will write to
    outputfile: name of .out file that SWMM will write binary outputs to
    exe:
    
    Outputs:
    project.rpt: SWMM report file
    project.out: SWMM output file
    '''
    
    p = sp.Popen([exe, inputfile, reportfile, outputfile], stdout=sp.PIPE, universal_newlines=True)          #run SWMM (and report process output)
    for line in p.stdout:          #for each line of process output
        if 'hour' not in line:     #if the line doesn't include the string 'hour' (to avoid a huge mass of text for each timestep)
            print(line) 
    