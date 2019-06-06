'''Mapping module for GemPy.
Chloé Fandel 2019.
Functions to visualize and export GemPy model outputs.
2D maps, custom diagonal cross-sections, export to GSLIB and VTK, import from GSLIB.
Based on original code from Elisa Heim.

See accompanying notebook for examples of how to use.'''

#Imports:
import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import gdal
import pyevtk
from copy import copy

#sys.path.append("../../..")   #optional: if gempy has been downloaded from GitHub rather than installed normally, look for it in the folders above the current folder
import gempy as gp


#############################################################################
def importDEM(filename, show=True):
    '''Import DEM from a tif file using gdal package.
    Return a dem object, and xyz extent and resolution.
    (this can be used to set the model extent)
    NOTE: vertical (z) resolution can't be extracted from the raster!
    
    filename: string indicating the filename (must be a rectangular tif)
    show:     option to show a plot of the DEM or not.
    
    Returns:     grid_info = [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz]  

    dem:      gdal dem objecct
    dema:     array of elevation values of dim: xres,yres
    xmin:     minimum x value (same for ymin, zmin)
    xmax:     maximum x value (same for ymax, zmax)
    xres:     x resolution, aka number of columns, aka number of cells along x axis (NOT pixel width)
    dx:       pixel width in x direction 
    etc.
    '''
    
    dem = gdal.Open(filename)    #DEM must be rectangular tif 
    dema = dem.ReadAsArray()     #copy of DEM as a numpy array (defaults to integers)
    dema = dema.astype(float)    #convert integer array to float array
    dema[dema==0] = np.nan       #replace zeros with NaNs (have to convert array to float first)

    ulx, pixelwidthx, xskew, uly, yskew, pixelheighty = dem.GetGeoTransform() #get resolution and coordinate info (for some reason the order of skew and pixel size is flipped for y axis?!)
    ncol = dem.RasterXSize            #number of columns (aka number of cells along x axis)
    nrow = dem.RasterYSize            #number of rows (aka number of cells along y axis)
    lrx = ulx + (ncol * pixelwidthx)  #lower right x coord = upper left x coord + (width of raster cells in x direction * number of raster cells in x direction)
    lry = uly + (nrow * pixelheighty)

    #Get min and max elevations (z):
    #note: gdal's built-in GetRasterBand and GetStatistics return an incorrect zmin (WHY?!)
    zmin = np.nanmin(dema)
    zmax = np.nanmax(dema)
    
    #Assign useful names:
    xmin = ulx
    xmax = lrx
    xres = ncol
    dx =   abs(pixelwidthx)
    ymin = lry
    ymax = uly
    dy =   abs(pixelheighty)
    yres = nrow
    zres = 'na'     #can't be extracted from raster

    #Print results & display raster:
    if show==True:
        print('Raster dimensions: \nxmin: {:<12} xmax: {:<12} xres: {} \nymin: {:<12} ymax: {:<12} yres: {} \nzmin: {:<12} zmax: {:<12} zres: {}'.format(
            xmin,xmax,xres,ymin,ymax,yres,zmin,zmax,zres))
        plt.imshow(dema, extent=(xmin,xmax,ymin,ymax), vmin=zmin, vmax=zmax) #plot raster as image
        #print(gdal.Info(dem))  #for more detailed file info, uncomment this line
        
    return dem,dema,xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres



#############################################################################
def get_surflith(dem, dema, interp_data, grid_info, output_filename='DEMxyz.csv'):
    '''Reshape DEM and use it to compute the lithology values of the GemPy model at the land surface.
    Returns an array of lith values at the surface z elevation for each xy point.
    
    dem:                dem object returned by importDEM() or dem = gdal.Open(filename)
    dema:               dem array returned by importDEM() or dema = dem.ReadAsArray()
    interp_data:        interpolated data returned by gempy.InterpolatorData()
    grid_info:          [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz] array of model grid and resolution info from importDEM()
    output_filename:    string to name gdal's output csv file (can be a throw-away - not used again)
    
    returns:
    surflith:           an array of lith values at the surface z elevation for each xy point, dim (yres,xres)'''
    
    #Get required grid info: 
    ##grid_info = [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz]
    xres = grid_info[2]
    yres = grid_info[6]
    
    #Get an array with xyz values from the DEM:
    #can this be streamlined to avoid having to export and re-import?
    translate_options = gdal.TranslateOptions(options = ['format'],format = "XYZ")  #set options for gdal.Translate()
    gdal.Translate(output_filename, dem, options=translate_options)  #convert dem to a csv with one column of points, each with an xyz value
    xyz = pn.read_csv(output_filename, header=None, sep = ' ')       #read xyz csv with pandas
    demlist = xyz.values  #convert to np array of (x,y,z) values with dim (ncol*nrow, 3)
    
    #Get and format the geologic data:
    surfall, fault2 = gp.compute_model_at(demlist, interp_data) #compute the model values at the locations specified (aka the land surface) (why is fault a required output?)
    surflith = surfall[0].reshape(dema.shape) #reshape lith block (a list) to an array with same dimensions as dem (yres,xres,zres) (note: xres*yres must equal length of lith)
    #now we have a discretized array with the same resolution as the dem, with a value for the lithology at the surface elevation for each xy point
    surflith[np.isnan(dema)] = np.nan       #crop to dem boundary (insert nan values everywhere where the dem has no data) 
    surflist = surflith.reshape(xres*yres)  #reshape cropped data to list format (with dimensions xres*yres)
    
    return surflith, surflist 



#############################################################################
def crop2elevation(lith, dema, grid_info):
    '''Discretizes lith block into an array matching the model dimensions (yres,xres,zres), i.e. (nrow,ncol,nlay), 
    then crops off values above the land surface and replaces them with np.nan.
    
    IMPORTANT: lith returned by gempy.compute_model() is two arrays (each with dim: model extent) with a formation number assigned to each cell, and the orientation of that formation at each cell
    The format of lith seems to change if fault is present or not, so need to index differently to get the right slice:
    fault present:    lith[0]
    no fault present: lith[0][0]
    
    lith:       array of lithological unit values of dimensions (slice of lith block array returned by gempy.compute_model() - either lith[0] if fault present, or lith[0][0] if no fault present)
    dema:       array of elevation values generated from DEM (use importDEM(), or dem = gdal.Open(filename) & dema = dem.ReadAsArray())
    grid_info:  [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz] array of model grid and resolution info from importDEM()
    
    returns:
    lithzcrop:  elevation-cropped array of lithologic unit indices of dimensions (nrow,ncol,nlay), i.e. (yres,xres,zres).'''

    #Get required grid info: 
    ##grid_info = [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz]
    xres = grid_info[2]
    yres = grid_info[6]
    zmin = grid_info[8]
    zmax = grid_info[9]
    zres = grid_info[10]
    dz =   grid_info[11]

    #Get lith array into a shape that matches the model dimensions:
    lith2 = lith.reshape([xres,yres,zres]) #reshape lith block (a list) to an array with same dimensions as model (xres,yres,zres) (note: length of lith mst equal mxres*myres*mzres) 
    lith3 = np.transpose(lith2,(1,0,2)) #transpose swaps order of axes (so here, I am flipping x and y (aka 0 and 1))
    lith4 = np.flip(lith3,axis=0)       #transpose didn't correctly map y (north-south) axis, so need to flip just that axis

    #Convert the DEM to an array of vertical cell indices:
    #i.e. how many cells (aka layers) up from the base of the model is each (x,y) location?
    zvals = np.linspace(zmin, zmax, zres)  #create linearly-spaced z values within model range
    zind = (dema - zmin) / dz              #calculate the cell index of each point in the dem array using the cell height (i.e. how many cells/layers up it is from the base)
    zind = zind.astype(int)                #convert to integers for use as vertical indices

    #Remove the model values above the land surface:
    lithzcrop = copy(lith4)                 #make a copy to avoid messing up original (is this necessary?)
    for row in range(yres):                 #loop over rows (y axis)
        for col in range(xres):             #loop over columns (x axis)
            z = zind[row,col]               #get z index at current row and col
            lithzcrop[row,col,z:] = np.nan  #assign nan to all cells greater than z of land surface
            
    return lithzcrop



#############################################################################
def crop2raster(lith, grid_info, rasterfilename, nanval=0):
    '''Crop the extent of geologic model to the extent of an irregularly-shaped imported raster with a set value indicating empty cells.
    
    lith:           array of lithologic unit indices of dimensions (yres,xres,zres) OR dimensions (yres,xres)
                    this can be the array of surface lith values, all lith values, or elevation-cropped lith values
    grid_info:  [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz] array of model grid and resolution info from importDEM()
    rasterfilename: string indicating the name of the raster file to use for cropping (can be bigger but not smaller than model extent) 
    nanval:         value indicating empty cells in raster file
    
    returns:
    lithxycrop:     array of same dimensions as input, but with empty cells filled with np.nan'''
    
    #Get grid info:
    xres = grid_info[2]
    yres = grid_info[6]
    
    #Import & format raster:
    geo = gdal.Open(rasterfilename)     #import raster file
    geoa = geo.ReadAsArray()            #read file as array (output will have an array for each color channel)
    geoa = geoa[-yres:,0:xres]          #if raster doesn't match DEM size, slice to size
    geoa = geoa.astype(float)           #convert integer array to float array
    geoa[geoa==0] = np.nan              #replace zeros with NaNs (have to convert array to float first)

    #Crop lith array:
    lithxycrop = lith.copy()                  #make a copy to avoid messing up original
    if lithxycrop.shape == geoa.shape:        #if arrays have same extent
        lithxycrop[np.isnan(geoa)] = np.nan   #crop lith to active cells in imported raster
    else:                                     #otherwise, assume lith is 3D
        lithxycrop[np.isnan(geoa),:] = np.nan #crop over all z cells
        lithxycrop = np.rot90(lithxycrop)     #rotate 90 degrees to give dimensions of (xres,yres,zres)
    
    return lithxycrop


########################################################################################################
def export2gslib(a, filename, grid):
    '''Exports a numpy array to a gslib file (GeoModeller-style format), using gempy grid objects to get the correct dimensions.
    
    Inputs: 
    a:        numpy array to be exported 
    filename: filename or path to save to
    grid:     gempy grid object to use for the dimensions of the array, obtained from geo_model.grid.gridnamehere. 
              A list of available grid names can be found with geo_model.grid.grid_types. 
              The active grids in this list can be seen with geo_model.grid.active_grids
    
    Output: 
    filename.gslib: a gslib file
    '''
    
    #Format data:
    a[np.isnan(a)] = 0              #assign zeros to nan values (for SKS)
    a = np.round(a)                 #round to integers
    a = a.astype(int)               #convert from floats to integers (so that gslib file will have integers)
    
    #Get grid info
    xres = grid.resolution[0]                          #get x resolution
    yres = grid.resolution[1]                          #get y resolution
    if len(grid.resolution)==3:                        #if 3D grid
        zres = grid.resolution[2]                      #get z resolution
       
    #Format array shape:
    try:                                                #try to reshape in 3D
        a = np.reshape(a,(xres,yres,zres))              #reshape 1D array to a 3D array with correct dimensions
        a = np.reshape(a, xres*yres*zres, order='F')    #reshape 3D array back to 1D array using Fortran indexing 
    except:                                             #if 3D doesn't work, try 2D 
        a = np.reshape(a,(xres,yres))                   #reshape 1D array to a 2D array with correct dimensions
        a = np.reshape(a, xres*yres, order='F')         #reshape 2D array back to 1D array using Fortran indexing 
        
    #Export:    
    df = pd.DataFrame(a)                            #store array in a pandas dataframe
    header = pd.DataFrame(['GemPy output',1,'lith'])    #set gslib file header
    df = header.append(df)                              #attach header and data
    df.to_csv(filename, header=False, index=False)      #write a text file in gslib format


#############################################################################
def importgslib(filename, grid):
    '''Imports a gslib file (GeoModeller-style) into a numpy array with dimensions taken from a GemPy grid object.
    
    Inputs:
    filename:  the name or path to the file to be imported
    grid:     gempy grid object to use for the dimensions of the array, obtained from geo_model.grid.gridnamehere. 
              A list of available grid names can be found with geo_model.grid.grid_types. 
              The active grids in this list can be seen with geo_model.grid.active_grids
              
    Returns:
    a:        2D or 3D numpy array representing the gslib data, with dimensions taken from the grid'''
    
    #Get grid info
    xres = grid.resolution[0]                          #get x resolution
    yres = grid.resolution[1]                          #get y resolution
    if len(grid.resolution)==3:                        #if 3D grid
        zres = grid.resolution[2]                      #get z resolution
    
    #Get data:
    a = pd.read_csv(filename, skiprows=2, dtype=float)  #read in gslib file to pandas df without header rows, as floats
    a = a.values                                        #get an array of the values
    a[a==0] = np.nan                                    #replace zeros with NaNs (array must be float first)

    #Reshape based on specified grid:
    try:                                             #for 3D files (blocks)
        a = np.reshape(a,(xres,yres,zres),order='F') #reshape to xyz grid using Fortran ordering
        a = np.rot90(a)                              #rotate 90 degrees CCW
    except:                                          #for 2D files (maps)
        a = np.reshape(a,(xres,yres),order='F')      #reshape to xyz grid using Fortran ordering
        #a = np.rot90(a)                              #rotate 90 degrees CCW
    
    return a



#############################################################################
def export2vtk(a, filename, grid):
    '''Exports a numpy array to a VTK file (for ParaView), using gempy grid objects to get the correct dimensions.
    
    Inputs: 
    a:        numpy array to be exported 
    filename: filename or path to save to
    grid:     gempy grid object to use for the dimensions of the array, obtained from geo_model.grid.gridnamehere. 
              A list of available grid names can be found with geo_model.grid.grid_types. 
              The active grids in this list can be seen with geo_model.grid.active_grids
    
    Output: 
    filename.vtr: a VTK file
    '''
    #Get resolution info
    xres = grid.resolution[0]  
    yres = grid.resolution[1]
    zres = grid.resolution[2]
    
    #Reshape if needed:
    try:                                                #try to reshape to 3D
        a = np.reshape(a,(xres,yres,zres))              #reshape 1D array to a 3D array with correct dimensions
    except:
        print('array must be 3D')                       #print error message
    
    #Export:
    pyevtk.hl.gridToVTK(filename, np.arange(xres+1), np.arange(yres+1), np.arange(zres+1), cellData={'data': a}) #export to VTK
    
    
#############################################################################
def plot_map(a, geo_model, plot_data=False, ref_points=None, ax=None):
    '''Plot 2D geologic map (generated by cropping lith_block with topography).
    Inputs:
    a:          Data representing the geologic map. Can be a 1D or 2D numpy array, or the path to a gslib file.
                To get 1D array: a = sol.geological_map. To get 2D array: a = np.reshape(sol.geological_map, (yres,xres)) (may need to also flip up-down)
                To import a gslib file: a = mapping.importgslib(filename, grid)
    plot_data:  Whether or not to display the input data on the map. Defaults to false.
    ref_points: Optional path to a csv file of reference points to be added to the map (could be springs, peaks, towns, etc.).
                CSV must have the following labeled columns: Name, X, Y, Z
    ax:         Optional axis object to plot to. If not specified, new axes will be created.
                
    Returns:
    f,ax:       The figure and axis objects being plotted on.'''
    
    #Get grid info
    xres = geo_model.grid.topography.resolution[0]       #get x resolution
    yres = geo_model.grid.topography.resolution[1]       #get y resolution

    #Load data:
    try:
        a = np.reshape(a, (xres,yres))                        #if 1D, try reshaping to 2D
    except: 
        pass
    try:                                                          #if gslib file, try to load
        a = mapping.importgslib(a, geo_model.grid.topography)     #import geologic map from gslib file
        a = np.flipud(a)                                          #flip map up-down to account for gempy v2 bug
    except:                             
        pass

    #Generate colormap:
    cmap = matplotlib.colors.ListedColormap(geo_model.surfaces.colors.colordict.values())     #set colormap
    norm = matplotlib.colors.Normalize(vmin=1, vmax=len(geo_model.surfaces.colors.colordict)) #set normalization

    #Create axes:
    if ax==None:
        fig = plt.figure(figsize=(10,10))           #create empty figure
        ax = fig.add_subplot(111)                   #add subplot axes
        
    #Plot map:
    ax.imshow(a, extent=geo_model.grid.topography.extent, cmap=cmap, norm=norm)              #plot map with coordinate extent and custom colormap

    #Plot data:
    if plot_data==True:
        pts = geo_model.orientations.df                                                                             #get orientations data points
        ax.scatter(pts.X,pts.Y, s=20, edgecolors='k', linewidths=1, c=pts.id, cmap=cmap, norm=norm, marker='<')    #plot as points with custom colormap
        pts = geo_model.surface_points.df                                                                           #get interface data points
        plt.scatter(pts.X,pts.Y, s=20, edgecolors='k', linewidths=1, c=pts.id, cmap=cmap, norm=norm, marker='o')    #plot as points with custom colormap
    if ref_points:
        pts = pd.read_csv(ref_points)                           #get reference points from file
        ax.scatter(pts.X,pts.Y, s=20, c='k')                   #plot as points 
        for i,pt in pts.iterrows():
            ax.annotate(pt.Name, (pt.X,pt.Y), fontsize=15)     #add labels to reference points

    #Add legend:
    patches = [matplotlib.patches.Patch(color=geo_model.surfaces.colors.colordict[key], label=key) for key in geo_model.surfaces.colors.colordict.keys()]  #loop to create a list of color patches & labels
    plt.legend(handles=patches, bbox_to_anchor=(1.35,1))  #add legend using color patches, at specific position



#############################################################################
def compare_geology(model_geo, true_geo , grid_info, colors):
    
    '''Plots the GemPy-generated geologic map next to the true geologic map for comparison. 
    model_geo:    gslib file of lithologic unit values at the land surface, of dimensions (yres,xres), i.e. (nrow,ncol)
    true_geo:     raster file of the actual geologic map, of same dimensions as model_geo
    gridinfo:     [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz] array of model grid and resolution info from importDEM()
    colors:       dictionary of strings indicating which colors to map to which lithologic unit name (must be in order from youngest to oldest)
    '''
    
    #Get grid info:
    xres = grid_info[2]
    yres = grid_info[6]
    zres = grid_info[10] 
    
    #Import data & reformat:
    #Model:
    model = pn.read_csv(model_geo, skiprows=2, dtype=float)     #read in gslib file to pandas df
    model = model.values                                        #get an array of the values
    model = np.reshape(model,(xres,yres),order='F')             #reshape to xyz grid using Fortran ordering
    model = np.rot90(model)                                     #rotate 90 degrees CCW
    model[model==0] = np.nan                                    #replace zeros with NaNs (array must be float first)
    
    #True:
    true = gdal.Open(true_geo)  #raster must be rectangular tif with zeros as NaN values and formation values corresponding to GemPy fm values 
    true = true.ReadAsArray()    #convert raster to a numpy array (defaults to integers)
    true = true.astype(float)    #convert integer array to float array (to be able to use NaNs)
    true[true==0] = np.nan       #replace zeros with NaNs (have to convert array to float first)

    #Plot:
    f = plt.figure(figsize=(20,20))                 #create figure
    cmap = matplotlib.colors.ListedColormap(colors.values()) #set colormap
    norm = matplotlib.colors.Normalize(vmin=1, vmax=len(colors)) #set normalization
    
    plt.subplot(121)                                #create subplot in position 1 (1 row, 2 col, position 1)
    plt.imshow(model, cmap=cmap)                    #plot model using colormap
    plt.title('model')
    plt.subplot(122)                                #create 2nd subplot
    plt.imshow(true, cmap=cmap)                     #plot true geologic map
    plt.title('data')

    #Add legend:
    patches = [matplotlib.patches.Patch(color=colors[key], label=key) for key in colors.keys()]  #loop to create a list of color patches & labels
    plt.legend(handles=patches, bbox_to_anchor=(1.35,1))  #add legend using color patches, at specific position

    return f


#############################################################################
def plotXsection(startpoints, endpoints, names, grid_info, lith, surflith, vscale=1, 
                 colors = 'gempy', unitnames = None):
    '''Plots an approximate cross-section between the two specified points, using an elevation-cropped array (cells above the land surface should have nan values). Does not work well for N-S or nearly N-S xsections - use gempy built-in for that.
    startpoints: [[x1,y1],[x2,y2],...] float or [[col1,row1],[col1,row2],...] integer array of coordinates of starting points A
    endpoints:   [[x1,y1],[x2,y2],...] float or [[col1,row1],[col1,row2],...] integer array of coordinates of ending points B
    names:       [['A','B'],['C','D'],...] string array of names for starting and ending points
    grid_info:   [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz] model grid info (can get these using importDEM() function if model grid is same as DEM grid)
    lith:        elevation-cropped array of lithologic unit indices of dimensions (nrow,ncol,nlay), i.e. (yres,xres,zres). Can use uncropped array, but will plot above land surface.
    surflith:    array of lithologic unit values at the land surface, of dimensions (yres,xres)
    vscale:      vertical exaggeration factor (y/x, defaults to 1)
    colors:      dictionary of color & unit names to use (in order from youngest to oldest) OR 'gempy' to use default gempy colormap
    '''
    
    #Get model grid and resolution info:
    xmin = grid_info[0]
    xres = grid_info[2]
    dx =   grid_info[3]
    ymin = grid_info[4]
    yres = grid_info[6]
    dy =   grid_info[7]
    zmin = grid_info[8]
    zmax = grid_info[9]
    zres = grid_info[10]
    dz =   grid_info[11]
    
    #Set colormap:
    if colors == 'gempy':                               #default gempy colors
        cmap=gp.plotting.colors.cmap                    #set colormap from GemPy
        norm=gp.plotting.colors.norm                    #set normalization
    else:
        cmap = matplotlib.colors.ListedColormap(colors.values()) #set colormap from list
        norm = matplotlib.colors.Normalize(vmin=1, vmax=len(colors)) #set normalization
    
    #Plot geologic map once for reference:
    f1,ax1 = plt.subplots(1,1,figsize=(10,10))                      #create empty figure
    plt.imshow(surflith, cmap=cmap, norm=norm)                  #plot geology (normalized to gempy color range)
    #Add legend:
    if colors != 'gempy':
        patches = [matplotlib.patches.Patch(color=colors[key], label=key) for key in colors.keys()]  #loop to create a list of color patches & labels
        plt.legend(handles=patches, bbox_to_anchor=(1.35,1))  #add legend using color patches, at specific position
        
    f2,ax2 = plt.subplots(len(startpoints),1,figsize=(15,20))   #create figure and axes objects for subplots (one per xsection)
    
    for i in range(len(startpoints)):   #loop over number of sections
        #Get starting coordinates:
        xA = startpoints[i][0]   #get starting x coordinate
        yA = startpoints[i][1]   #get starting y coordinate
        xB = endpoints[i][0]     #get ending x coordinate
        yB = endpoints[i][1]     #get ending y coordinate

        #Calculate corresponding row,col
        if type(xA) != int:                     #if coordinates are NOT integers (i.e. not row,col numbers), convert them
            colA = (xA - xmin)//dx              #col:x calculate column index  c = (x1-x0)/dx 
            rowA = yres - ((yA - ymin)//dy)     #row:y calculate row index     r = ymax - (y1-y0)/dy
            colB = (xB - xmin)//dx                 
            rowB = yres - ((yB - ymin)//dy) 
        else:                                  #if coordinates are already in row,col format
            colA = xA
            rowA = yA
            colB = xB
            rowB = yB

        #Calculate line equation between points A and B:
        m = (rowB - rowA) / (colB - colA)   #calculate slope     m = (y2-y1)/(x2-x1)
        b = -m*colA + rowA                  #calculate intercept b = m*x1 + y1 (slope is neg here bc y axis is flipped)
        
        #Calculate true distance (not # of cells) between points A and B:
        #distance = ((xB-xA)**2 + (yB-yA)**2)**.5
        #xsizes.append(distance*.001)

        #Get xy indices for cells intersected by the x-sec line, then get z values for those xy points:
        xvals = np.arange(colA,colB)    #generate array of x values between the two points
        xvals = xvals.astype(int)       #convert to integer
        yvals = m*xvals + b             #calculate corresponding y values  y = mx + b 
        yvals = yvals.astype(int)       #convert to integers to be able to slice

        xsec = lith[yvals,xvals,:].T    #select x-sec to plot and transpose to make it plot horizontally
        
        #Plotting:
        #Add xsection lines to geologic map:
        plt.figure(f1.number)                           #make the map the active figure
        plt.plot([colA,colB],[rowA,rowB],'k')           #plot x-sec location line
        plt.annotate(names[i][0],xy=(colA,rowA),xytext=(colA-4,rowA+4)) #annotate start point
        plt.annotate(names[i][1],xy=(colB,rowB),xytext=(colB+1,rowB-1)) #annotate start point
        plt.ylim(bottom=yres, top=0) 
        plt.xlim(left=0, right=xres)

        #Plot cross-sections in a new figure:
        #Set and get correct subplot axes:
        if len(startpoints) == 1:               #check if there are more than 1 subplots (for indexing purposes)
            plt.sca(ax2)                        #make current subplot axes active (automatically makes fig active too)
            cax = plt.gca()                     #get current axes object
        else:          
            plt.sca(ax2[i])                     
            cax = plt.gca()
        cax.imshow(xsec, origin="lower", cmap=cmap, norm=norm)   #plot (with down=lower z indices)
        #cax.imshow(xsec, origin="lower", cmap=cmap)   #plot (with down=lower z indices)
        cax.set_aspect(vscale*dz/dx)                             #apply vertical exaggeration
        cax.set_ylim(bottom=0, top=zres)                         #set y limit to zres
        cax.set_title(names[i][0]+names[i][1])
        cax.set_anchor('W')                                      #align left (West)

        #Set ticks to accurately reflect elevation (masl):
        locs = cax.get_yticks()                               #get tick locations
        nlabels = len(cax.get_yticklabels())                  #get number of initial ticks 
        labels = np.linspace(zmin, zmax, nlabels)             #generate list of tick labels
        ticks = cax.set(yticks=locs,yticklabels=labels)       #set tick locations and labels
        
    return f1,ax1,f2,ax2


    
##############################################################################################
def xy2rowcol(X, Y, grid_info, flip=False):

    '''Converts between X,Y coordinates and row,col coordinates.
    Inputs:
    X: list of x coordinates to convert to columns
    Y: list of y coordinates to convert to rows
    grid_info:  [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz] model grid info (can get these using importDEM() function if model grid is same as DEM grid)
    flip: if False (default), converts X,Y to row,col. If True, converts row,col to X,Y'''
    
    #Get grid info:
    xmin = grid_info[0]
    dx =   grid_info[3]
    ymin = grid_info[4]
    yres = grid_info[6]
    dy =   grid_info[7]
    
    #Convert X,Y to row, col
    if flip == False:
        cols = []
        for x in X:
            col = (x - xmin)//dx              #col:x calculate column index  c = (x1-x0)/dx 
            cols.append(col)

        rows = []
        for y in Y:
            row = yres - ((y - ymin)//dy)     #row:y calculate row index     r = ymax - (y1-y0)/dy
            rows.append(row)
            
        return rows, cols