#Generating plots module

def prettyplot(x,y,title,xlabel,ylabel,color='steelblue'):
    
    '''Generate line plot with specified title, axis labels, and color.
    Inputs:
    x: list of x values
    y: list of y values (must be same size as x)
    title: string for plot title
    xlabel: string for x axis label
    ylabel: string for y axis label
    color: string for color (see options here: https://matplotlib.org/users/colors.html)'''
    
    from matplotlib import pyplot as plt
    plt.plot(x,y,color)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
