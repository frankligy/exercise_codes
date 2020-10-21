'''
This script is to play how to use python to generate scientific figures
'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Avenir'


if __name__ == '__main__':
    fig = plt.figure()   # create a new figure, default figsize is (6.4,4.8), unit is inch
    axes1 = plt.axes([0.1,0.1,0.65,0.65])   # subplot with [left,bottom,width,height]
    axes2 = plt.axes([0.1,0.78,0.65,0.2])
    fig,axes = plt.subplots(3,3)  # nrows=3, ncols=3, each subplot will be accessed by axes[2,2]
    axes1 = plt.subplot(221)
    axes2 = plt.subplot(222)


    ax = fig.add_axes([0,0,1,1])   # add an ax to the figure
    plt.show()    # looks for all current active figure objects, opens interactive windows and display for figures

    import matplotlib.font_manager as fm
    font_names = [f.name for f in fm.fontManager.ttflist]
    print(font_names)  # check all available font in my computer

    colors = cm.get_cmap('tab10',2)


    # remove spines
    axes1.spines['right'].set_visible(False)
    axes1.spines['top'].set_visible(False)

    # ticks
    axes1.xaxis.set_tick_params(which='major',size=10,width=2,direction='in',top='on')
    '''
    which: major, minor, both
    size: length of ticks
    width: line width of ticks
    direction: in, out, inout, which direction ticks will face
    top/right: 'on', will ticks show on the secondary axes (top and right spines)  
    '''

    ax.plot(data,color=colors[0])

    ax.set_xlim(10,100)

    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(100))  # set the major ticks every multiple of 100, FixedLocator for mannually add ticks


    ax.set_xlabel(r'$\mathregular{\lambda}$',labelpad=10)   # labelpad: extra padding between the tick labels and the axis label

    ax.legend(bbox_to_anchor=(1,1),loc=1,frameon=False,fontsize=16)
    # a anchor point of bound box, loc determine which anchor point, like lower left or upper right,
    # bbox_to_anchor is the coordinate of this anchor point, remember the axes object is 0 to 1, not the x,y label itself
    # frameon determine whether there will be a frame around legend























