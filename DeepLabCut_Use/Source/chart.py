import numpy as np
from matplotlib import pyplot


def display_chart(dist,time,name):

    pyplot.plot(time, dist)
    pyplot.xlabel('Time ')
    pyplot.ylabel('Distance between each frame')

    pyplot.savefig(name + ".png")
    pyplot.close('all')


def display_double_chart(dist, movingness,temps, name):


    # Create some mock data
    fig, ax1 = pyplot.subplots()
    ax1.title.set_text(name)

    color = 'tab:grey'
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Is moving', color=color)
    ax1.plot(temps, movingness, color=color)
    ax1.tick_params(axis='y', labelcolor=color)


    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('Distance between each frame', color=color)  # we already handled the x-label with ax1
    ax2.plot(temps, dist, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    pyplot.savefig(str(name) + ".png")
    pyplot.close('all')


def display_threshold_chart(dist,time,threshold, name):

    # line 1 points
    x1 = time
    y1 = dist[:,0]
    # plotting the line 1 points
    pyplot.plot(x1, y1, label="dist")
    # line 2 points
    x2 = time
    y2 = np.repeat(threshold,np.size(dist[:,0],0))
    # plotting the line 2 points
    pyplot.plot(x2, y2, label="threshold")
    pyplot.xlabel('time (s)')
    # Set the y axis label of the current axis.
    pyplot.ylabel('distance between frame (pixels)')
    # Set a title of the current axes.
    pyplot.title('Two or more lines on same plot with suitable legends ')
    # show a legend on the plot
    pyplot.legend()

    pyplot.savefig(name + ".png")
    pyplot.close('all')