import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from numpy.random import default_rng



def plot_trajectory(trajectory):
    """
    Plots a trajectory.
    
    The trajectory is drawn as a line and
    the current position is drawn as a point.
        
    Args:
        trajectory (array): list of coordinate pairs
    """
    
    # separate x and y coordinates by slicing the array
    xcoordinates = trajectory[:,0]
    ycoordinates = trajectory[:,1]
    
    plt.clf()
    ax = plt.axes()
    ax.set_xlim([-200, 200])
    ax.set_ylim([-200, 200])
    ax.set_aspect('equal')
    
    # plot the path as a line
    plt.plot( xcoordinates, ycoordinates, 'k-' )

    # plot the current position as a point
    plt.plot( xcoordinates[-1], ycoordinates[-1], 'ro' )
    
    plt.show()


def draw(frame):
    """
    Draws one frame for animation.
    
    A trajectory of the shifting coordinates are saved
    in the global variable trajectory as a list of coordinates. 
    This function draws the system as defined in trajectory[frame].
    
    Args:
        frame: index of the frame to be drawn.
    """

    plt.clf()
    ax = plt.axes()
    ax.set_xlim([-200, 200])
    ax.set_ylim([-200, 200])
    ax.set_aspect('equal')
    
    # plot the path as a line
    plt.plot( trajectory[:frame+1,0], trajectory[:frame+1,1], 'k-' )

    # plot the current position as a point
    plt.plot( trajectory[frame,0], trajectory[frame,1], 'ro' )


def animate():
    """
    Animates the trajectory.
    
    A trajectory of the shifting coordinates are saved
    in the global variable trajectory as a list of coordinates. 
    This function animates the system up to the specified
    frame.
    """

    last_frame = len(trajectory)
    fig = plt.figure()
    motion = ani.FuncAnimation(fig, draw, last_frame)
    plt.show()


def write_file(trajectory, filename):
    """
    Writes the trajectory to file.
    
    Each line in the file will contain a pair x and y coordinates
    separated by whitespace:
    
    .. code-block::
    
        x0   y0
        x1   y1
        x2   y2
        ...
    
    Args:
        trajectory (array): list of coordinate pairs
        filename (str): name of file to write
    """
    f = open(filename, 'w')
    
    for point in trajectory:
        x = point[0]
        y = point[1]
        line = str(x)+"   "+str(y)+"\n"
        f.write(line)

    f.close()
    
    
def read_file(filename):
    """
    Reads a trajectory from file.
    
    The file must obey the format specified in :meth:`write_file`.
    
    Args:
        filename (str): name of the file to read.
        
    Returns:
        array: trajectory as an array of coordinate pairs
    """
    
    traj = np.genfromtxt(filename)
    
    return traj


def move(x,y,step,angle):
    """
    Calculates new coordinates by taking a step from a given starting point. 
    
    The function starts from position :math:`(x, y)` 
    and moves the distance :math:`L`
    in the direction defined by the angle :math:`\\theta`,
    where :math:`\\theta = 0` means moving in positive x direction
    and :math:`\\theta = \\pi/2` means moving in positive y direction.
    
    The function returns the coordinates of the final position.
    
    Args:
        x (float): initial x coordinate
        y (float): initial y coordinate
        step (float): step length :math:`L`
        angle (float): step direction :math:`\\theta` in radians
        
    Returns:
        float, float: final coordinates (x, y)
    """

    nx = x + step*np.cos(angle)
    ny = y + step*np.sin(angle)
    
    return nx, ny


def move_randomly(x,y,std):
    """
    Calculates new coordinates by taking a step from a given starting point. 
    
    The function starts from position :math:`(x, y)` 
    and takes a random step so that both coordinates are
    changed by a normally distributed shift.
    
    The function returns the coordinates of the final position.
    
    Args:
        x (float): initial x coordinate
        y (float): initial y coordinate
        std (float): standard deviation for the coordinate shifts.
        
    Returns:
        float, float: final coordinates (x, y)
    """

    nx = x + std * random.standard_normal()
    ny = y + std * random.standard_normal()
    
    return nx, ny


def main(x_start, y_start, angle_start, simulation_length = 500):
    """
    Moves a point step by step and draws the trajectory.
    
    The point is moved using :meth:`move`.
    The path is visualized using :meth:`plot_trajectory` and :meth:`animate`.
    
    Args:
        x_start (float): starting point x coordinate
        y_start (float): starting point y coordinate
        angle_start (float): angle defining the starting direction
        simulation_length (int): total number of steps to take
    """
    
    x = x_start
    y = y_start
    
    # Declare trajectory as a global variable for storing the trajectory.
    # Global variables are usually not good style, but it makes animation simpler.
    # We will later see how to do this without global variables.
    global trajectory
    trajectory = [[x,y]]

    delta_angle = np.pi / 50

    for i in range(simulation_length):
        x, y = move_randomly(x,y, 4)
        trajectory.append([x,y])

    # finally, turn trajectory into an array
    trajectory = np.array(trajectory)
    
    plot_trajectory(trajectory)

    animate()


if __name__ == "__main__":
    random = default_rng()
    main(0, 0, 0)
    