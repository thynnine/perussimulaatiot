# /usr/bin/env python
from numpy.random import default_rng
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
from math import *
from scipy.optimize import curve_fit


def print_progress(step, total):
    """
    Prints a progress bar.
    
    Args:
        step (int): progress counter
        total (int): counter at completion
    """

    message = "simulation progress ["
    total_bar_length = 60
    percentage = int(step / total * 100)
    bar_fill = int(step / total * total_bar_length)
    for i in range(total_bar_length):
        if i < bar_fill:
            message += "|"
        else:
            message += " "
    
    message += "] "+str(percentage)+" %"
    if step < total:
        print(message, end="\r")     
    else:
        print(message) 
        


def generate_empty_lattice(lattice_size=100):
    """
    Creates an empty lattice as an :math:`N \\times N` array of zeros.
    
    Args:
        lattice size (int): system size :math:`N`
    
    Returns:
        array: the lattice
    """

    grid = np.zeros( [ lattice_size, lattice_size ] )
    return grid


def generate_lattice(particles, lattice_size, label=1, additive=True):
    """
    Creates a lattice with particles as an :math:`N \\times N` array.
    
    The lattice sites with no particles are given the value of zero.
    The sites with 1 or more particle are given a non-zero value.
    
    Args:
        particles (list): list of :class:`Walker` objects
        lattice size (int): system size :math:`N`
        label (int): value given to sites with particles
        additive (bool): If True, the value given to sites is proportional
           to the number of particles on that site. If False, all sites
           with however many particles are given the same value.
    Returns:
        array: the lattice
    """
       
    grid = generate_empty_lattice(lattice_size)
    
    for part in particles:
        try:
            if additive:
                grid[part.x, part.y] += label
            else:
                grid[part.x, part.y] = label
        except:
            pass
            
    return grid


def show_system(particles, lattice_size, scale=5):
    """
    Draws the system as a :math:`N \\times N` pixel image.
    
    Args:
        particles (list): list of :class:`Walker` objects
        lattice size (int): system size :math:`N`
        scale (int): values above this will be shown as black
    """
    plt.clf()
    ax = plt.axes()
    ax.set_aspect('equal') 
    plt.pcolormesh(generate_lattice(particles, lattice_size), cmap='Greys', vmax = scale)
    plt.show()
    
    
def draw(frame, history, scale=5):
    """
    Draws the system for animation.
    
    Args:
        frame (int): index of the frame to draw
        history (list): list of systems at different times
        scale (int): values above this will be shown as black
    """
    plt.clf()
    ax = plt.axes()
    ax.set_aspect('equal') 
    plt.pcolormesh( history[frame], cmap='Greys', vmax=scale)


def animate(history):
    """
    Animate the simulation.
    
    Args:
        history (list): list of systems at different times
        lattice size (int): system size :math:`N`
    """

    nframes = len(history)
    print("animating "+str(nframes)+" frames")
    fig = plt.figure()
    motion = ani.FuncAnimation(fig, draw, nframes, fargs=( history, ) )
    plt.show()


class Walker:
    """
    A particle moving on a lattice.
    
    Since the particle can only move on a lattice,
    its coordinates are always integers.
    
    Args:
        x (int): initial x coordinate
        y (int): initial y coordinate
    """

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.initial_x = x
        self.initial_y = y
        
        
    def move_randomly(self):
        """
        Makes the particle move.
        
        The particle takes a step of length 1 in one of the
        cardinal directions: up, down, left or right.
        
        There is a 20 % chance for each direction.
        In addition, there is a 20 % chance the particle
        stays in place.        
        """

        randomizer = random.random()        
        
        if randomizer < 0.20:
            self.x += 1
        elif randomizer < 0.40:
            self.y += 1
        elif randomizer < 0.60:
            self.x -= 1
        elif randomizer < 0.80:
            self.y -= 1


    def get_distance_sq(self):
        """
        Calculates the squared distance between the initial
        and current coordinates of this particle.
        
        Returns:
            float: squared distance :math:`r^2 = (x-x_0)^2 + (y-y_0)^2`
        """
        
        dx = self.x - self.initial_x
        dy = self.y - self.initial_y        
        return ( dx*dx + dy*dy )




def move(particles):
    """
    Makes all particles move.
    
    All particles take a random step using :meth:`Walker.move_randomly`.
    
    Args:
        particles (list): list of :class:`Walker` objects
    """
    for part in particles:
        part.move_randomly()


def calculate_distance_statistics(particles):
    """
    Calculates the average and error of mean of all particles squared displacement. 
    
    The squared displacement :math:`r_i^2` is calculated with 
    :meth:`Walker.get_distance_sq` for every particle :math:`i`. 
    The function then calculates the mean 
    
    .. math::
        \\langle r^2 \\rangle = \\frac{1}{N} \\sum_i r_i^2
    
    and error of mean

    .. math::
        \\Delta (r^2) = \\frac{s_{r^2}}{\\sqrt{N}},
        
    where :math:`s_{r^2}` is the sample standard deviation
    
    .. math::
        s_{r^2} = \\sqrt{ \\frac{1}{N-1} \\sum_i ( r_i^2 - \\langle r^2 \\rangle )^2 }.
    
    Args:
        particles (list): list of :class:`Walker` objects
    
    Returns:
        float, float: average :math:`\\langle r^2 \\rangle`, error :math:`\\Delta (r^2)`
    """
    n_particles = len(particles)
    sq_average = 0.0
    sq_deviation = 0.0
    sq_distances = np.zeros(n_particles)
    for i in range(n_particles):
        sq_distances[i] = particles[i].get_distance_sq()
    
    for dist in sq_distances:
        sq_average += dist

    sq_average /= n_particles

    for dist in sq_distances:
        sq_deviation += (dist - sq_average)*(dist - sq_average)

    sq_deviation = sqrt( sq_deviation / (n_particles - 1) )
    sq_error = sq_deviation / sqrt(n_particles)
    
    return sq_average, sq_error
    

def calculate_density_histogram(particles, range_steps = 10, max_range = 50):
    """
    Calculates a histogram of particle density as a function of displacement.
    
    The, function only takes into account the particles for which
    the displacement :math:`r_i` is not greater than a given maximum :math:`R`.
    The range of possible values is split in
    :math:`k` separate subranges or *bins*, :math:`[0,R/k), [R/k, 2R/k)` etc.
    The midpoints of each bin, :math:`\\frac{R}{2k}, \\frac{3R}{2k}, \\frac{5R}{2k}` etc.
    are also saved.
    
    Having defined these bins, the function calculates the displacement :math:`r_i`
    for each particle using :meth:`Walker.get_distance_sq` and checks
    which bin the displacement belongs to.
    The function then counts how many particles are placed in each bin, :math:`n(r)`.
    
    Finally, the function calculates the surface density of particles at each
    displacement range. That is, the number of particles in a given region
    is divided by the area of that region,
    
    .. math::
    
        \\sigma(r) = \\frac{n(r)}{A(r)}.
        
    As this is a 2D simulation, the particles whose displacement hit the range
    :math:`[r, r+\\Delta r]` are inside a circular edge with inner radius
    :math:`r` and thickness :math:`\\Delta r`.
    The area of this edge is :math:`A(r) = \\pi (r + \\frac{1}{2} \\Delta r) \\Delta r`.
    
    Args:
        particles (list): list of :class:`Walker` objects
        range_steps (int): the number of bins, :math:`k`
        max_range (float): maximum diaplacement, :math:`R`
        
    Returns:
        array, array, array: average :math:`r`, particle count :math:`n(r)`, particle density :math:`\\sigma(r)` in each bin
    """

    # array of counters for how many particles hit each bin
    range_bins = np.zeros(range_steps)
    
    # array of density scaled particle counters
    density_bins = np.zeros(range_steps)
    
    # the middle points for each subrange
    avr_ranges = (np.array( range(range_steps) ) + 0.5 ) * max_range / (range_steps)
    dr = avr_ranges[1]-avr_ranges[0]
    
    for part in particles:
        distance = sqrt( part.get_distance_sq() )
        #range_index = min( floor( distance / max_range * range_steps  ), range_steps - 1 )
        range_index = floor( distance / max_range * range_steps  )
        if range_index < range_steps:
            range_bins[range_index] += 1
    
    # range_bins stores the number of particles at a distance [r, r+dr] from the starting point
    # to get the density of particles at this slice, divide by the area dA = 2 pi r dr
    return avr_ranges, range_bins, range_bins/(2*pi*avr_ranges*dr)
    
  
     
def linear(x, a):
    """Calculates the linear function :math:`y=ax` and returns the result.
    
    Args:
        x (float): the variable
        a (float): the slope
    
    Returns: 
        float: the result 
    """
    return a*x 


def linear_fit(xdata, ydata, name):
    """Fit a linear function :math:`y = ax` to the given xy data.
    Also return the optimal fit values obtained for slope a.
    
    To identify this information, you may also name the operation. 

    Args:
        xdata (array): x values
        ydata (array): y values
        name (str): identifier for printing the result 
    
    Returns: 
        float: slope
    """
    params, covariance = curve_fit(linear, 
                            np.array(xdata), 
                            np.array(ydata))
    print( "" )
    print( "slope for "+name+" fit: "+str(params[0]) )
    print( "" )

    return params[0]
    





def main(n_particles, n_steps, n_plots):
    """
    The main program.
    
    Creates particles at the center of the system and then allows them to diffuse
    by performing a random walk on a lattice.
    
    Animates the motion and if the number of particles is large enough,
    also calculates statistics for the particle distribution.
    
    Args:
        n_particles (int): number of particles
        n_steps (int): number of simulation steps
        n_plots (int): number of times data is recorded
    """

    lattice_size = int( sqrt(20*n_steps) )
    
    # if there are at least this many particles,
    # calculate statistics
    limit = 1000
    
    particles = []
    for i in range(n_particles):
        x = lattice_size//2
        y = lattice_size//2
        particles.append( Walker(x,y) )
 
    averages = []
    errors = []
    steps = []
    bins = []
    histograms = []
    ranges = []
    history = []

    for i in range(n_steps+1):
        
        # move particles
        move(particles)
        print_progress( i, n_steps )
        
        # record the system every now and then
        if i%(n_steps//n_plots) == 0:            
            history.append( generate_lattice(particles, lattice_size) )
            
            # calculate stats if we have a lot of particles
            if n_particles > limit:
                ranges, bins, density_histogram = calculate_density_histogram(particles, int(sqrt(n_particles)/10), max(50, lattice_size//2) )    
                histograms.append( density_histogram ) 
        
                avr, err = calculate_distance_statistics(particles)
                averages.append(avr)
                errors.append(2.0*err)
                steps.append(i)
        

    # show simulation
    animate( history )

    # show statistics if we have a lot of particles
    if n_particles > limit:
    
        steps = np.array(steps)
        averages = np.array(averages)
        errors = np.array(errors)

        # do a linear fit to the data
        slope = linear_fit(steps, averages, "diffusion")
        plt.errorbar(steps, averages, yerr=errors, fmt='-o')
        plt.plot(steps, linear(steps, slope), '-')
        plt.xlabel("t")
        plt.ylabel("$<\\Delta r^2>$")
        plt.show()

        for hist in histograms[1:]:
            plt.plot(ranges, hist, '-')
            
        plt.xlabel("r")
        plt.ylabel("$\\sigma (r)$")
        plt.show()
        
    



if __name__ == "__main__":

    random = default_rng()
    n_particles = 10000
    n_steps = 2000
    n_plots = 20
    
    main(n_particles, n_steps, n_plots)