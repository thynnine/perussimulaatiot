import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
from numpy.random import default_rng
import copy



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


def generate_system(lattice_size, cluster=False):
    """
    Creates the system.
    
    The system is represented by a 2D :math:`N \\times N` lattice 
    where each lattice site
    may be either occupied (1) or unoccupied (0).
    
    The system may be a surface, in which case one of its sides
    has a row of empty sites, or a cluster, in which case all
    sides have empty sites.
    
    Args:
        lattice_size (int): lattice size :math:`N`
        cluster (bool): if True, the system is a cluster
        
    Returns:
        array: the system as an integer array
    """

    grid = np.array( [ [1]*lattice_size ]*lattice_size )   

    if cluster: # create a cluster by clearing sites on all sides
        grid[:,0] = 0
        grid[:,-1] = 0
        grid[0,:] = 0
        grid[-1,:] = 0
                
    else: # create a surface by clearing sites on one side
        grid[-1,:] = 0

    return grid
          

def show_system(grid):
    """
    Draws the system.
    
    Args:
        grid (array): the system as an integer array
    """
    plt.clf()
    ax = plt.axes()
    ax.set_aspect('equal') 
    plt.pcolormesh(grid, cmap='Greys')
    plt.show()
    
    
def draw(frame, history):
    """
    Draws the system for animation.
    
    Args:
        frame (int): index of the frame to draw
        history (list): list of systems at different times
    """
    plt.clf()
    ax = plt.axes()
    ax.set_aspect('equal') 
    plt.pcolormesh(history[frame], cmap='Greys')


def animate(history):
    """
    Animate the simulation.
    
    Args:
        history (list): list of systems at different times
    """

    nframes = len(history)
    print("animating "+str(nframes)+" frames")
    fig = plt.figure()
    motion = ani.FuncAnimation(fig, draw, nframes, fargs=( history, ) )
    plt.show()


def count_neighbors(i,j,grid):
    """
    Count the number of occupied neighbors for the lattice site (i,j).
    
    The system is represented as a square lattice so the first neighbors of
    a given site are the 4 sites directly next to it
    left, right, up and down.
    
    The second neighbors of a site are the 4 sites diagonally next to it:
    up-left, up-right, down-left, down-right.
    
    This function counts how many of those sites are occupied (have the value 1).
    
    The edges of the system wrap around.
    It means that for a site that is, say, on the left edge of the system,
    the left-side neighbor is a site on the right edge of the system.
    
    Args:
        i (int): row index
        j (int): column index
        grid (array): the system as an integer array
        
    Returns:
        int, int: number of first neighbors, number of second neighbors
    """
            
    lattice_size = len(grid)
    
    left = (j-1)%lattice_size # i-1 (mod n). This is n-1, if i = 0.
    right = (j+1)%lattice_size # i+1 (mod n). This is 0, if i = n-1.
    up = (i+1)%lattice_size
    
    # The bottom layer is supposed to be bulk matter, so we count the atoms there 
    # to be their own neighbors. This makes them harder to remove so that surfaces 
    # are really etched from the top.
    down = max( (i-1), 0 )  
    
    firsts = grid[i,left] + grid[i,right] + grid[up,j] + grid[down,j]
    seconds = grid[up,left] + grid[up,right] + grid[down,left] + grid[down,right] 

    return firsts, seconds


def get_removal_rate( i, j, grid, unremovables = [] ):       
    """
    Calculates the removal rate for a site.
    
    For normal sites, the removal rate is evaluated based on
    how many 1st and 2nd neighbors the site has.
    These are calculated with :meth:`count_neighbors`.
    
    Certain sites may be listed as protected sites that may never be
    removed. Such sites mimic the use of protective masks.
    
    .. note ::
        Edit to change the physics!
    
    Args:
        i (int): row index
        j (int): column index
        grid (array): the system as an integer array
        unremovables (list): indices of sites that are never removed
    
    Return:
        float: the removal rate
    """
    
    for protected in unremovables:
        if i == protected[0] and j == protected[1]:
            return 0
            
    # If this is a normal site, we evaluate its rate
    # of being removed according to how many neighbors it has.
    n_1st, n_2nd = count_neighbors(i,j,grid)
    
    # There are a few possible sets of removal rates below.

    # kinks are much more soluble than planes,
    # horizontal and diagonal planes are equally soluble
    rates_A = np.array( 
        [ [10**20, 10**15, 10**10, 1],
          [10**15, 10**10, 10**5,  0],
          [10**10, 10**5,  100,    0],
          [1,      0,      0,      0] ] )
          
    # kinks are much more soluble than planes,
    # horizontal planes are more soluble than diagonal
    rates_B = np.array( 
        [ [10**20, 10**15, 10**10, 1],
          [10**15, 10**10, 1000,   0],
          [10**10, 10**5,  100,    0],
          [1,      0,      0,      0] ] )
          
    # kinks are much more soluble than planes,
    # horizontal planes are less soluble than diagonal
    rates_C = np.array( 
        [ [10**20, 10**15, 10**10, 1],
          [10**15, 10**10, 10**5,  0],
          [10**10, 1000,   100,    0],
          [1,      0,      0,      0] ] )

    # kinks and planes are equally soluble
    rates_D = np.array( 
        [ [10**20, 10**15, 10**10, 1],
          [10**15, 10**5,  10**5,  0],
          [10**10, 10**5,  100,    0],
          [1,      0,      0,      0] ] )

    # Choose which set of rates to use:
    rates = rates_A
    
    return rates[n_1st-1, n_2nd-1]
    


def count_probabilities(grid):
    """
    Calculates, for each site, the probability that this particular site is
    the next to be removed.
    
    The function loops over all occupied sites :math:`(i,j)` and calculates 
    the removal rate :math:`r_{i,j}` for each site.
    
    The total removal rate is defined as
    
    .. math ::
    
        R = \\sum_{i,j} r_{i,j}
        
    and the probability that site :math:`(i,j)` is the next to be removed, is

    .. math ::
    
        p(i,j) = \\frac{1}{R} r_{i,j}.
        
    Finally, the function lists the coordinates of all occupied sites in a
    single site vector :math:`S` and constructs the corresponding vectors 
    containing the site specific probabilities :math:`p` and accumulated 
    probabilities :math:`P`. 
    
    The elements of vector :math:`S` are pairs of integers, so you may have, 
    e.g., :math:`S_0 = (0,0), S_1 = (1,0), S_2 = (2,0)` etc.
    
    The elements of vector :math:`p` are probabilities, so you may have,
    e.g., :math:`p_0 = p(0,0), p_1 = p(1,0), p_2 = p(2,0)` etc.
    
    The elements of vector :math:`P` are sums of probabilities,
    
    .. math::
    
        P_n = \\sum_{k = 0}^n p_k.
        
    So :math:`P_0 = p_0, P_1 = p_0 + p_1, P_2 = p_0 + p_1 + p_2` etc.
    Especially the probabilities must sum up to 1, so the last
    element of vector :math:`P` is always 1.
     
    Args:
        grid (array): the system as an integer array
        
    Returns:
        array, array, array: :math:`p, P, S`
    """
    
    accumulated_probabilities = []
    probabilities = []
    coordinates = []
    total_rate = 0.0

    # Loop over all sites except the bottom row.
    # We assume the bottom row is connected to bulk
    # and will not be removed.
    for i in range(lattice_size-1):
        for j in range(lattice_size):
        
            # only consider occupied sites
            if grid[i,j] > 0:
            
                # calculate the rate and sum of rates
                removal_rate = get_removal_rate( i, j, grid )
                total_rate += removal_rate
                
                # only consider sites that can be removed
                if removal_rate > 0:
                
                    # save the rate, sum of rates and current coordinates
                    probabilities.append( removal_rate )
                    accumulated_probabilities.append( total_rate )
                    coordinates.append( [i,j] )
                
    # Divide rates by the final total rate to get probabilities.    
    accumulated_probabilities = np.array( accumulated_probabilities ) / total_rate
    probabilities = np.array( probabilities ) / total_rate
    coordinates = np.array( coordinates )
    
    return probabilities, accumulated_probabilities, coordinates



def pick_random_event(grid):
    """
    Randomly choose the next event (removal of a site) to happen.
    
    The random event is chosen as follows:
    
    * The probability :math:`p` for choosing each site is calculated 
      with :meth:`count_probabilities`.
    * Also the accumulated probabilities :math:`P` are calculated.
    * A random number :math:`R \\in [0,1]` is drawn.
    * The event :math:`n` for which :math:`P_{n-1} < R \le P_{n}` is chosen.
      (For :math:`n = 0`, it is enough that :math:`R \le P_{0}`.)
    * The function returns the coordinates of the site for this event, :math:`S_n`.
    
    Args:
        grid (array): the system as an integer array
        
    Returns:
        array: coordinates :math:`(i,j)` of the chosen site
    """
    
    randomizer = random.random()
    
    prob, acc_prob, coords = count_probabilities(grid)
    
    for i in range(len(acc_prob)):
        if acc_prob[i] >= randomizer:
            return coords[i]       


def randomly_remove_atom(grid):
    """
    Randomly choose one site for removal and make it unoccupied.
    
    The site is chosen using :meth:`pick_random_event`.
    The chosen site is marked as unoccupied by
    setting the grid value to 0.
    
    Args:
        grid (array): the system as an integer array
    """

    atom_coordinates = pick_random_event(grid)
    grid[atom_coordinates[0], atom_coordinates[1]] = 0




def main(lattice_size, n_steps, n_plots, show_as_you_go = False, cluster = False):
    """
    Main program.
    
    Creates a system as a :math:`N \\times N` square lattice, 
    then removes atoms from it.
    
    Args:
        lattice_size (int): lattice size :math:`N`
        n_steps (int): the total number of atoms to remove
        n_plots (int): the number of images to create
        show_as_you_go (bool): If True, will show images of the
            system during simulation. If False, an animation is
            created at the end. Only set to True for short tests.
        cluster (bool): if True, the system is a cluster
    """

    grid = generate_system(lattice_size, cluster)
    history = [grid]

    for i in range(n_steps):
        randomly_remove_atom(grid)
        print_progress(i+1, n_steps)
    
        # take a look at the system every now and then
        if i%(n_steps/n_plots) == 0:
            
            # show progress as simulation proceeds
            if show_as_you_go:
                show_system(grid)
            else:
                history.append(copy.copy(grid))
                
    if not show_as_you_go:
        animate(history)


if __name__ == "__main__":

    lattice_size = 40 # system size
    n_steps = 1000    # simulation length
    n_plots = 200     # number of images to produce
    show_as_you_go = False # draw the system during simulation or at the end?
  
    random = default_rng()
    main(lattice_size, n_steps, n_plots, show_as_you_go, False)

