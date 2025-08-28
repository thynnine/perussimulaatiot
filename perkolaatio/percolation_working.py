# /usr/bin/env python
import random as rng
import matplotlib.pyplot as plt
import numpy as np


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

def generate_network(lattice_size, node_probability):
    """
    Creates a random :math:`N \\times N` square lattice.
    
    Each lattice site is given a random value of 0 or 1.
    A value of 0 represents an unoccupied site.
    A value of 1 represents an occupied site.
    
    Args:
        lattice_size (int): system length :math:`N`
        node_probability (float): probability for each node to get value 1
        
    Returns:
        array: the network
    """

    grid = np.array( [ [0]*lattice_size ]*lattice_size )   

    for i in range(lattice_size):
        for j in range(lattice_size):
            if rng.random() < node_probability:
                grid[i,j] = 1

    return grid
          

def show_network(grid):
    """
    Plots the system.
    
    Args:
        grid (array): the system as an integer lattice
    """
    plt.imshow(grid, cmap='Greys', interpolation='nearest')
    plt.show()


def find_top_cluster(grid):
    """
    Find and mark the lattices sites with an unbroken connection to the
    top of the system.
    
    Two occupied sites (value 1) are connected, if they are directly adjacent,
    i.e., connected by a common side left-right or up-down.
    
    Two sites are connected by an unbroken connection, if there is
    and unbroken chain of connected sites between the two.

    This algorithm starts from the top of the system and systematically
    searches for all the sites with an unbroken connection to the top.
        
    All the connected sites will be given the value of 2
    for later recognition.
    
    Args:
        grid (array): the system as an integer lattice
    """
    lattice_size = len(grid[0,:])

    # Find the occupied sites in top row.
    #
    # Use the temporary value of 10 to mark sites that
    # are connected to the top of the lattice.
    # Since this is the top row, all sites are.
    for j in range(lattice_size):
        if grid[0,j] == 1:
            grid[0,j] = 10


    # Repeatedly loop over the system searching
    # for occupied sites (value 1) next to sites
    # connected to the top of the system (value 10).
    # If any new sites are found, the variable 
    # 'changed' is made True to force another check.
    changed = True
    while changed:

        # Starting a new loop, we have not yet found any
        # new connected sites, so changed is False.
        # If any new sites are found, changed becomes False
        # and the while loop runs again.
        changed = False

        # Loop over all lattice sites.
        for i in range(lattice_size):
            for j in range(lattice_size):

                # Value 10 means the site is connected to
                # the top of the lattice and its neighbors
                # have not yet been checked. Change the value
                # of the site to 2 and check if any of the neighbors
                # is occupied.
                if grid[i,j] == 10:

                    # The value of 2 means this site is connected
                    # to the top of the lattice AND its neighbors
                    # have been checked for occupation.
                    grid[i,j] = 2                     

                    # Check the neighbors of the site.
                    # If any of them is occupied (value 1)
                    # but not yet marked as being connected
                    # (value 2 or 10), mark them as being connected
                    # (value 10).
                    
                    # Check the site (i, j+1).
                    # Note that the maximum value of j is
                    # lattice_size-1, so this site exists only
                    # if j is smaller than this. If not,
                    # we are at the edge of the system. 
                    if j<lattice_size-1:
                        if grid[i,j+1] == 1:
                            grid[i,j+1] = 10
                            changed = True
                            
                    # Check (i, j-1)
                    if j>0:
                        if grid[i,j-1] == 1:
                            grid[i,j-1] = 10
                            changed = True

                    # Check (i+1, j)
                    if i<lattice_size-1:
                        if grid[i+1,j] == 1:
                            grid[i+1,j] = 10
                            changed = True

                    # Check (i-1, j)
                    if i>0:
                        if grid[i-1,j] == 1:
                            grid[i-1,j] = 10
                            changed = True
                    
                

  
def is_vertically_connected(grid):
    """
    Checks if there is an unbroken connection between
    the top and bottom sides of the system.
    
    The function first calls the function :meth:`find_top_cluster`. 
    to locate the percolation cluster. 
    Then, this routine checks if this cluster reaches 
    the bottom of the lattice.
    
    Args:
        grid (array): the system as an integer lattice
        
    Returns:
        bool: True if a connection exists
    """
    lattice_size = len(grid[0,:])
    find_top_cluster(grid)
    
    connected = False
    
    # Loop over the bottom row of the lattice.
    # If any of the bottom sites is connected to
    # the top, there is a percolation path.
    for j in range(lattice_size):
        if grid[lattice_size-1, j] == 2:
            connected = True
            break
            
    return connected
      


def percolation_mean_and_error(successes, repeats):
    """
    Calculates the average and error of mean for percolation probability.
    
    Assume that the percolation probability for the system is :math:`P`.
    If you generate :math:`n` systems, on average, you will find percolation
    :math:`nP` times. The precise number is, however, a random variable.
    
    Let us define :math:`p` as the percolation result in a *single simulation*.
    If there is percolation, :math:`p = 1`, and if not, :math:`p = 0`.
    Let us also denote the total number of systems with :math:`p = 1` as 
    :math:`k = \\sum_p p`.
    The average of :math:`p` is then
    
    .. math ::
    
        \\mu_p = \\frac{1}{n} \\sum_p p = \\frac{k}{n}.
    
    Since on average :math:`k = nP`, we have :math:`P = k/n`, and so
    the average of :math:`p` is an unbiased estimate for the true probability :math:`P`.
    
    Similarly, the sample variance of :math:`p` is, by definition,
    
    .. math ::
    
        s_p^2 & = \\frac{1}{n-1} \\sum_{p} (p - \\mu_p)^2 \\\\
              & = \\frac{1}{n-1} \\sum_{p} \\left( p^2 - 2 p \\frac{k}{n} + \\frac{k^2}{n^2} \\right) \\\\
              & = \\frac{1}{n-1} \\left( k - 2 \\frac{k^2}{n} + \\frac{k^2}{n} \\right) \\\\
              & = \\frac{k}{n-1} \\left( 1 - \\frac{k}{n} \\right) \\\\
              & = \\frac{k (n-k) }{n (n-1) } .
          
    (Note :math:`\\sum_p p^2 = \\sum_p p = k`, and :math:`\\sum_p 1 = n`.)
          
    The standard error of mean for our estimated :math:`\\mu_p` is then
    
    .. math ::
    
        \\Delta p = \\frac{s_p}{\\sqrt{n}} = \\frac{1}{n} \\sqrt{ \\frac{ k (n-k) }{ (n-1) } }.
    
    Args:
        successes (int): number of systems with percolation, :math:`k`
        repeats (int): total number of systems, :math:`n`
        
    Returns:
        float, float: sample mean :math:`\\mu_p`, error of mean :math:`\\Delta p`
    """
    
    
    return successes/repeats, np.sqrt(successes * (repeats - successes) / (repeats - 1) ) / repeats



def run_simulation(node_probability, lattice_size):
    """
    Creates and analyses a single system.
    
    The created system will be a :math:`N \\times N` square lattice where
    each site is occupied with probability :math:`p`.
    
    The function reports if there is a connection between the top and bottom 
    of the system and shows a picture of the system.
    
    Args:
        node_probability (float): node occupation probability :math:`p`
        lattice_size (int): system size :math:`N`
        
    Returns:
        array: the created lattice
    """

    
    # create a lattice
    grid = generate_network(lattice_size, node_probability)

    # find the lattice sites connected to the top of the lattice
    print("percolation cluster exists:", is_vertically_connected(grid))

    # plot the system
    show_network(grid)
    
    return grid


def run_series(p_min, p_max, p_steps, lattice_size, n_repeats):
    """
    Runs a series of simulations and calculates the probability
    that a connection between the top and bottom of the system exists.
    
    The method varies node occupation probability :math:`p`
    and creates, for each probability, :math:`M` independent systems.
    These systems are :math:`N \\times N` square lattices.
    For each system, the algorithm checks if a connection exists.
    Based on the results, the function estimates the percolation probability
    :math:`P` and its error estimate using :meth:`percolation_mean_and_error`.
    
    Finally, the function plots the calculated percolation probabilities :math:`P`
    as a function of node probability :math:`p`.
    
    Args:
        p_min (float): smallest node probability :math:`p` to check
        p_max (float): highest node probability :math:`p` to check
        p_steps (int): number of different probabilities to check
        lattice_size (int): system size :math:`N`
        n_repeats (int): number of random systems to create for each probability, :math:`M`
        
    Returns:
        list, list: percolation probabilities :math:`P`, error estimates :math:`\\Delta P`
    """

    percolation_ps = []
    percolation_errors = []
    node_ps = np.linspace(p_min, p_max, p_steps)
    
    case = 0
    
    for p in node_ps:
        connected = 0
        
        for i in range(n_repeats):
        
            # create a lattice
            grid = generate_network(lattice_size, p)
            if is_vertically_connected(grid):
                connected += 1
            
            case += 1
            print_progress(case, p_steps*n_repeats)
                
        pp, dp = percolation_mean_and_error(connected, n_repeats)
        percolation_ps.append(pp)
        percolation_errors.append(2*dp) # use 2 x error of mean for a 95 % confidence interval
        
    plt.errorbar(node_ps, percolation_ps, yerr=percolation_errors, fmt='-o', ms=4, lw=2)
    plt.xlabel("node occupation probability, $p$")
    plt.ylabel("percolation probability, $P$")
    plt.ylim([0,1])
    plt.grid()
    plt.show()
        
    return percolation_ps, percolation_errors
    



if __name__ == "__main__":

    # simulation parameters
    lattice_size = 100
    node_probability = 0.6
  
    #run_simulation(node_probability, lattice_size)
    
    lattices = [25, 50, 100, 150]
    reps = 400
    pps = []
    pes = []
    p_min = 0.55
    p_max = 0.65
    nps = 11
    ps = np.linspace(p_min, p_max, nps)
    
    for n in lattices:
        p, e = run_series(p_min, p_max, nps, n, reps)
        pps += [p]
        pes += [e]
    
    for n in range(len(lattices)):
        plt.errorbar(ps, pps[n], yerr=pes[n], label="N = "+str(lattices[n]), fmt='-o', ms=4, lw=2)
        #print(ps)
        #print(pps[n])
        #print(pes[n])
        
    plt.xlabel("node occupation probability, $p$")
    plt.ylabel("percolation probability, $P$")
    plt.legend()
    plt.grid()
    plt.ylim([0,1])
    plt.xlim([p_min, p_max])
    plt.show()
