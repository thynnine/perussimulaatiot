import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
from math import *
from scipy.optimize import curve_fit
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
        

def show_grid(grid):
    """
    Draws the system as a :math:`N \\times N` pixel image.
    
    Args:
        grid (array): the system as an integer array
    """
    plt.clf()
    ax = plt.axes()
    ax.set_aspect('equal') 
    plt.pcolormesh(grid, cmap='Greys', vmin=-1, vmax=1)
    plt.show()
        
        
        
def save_grid_image(grid, index):
    """
    Saves an image of the system as a :math:`N \\times N` pixel image.
    
    The image is saved in a file named "ising_X.png" where
    X is the given index.
    
    Args:
        grid (array): the system as an integer array
        index (int): index for naming the output file
    """
    plt.clf()
    ax = plt.axes()
    ax.set_aspect('equal') 
    plt.pcolormesh(grid, cmap='Greys')
    plt.savefig('ising_'+str(i)+'.png', bbox_inches='tight', vmin=-1, vmax=1)
    plt.close()

    
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
    plt.pcolormesh(history[frame], cmap='Greys', vmin=-1, vmax=1)


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


    #print("calling ffmpeg")
    #filename = "L2_ising_T.mp4"
    #motion.save(filename,fps=30,codec='mpeg4')
    #print("animation saved as "+filename)
    #plt.clf()


def generate_grid(lattice_size, randomly = False):
    """
    Creates the system as :math:`N \\times N` array.
    
    Each value in the array represents a spin with two possible values, up and down.
    The value up is represented by the value 1
    and the value down by -1.
    
    Args:
        lattice_size (int): system size :math:`N`
        randomly (bool): If True, spins will be randomly oriented. Otherwise all spins point up.
    """

    if randomly:
        grid = 2*random.integers(low=0,high=2,size=(lattice_size, lattice_size))-1
    else:
        grid = np.ones( [ lattice_size, lattice_size] )

    return grid





def calculate_flip_energy(grid, i,j):
    """
    Calculates how much the total energy of the system would change,
    if the spin :math:`(i,j)` was flipped.
    
    The system is not changed in any way. The function only calculates
    how much the energy would change, if a spin was flipped.
    
    The total energy of the 2D Ising model is
    
    .. math::
    
        E = - \\sum_{a,b} J_{a,b} s_{a} s_{b},
        
    where :math:`a = (i,j)` and :math:`b = (k,l)` are spin
    coordinates, :math:`s_a, s_b` are the spin values at these
    coordinates (either 1 or -1), and :math:`J_{a,b}` is the
    coupling parameter between these two spins.
    
    * If :math:`J > 0`, the spins interact ferromagnetically, i.e., 
      they prefer to point in the same direction.
    * If :math:`J < 0`, the spins interact antiferromagnetically
      and prefer to point in opposite directions.
    * If :math:`J = 0`, there is no interaction between the spins.
    
    In our model, :math:`J = 1` for all spins which are directly
    next to each other (vertical and horizontal neighbors), 
    and :math:`J = 0` for all other spins (in some suitable units).
    This means that if only one spin, :math:`s_(i,j)` changes its sign,
    there are exactly four terms in the energy sum that change: 
    the terms with :math:`s_{i,j}` and its neighbors up, down, right and left,
    i.e., :math:`s_{i+1,j}, s_{i-1,j}, s_{i,j+1}, s_{i,j-1}`.
    Furthermore they change by changing their sign.
    For instance, if :math:`s_{i,j,\\text{before}} = 1` and :math:`s_{i+1,j} = -1`,
    initially the energy of this pair is        
    :math:`-J_{(i,j),(i+1,j)} s_{i,j,\\text{before}} s_{i+1,j} = 1`.
    After the flip, we would have :math:`s_{i,j,\\text{after}} = -1`, 
    and the energy of the pair would become -1.
    The *change* in energy would therefore be
    
    .. math ::
        \\Delta E_1 = 2 J_{(i,j),(i+1,j)} s_{i,j,\\text{before}} s_{i+1,j} = -2.
    
    The energies for the three other affected pairs are calculated similarly,
    so substituting :math:`J = 1` we get the total energy change
    
    .. math::
    
        \\Delta E = 2s_{i,j,\\text{before}} ( s_{i+1,j} + s_{i-1,j} + s_{i,j+1} + s_{i,j-1}).
        
    Args:
        grid (array): the system as an integer array
        i (int): row of the spin to flip
        j (int): column of the spin to flip
    
    Returns:
        int: energy change :math:`\\Delta E`
    """
    energy = 0.0
    lattice_size = len(grid[0])
    energy += 2*(grid[i,j]*grid[(i-1)%lattice_size,j])
    energy += 2*(grid[i,j]*grid[(i+1)%lattice_size,j]) 
    energy += 2*(grid[i,j]*grid[i,(j-1)%lattice_size]) 
    energy += 2*(grid[i,j]*grid[i,(j+1)%lattice_size]) 

    return energy


def flip(grid, i,j):
    """
    Flips the spin :math:`s_{i,j}` in the system.
    
    Since spin values up and down are represented by values of +1 and -1,
    respectively, the function simply changes the sign of the given
    element in the array.
    
    Args:
        grid (array): the system as an integer array
        i (int): row of the spin to flip
        j (int): column of the spin to flip
    """
    grid[i,j] = -grid[i,j]
    

def metropolis_step(grid, temperature):
    """
    Tries to flip one randomly chosen spin using the Metropolis algorithm.
    
    The flip is either accepted or discarded according to a Boltzmann-like
    probability:
    
        * If the spin flip leads to *lower* total energy, it is always accepted.
        * If the spin flip leads to *higher* total energy, the new configuration
          is accepted with probability :math:`P = e^{- \\frac{ \Delta E }{k_B T} }`,
          where :math:`\\Delta E` is the energy change, :math:`k_B` is the Boltzmann
          constant and :math:`T` is the temperature.
    
    If the spin flip is accepted, the change is applied to the given system.
    
    The function uses a unitless temperature scale with :math:`k_B = 1`.
    
    Args:
        grid (array): the system as an integer array
        temperature (float): temperature
    """
    # pick a random spin
    lattice_size = len(grid[0])
    coords = random.integers(low=0, high=lattice_size, size=2)
    i = coords[0]
    j = coords[1]
    delta_ene = calculate_flip_energy(grid,i,j)
    
    # if the Metropolis test succeeds, the spin is flipped
    randomizer = random.random()
    if randomizer < exp( - delta_ene / temperature ):
        flip(grid,i,j)
       
           
def calculate_magnetic_moment(grid):
    """
    Calculates the total magnetization of the system.
    
    A magnetic moment is associated with electron spins,
    :math:`\\vec{\\mu} = -g\\frac{e}{2 m_e} \\vec{s}`.
    In the ising model, spins only point either up or down
    and get values of +1 or -1, and so
    also the magnetic moment points only up or down.
    Furthermore, by choosing a unitless scale :math:`g\\frac{e}{2 m_e} = 1`,
    the magnetic moment becomes the same as the spin (only reversed), and
    we can calculate the total magnetization as a sum over
    all spins in the system,
    
    .. math::
    
        M = -\\frac{1}{N^2} \\sum_{i,j} s_{i,j}.
        
    Here :math:`N \\times N` is the number of spins in the system.
    Thus :math:`M` is scaled so that is only gets values between
    -1 and 1. If :math:`|M| = 1`, the magnetization is complete and
    all spins point in the same direction. If :math:`M = 0`, there
    are equally many spins pointing up and down, and there is no
    macroscopic magnetic moment.
    
    Args:
        grid (array): the system as an integer array
        
    Returns:
        float: magnetization absolute value :math:`|M|`
    """
    magnetic_moment = 0
    lattice_size = len(grid[0])
    for i in range(lattice_size):
        for j in range(lattice_size):
            magnetic_moment += grid[i,j]
            
    return np.abs(magnetic_moment)/(lattice_size**2)
    
    
def calculate_Binder_cumulant(magnetizations):
    """
    Calculates the Binder cumulant.
    
    The transition from an ordered, ferromagnetic state to a disordered,
    paramegnetic state should be discontinuous at the thermodynamic limit,
    but in a finite system, the transition is gradual. This makes it difficult
    to pinpoint the critical temperature :math:`T_C` from a finite simulation.
    
    This can be done very precisely using the cumulant method.
    The quantity
    
    .. math ::
    
        U_4 = 1 - \\frac{ \\langle M^4 \\rangle }{ 3 \\langle M^2 \\rangle^2 },
        
    known as the Binder cumulant, is also dependent on system size.
    At temperatures below :math:`T_C`, it approaches 0 as the system size grows.
    At temperatures above :math:`T_C`, it approaches 2/3 for the Ising model.
    At :math:`T_C`, it has a non-trivial system dependent value, which is approximately
    independent of system size.
    
    Therefore one can pinpoint the transition temperature by calculating the 
    cumulant as a function of temperature at different system sizes
    and plotting the results. The cumulants will cross approximately at :math:`T_C`.
    
    Args:
        magnetizations (array): sequence of magnetization values
        
    Returns:
        float: the cumulant :math:`U_4` 
    """

    m4 = 0.0
    m2 = 0.0
    n = 0
    
    for m in magnetizations:
        m4 += m**4
        m2 += m**2
        n += 1
        
    m4 /= n
    m2 /= n
    
    return 1 - m4/(3*m2**2)
    
    
def write_data_to_file(lattice_size, temps, mags):
    """
    Writes a data series in a file.
    
    The data is written in file "ising_data_N.txt" where
    N is the lattice size.
    
    The data is written so that each column holds
    a temperature and the calculated magnetization at that temperature.
    
    Args:
        lattice_size (int): system size :math:`N`
        temps (list): a list of temperatures
        mags (list): a list of magnetizations
    """

    writelines = ""
    for t, m in zip(temps, mags):
        writelines += str(t)+" "+str(m)+" \n"

    f = file('ising_data_'+str(lattice_size)+".txt",'w')
    f.write(writelines)
    f.close()


    

def calculate_statistics(samples):
    """
    Calculates and returns the sample mean, variance, 
    standard deviation and error of mean for the given set of samples.
    
    For a sequence of :math:`N` values for a random variable :math:`X`,
    the mean is estimated by the average
    
    .. math ::
        
        \\mu_X = \\langle X \\rangle = \\frac{1}{N} \\sum_i X_i
        
    and the variance by the sample variance
    
    .. math ::
    
        \\sigma^2_X = \\frac{1}{N-1} \\sum_i (X_i - \\mu_X)^2.
        
    Standard deviation is the square root of variance
    
    .. math ::
    
        \\sigma_X = \\sqrt{ \\sigma^2_X }
        
    and the standard error of mean is
    
    .. math ::
    
        \\Delta X = \\frac{ \\sigma_X }{ \\sqrt{ N } }.
    
    Args:
        samples (list): sequence of random results, :math:`X_i`
        
    Returns:
        float, float, float, float: mean, variance, standard deviation, error of mean
    """
    expect = 0.0
    variance = 0.0
    
    expect = 0.0
    variance = 0.0
    n_repeats = len(samples)
    
    for i in samples:
        expect += i
    expect /= float(n_repeats)
    
    for i in samples:
        variance += (i-expect)*(i-expect)
    variance /= float(n_repeats-1)
    
    std = np.array( sqrt(variance) )
    error = std / sqrt(n_repeats)
    
    return expect, variance, std, error
    
    

def run_cumulant_series(n_steps = 10000, thermalization_steps = 1000, lattice_sizes = [5, 10, 15],
                        min_temperature = 1.0, max_temperature = 3.0, t_steps = 11):
    """
    Runs a series of simulations at different temperatures and lattice sizes.
    
    The function calculates Binder cumulants for different systems and plots these.
    If the statistics are good enough, these cumulants should intersect approximately
    at the critical point :math:`T_C`.
    
    
    Args:
        n_steps (int): total simulation cycles (each cycle is :math:`N^2` spin flip attempts)
        thermalization_steps (int): thermalization cycles in the beginning of the simulation
        lattice_sizes (list): list of lattice sizes :math:`N`
        min_temperature (float): the inital temperature
        max_temperature (float): the final temperature (assumed larger than min_temperature)
        t_steps (int): the number of different temperatures to simulate
    """

    cumulants = []
    
    for lattice_size in lattice_sizes:
        print("simulating system N = "+str(lattice_size))
        grid = generate_grid(lattice_size)
        ts, ms, es, cs = run_temperature_series(grid, 
                                            min_temperature,
                                            max_temperature,
                                            t_steps,
                                            n_steps,
                                            thermalization_steps)

        cumulants.append(cs)

    for i in range( len(sizes) ):
        plt.plot(ts,cumulants[i],'-',label="N = "+str(sizes[i]))
    plt.legend()
    plt.show()  



def run_temperature_series(grid, min_temperature, max_temperature, 
                           temp_steps, run_steps, thermalization_steps):
                           
    """
    Runs a series of simulations at different temperatures.
    
    The function calculates the average magnetization and its error estimate (as
    the standard error of mean) at each temperature and returns these as arrays.
    
    Args:
        grid (array): the system as an integer array
        min_temperature (float): the inital temperature
        max_temperature (float): the final temperature (assumed larger than min_temperature)
        temp_steps (int): the number of different temperatures to simulate
        run_steps (int): total simulation cycles (each cycle is :math:`N^2` spin flip attempts)
        thermalization_steps (int): thermalization cycles in the beginning of the simulation
        
    Returns:
        array, array, array, array: temperatures, magnetizations, error estimates, cumulants
    """
     
    # lists for storing measurements
    m_averages = []
    m_errors = []
    cumulants = []     
                           
    # list of temperatures to simulate
    temperatures = np.linspace(min_temperature, max_temperature, temp_steps)
          
    # Run a separate simulation at each temperature.
    # Use the final state of the previous simulation
    # as a starting point for the next one.
    for T in temperatures:
    
        # run simulation
        print("starting simulation at T = "+str(T))
        magnetization = run_metropolis_algorithm(grid, T, run_steps, thermalization_steps)
        
        # calculate average magnatization and its error
        avr, var, std, err = calculate_statistics(magnetization)
        m_averages.append(avr)
        m_errors.append(err)
        
        # calculate cumulant
        c = calculate_Binder_cumulant(magnetization)
        cumulants.append(c)
     
    return temperatures, np.array(m_averages), np.array(m_errors), np.array(cumulants)
    
    
def run_metropolis_algorithm(grid, temperature, run_steps, thermalization_steps, recording_interval = 10, animated=False):
    """
    Runs the Metropolis algorithm at a contant temperature.
    
    The algorithm repeatedly chooses random spins in the system
    and tries to flip them using :meth:`metropolis_step`.
    The changes are applied directly on the given grid.
    
    After a thermalization period, the algorithm also calculates
    the magnetization of the system. 
    By default, the recording is not done at every step because too
    closely measured values will be correlated.
    
    Args:
        grid (array): the system as an integer array
        temperature (float): the temperature
        run_steps (int): total simulation cycles (each cycle is :math:`N^2` spin flip attempts)
        thermalization_steps (int): thermalization cycles in the beginning of the simulation
        recording_interval (int): number of cycles between data recording
        animated (bool): if True, shows an animation of the simulation in the end
        
    Returns:
        array: sequence of magnetizations
    """

    magmoms = []
    history = []
    n_spins = len(grid)**2 # total number of spins
    
    # run the simulation for the required time
    for i in range(run_steps):
        print_progress(i,run_steps)

        # for each cycle, try to flip spins as many times
        # as there are spins in the system
        for j in range(n_spins):
            metropolis_step(grid, temperature)
        
        # only record data after thermalization is complete
        if i >= thermalization_steps and i%recording_interval == 0:
            moment = calculate_magnetic_moment(grid)
            magmoms.append( abs(moment) )
            
            if animated:
                history.append(copy.deepcopy(grid))
        
    print_progress(run_steps, run_steps)
    
    if animated:
        animate(history)
    
    return np.array(magmoms)


def main():
    """
    Main program.
    """

    lattice_size = 100
    n_steps = 1000
    thermalization_steps = 0
    T = 2.20
    
    grid = generate_grid(lattice_size, randomly=False)
    mag = run_metropolis_algorithm(grid, T, n_steps, thermalization_steps, recording_interval=20, animated=True)

    # calculate statistics
    m_avr, m_var, m_std, m_err = calculate_statistics(mag)

    # plot magnetization and the calculated average +- confidence interval
    plt.plot(mag)
    plt.plot(np.linspace(0, len(mag), len(mag)), m_avr*np.ones(len(mag)), 'b:')
    plt.fill_between(np.linspace(0, len(mag), len(mag)), (m_avr-2*m_err)*np.ones(len(mag)), (m_avr+2*m_err)*np.ones(len(mag)), color = 'b', alpha=0.1 )
    plt.xlabel("simulation step")
    plt.ylabel("magnetization")
    plt.show()
 
    quit()
    
    min_temperature = 2.1
    max_temperature = 2.5
    temp_steps = 5
    ts, ms, mes, cs = run_temperature_series(grid, min_temperature, max_temperature, 
                           temp_steps, n_steps, thermalization_steps)
    plt.errorbar(ts, ms, yerr=mes)
    plt.show()
    
    


if __name__ == "__main__":
    random = default_rng()
    main()
else:
    random = default_rng()
    
