# /usr/bin/env python
import sys
import copy
from math import *
import cmath
import numpy as np
import numpy.fft as ft
import matplotlib.pyplot as plt
import matplotlib.animation as ani


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


def draw(frame, wave, x):
    """
    Draws the wave. 
    
    Used for animation.
    
    Args:
        frame (int): index of the frame to be drawn
        wave (list): list of waves as a time series
        x (array): x coordinates
    """
    plt.clf()
    ax = plt.axes()
    ax.set_ylim([-2,2])
    plt.xlabel("position, x")
    plt.ylabel("wavefunction, u")
    plt.plot(x, wave[frame])
    
    
def animate(wave, L):
    """
    Animates the simulation.
    
    Args:
        wave (list): list of waves as a time series
        L (float): length of simulated region :math:`L`
    """
    
    nframes = len(wave)
    n_nodes = len(wave[0])
    x = np.linspace(0,L,n_nodes)
    fig = plt.figure()
    motion = ani.FuncAnimation(fig, draw, nframes, interval=10, fargs=( wave, x ) )
    plt.show()
    
      
    
def add_gaussian(string, a = 1, center = 0.5):
    """
    Initialize the wave to a gaussian initial shape.
    
    The function adjusts the shape of the wave by adding
    a narrow gaussian function to the wave function. 
    Since this operation is additive,
    it can be called several times to add several pulses.
    
    Args:
        string (array): discretized wave function to be modified
        a (float): pulse height
        center (float): fractional position of the peak, should be between 0 and 1
    """
    n_nodes = len(string)
    for i in range(n_nodes):
        string[i] += a * exp( - ( 5.0/float(n_nodes) ) * (i-n_nodes*center+0.5)**2 )


def add_sinusoidal(string, n, a = 1):
    """
    Initialize the wave to a sinusoidal initial shape.
    
    The function adjusts the shape of the wave by adding
    a sinus function to the wave function. 
    The sinus function will go to zero at the edges of the system.
    Since this operation is additive,
    it can be called several times to build the wave as a
    superposition of sinus functions.
    
    Args:
        string (array): discretized wave function to be modified
        n (int): number of antinodes
        a (float): amplitude
    """
    n_nodes = len(string)
    for i in range(n_nodes):
        string[i] += sin(n*pi*(i+1)/(n_nodes+1))
    

def add_triangle(string, a = 1, peak = 0.5):
    """
    Initialize the wave to a triangular initial shape.
    
    The function adjusts the shape of the wave by adding
    a triangular function to the wave function.
    The triangle will go to zero at the edges of the system.
    Since this operation is additive,
    it can be called several times to add several triangles.
    
    Args:
        string (array): discretized wave function to be modified
        a (float): amplitude
        peak (float): fractional position of the peak, should be between 0 and 1
    """

    n_nodes = len(string)
    top = peak*n_nodes-0.5
    for i in range(n_nodes):
        if i < top:
            string[i] = a*(i+1)/(top+1)
        else:
            string[i] = a*(i-n_nodes)/(top-n_nodes)
        


def fourier_analysis(wave):
    """
    Calculates coefficients of the discrete complex Fourier transformation.
    
    A wave function defined in :math:`[0, L]` can be represented as a Fourier series
    
    .. math ::
    
        u(x) = \\sum_{n = 0}^N c_n e^{i \\frac{ 2 \\pi n }{ L } x }.
        
    This function calculates the coefficients :math:`c_n` using numpy's
    fast Fourier transform functionality.
    
    Args:
        wave (array): discretized wave function
        
    Returns:
        array: Fourier coefficients :math:`c_n`
    """
    
    return ft.rfft(wave)


def fourier_synthesis(c, v, L, t, n_nodes):
    """
    Calculates the wave function from Fourier coefficients.
    
    A wave function defined in :math:`[0, L]` can be represented as a time
    dependent Fourier series
    
    .. math ::
    
        u(x,t) = \\sum_{n = 0}^N c_n e^{i \\left( \\frac{ 2 \\pi n }{ L } x - \\omega_n t \\right) }
        
    where the angular frequency of component :math:`n` is 
    :math:`\\omega_n = 2 \\pi n v / L`.
        
    This function reconstructs a discretized approximation for the wave function 
    at a given time from the Fourier coefficients.
    
    Args:
        c (array): Fourier coefficients :math:`c_n`
        v (float): wave speed :math:`v`
        L (float): length of simulated region :math:`L`
        t (float): time :math:`t`
        n_nodes (int): number of grid points in the original wave function
        
    Returns:
        array: discretized wave function :math:`u`
    """

    ct = [0.0] * n_nodes
    
    I = complex(0,1)
        
    cutoff = len(c)
    for n in range(cutoff):
        ct[n] = c[n]*cmath.exp( -2*pi*I*n*v/L*t )
    
    wave = ft.irfft(ct, n_nodes)
    
    return wave



def cos_analysis(wave):
    """
    Calculates coefficients of the discrete cosine series.
    
    A wave function defined in :math:`[0, L]` can be represented as a series
    
    .. math ::
    
        u(x) = \\sum_{n = 0}^N a_n \\cos\\left( \\frac{ \\pi n }{ L } x \\right),
        
    if slope of the function is zero at both boundaries, :math:`u'(0) = u'(L) = 0`.
        
    This function calculates the coefficients :math:`a_n` using numpy's
    fast Fourier transform functionality. This is done by taking the original function and
    extending it to :math:`[0, 2L]` so that the result is symmetric with respect to
    the point :math:`x = L`. Due to symmetry, the complex Fourier transform of this function
    will yield real valued coefficients which are also the coefficients for the
    cosine series of the original function.
    
    Args:
        wave (array): discretized wave function
        
    Returns:
        array: Expansion coefficients :math:`a_n`
    """

    extension = copy.copy(wave)
    extension = np.flip(extension)
    extended_function = np.concatenate( ( wave, extension ) )
    
    return np.real( ft.rfft(extended_function) )



def cos_synthesis(a, v, L, t, n_nodes, visual = False):
    """
    Calculates the wave function from cosine series coefficients.
    
    A wave function defined in :math:`[0, L]` obeying :math:`u'(0) = u'(L) = 0` 
    can be represented as a time dependent cosine series
    
    .. math ::
    
        u(x,t) = \\sum_{n = 0}^N a_n \\cos\\left( \\frac{ \\pi n }{ L } x \\right) \\cos (\\omega_n t)
        
    where the angular frequency of component :math:`n` is 
    :math:`\\omega_n = \\pi n v / L`.
    This assumes the wave was stationary at :math:`t = 0`.
        
    This function reconstructs a discretized approximation for the wave function 
    at a given time from the expansion coefficients.
    
    Note that if the original discretized wave function was known at :math:`N` grid points,
    this function reconstruct a wave function at :math:`N+2` points.
    This is because the original function can be almost anything but the end result must
    obey the boundary conditions :math:`u'(0) = u'(L) = 0`. Therefore the reconstruction
    adds one grid point at both ends of the function so that it
    explicitly satisfies this condition.
    
    Args:
        a (array): cosine series coefficients :math:`a_n`
        v (float): wave speed :math:`v`
        L (float): length of simulated region :math:`L`
        t (float): time :math:`t`
        n_nodes (int): number of grid points in the original wave function, :math:`N`
        visual (bool): if True, the harmonic components and their superposition is plotted
        
    Returns:
        array: discretized wave function :math:`u`
    """

    wave = np.zeros(n_nodes+2)
    
    if visual:
        plt.clf()
        ax = plt.axes()
        ax.set_ylim([-2,2])
    
    for k in range(len(a)):
        # the harmonic component
        pwave = 1/(n_nodes+1) * np.cos( pi*k* np.linspace( 0, 1 , n_nodes+2 ) )

        # time evolution
        pwave *= np.cos( pi*k*v*t/L )
                
        # expansion coefficient
        pwave *= a[k]
        
        # series expansion
        wave += pwave

        if visual:
            plt.plot(pwave)
        
    if visual:
        plt.plot(wave)
        plt.show()
    
    return wave


def sin_analysis(wave):
    """
    Calculates coefficients of the discrete sine series.
    
    A wave function defined in :math:`[0, L]` can be represented as a series
    
    .. math ::
    
        u(x) = \\sum_{n = 0}^N b_n \\sin \\left( \\frac{ \\pi n }{ L } x \\right),
        
    if the function is zero at both boundaries, :math:`u(0) = u(L) = 0`.
        
    This function calculates the coefficients :math:`a_n` using numpy's
    fast Fourier transform functionality. This is done by taking the original function and
    extending it to :math:`[0, 2L]` so that the result is antisymmetric with respect to
    the point :math:`x = L`. Due to symmetry, the complex Fourier transform of this function
    will yield imaginary valued coefficients the absolute values of which are also the 
    coefficients for the sine series of the original function.
    
    Args:
        wave (array): discretized wave function
        
    Returns:
        array: Expansion coefficients :math:`b_n`
    """

    extension = copy.copy(wave)
    extension = np.flip(-extension)
    extended_function = np.concatenate(([0], wave, [0], extension ))
    
    return np.imag( ft.rfft(-extended_function) )



def sin_synthesis(b, v, L, t, n_nodes, visual = False):
    """
    Calculates the wave function from sine series coefficients.
    
    A wave function defined in :math:`[0, L]` obeying :math:`u(0) = u(L) = 0` 
    can be represented as a time dependent sine series
    
    .. math ::
    
        u(x,t) = \\sum_{n = 0}^N b_n \\sin \\left( \\frac{ \\pi n }{ L } x \\right) \\cos (\\omega_n t)
        
    where the angular frequency of component :math:`n` is 
    :math:`\\omega_n = \\pi n v / L`.
    This assumes the wave was stationary at :math:`t = 0`.
        
    This function reconstructs a discretized approximation for the wave function 
    at a given time from the expansion coefficients.
    
    Note that if the original discretized wave function was known at :math:`N` grid points,
    this function reconstruct a wave function at :math:`N+2` points.
    This is because the original function can be almost anything but the end result must
    obey the boundary conditions :math:`u(0) = u(L) = 0`. Therefore the reconstruction
    adds one grid point at both ends of the function so that it
    explicitly satisfies this condition.
    
    Args:
        b (array): sine series coefficients :math:`b_n`
        v (float): wave speed :math:`v`
        L (float): length of simulated region :math:`L`
        t (float): time :math:`t`
        n_nodes (int): number of grid points in the original wave function, :math:`N`
        visual (bool): if True, the harmonic components and their superposition is plotted
        
    Returns:
        array: discretized wave function :math:`u`
    """

    wave = np.zeros(n_nodes+2)
    
    if visual:
        plt.clf()
        ax = plt.axes()
        ax.set_ylim([-2,2])
    
    for k in range(len(b)):
        # the harmonic component
        pwave = 1/(n_nodes+1) * np.sin( pi*k* np.linspace( 0, 1 , n_nodes+2 ) )

        # time evolution
        pwave *= np.cos( pi*k*v*t/L )
        
        # expansion coefficient
        pwave *= b[k]
        
        # series expansion
        wave += pwave

        if visual:
            plt.plot(pwave)
        
    if visual:
        plt.plot(wave)
        plt.show()
        
    return wave
    
    
def calculate_transformation(type, wave):
    """
    Calculates expansion coefficients.
    
    The function calls one of the Fourier analysis methods:
    
        * :meth:`fourier_analysis` if type is "exp"
        * :meth:`sin_analysis` if type is "sin"
        * :meth:`cos_analysis` if type is "cos"
        
    Args:
        type (str): one of "exp", "sin" or "cos"
        wave (array): discretized wave function
        
    Returns:
        array: Expansion coefficients
    """
    
    if type == "sin":
        coefficients = sin_analysis( wave )
    elif type == "cos":
        coefficients = cos_analysis( wave )      
    elif type == "exp":
        coefficients = fourier_analysis( wave )

    return coefficients
        

def calculate_inverse_transformation(type, coefficients, v, L, t, n_nodes):
    """
    Reconstructs wave function from expansion coefficients.
    
    The function calls one of the Fourier analysis methods:
    
        * :meth:`fourier_synthesis` if type is "exp"
        * :meth:`sin_synthesis` if type is "sin"
        * :meth:`cos_synthesis` if type is "cos"
        
    Args:
        type (str): one of "exp", "sin" or "cos"
        coefficients (array): series coefficients :math:`b_n`
        v (float): wave speed :math:`v`
        L (float): length of simulated region :math:`L`
        t (float): time :math:`t`
        n_nodes (int): number of grid points in the original wave function, :math:`N`
        
    Returns:
        array: discretized wave function :math:`u`
    """
    
    if type == "sin":
        wave = sin_synthesis( coefficients, v, L, t, n_nodes )
    elif type == "cos":
        wave = cos_synthesis( coefficients, v, L, t, n_nodes )        
    elif type == "exp":
        wave = fourier_synthesis( coefficients, v, L, t, n_nodes )

    return wave
    
    
def run_simulation(type, initial_wave, v, L, time, dt, cutoff, plot_coefficients = False):
    """
    Calculates the wave function at different time steps.
    
    The wave function is solved at equally spaced time steps
    by first calculating the harmonic series expansion coefficients
    of the wave function at time :math:`t = 0`.
    Once the series expansion is known, the wave function
    can be calculated at any time by calculating the phases of all the
    harmonic components at the given time and reconstructing the
    wave function from the exansion coefficients and the time-dependent
    harmonic functions.
    
    It is also possible to calculate approximate solutions by only including
    some of the first expansion coefficients in the calculation.
    This saves both calculation time and memory.
    
    Args:
        type (str): one of "exp", "sin" or "cos"
        initial_wave (array): discretized wave function at the beginning
        v (float): wave speed :math:`v/L`
        L (float): length of simulated region :math:`L`
        time (float): total simulation time
        dt (float): time between recorded wave functions
        cutoff (int): the number of expansion coefficients to include
        plot_coefficients (bool): If True, expansion coefficients are plotted.
            All coefficients are shown as black dots. The ones included in
            the simulation are shown with a white center.
        
    Returns:
        list: time evolution of the wave function
    """
    
    n_nodes = len(initial_wave)
    n_steps = int(time/dt)+1
    
    coefficients = calculate_transformation(type, initial_wave)
    if cutoff < len(coefficients):
        c = coefficients[:cutoff]
    else:
        print("cutoff is too large, using full accuracy")
    
    history = []
    
    for i in range(n_steps):

        t = i*dt
        wave = calculate_inverse_transformation(type, c, v, L, t, n_nodes )
        history += [wave]
        
        print_progress(i+1, n_steps)
    
    if plot_coefficients:
        plt.plot(np.real(coefficients), 'ko', markersize=3, label="all coefficients")
        plt.plot(np.real(c), 'wo', markersize=2, label="included in simulation")
        plt.xlabel('k')
        plt.ylabel('Re[$c_k$]')
        plt.legend()
        plt.show()
    
    return history
    

def main():
    """
    Main function. Simulates a 1D wave using harmonic basis functions.
    """

    # simulation parameters
    n_nodes = 100
    length = 2.0
    cutoff = 30
    type = "sin"

    # wave speed
    v = 1.0    
    
    # timing
    simulation_time = 4.0
    sample_dt = 0.05
    
    # initial shape of the wave
    wave_initial = np.zeros(n_nodes)
    add_triangle(wave_initial, a=1, peak=0.2)
    #add_gaussian(wave_initial, center=0.2)
    
    # check the expansion at starting time
    ft = sin_analysis(wave_initial)
    sin_synthesis(ft, v, length, 0, n_nodes, visual=True)

    # run simulation
    wavefunction = run_simulation(type, wave_initial, v, length, simulation_time, sample_dt, cutoff, True)

    # show the result
    animate(wavefunction, length)
    
    
    
if __name__ == "__main__":
    main()
