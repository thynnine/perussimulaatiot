import sys
import copy
from math import *
import numpy as np
from numpy.linalg import inv, eig
import matplotlib.pyplot as plt
import matplotlib.animation as ani



def draw(frame, mesh, T_trajectory, Tlims):
    """
    Draws the temperature profile.
    
    Args:
        frame (int): index of the frame to be drawn
        mesh (array): x-coordinates
        T_trajectory (list): list of arrays of temperatures (at different times)
        Tlims (list): the bottom and top values of the temperature axis: [Tmin, Tmax]
    """

    plt.clf()
    ax = plt.axes()
    ax.set_xlim([mesh[0],mesh[-1]])
    ax.set_ylim(Tlims)
    plt.xlabel("position, x")
    plt.ylabel("temperature, T")
    plt.plot(mesh, T_trajectory[frame], 'o-', markersize=2 )


def animate(mesh, T_trajectory):
    """
    Animates the development of the temperature profile using :meth:`draw`.
    
    Args:
        mesh (array): x-coordinates
        T_trajectory (list): list of arrays of temperatures (at different times)
    """

    nframes = len(T_trajectory)    
    
    print("animating "+str(nframes)+" frames")


    Tmax = T_trajectory[0][0]
    Tmin = T_trajectory[0][0]
    for Ts in T_trajectory:
        for t in Ts:
            if t > Tmax:
                Tmax = t
            if t < Tmin:
                Tmin = t

    Tmax = round(Tmax * 0.12, 0)*10 
    Tmin = round(Tmin * 0.08, 0)*10 
    
    fig = plt.figure()
    motion = ani.FuncAnimation(fig, draw, nframes, interval=100, fargs=(mesh, T_trajectory, [Tmin,Tmax]) )
    plt.show()


def show_history(mesh, T_trajectory, last=0):
    """
    Plots several temperature profiles, from different simulation times, in the same plot.
    
    Args:
        mesh (array): x-coordinates
        T_trajectory (list): list of arrays of temperatures (at different times)
        last (int): only plot this many temperature profiles from the end of the simulation
                    (by default, all data is plotted)
    """
    for values in T_trajectory[-last:]:
        plt.plot(mesh, values, 'o-', markersize=2)
    plt.xlabel("position, x")
    plt.ylabel("temperature, T")
    plt.show()
   

def show_temperature(mesh, temperature):
    """
    Plots the temperature profile.
    
    Args:
        mesh (array): x-coordinates
        temperature (array): temperature at mesh points
    """
    plt.plot(mesh, temperature, 'o-', markersize=2)
    plt.xlabel("position, x")
    plt.ylabel("temperature, T")
    plt.show()


def calculate_stiffness_matrix(dx, n_mesh, cond, printout=False):
    """
    Calculates and returns the stiffness matrix S.
    
    S is a matrix whose elements
    are the overlap integrals of the derivatives of the FEM basis functions.
    
    We define :math:`\\varphi_n(x)` to be the basis function that is
    1 at mesh point :math:`x_n` and zero everywhere else.
    We also denote by :math:`k_{n,m}` the conductivity between mesh points
    :math:`x_n` and :math:`x_m`.    
    The elements of the stiffness matrix are then defined as
    
    .. math ::
    
        S_{i,j} = k_{i,j} \\int_0^L \\varphi'_i(x) \\varphi'_j(x) d x
        
    except where boundary conditions override this definition.
    
    Args:
        dx (float): distance between mesh points
        n_mesh (int): number of mesh points
        cond (float): heat diffusivity :math:`k`, assumed constant
        printout (bool) : If True, the matrix is printed on screen.
                          
    Returns: 
        array: the S matrix
    """
    S = np.zeros( [n_mesh, n_mesh] )

    for i in range(0,n_mesh-1):
        if i > 0:
            S[i,i] = 2.0 / dx * cond
        S[i,i+1] = -1.0 / dx * cond
        S[i+1,i] = -1.0 / dx * cond
    
    S[0,0] = 1.0
    S[0,1] = 0.0
    S[n_mesh-1,n_mesh-1] = 1.0
    S[n_mesh-1,n_mesh-2] = 0.0
    
    if printout:
        # print S with rounded values
        # choose rounding precision according to the smallest non-zero element
        maximum = np.amax( np.abs(S) )
        minimum = np.amin( np.where( np.abs(S)>0, np.abs(S), maximum ) )
        decimals = int( np.ceil( -np.log10( minimum ) ) )        
        print( "" )
        print( "stiffness matrix" )
        print( np.round( S, decimals ) )
        print( "" )

    return S
    
 

def calculate_mass_matrix(dx, n_mesh, printout=False):
    """
    Calculates and returns the mass matrix M.
    
    M is a matrix whose elements
    are the overlap integrals of the FEM basis functions.
    
    We define :math:`\\varphi_n(x)` to be the basis function that is
    1 at mesh point :math:`x_n` and zero everywhere else.    
    The elements of the stiffness matrix are then defined as
    
    .. math ::
    
        M_{i,j} = \\int_0^L \\varphi_i(x) \\varphi_j(x) d x
        
    except where boundary conditions override this definition.
    
    Args:
        dx (float): distance between mesh points
        n_mesh (int): number of mesh points
        printout (bool) : If True, the matrix is printed on screen.
                          
    Returns: 
        array: the M matrix :math:`M`
    """

    M = np.zeros( [n_mesh, n_mesh] )

    for i in range(0,n_mesh-1):
        if i > 0:
            M[i,i] = 2.0 / 3.0 * dx
        M[i,i+1] = dx / 6.0
        M[i+1,i] = dx / 6.0

    M[0,0] = 1.0
    M[0,1] = 0.0
    M[1,0] = 0.0
    M[n_mesh-1,n_mesh-1] = 1.0
    M[n_mesh-1,n_mesh-2] = 0.0
    M[n_mesh-2,n_mesh-1] = 0.0

    if printout:
        # print M with rounded values
        # choose rounding precision according to the smallest non-zero element
        maximum = np.amax( np.abs(M) )
        minimum = np.amin( np.where( np.abs(M)>0, np.abs(M), maximum ) )
        decimals = int( np.ceil( -np.log10( minimum ) ) )
        print( "" )
        print( "mass matrix" )
        print( np.round( M, decimals ) )
        print( "" )
        
    return M    
    
    
def calculate_power_vector(T_left, T_right, heat, M,
                           printout = False):
    """
    Calculates the power vector :math:`p`.
    
    For a fixed boundary (temperature is kept constant):

        * Replace the boundary node values of the heating power density vector
          :math:`\\rho` with the temperature values at the boundary.
          You get :math:`\\rho = [T(x_0), \\rho(x_1), \\rho(x_2), \\ldots, \\rho(x_{N}), T(x_{N+1})]`.
        * Calculate the power vector as :math:`p = M \\rho`.
    
    Args:
        T_left (float): temperature at the left boundary
        T_right (float): temperature at the right boundary
        heat (array): heating power densities :math:`\\rho`
        M (array): the mass matrix :math:`M`
        printout (bool) : If True, the matrix is printed on screen.
                          
    Returns: 
        array: power vector :math:`p`
    """

    # When the temperature is fixed, we set this value in
    # the power vector.
    heat[0] = T_left
    heat[-1] = T_right

    power = M @ heat
    
    if printout:
        print( "" )
        print( "power vector" )
        print( np.round( power, 3 ) )
        print( "" )

    return power
    

def calculate_propagator(S, M, dt):
    """
    Calculates the propagator matrix for a dynamic simulation.
    
    Time is advanced according to the equation
    
    .. math ::
    
        T(t + \\Delta t) = P T(t) + L.
        
    This function calculates the matrix :math:`P` as
    
    .. math ::
    
        P = \\left( M + \\frac{1}{2} \\Delta t S \\right)^{-1} \\left( M - \\frac{1}{2} \\Delta t S \\right).
    
    Args:
        S (array): the stiffness matrix :math:`S`
        M (array): the mass matrix :math:`M`
        dt (float): time step :math:`\\Delta t`
    
    Returns: 
        array: the propagator matrix :math:`P`
    """
    return inv(M+0.5*dt*S) @ (M-0.5*dt*S) 


def calculate_load(S, M, p, dt):
    """
    Calculates the load vector for a dynamic simulation.
    
    Time is advanced according to the equation
    
    .. math ::
    
        T(t + \\Delta t) = P T(t) + L.
        
    This function calculates the vector :math:`L` as

    .. math ::
    
        L = \\left( M + \\frac{1}{2} \\Delta t S \\right)^{-1}  p \\Delta t .
    
    Args:
        S (array): the stiffness matrix :math:`S`
        M (array): the mass matrix :math:`M`
        p (array): the power vector :math:`p`
        dt (float): time step :math:`\\Delta t`
    
    Returns: 
        array: the load vector :math:`L`
    """
    return inv(M+0.5*dt*S) @ (dt*p)


def solve_equilibrium(S, p):
    """
    Solves the equilibrium temperature profile.
    
    The equilibrium satisfies :math:`ST = p` so
    the temperature is given by
    
    .. math ::
    
        T = S^{-1} p.
        
    For fixed boundaries, the matrix :math:`S` should always
    be non-singular and a solution can be found.
    For arbitrary boundary conditions this may not be true.
    
    Args:
        S (array): the stiffness matrix :math:`S`
        p (array): the power vector :math:`p`
        
    Returns:
        array: equilibrium temperatures
    """

    return inv(S) @ p


def run_simulation(temperatures, propagator, load, n_steps):
    """
    Runs a time dependent heat simulation.
    
    Time is advanced according to the equation
    
    .. math ::
    
        T(t + \\Delta t) = P T(t) + L.

    The function does not change the original temperatures array.
    Instead, the final temperatures are returned as a new array.
    
    Args:
        temperatures (array): temperatures
        propagator (array): the propagator matrix :math:`P`
        load (array): the load vector :math:`L`
        n_steps (int): simulation length in number of time steps
    
    Returns: 
        array: temperatures at the end of the simulation
    """

    for i in range(n_steps):
        temperatures = propagator @ temperatures + load
    
    return temperatures
        
    
def main(dx=0.1, n_mesh=11,
         T_left=300, T_right=400,
         heating_power = None,
         dt=0.001, recording_dt=0.01, simulation_time=0):
    """
    The main program.
    
    If simulation_time is 0, only the equilibrium temperature
    distribution is calculated.
    
    If simulation_time > 0, a dynamic simulation is run.
    
    Args:
        dx (float): mesh point separation
        n_mesh (int): number of mesh points
        T_left (float): temperature at the left boundary
        T_right (float): temperateru at the right boundary
        heating_power (array): heating power density at mesh points
        dt (float): time step :math:`\\Delta t`
        recording_dt (float): time between recorded temperature profiles
        simulation_time (float): total time to simulate
    """

    # create a uniform mesh
    mesh = np.linspace(0, (n_mesh-1)*dx, n_mesh)

    # thermal conductivity (heat diffusivity)
    conductance = 1.0   

    if heating_power is None:
        heating_power = np.zeros(n_mesh)

    # calculate matrices
    S = calculate_stiffness_matrix(dx, n_mesh, conductance, printout=True)
    M = calculate_mass_matrix(dx, n_mesh, printout=True)
    p = calculate_power_vector(T_left, T_right, heating_power, M, printout=False)
        
    if simulation_time == 0: # just calculate the equilibrium
    
        # solve the equilibrium temperature profile
        temperatures = solve_equilibrium(S, p)

        print( "equilibrium temperatures ", np.round(temperatures,2) )
    
        # plot the temperature profile
        show_temperature(mesh, temperatures)
    
    else: # run a dynamic simulation    
    
        # start from a linear temperature profile
        # (or if you like, specify the initial profile here)
        temperatures = np.linspace(T_left, T_right, n_mesh)
        
        # simulation parameters
        interval_steps = int(recording_dt / dt)
        records = int(simulation_time / recording_dt)
        
        # calculate propagator and load
        propagator = calculate_propagator(S, M, dt)
        load = calculate_load(S, M, p, dt)
        
        # list for storing the temperature profiles at different times    
        history = [temperatures]
    
        # run the dynamic simulation
        #
        # the simulation is run in short sequences
        # and temperatures are saved between each sequence
        for i in range(records):
            temperatures = run_simulation(temperatures, propagator, load, interval_steps)
            history.append(temperatures)

        
        print( "final temperatures ", np.round(temperatures,2) )

        # Plot the development of the temperature profile.
        show_history(mesh, history)
        animate(mesh, history)
        
    
if __name__ == "__main__":

    
    dx = 0.1 # mesh spacing
    n_mesh = 11 # number of mesh points 

    # bounding temperatures
    T_left = 370
    T_right = 270
    
    # heating power density
    heat = np.array( [0,0,0,0,1000,1000,1000,0,0,0,0] )

    main( dx, n_mesh, T_left, T_right, heat, simulation_time = 0 )
