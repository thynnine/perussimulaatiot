# /usr/bin/env python
import sys
import copy
from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from numpy.random import default_rng



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


def draw(frame, xtraj, ytraj, ztraj, bounds):
    """
    Draws a representation of the particle system as a scatter plot. 
    
    Used for animation.
    
    Args:
        frame (int): index of the frame to be drawn
        xtraj (array): x-coordinates of all particles at different animation frames
        ytraj (array): y-coordinates at all particles at different animation frames
        ztraj (array): z-coordinates at all particles at different animation frames
        bounds (array): list of lower and upper bounds for the plot as [[xmin, xmax], [ymin, ymax]]
    """

    plt.clf()
    ax = plt.axes()
    ax.set_xlim(bounds[0])
    ax.set_ylim(bounds[1])
    ax.set_aspect('equal')
    size = (5*ztraj[frame])*2
    size = np.where( size > 1, size, np.ones( len( ztraj[frame] ) ) )
    plt.scatter(xtraj[frame], ytraj[frame], marker='o', s=size )
        


def animate( particles, box, multiply = [3,3] ):
    """
    Animates the simulation.
    
    Args:
        particles (list): list of :class:`temp.Atom` objects
        box (temp.PeriodicBox): supercell
        multiply (array): number of periodic images to draw in x and y directions
    """

    nframes = len(particles[0].trajectory)    
    
    print("animating "+str(nframes)+" frames")

    xtraj = []
    ytraj = []
    ztraj = []
    
    Lx = box.lattice[0]
    Ly = box.lattice[1]

    # number of periodic images in x and y directions
    multix = multiply[0]
    multiy = multiply[1]
    margin = 0.1

    bounds = np.zeros([3,2])
    bounds[0,1] = Lx*multix
    bounds[1,1] = Ly*multiy
    bounds[:,0] -= margin
    bounds[:,1] += margin

    for i in range(nframes):
    
        xtraj.append([])
        ytraj.append([])
        ztraj.append([])
        
        for p in particles:
            box.shift_inside_box(p.trajectory[i])
            for ix in range(multix):
                for iy in range(multiy):
                    xtraj[-1].append(p.trajectory[i][0] + ix*Lx)
                    ytraj[-1].append(p.trajectory[i][1] + iy*Ly)
                    ztraj[-1].append(p.trajectory[i][2])
            
    xtraj=np.array(xtraj)
    ytraj=np.array(ytraj)   
    ztraj=np.array(ztraj)   
    
    fig = plt.figure()
    motion = ani.FuncAnimation(fig, draw, nframes, interval=10, fargs=(xtraj, ytraj, ztraj, bounds) )
    plt.show()



def show_particles(particles):
    """
    Plot a 2D-projection (xy-coordinates) of the system.
    
    Args:
        particles (list): list of :class:`temp.Atom` objects
    """

    N = len(particles)
    xcoordinates = np.array( [0.0]*N )
    ycoordinates = np.array( [0.0]*N )
    zcoordinates = np.array( [0.0]*N )

    i = 0
    for part in particles:
        xcoordinates[i] = part.position[0]
        ycoordinates[i] = part.position[1]
        zcoordinates[i] = part.position[2]
        i += 1

    plt.clf()
    ax = plt.axes()
    ax.set_aspect('equal','box')
    plt.plot(xcoordinates, ycoordinates, 'o')
    plt.show()

    


def show_trajectories(particles,box,tail=10,skip=10,multiply=[3,3]):
    """
    Plot a 2D-projection of the trajectory of the system.
    
    The function creates a plot showing the current and past
    positions of particles. 
    
    Args:
        particles (list): list of :class:`Planet` objects
        box (temp.PeriodicBox): supercell
        tail (int): the number of past positions to include in the plot
        skip (int): only every nth past position is plotted - skip is the number
            n, specifying how far apart the plotted positions are in time
        multiply (array): number of periodic images to draw in x and y directions
    """
    N = len(particles)
    nx = multiply[0]
    ny = multiply[1]

    plt.clf()
    ax = plt.axes()
    ax.set_aspect('equal','box')

    for part in particles:
        xcoordinates = []
        ycoordinates = []

        for pos in reversed(part.trajectory[-tail*skip::skip]):
            box.shift_inside_box(pos)
            for ix in range(nx):
                for iy in range(ny):
                    xcoordinates.append( pos[0] + ix*box.lattice[0] )
                    ycoordinates.append( pos[1] + iy*box.lattice[1] )

            plt.plot( np.array(xcoordinates) , np.array(ycoordinates) , 'o')

    plt.show()
    


def expand_supercell(particles, box, multiplier):
    """
    Expands a periodic system.
    
    The periodic system is represented by the particles and the box.
    This method creates a new, similar system, which is larger by the
    factor 'multiplier' in x, y, and z directions. That is, the list
    'particles' will be expanded by a factor of multiplier^3.
    
    For example, if particles contains Atoms at positions
    [0,0,0] and [1,1,1], and box is a cube with edge length 2,
    calling this function with multiplier = 2 will change particles
    to contain the Atoms at positions
    [0,0,0], [1,1,1], 
    [2,0,0], [3,1,1],
    [0,2,0], [1,3,1], 
    [0,0,2], [1,1,3], 
    [2,2,0], [3,3,1], 
    [2,0,2], [3,1,3], 
    [0,2,2], [1,3,3], 
    [2,2,2], [3,3,3] 
    and box will be expanded to a square with edge lengths 4.
    """
    n_particles = len(particles)
    for i in range(n_particles):
        for ix in range(multiplier):
            for iy in range(multiplier):
                for iz in range(multiplier):
                    if ix > 0 or iy > 0 or iz > 0:

                        position = particles[i].position + \
                            np.array( [ix * box.lattice[0], \
                            iy * box.lattice[1], \
                            iz * box.lattice[2] ])
                            
                        velocity = particles[i].velocity
                        mass = particles[i].mass
                        particles.append( Atom(position, velocity, mass) )
            
    box.lattice *= multiplier
    

class PeriodicBox:
    """
    Class representing a simulation box with periodic boundaries.
    
    The box is orthogonal, i.e., a rectangular volume. As such,
    it is specified by the lengths of its edges (lattice constants).
    
    Args:
        lattice (array): lattice constants
    """

    def __init__(self, lattice):
        self.lattice = lattice
        
        
    def shift_inside_box(self, position):
        """
        If the given position (3-vector) is outside the
        box, it is shifted by multiple of lattice vectors until
        the new position is inside the box. That is, the function
        transforms the position vector to an equivalen position 
        inside the box.
        
        Args:
            position (array): the position to be shifted
        """
        
        # go over x, y and z coordinates
        for i in range(3):
        
            while position[i] < 0:
                position[i] += self.lattice[i]
                
            while position[i] > self.lattice[i]:
                position[i] -= self.lattice[i]



    def distance_squared(self, particle1, particle2):
        """
        Calculates and returns the square of the 
        distance between two particles.
        
        .. math ::
            r^2_{ij} = (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2.
        
        In a periodic system, each particle has an infinite number of
        periodic copies. Therefore the distance between two particles is
        not unique. The function returns the shortest such distance,
        that is, the distance between the the periodic copies which are
        closest ot each other.
        
        Args:
            particle1 (Molecule): the first body
            particle2 (Molecule): the second body
    
        Returns:
            float: the squared distance :math:`r^2_{ij}`
        """
        pass
        
        v = self.vector(particle1,particle2)
        return v @ v
        #############################
        # implement the calculation #
        #############################



    def vector(self, particle1, particle2):
        """
        Returns the vector pointing from the position of
        particle1 to the position of particle2.
        
        .. math ::
            \\vec{r}_{i \\to j} = \\vec{r}_j - \\vec{r}_i
        
        In a periodic system, each particle has an infinite number of
        periodic copies. Therefore the displacement between two particles is
        not unique. The function returns the shortest such displacement
        vector.
        
        Args:
            particle1 (Molecule): the first body
            particle2 (Molecule): the second body
            
        Returns:
            array: components of :math:`\\vec{r}_{i \\to j}`, :math:`[x_{i \\to j}, y_{i \\to j}, z_{i \\to j}]`
        """
        
        vector_1to2 = particle2.position - particle1.position

        # loop over x, y, and z directions
        for i in range(3):
        
            # If the absolute value of the separation 
            # in this direction is larger than half the length
            # of the simulation box, there must be another
            # periodic image of particle2 closer to particle 1.
            #
            # We can find it by translating the coordinates by
            # an integer multiple of lattice vectors.
            #
            # Note: there is a more efficient way to calculate this.
            
            while vector_1to2[i] < -self.lattice[i]/2:
                vector_1to2[i] += self.lattice[i]
                
            while vector_1to2[i] > self.lattice[i]/2:
                vector_1to2[i] -= self.lattice[i]      

        return vector_1to2



class Atom:
    """
    A point like object.
    
    An atom has a position (a 3-vector), a velocity (3-vector)
    and a mass (a scalar).
    
    Args:
        position (array): coordinates :math:`[x, y, z]`
        velocity (array): velocity components :math:`[v_x, v_y, v_z]`
        mass (float): mass :math:`m`
    """

    def __init__(self, position, velocity, mass = 1.0):
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.mass = mass
        self.trajectory = []
        
        
    def move(self, force, dt):
        """
        Move the atom.
        
        Args:
            shift (array): coordinate change :math:`[\Delta x, \Delta y, \Delta z]`
        """
        self.position += self.velocity * dt + 0.5 * force/self.mass * dt*dt        
        
    def accelerate(self, force, dt):
        """
        Set a new velocity for the particle as
        
        .. math::

            \\vec{v}(t+\\Delta t) = \\vec{v}(t) + \\frac{1}{2m}\\vec{F} \Delta t
           
        Args:
            force (array): force acting on the planet :math:`[F_x, F_y, F_z]`
            dt (float): time step :math:`\\Delta t`
        """
        self.velocity += force * dt/self.mass
        
    def save_position(self):
        """
        Save the current position of the particle
        in the list 'trajectory'.
        
        Note: in a real large-scale simulation one would
        never save trajectories in memory. Instead, these
        would be written to a file for later analysis.
        """
        self.trajectory.append( [ self.position[0], self.position[1], self.position[2], ] )
        

def read_particles_from_file(filename):
    """
    Reads the properties of planets from a file.
    
    Each line should define a single :class:`Planet` listing its
    position in cartesian coordinates, velocity components and mass,
    separated by whitespace:
    
    .. code-block::
    
        x0 y0 z0 vx0 vy0 vz0 m0
        x1 y1 z1 vx1 vy1 vz1 m1
        x2 y2 z2 vx2 vy2 vz2 m2
        x3 y3 z3 vx3 vy3 vz3 m3
        ...
    
    Args:
        filename (str): name of the file to read
        
    Returns:
        list: list of :class:`temp.Atom` objects   
    """
    f = open(filename)
    lines = f.readlines()
    f.close()
    
    particles = []
    lattice = []
    
    parts = lines[0].split()
    for i in range(3):
        lattice.append(float(parts[i]))
        
    box = PeriodicBox(np.array(lattice))
    
    for line in lines[1:]:
        parts = line.split()
        if len(parts) > 0:
            position = np.array( [0.0]*3 )
            velocity = np.array( [0.0]*3 )
            
            position[0] = float(parts[0])
            position[1] = float(parts[1])
            position[2] = float(parts[2])

            velocity[0] = float(parts[3])
            velocity[1] = float(parts[4])
            velocity[2] = float(parts[5])
        
            mass = float(parts[6])
        
            particles.append( Atom( position, velocity, mass ) )

    return particles, box




def write_particles_to_file(particles, box, filename):
    """
    Write the configuration of the system in a file.
    
    The format is the same as that specified in :meth:`read_particles_from_file`.
    
    Args:
        particles (list): list of :class:`temp.Atom` objects
        box (temp.PeriodicBox): supercell
        filename (str): name of the file to write
    """
    writelines = "" 
    writelines += str(box.lattice[0])+" "+str(box.lattice[1])+" "+str(box.lattice[2])+"\n"
    for part in particles:
        writelines += str(part.position[0])+" "+str(part.position[1])+" "+str(part.position[2])+" "
        writelines += str(part.velocity[0])+" "+str(part.velocity[1])+" "+str(part.velocity[2])+" "
        writelines += str(part.mass)+"\n"
    f = open(filename,'w')
    f.write(writelines)
    f.close()    


def write_xyz_file(particles, filename):
    """
    Write the configuration of the system in a file.
    
    The information is written in so called xyz format, which many
    programs can parse.
    
    Args:
        particles (list): list of :class:`Planet` objects
        filename (str): name of the file to write
    """
    writelines = ""
    N = len(particles)
    steps = len(particles[0].trajectory) 
    
    for i in range(steps):
        writelines += str(N)+"\n \n"

        for part in particles:
            writelines += "H "+str(part.trajectory[i][0])+" "+str(part.trajectory[i][1])+" "+str(part.trajectory[i][2])+"\n"


    f = open(filename,'w')
    f.write(writelines)
    f.close() 



    
def calculate_forces(particles, box, sigma=1.0, epsilon=0.1, cutoff=1.5):
    """
    Calculates the total force applied on each atom.
    
    The forces are returned as a numpy array where
    each row contains the force on an atom
    and the columns contain the x, y and z components.
    
    .. code-block::
    
        [ [ fx0, fy0, fz0 ],
          [ fx1, fy1, fz1 ],
          [ fx2, fy2, fz2 ],
          ...               ]
          
    The function also calculates the virial, 
    
    .. math::
    
        \\sum_{i < j} U'(r_{ij}) r_{ij},
        
    which is needed for pressure calculation.
    
    Args:
        atoms (list): a list of :class:`temp.Atom` objects
        box (temp.PeriodicBox): supercell
        sigma (float): Lennard-Jones parameter :math:`\\sigma`
        epsilon (float): Lennard-Jones parameter :math:`\\epsilon`
        cutoff (float): maximum distance for interactions
        
    Returns:
        array, float: forces, virial
    """

    sigma_six = sigma*sigma*sigma*sigma*sigma*sigma
    
    virial = 0.0            
    forces = np.array( [[0.0, 0.0, 0.0]]*len(particles) )
    for i in range(1,len(particles)):
        for j in range(0,i):
        
            part_i = particles[i]
            part_j = particles[j]
        
            dist_sq = box.distance_squared(part_i, part_j)
            if dist_sq < cutoff*cutoff:
                
                dist_six = dist_sq * dist_sq * dist_sq
                dist_12 = dist_six * dist_six
            
                force_j_to_i = - 4.0 * epsilon * (12.0 * sigma_six*sigma_six / dist_12 - \
                        6.0 * sigma_six / dist_six ) * \
                        box.vector(part_i, part_j) / dist_sq
                        
                
                virial += - 4.0 * epsilon * (12.0 * sigma_six*sigma_six / dist_12 - \
                    6.0 * sigma_six / dist_six ) * dist_sq
                    
                forces[i, :] += force_j_to_i[:]
                forces[j, :] -= force_j_to_i[:]
                
    return forces, virial
            

def update_positions(particles, forces, dt):
    """
    Update the positions of all particles 
    using :meth:`temp.Atom.move` according to
    
    .. math::

        \\vec{r}(t+\\Delta t) = \\vec{r}(t) + \\vec{v} \Delta t + \\frac{1}{2m}\\vec{F} (\\Delta t)^2
             
    Args:
        particles (list): a list of :class:`Planet` objects
        force (array): array of forces on all bodies
        dt (float): time step :math:`\\Delta t`
    """
    for i in range(len(particles)):
        part = particles[i]
        force = forces[i,:]
        part.move(force, dt)

def update_positions_no_force(particles, dt):
    """
    Update the positions of all particles 
    using :meth:`temp.Atom.move` according to
    
    .. math::

        \\vec{r}(t+\\Delta t) = \\vec{r}(t) + \\vec{v} \Delta t
             
    Args:
        particles (list): a list of :class:`Planet` objects
        dt (float): time step :math:`\\Delta t`
    """
    for i in range(len(particles)):
        part = particles[i]
        part.move(0.0, dt)
    

def update_velocities(particles, forces, dt): 
    """
    Update the positions of all particles 
    using :meth:`temp.Atom.accelerate` according to
    
    .. math::

        \\vec{v}(t+\\Delta t) = \\vec{v}(t) + \\frac{1}{m}\\vec{F} \Delta t 
             
    Args:
        particles (list): a list of :class:`Planet` objects
        force (array): array of forces on all bodies
        dt (float): time step :math:`\\Delta t`
    """ 
    for i in range(len(particles)):
        part = particles[i]
        force = forces[i,:]
        part.accelerate(force, dt)



def velocity_verlet(particles, box, dt, time, trajectory_dt = 1.0, temperature = 0, tau = 0):
    """
    Verlet algorithm for integrating the equations of motion,
    i.e., advancing time.
    
    There are a few ways to implement Verlet. The leapfrog
    version works as follows: First, forces are calculated
    for all particles and velocities are updated by half a
    time step, 
    :math:`\\vec{v}(t+\\frac{1}{2}\\Delta t) = \\vec{v}(t) + \\frac{1}{2m}\\vec{F} \Delta t`.
    Then, these steps are repeated:
    
        * Positions are updated by a full time step using velocities but not forces,
            .. math ::
                \\vec{r}(t+\\Delta t) = \\vec{r}(t) + \\vec{v}(t+\\frac{1}{2}\\Delta t) \Delta t.
        * Forces are calculated at the new positions, :math:`\\vec{F}(t + \\Delta t)`.
        * Velocities are updated by a full time step using the forces
            .. math ::
                \\vec{v}(t+\\frac{3}{2}\\Delta t) = \\vec{v}(t+\\frac{1}{2}\\Delta t) + \\frac{1}{m}\\vec{F}(t+\\Delta t) \Delta t
    
    
    These operations are done using the methods
    :meth:`calculate_forces`,
    :meth:`update_velocities` and
    :meth:`update_positions_no_force`.
    
    Because velocities were updated by half a time step in the beginning of the
    simulation, positions and velocities are always offset by half a timestep.
    You always use the one that has advanced further to update the other and this
    results in a stable algorithm.
    
    Args:
        particles (list): a list of :class:`Planet` objects
        box (temp.PeriodicBox): supercell
        dt (float): time step :math:`\\Delta t`
        time (float): the total system time to be simulated
        trajectory_dt (float): the positions of particles are
            saved at these these time intervals - does not affect
            the dynamics in any way
    """
    steps = int(time/dt)
    trajectory_wait = int(trajectory_dt / dt)
    forces, virial = calculate_forces(particles, box)
    
    # get velocities at half time-step
    update_velocities(particles, forces, 0.5*dt)
    
    average_t = 0
    average_p = 0
    
    for i in range(steps):
        update_positions_no_force(particles, dt) 
        forces, virial = calculate_forces(particles, box)            
        update_velocities(particles, forces, dt)
        
        if tau > 0:
            # apply Berendsen thermostat velocity scaling
            berendsen_thermostat(particles, dt, tau, temperature) 
        
        # calculate temperature
        t = calculate_temperature(particles)
        average_t += t
                  
        # calculate the pressure
        volume = box.lattice[0]*box.lattice[1]*box.lattice[2]
        average_p += calculate_pressure(particles, box, virial, t) 
                
        if i%trajectory_wait == 0:
            for part in particles:
                part.save_position()
            
    # to calculate proper kinetic energy, return to proper timestep
    update_velocities(particles, forces, -0.5*dt)
    
    return average_t/steps, average_p/steps
        
    
def calculate_momentum(particles):
    """
    Calculates the total momentum of the system.
    
    .. math ::
        \\vec{p}_\\text{total} = \\sum_i \\vec{p}_i =  \\sum_i m_i \\vec{v}_i
    
    Args:
        particles (list): a list of :class:`Planet` objects
        
    Returns:
        array: momentum components :math:`[p_x, p_y, p_z]`
    """
    p = np.array( [0.0, 0.0, 0.0] )
    for part in particles:
        p += part.mass * part.velocity
        
    return p
    
def calculate_kinetic_energy(particles):
    """
    Calculates the total kinetic energy of the system.
    
    .. math ::
        K_\\text{total} = \\sum_i \\frac{1}{2} m_i v_i^2.
    
    Args:
        particles (list): a list of :class:`temp.Atom` objects
        
    Returns:
        float: kinetic energy :math:`K`
    """
    k = 0.0
    for part in particles:
        vx = part.velocity[0]
        vy = part.velocity[1]
        vz = part.velocity[2]

        k += 0.5 * part.mass * (vx*vx + vy*vy + vz*vz)
    
    return k
    
def calculate_potential_energy(particles, box, sigma=1.0, epsilon=0.1, cutoff=1.5):
    """
    Calculates the total potential energy of the system.
    
    The potential energy is calculated using the Lennard-Jones model
    
    .. math ::
        U = \\sum_{i \\ne j} 4 \\epsilon \\left( \\frac{ \\sigma^{12} }{ r^{12}_{ij} } 
        - \\frac{ \\sigma^6 }{ r^6_{ij} } \\right).
    
    Args:
        particles (list): a list of :class:`temp.Atom` objects
        box (temp.PeriodicBox): supercell
        sigma (float): Lennard-Jones parameter :math:`\\sigma`
        epsilon (float): Lennard-Jones parameter :math:`\\epsilon`
        cutoff (float): maximum distance for interactions
        
    Returns:
        float: potential energy :math:`U`
    """

    sigma_six = sigma*sigma*sigma*sigma*sigma*sigma
    
    u = 0.0
    for i in range(len(particles)):
        for j in range(0,i):
        
            part_i = particles[i]
            part_j = particles[j]
        
            dist_sq = box.distance_squared(part_i, part_j)
            
            # if we are beyond the cutoff, make U constant
            if dist_sq > cutoff*cutoff:
                dist_sq = cutoff*cutoff
                          
            dist_six = dist_sq * dist_sq * dist_sq
            dist_12 = dist_six * dist_six
            u += 4.0 * epsilon * \
                ( sigma_six*sigma_six/dist_12 - sigma_six/dist_six)                
                
    return u
    
    


def scale_velocities(particles, scale_factor):
    """
    Scale the velocities of all particles by the scaling factor.
    
    Args:
        particles (list): a list of :class:`temp.Atom` objects
        scale_factor (float): scaling factor
    """
    for part in particles:
        part.velocity *= scale_factor


def berendsen_thermostat(particles, dt, tau = 5.0, t0 = 0.1):
    """
    Implements the velocity scaling of the Berendsen thermostat.

    A thermostat is an algorithm which couples the simulated system
    to an external, fictious heat bath at some constant temperature :math:`T_0`.
    If the system is hotter than this, the thermostat removes energy
    from the system. And vice versa, if the system is cooler than the
    heat bath, energy is brought in the system.

    The Berendsen thermostat aims at scaling the temperature :math:`T` 
    of the system according to
    
    .. math::
    
        \\frac{d T}{d t} = \\frac{1}{\\tau}(T_0 - T),
        
    where :math:`t` is time and :math:`\\tau` is a time constant.
    This makes :math:`T` approach :math:`T_0` exponentially.
    
    In practice, the temperature is changed by scaling all velocities
    with a scaling factor :math:`\\lambda` using :meth:`scale_velocities`.
    In time step :math:`\\Delta t` one expects approximately
    
    .. math::
    
        \\Delta T = \\frac{\\Delta t}{\\tau}(T_0 - T)
        
    and solving for :math:`\\lambda` from the definition of kinetic 
    temperature yields
    
    .. math::
    
        \\lambda = \\sqrt{ 1 + \\frac{\\Delta t}{\\tau} \\left[ \\frac{T_0}{T} - 1 \\right] }.

    .. note ::
        This function is incomplete!

    Args:
        dt (float): timestep :math:`\\Delta t`
        tau (float): time constant :math:`\\tau`
        t0 (float): external temperature :math:`T_0`
    """

    pass
    
    # todo
    

    

def calculate_temperature(particles):
    """
    Calculate and return the current temperature.
    """
    
    t = 2.0/(3.0*len(particles))*calculate_kinetic_energy(particles)
    
    return t
        
    

    
def assign_maxwell_boltzmann_velocities(particles, temperature):
    """
    Randomly pick velocities for all particles. The velocities are
    chosen according to the Maxwell-Boltzmann distribution.
    """

    random = default_rng()
    total_mass = 0.0
    total_momentum = np.array( [0.0, 0.0, 0.0] )
    for part in particles:
        m = part.mass
        total_mass += m
        for i in range(3):
            v = random.standard_normal()*sqrt(temperature/m)
            part.velocity[i] = v
            total_momentum[i] += m*v
            
    deltav = -total_momentum / total_mass
    for part in particles:
        part.velocity += deltav
    

    
def calculate_average_and_error(values, start=0):
    """
    Calculates the average and standard error of mean of a sequence.
    
    The values in the sequence are assumed to be uncorrelated.
    
    If the beginning of the sequence cannot be used in the analysis (equilibrium
    has not yet been reached), one can ignore the early values by specifying a
    starting index.
    
    Args:
        values (array): values to analyse
        start (int): index of the first value to be included in the analysis
    """

    avr_x = 0.0
    avr_sq = 0.0
    for x in values[start:]:
        avr_x += x
        avr_sq += x*x
        
    n = float(len(values)-start)
    avr_x /= n
    avr_sq /= n
    variance = (avr_sq - avr_x*avr_x)*n/(n-1)
    error = sqrt(variance/n)
    
    return avr_x, error
    
    
def calculate_pressure(particles, box, virial, temperature):
    """
    Calculate the current pressure.
    
    For a molecular simulation with constant pressure, volume and temperature, 
    one can derive the relation
    
    .. math::
    
       pV = Nk_B T + \\frac{1}{d} \\langle \\sum_i \\vec{F}_i \\cdot \\vec{r}_i \\rangle,
       
    where :math:`p, V, N, k_B, T, d, \\vec{F}_i, \\vec{r}_i` are, respectively,
    pressure, volume, number of atoms, Boltzmann constant, temperature,
    number of dimensions, force acting on atom :math:`i` and position of atom :math:`i`.
    
    This function uses this relation to solve the effective instantaneous pressure as
    
    .. math ::

       p = \\frac{1}{V} Nk_B T + \\frac{1}{dV} \\sum_i \\vec{F}_i \\cdot \\vec{r}_i,

    where the sum is called the virial.
    
    This is not necessarily the true instantaneous pressure, but calculating
    the average of this quantity over an extended simulation should converge
    towards the true pressure.
    
    Args:
        particles (list): list of :class:`temp.Atom` objects
        box (temp.PeriodicBox): supercell
        virial (float): virial
        temperature (float): temperature
    
    Returns: 
        float: pressure
    """
    
    n = len(particles)
    volume = box.lattice[0]*box.lattice[1]*box.lattice[2]
    return (n*temperature - virial/3) / volume   
    
    
def main(filename, external_temperature, tau, dt, sample_interval, simulation_time, thermalization_time):
    """
    The main program.
    
    The program reads the system from a file, runs the simulation.
    
    Atomic velocities are initialized according to the given
    temperature. In addition, if a time constant is given,
    a thermostat is applied to drive temperature towards
    this temperature.
    
    Args:
        filename (str): name of the file to read
        external_temperature (float): temperature
        tau (float): thermostat time constant
    """

    # Create the system.
    particles, box = read_particles_from_file(filename)
    # Duplicate the system in all directions to make a larger system.
    expand_supercell(particles, box, 2)
    # Initialize the velocities of particles.
    assign_maxwell_boltzmann_velocities(particles, external_temperature)

    # Print initial properties.
    print( "starting simulation" )
    print( "number of particles ", len(particles) )
    
    #t = calculate_temperature(particles)
    #forces, virial = calculate_forces(particles, box)
    #p = calculate_pressure(particles, box, virial, t) 
    #k = calculate_kinetic_energy(particles)
    u0 = calculate_potential_energy(particles, box)

    # Set up lists to gather data during simulation.
    temperatures = []
    pressures = []
    kinetic_energies = []
    potential_energies = [] 
    times = []
    
    n_thermalization = int(thermalization_time / sample_interval)
    n_samples = int(simulation_time / sample_interval)

    # Run the simulation.
    #
    # The run is split in short stints, and properties
    # are recorded between stints.
    for i in range(n_samples):
    
        # Run a short simulation.
        t, p = velocity_verlet(particles, box, dt, sample_interval, 1.0, external_temperature, tau)
        
        print_progress(i+1, n_samples)
        
        # Evaluate properties.
        k = calculate_kinetic_energy(particles)
        u = calculate_potential_energy(particles, box)

        # Save statistics
        times.append( (i+0.5)*sample_interval )
        temperatures.append( t )
        pressures.append( p )
        kinetic_energies.append( k )
        potential_energies.append( u-u0 ) # measure U with respect to starting configuration
        
        
    t_avr, t_err = calculate_average_and_error(temperatures, n_thermalization)
    p_avr, p_err = calculate_average_and_error(pressures, n_thermalization)
        
    # error of mean gives < 70 % confidence interval
    # 2*error is about 95 % confidence interval
    print("temperature: "+str(t_avr)+" +- "+str(2*t_err))
    print("pressure:    "+str(p_avr)+" +- "+str(2*p_err))
        
    # Plot the temperature 
    plt.plot(times, t_avr*np.ones(len(times)), 'b:')
    plt.fill_between(times, (t_avr-2*t_err)*np.ones(len(times)), (t_avr+2*t_err)*np.ones(len(times)), color = 'b', alpha=0.1 )
    plt.plot(times, temperatures, 'r')
    plt.xlabel("time")
    plt.ylabel("temperature")
    plt.show()
    
    # Plot the pressure    
    plt.plot(times, p_avr*np.ones(len(times)), 'b:')
    plt.fill_between(times, (p_avr-2*p_err)*np.ones(len(times)), (p_avr+2*p_err)*np.ones(len(times)), color = 'b', alpha=0.1 )
    plt.plot(times, pressures, 'r')
    plt.xlabel("time")
    plt.ylabel("pressure")
    plt.show()
    
    # Plot the kinetic, potential and total energy
    kinetic_energies = np.array( kinetic_energies )
    potential_energies = np.array( potential_energies )    
    plt.plot(times, kinetic_energies, label="kinetic energy")
    plt.plot(times, potential_energies, label="potential energy")
    plt.plot(times, potential_energies+kinetic_energies, label="total energy")
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("energy")
    plt.show()

    # Plot the system
    animate(particles, box)


if __name__ == "__main__":

    filename = 'fcc.txt'
    T = 0.1
    tau = 5
    dt = 0.02
    sample_time = 1.0
    simulation_time = 100.0
    thermalization_time = 50.0
    main(filename, T, tau, dt, sample_time, simulation_time, thermalization_time)
