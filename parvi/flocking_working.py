import sys
import copy
from math import *
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
import matplotlib.animation as ani

#
# input file keywords
#
# note: global variables like these are
# often written in ALL CAPS
#
COHESION_TAG = "cohesion"
SEPARATION_TAG = "separation"
ALIGNMENT_TAG = "alignment"
AVOIDANCE_TAG = "avoidance"
CHASE_TAG = "chase"
MAX_TAG = "max"
RADIUS_TAG = "radius"
STRENGTH_TAG = "strength"

PREY_TAG = "prey"
PREDATOR_TAG = "predator"
SKY_TAG = "sky"
TIME_TAG = "time"
N_TAG = "count"
SPEED_TAG = "speed"
D_TAG = "dimension"
SIZE_TAG = "size"
TOTAL_TAG = "total"
DT_TAG = "dt"
RDT_TAG = "recording-dt"

STRENGTH_INDEX = 0
MAX_INDEX = 1
RADIUS_INDEX = 2
N_INDEX = 0
SPEED_INDEX = 1
D_INDEX = 0
SIZE_INDEX = 1
DT_INDEX = 0
TOTAL_INDEX = 1
RDT_INDEX = 2


def find_info(lines, tag):
    """
    Searches for the information wrapped in the given tag
    among the given lines of text.
    
    If tag is, e.g., "foo", the function searches for the start tag
    <foo> and the end tag </foo> and returns the lines of information
    between them.
    
    The function only finds the first instance of the given tag.
    However, in order to catch multiple instances of the tag, the
    function also returns all the information in lines following
    the end tag.
    
    For instance, if lines contains the strings:

    .. code-block ::
    
        aa
    
        <foo>
        bb
        cc
        </foo>
    
        dd
        ee
        ff
    
    the function will return two lists: ["bb", "cc"], ["", "dd", "ee", "ff"].
    
    Args:
        lines (list): the information as a list of strings
        tag (str): the tag to search
    
    Returns: 
        list, list: the lines between start and end tags, the lines following the end tag
    """
    info = []
    is_relevant = False
    line_number = 0
        
    # go through the data
    for i in range(len(lines)):
        line = lines[i]
        
        if is_relevant: # if we have found the starting tag, record information 
            info.append(line)
            
        contents = line.strip() # remove whitespace at the start and end of the line
        
        if len(contents) > 0: # skip empty lines
        
            if contents[0] == "<" and contents[-1] == ">": # is this a tag?
            
                if contents[1:-1] == tag: # found the starting tag

                    if not is_relevant: # we had not yet found the tag
                        is_relevant = True # the following lines are relevant
                        line_number = i
                        
                    else: # we had already started this tag
                        print("Found tag <"+tag+"> while already reading <"+tag+">")
                        raise Exception("parsing error")
                        
                if contents[1:-1] == "/"+tag: # found the end tag
                    return info, lines[i+1:]
    
        
    # we end up here, if we reach the end of the file
    
    if is_relevant: # the file ends while reading info (start tag was found, but no end tag)
        print("Reached the end of file while parsing <"+tag+"> from line "+str(line_number+1))
        raise Exception("parsing error")
        
    elif info == []: # the tag was not found
        print("Tag <"+tag+"> was not found")
        return [], lines
        

        
def parse_line(line):
    """
    Separates tag and info on a line of text.
    
    The function also removes extra whitespace and comments separated with #.
    
    For instance if line is " x :  1.23  # the x coordinate",
    the function returns ("x", "1.23").
    
    Args:
        line (str): a string of information
    
    Returns: 
        str, str: tag, info
    """

    parts = line.split(":")
    tag = ""
    info = ""
    
    if len(parts) > 1:
        tag = parts[0].strip()
        info = parts[1].split("#")[0].strip()
        
    return tag, info
    

def read_force_info(lines, default=1.0):
    """
    Reads parameters from given lines.
    
    Args:
        lines (list): information as a list of strings
        default (float): the default lattice parameter in all directions
    
    Returns: 
        list: parameters
    """
    parameters = [default]*3
    
    for line in lines:
        tag, info = parse_line(line)
        if tag == STRENGTH_TAG:
            parameters[STRENGTH_INDEX] = float(info)
        elif tag == RADIUS_TAG:
            parameters[RADIUS_INDEX] = float(info)
        elif tag == MAX_TAG:
            parameters[MAX_INDEX] = float(info)

    return parameters
    
    
def read_animal_info(lines, default=1):
    """
    Reads parameters from given lines.
    
    Args:
        lines (list): information as a list of strings
        default (float): the default lattice parameter in all directions
    
    Returns: 
        list: parameters
    """
    parameters = [default]*2
    
    for line in lines:
        tag, info = parse_line(line)
        if tag == N_TAG:
            parameters[N_INDEX] = float(info)
        elif tag == SPEED_TAG:
            parameters[SPEED_INDEX] = float(info)

    return parameters
    
    
def read_box_info(lines, default=2):
    """
    Reads parameters from given lines.
    
    Args:
        lines (list): information as a list of strings
        default (float): the default lattice parameter in all directions
    
    Returns: 
        list: parameters
    """
    parameters = [default]*2
    
    for line in lines:
        tag, info = parse_line(line)
        if tag == D_TAG:
            parameters[D_INDEX] = float(info)
        elif tag == SIZE_TAG:
            parameters[SIZE_INDEX] = float(info)

    return parameters
    
    
def read_time_info(lines, default=0.1):
    """
    Reads parameters from given lines.
    
    Args:
        lines (list): information as a list of strings
        default (float): the default lattice parameter in all directions
    
    Returns: 
        list: parameters
    """
    parameters = [default]*3
    
    for line in lines:
        tag, info = parse_line(line)
        if tag == DT_TAG:
            parameters[DT_INDEX] = float(info)
        elif tag == TOTAL_TAG:
            parameters[TOTAL_INDEX] = float(info)
        elif tag == RDT_TAG:
            parameters[RDT_INDEX] = float(info)

    return parameters
    
    
def read_dynamics_file(filename):
    """
    Reads the given file for force parameters.
    
    Args:
        filename (str): the file to read
        
    Returns:
        lists: cohesion, alignment, separation, avoidance, chase parameters
    """

    f = open(filename)
    lines = f.readlines()
    f.close()
        
    info, dummy = find_info( lines, COHESION_TAG )
    cohesion_params = read_force_info( info )
    
    info, dummy = find_info( lines, SEPARATION_TAG )
    separation_params = read_force_info( info )
    
    info, dummy = find_info( lines, ALIGNMENT_TAG )
    alignment_params = read_force_info( info )

    info, dummy = find_info( lines, AVOIDANCE_TAG )
    avoidance_params = read_force_info( info )
    
    info, dummy = find_info( lines, CHASE_TAG )
    chase_params = read_force_info( info )

    return cohesion_params, alignment_params, separation_params, avoidance_params, chase_params


def read_system_file(filename):
    """
    Reads the given file for system parameters.
    
    Args:
        filename (str): the file to read
        
    Returns:
        lists: prey, predator, box, time parameters
    """

    f = open(filename)
    lines = f.readlines()
    f.close()
        
    info, dummy = find_info( lines, PREY_TAG )
    prey_params = read_animal_info( info )
    
    info, dummy = find_info( lines, PREDATOR_TAG )
    predator_params = read_animal_info( info )
    
    info, dummy = find_info( lines, SKY_TAG )
    box_params = read_box_info( info )
    
    info, dummy = find_info( lines, TIME_TAG )
    time_params = read_time_info( info )

    return prey_params, predator_params, box_params, time_params



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
    Draws a representation of the system as a scatter plot. 
    
    Used for animation.
    
    Args:
        frame (int): index of the frame to be drawn
        xtraj (array): x-coordinates of all animals at different animation frames
        ytraj (array): y-coordinates at all animals at different animation frames
        ztraj (array): z-coordinates at all animals at different animation frames
        bounds (array): list of lower and upper bounds for the plot as [[xmin, xmax], [ymin, ymax]]
    """

    plt.clf()
    ax = plt.axes()
    ax.set_xlim(bounds[0])
    ax.set_ylim(bounds[1])
    ax.set_aspect('equal')
    size = (ztraj[frame]+1)*2
    plt.scatter(xtraj[frame], ytraj[frame], marker='o', s=size )
        


def animate( prey, predators, box, multiply = [3,3] ):
    """
    Animates the simulation.
    
    Args:
        prey (list): list of :class:`Prey` objects
        predators (list): list of :class:`Predator` objects
        box (flocking.PeriodicBox): supercell
        multiply (array): number of periodic images to draw in x and y directions
    """

    nframes = len(prey[0].trajectory)    
    
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
        
        for p in prey:            
            box.shift_inside_box(p.trajectory[i])
            for ix in range(multix):
                for iy in range(multiy):
                    xtraj[-1].append(p.trajectory[i][0] + ix*Lx)
                    ytraj[-1].append(p.trajectory[i][1] + iy*Ly)
                    ztraj[-1].append(1)
                    
        
        for p in predators:            
            box.shift_inside_box(p.trajectory[i])
            for ix in range(multix):
                for iy in range(multiy):
                    xtraj[-1].append(p.trajectory[i][0] + ix*Lx)
                    ytraj[-1].append(p.trajectory[i][1] + iy*Ly)
                    ztraj[-1].append(10)
                    
                    
            
    xtraj=np.array(xtraj)
    ytraj=np.array(ytraj)   
    ztraj=np.array(ztraj)   
    
    fig = plt.figure()
    motion = ani.FuncAnimation(fig, draw, nframes, interval=10, fargs=(xtraj, ytraj, ztraj, bounds) )
    plt.show()




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



    def distance_squared(self, position1, position2):
        """
        Calculates and returns the square of the 
        distance between two points,
        
        .. math ::
            r^2_{ij} = (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2.
        
        In a periodic system, each boid has an infinite number of
        periodic copies. Therefore the distance between two points is
        not unique. The function returns the shortest such distance,
        that is, the distance between the the periodic copies which are
        closest ot each other.
        
        Args:
            position1 (array): the first point
            position2 (array): the second point
    
        Returns:
            float: the squared distance :math:`r^2_{ij}`
        """
        pass
        
        v = self.vector(position1,position2)
        return v @ v



    def vector(self, position1, position2):
        """
        Returns the vector pointing from position1 to position2,
        
        .. math ::
            \\vec{r}_{i \\to j} = \\vec{r}_j - \\vec{r}_i
        
        In a periodic system, each boid has an infinite number of
        periodic copies. Therefore the displacement between two points is
        not unique. The function returns the shortest such displacement
        vector.
        
        Args:
            position1 (array): the first point
            position2 (array): the second point
            
        Returns:
            array: components of :math:`\\vec{r}_{i \\to j}`, :math:`[x_{i \\to j}, y_{i \\to j}, z_{i \\to j}]`
        """
        
        vector_1to2 = position2 - position1

        # loop over x, y, and z directions
        for i in range(3):
        
            # If the absolute value of the separation 
            # in this direction is larger than half the length
            # of the simulation box, there must be another
            # periodic image of boid2 closer to boid 1.
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



class Animal:
    """
    A class representing an animal.
    
    All animals are treated as point-like agents in this model.
    
    Args:
        position (array): coordinates :math:`[x, y, z]`
        velocity (array): velocity components :math:`[v_x, v_y, v_z]`
        top_speed (fload): the maximum speed allowed for this animal
        mass (float): mass :math:`m`    
    """

    def __init__(self, position, velocity, top_speed, mass = 1.0):
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.mass = mass
        self.force = np.zeros(3)
        self.trajectory = []
        self.neighbors = []
        self.top_speed = top_speed


    def move(self, dt):
        """
        Move the animal.
        
        Args:
            dt (float): time step :math:`\\Delta t`
        """
        self.position += self.velocity * dt        
        
        
    def accelerate(self, dt):
        """
        Set a new velocity as
        
        .. math::

            \\vec{v}(t+\\Delta t) = \\vec{v}(t) + \\frac{1}{2m}\\vec{F} \Delta t
            
        If the speed would exceed the maximum speed of the animal, it is
        scaled down to the max value. Similarly, the animal must move at least
        with a speed that is half of the maximum speed.
           
        Args:
            dt (float): time step :math:`\\Delta t`
        """
        self.velocity += self.force * dt/self.mass
        v2 = self.velocity @ self.velocity
        
        if v2 > self.top_speed**2:
            self.velocity *= self.top_speed/np.sqrt(v2)
        elif v2 < self.top_speed**2/4:
            self.velocity *= self.top_speed/np.sqrt(v2)/2
            
        
    def save_position(self):
        """
        Save the current position.
        
        Note: in a real large-scale simulation one would
        never save trajectories in memory. Instead, these
        would be written to a file for later analysis.
        """
        self.trajectory.append( [ self.position[0], self.position[1], self.position[2] ] )


class Predator(Animal):
    """
    A point-like model for a predator :class:`Animal`.
    
    Likes to harrass :class:`Prey` animals.
    
    Args:
        position (array): coordinates :math:`[x, y, z]`
        velocity (array): velocity components :math:`[v_x, v_y, v_z]`
        top_speed (fload): the maximum speed allowed for this animal
        mass (float): mass :math:`m`
    """

    def __init__(self, position, velocity, prey, top_speed, mass = 1.0):
        Animal.__init__(self, position, velocity, top_speed, mass)
        
        # pick a random prey as the initial target    
        self.target = random.choice( prey )
        


    def pick_target(self, prey, box, radius):
        """
        Randomly picks one :class:`Prey` within the given radius as the new target to chase.
        
        If there are no prey within range, nothing happens.
        
        Args:        
            prey (list): list of :class:`Prey` objects
            box (flocking.PeriodicBox): supercell
            radius (float): maximum range for choosing a target
        """
        targets = []

        for boid in prey:
            dist = box.distance_squared( self.position, boid.position )
            if dist < radius*radius:
                targets.append( boid )
            
        if len(targets) > 0:
            self.target = random.choice( targets )        

        

    def add_chase_force(self, box, parameters):
        """
        Adds a force driving the :class:`Predator` towards its target :class:`Prey`.
        
        Args:
            box (flocking.PeriodicBox): supercell
            parameters (list): force parameters
        """
    
        chase_strength = parameters[STRENGTH_INDEX]
        max_force = parameters[MAX_INDEX]

        to_target = box.vector(self.position, self.target.position)
        force = chase_strength*to_target
        
        f2 = force @ force + 0.001
        if f2 > max_force*max_force:
            force *= max_force/np.sqrt(f2)
        
        self.force += force
        
    
    def add_separation_force(self, box, parameters):
        """
        Adds a force to make predators avoid collisions.
        
        This is a short ranged repulsive force applied only if the
        animals are very close to each other.
        
        Args:
            box (flocking.PeriodicBox): supercell
            parameters (list): force parameters
        """
    
        separation_strength = parameters[STRENGTH_INDEX]
        min_distance = parameters[RADIUS_INDEX]
        max_force = parameters[MAX_INDEX]
    
        force = np.zeros(3)
        for boid in self.neighbors:
            to_boid = box.vector(self.position, boid.position)            
            
            d2 = to_boid @ to_boid + 0.001
            if d2 < min_distance*min_distance:
                force -= separation_strength*to_boid/d2
        
        f2 = force @ force + 0.001
        if f2 > max_force*max_force:
            force *= max_force/np.sqrt(f2)
            
        self.force += force
        

    def calculate_forces(self, box, force_parameters):
        """
        Calculates all forces acting on this animal.
        
        The force sum is saved in Predator.force.
        
        Args:
            box (flocking.PeriodicBox): supercell
            force_parameters (list): force parameters
        """
        
        separation_parameters = force_parameters[2]
        chase_parameters = force_parameters[4]

        self.force = np.zeros(3)
        self.add_chase_force(box, chase_parameters)
        self.add_separation_force(box, separation_parameters)





class Prey(Animal):
    """
    A point-like model for a flock-forming :class:`Animal`.
    
    Can be harrassed by a :class:`Predator`.
    
    Args:
        position (array): coordinates :math:`[x, y, z]`
        velocity (array): velocity components :math:`[v_x, v_y, v_z]`
        top_speed (fload): the maximum speed allowed for this animal
        mass (float): mass :math:`m`
    """

    def __init__(self, position, velocity, top_speed, mass = 1.0):
        Animal.__init__(self, position, velocity, top_speed, mass)   

        
    def add_cohesion_force(self, box, parameters):    
        """
        Adds a force to make :class:`Prey` seek each other.
        
        Cohesion force is a force pointing to the center
        (average coordinate) of all the neighbors of the animal.
        
        Args:
            box (flocking.PeriodicBox): supercell
            parameters (list): force parameters
        """
        
        cohesion_strength = parameters[STRENGTH_INDEX]
        max_force = parameters[MAX_INDEX]
        
        to_center = np.zeros(3)
        n_neighbors = len(self.neighbors)
        if n_neighbors == 0:
            return
        
        # calculate the vector from the animal to the center of its neighbors
        for boid in self.neighbors:
            to_center += box.vector(self.position, boid.position)
            
        to_center /= n_neighbors
        
        # normalize to_center and multiply by the strength parameter
        force = cohesion_strength*to_center / np.sqrt(to_center @ to_center + 0.001)
        
        # check that the force is not too large
        f2 = force @ force + 0.001
        if f2 > max_force*max_force:
            force *= max_force/np.sqrt(f2)
        
        self.force += force
        
        
    def add_alignment_force(self, box, parameters):
        
        """
        Adds a force to make :class:`Prey` move in the same direction
        
        Alignment force is calculated by first determining the
        average velocity vector of the neighbors. The force is then
        proportional to the difference between the current velocity
        of the animal and the calculated average.
        
        Args:
            box (flocking.PeriodicBox): supercell
            parameters (list): force parameters
        """
        
        alignment_strength = parameters[STRENGTH_INDEX]
        max_distance = parameters[RADIUS_INDEX]
        max_force = parameters[MAX_INDEX]
    
        velocity = np.zeros(3)
    
        n_neighbors = len(self.neighbors)
        if n_neighbors == 0:
            return
        
        # calculate the average of the neighbor velocities
        for boid in self.neighbors:
            d2 = box.distance_squared(self.position, boid.position)
            if d2 < max_distance*max_distance:
                velocity += boid.velocity

        velocity /= n_neighbors
             
        # calculate difference between the velocity of the animal and the average
        # also multiply by the strength parameter
        to_velo = alignment_strength*(velocity-self.velocity)
        
        # normalize to_velo
        to_velo = to_velo / np.sqrt(to_velo @ to_velo + 0.001)
        
        # check that the force is not too large
        force = alignment_strength*to_velo
        f2 = force @ force + 0.001
        if f2 > max_force*max_force:
            force *= max_force/np.sqrt(f2)
                    
        self.force += force
        
        
    def add_separation_force(self, box, parameters):
        """
        Adds a force to make prey avoid collisions.
        
        This is a short ranged repulsive force applied only if the
        animals are very close to each other.
        
        Args:
            box (flocking.PeriodicBox): supercell
            parameters (list): force parameters
        """
    
        separation_strength = parameters[STRENGTH_INDEX]
        min_distance = parameters[RADIUS_INDEX]
        max_force = parameters[MAX_INDEX]
        
        force = np.zeros(3)
        for boid in self.neighbors:
            to_boid = box.vector(self.position, boid.position)            
            
            d2 = to_boid @ to_boid + 0.001
            
            # apply repulsive force for too close neighbors
            # this implementation uses 1/r force
            if d2 < min_distance*min_distance:
                force -= separation_strength*to_boid/d2
        
        # check that the force is not too large
        f2 = force @ force + 0.001
        if f2 > max_force*max_force:
            force *= max_force/np.sqrt(f2)
            
        self.force += force
        

    def add_avoidance_force(self, predators, box, parameters):        
        """
        Adds a force to make :class:`Prey` flee from :class:`Predator` animals.
        
        This force is calculated similarly to :meth:`Prey.add_separation_force`,
        but it is triggered by predators instead of other prey.
        Typically, this force should also be stronger and have a longer range
        than the separation force.
        
        Args:
            box (flocking.PeriodicBox): supercell
            parameters (list): force parameters
        """
        
        avoid_strength = parameters[STRENGTH_INDEX]
        min_distance = parameters[RADIUS_INDEX]
        max_force = parameters[MAX_INDEX]

        force = np.zeros(3)
        
        for predator in predators:
            to_predator = box.vector(self.position, predator.position)

            d2 = to_predator @ to_predator
            
            # apply repulsive force for too close predators
            # this implementation uses 1/r force
            if d2 < min_distance*min_distance:
                force -= avoid_strength*to_predator/d2
        
        f2 = force @ force + 0.001

        # check that the force is not too large
        if f2 > max_force*max_force:
            force *= max_force/np.sqrt(f2)
            
        self.force += force
        
        
    def calculate_forces(self, predators, box, force_parameters):
        """
        Calculates all forces acting on this animal.
        
        The force sum is saved in Prey.force.
        
        Args:
            predators (list): list of :class:`Predator` objects
            box (flocking.PeriodicBox): supercell
            force_parameters (list): force parameters
        """

        cohesion_parameters = force_parameters[0]
        alignment_parameters = force_parameters[1]
        separation_parameters = force_parameters[2]
        avoid_parameters = force_parameters[3]

        self.force = np.zeros(3)
        self.add_cohesion_force(box, cohesion_parameters)
        self.add_alignment_force(box, alignment_parameters)
        self.add_separation_force(box, separation_parameters)
        self.add_avoidance_force(predators, box, avoid_parameters)



def update_neighbors(animals, box, radius):
    """
    Finds for each animal the other animals of the same type within the given radius.
    
    The animals can be both :class:`Prey` or :class:`Predator`.
    
    The neighbors are saved in the animal object as a list of animal objects.
    
    Args:
        animals (list): list of animals
        box (flocking.PeriodicBox): supercell
        radius (float): maximum distance to neighbors
    """
    
    neighbors = []
    r2 = radius*radius
    n_animals = len(animals)
    
    # empty the neighbor lists
    for i in range(n_animals):
        boid_i = animals[i]
        boid_i.neighbors = []
        
    # loop over all pairs of animals and add new neighbors
    #
    # note: This is the most computationally expensive step in the whole
    # simulation. There are ways to make it faster, but this simple implementation
    # is sufficient for us.
    for i in range(n_animals):
        for j in range(i+1,n_animals):
            boid_i = animals[i]
            boid_j = animals[j]
            
            # if the animals are close to each other, add them as neighbors to both
            if box.distance_squared(boid_i.position, boid_j.position) <= r2:
                boid_i.neighbors.append(boid_j)
                boid_j.neighbors.append(boid_i)
        
            
def update_targets(prey, predators, box, radius):
    """
    Possible updates the target for every predator.
    
    For each :class:`Predator`, there is a 10 % chance it seeks a new target :class:`Prey`
    via :meth:`Predator.pick_target`.
    
    Args:
        prey (list): list of :class:`Prey` objects
        predators (list): list of :class:`Predator` objects
        box (flocking.PeriodicBox): supercell
        radius (float): maximum distance to new targets
    """
    for p in predators:
        r = random.random()
        if r < 0.1 or p.target == None:
            p.pick_target(prey, box, radius)



def calculate_forces(prey, predators, box, force_parameters):
    """
    Calculates forces on all animals.
    
    Args:
        prey (list): list of :class:`Prey` objects
        predators (list): list of :class:`Predator` objects
        box (flocking.PeriodicBox): supercell
        force_parameters (list): force parameters
    """

    radius = force_parameters[0][RADIUS_INDEX]

    update_neighbors(prey, box, radius)
    update_neighbors(predators, box, radius)
    update_targets(prey, predators, box, radius)
    
    for p in prey:
        p.calculate_forces(predators, box, force_parameters)
        
    for p in predators:
        p.calculate_forces(box, force_parameters)
        


def create_random_flock(prey_parameters, predator_parameters, box):
    """
    Creates a flock of :class:`Prey` and possible also :class:`Predator` animals.
    
    The animals are placed randomly and they are given random velocities.
        
    Args:
        prey_parameters (list): amount and top speed of prey
        predator_parameters (list): amount and top speed of predators
        box (flocking.PeriodicBox): supercell
    """
    
    n_prey = prey_parameters[N_INDEX]
    n_predator = predator_parameters[N_INDEX]
    s_prey = prey_parameters[SPEED_INDEX]
    s_predator = predator_parameters[SPEED_INDEX]
    
    # check if this is a 2D or 3D simulation
    if box.lattice[2] < 1:
        dimension = 2
    else:
        dimension = 3
    
    flock = []
    predators = []
    
    # create prey
    while len(flock) < n_prey:
        
        r = np.zeros(3)
        v = np.zeros(3)
    
        for i in range(dimension):
            r[i] = random.random()*box.lattice[i]
            v[i] = (2.5*random.random()-1)*s_prey/2
            
        # don't put the prey too close to each other
        ok = True        
        for boid in flock:
            if box.distance_squared(r, boid.position) < 0.5:
                ok = False
                break
                
        if ok:
            flock.append( Prey( r, v, s_prey ) )
            
    # create predators
    while len(predators) < n_predator:
        
        r = np.zeros(3)
        v = np.zeros(3)
    
        for i in range(dimension):
            r[i] = random.random()*box.lattice[i]
            v[i] = (2.5*random.random()-1)*s_predator/2
            
        ok = True
        for p in predators:
            if box.distance_squared(r, p.position) < 0.5:
                ok = False
                break
                
        if ok:
            predators.append( Predator( r, v, flock, s_predator ) )
                
    return flock, predators


        


def update_positions_no_force(prey, predators, dt):
    """
    Update the positions of all animals.
    
    .. math::

        \\vec{r}(t+\\Delta t) = \\vec{r}(t) + \\vec{v} \Delta t
             
    Args:
        prey (list): a list of :class:`Prey` objects
        predators (list): a list of :class:`Predator` objects
        dt (float): time step :math:`\\Delta t`
    """
    for p in prey:
        p.move(dt)

    for p in predators:
        p.move(dt)
            

def update_velocities(prey, predators, dt): 
    """
    Update the positions of all animals according to
    
    .. math::

        \\vec{v}(t+\\Delta t) = \\vec{v}(t) + \\frac{1}{m}\\vec{F} \Delta t 
             
    Args:
        prey (list): a list of :class:`Prey` objects
        predators (list): a list of :class:`Predator` objects
        force (array): array of forces on all bodies
        dt (float): time step :math:`\\Delta t`
    """ 
    for p in prey:
        p.accelerate(dt)

    for p in predators:
        p.accelerate(dt)



def velocity_verlet(prey, predators, box, force_parameters, dt, time, trajectory_dt = 1.0):
    """
    Verlet algorithm for integrating the equations of motion,
    i.e., advancing time.
    
    There are a few ways to implement Verlet. The leapfrog
    version works as follows: First, forces are calculated
    for all animals and velocities are updated by half a
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
        prey (list): a list of :class:`Prey` objects
        predators (list): a list of :class:`Predator` objects
        box (flocking.PeriodicBox): supercell
        force_parameters (list): force parameters
        dt (float): time step :math:`\\Delta t`
        time (float): the total system time to be simulated
        trajectory_dt (float): the positions of prey are
            saved at these these time intervals - does not affect
            the dynamics in any way
    """
    steps = int(time/dt)
    trajectory_wait = int(trajectory_dt / dt)
    calculate_forces(prey, predators, box, force_parameters)
    
    # get velocities at half time-step
    update_velocities(prey, predators, 0.5*dt)
    
    for i in range(steps):
        update_positions_no_force(prey, predators, dt) 
        calculate_forces(prey, predators, box, force_parameters)    
        update_velocities(prey, predators, dt) 
                
        print_progress(i+1,steps)
        if i%trajectory_wait == 0:
            for part in prey:
                part.save_position()
            for part in predators:
                part.save_position()
            
    # to calculate proper kinetic energy, return to proper timestep
    update_velocities(prey, predators, -0.5*dt)
        

    
    
    
def main(dynamic_file = "flock_dynamics.txt", system_file = "flock_system.txt"):
    """
    The main program.
    
    The program reads info from files, runs the simulation,
    and plots the trajectory.
    
    Args:
        dynamic_file (str): file with force parameters
        system_file (str): file with simulation parameters
    """

    # read force parameters
    cohesion_p, alignment_p, separation_p, avoid_p, chase_p = read_dynamics_file(dynamic_file)
    force_p = (cohesion_p, alignment_p, separation_p, avoid_p, chase_p)
    
    # read simulation parameters
    prey_p, predator_p, sky_p, time_p = read_system_file(system_file)
    
    # create simulation box
    dimension = sky_p[D_INDEX]
    L = sky_p[SIZE_INDEX]
    if dimension == 2:
        box = PeriodicBox( np.array([L, L, 0.1]) )
    else:
        box = PeriodicBox( np.array([L, L, L]) )
    
    # create a random set of animals
    birds, predators = create_random_flock(prey_p, predator_p, box)
    
    dt = time_p[DT_INDEX]
    simulation_time = time_p[TOTAL_INDEX]
    recording_dt = time_p[RDT_INDEX]
    
    # run simulaiton
    velocity_verlet(birds, predators, box, force_p, dt, simulation_time, recording_dt)

    # animate simulation
    animate( birds, predators, box, [1,1] )
        
    
if __name__ == "__main__":
    random = random.default_rng()
    main()
