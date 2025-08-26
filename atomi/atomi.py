#
# Piirtää aaltofunktioita ja todennäköisyystiheyksiä
# elektronille vetyatomin ominaistiloilla.
#
# Teemu Hynninen 2023
#

import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
import matplotlib.animation as ani

import parametrit as par

n = par.N
l = par.L
m = par.M

if l >= n:
    print("Kvanttiluvun l on oltava pienempi kuin n.")
    quit()
if np.abs(m) > l:
    print("Kvanttiluvun l on oltava suurempi kuin |m|.")
    quit()

R_MAX = par.R_MAKSIMI
N_PTS = par.VOIMAKKUUS * 1000
RHO_MAX = 10**(-(0.7*np.log(par.VOIMAKKUUS) + 0.3*np.sqrt(par.R_MAKSIMI)))
KUVAFORMAATTI = par.KUVAFORMAATTI

BOHR_RADIUS = 0.529177210903 # Ångströms

# ket-esitys
def ket(n,l,m):

    return "|"+str(n)+","+str(l)+","+str(m)+")"


# pallokoordinaatit karteesisiksi
def to_cartesian(spherical):

    r = spherical[:,0]
    theta = spherical[:,1]
    phi = spherical[:,2]
    n = len(r)

    xyz = np.zeros([n,3])

    xyz[:,0] = r * np.sin(theta) * np.cos(phi) # x
    xyz[:,1] = r * np.sin(theta) * np.sin(phi) # y
    xyz[:,2] = r * np.cos(theta) # z
    
    return xyz


# karteesiset pallokoordinaateiksi
def to_spherical(xyz):

    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    n = len(x)

    spherical = np.zeros([n,3])

    xy = xyz[:,0]**2 + xyz[:,1]**2
    spherical[:,0] = np.sqrt(xy + xyz[:,2]**2) # r
    spherical[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2]) # theta
    spherical[:,2] = np.arctan2(xyz[:,1], xyz[:,0]) # phi
    
    return spherical


# aaltofunktion kulmariippuva osuuus
def angular_function(l, m, angle_from_z, angle_around_z):

    return sp.special.sph_harm(m, l, angle_around_z, angle_from_z)


# aaltofunktion säteestä riippuva osuus
def radial_function(n, l, distance):

    x = 2*distance/n
    laguerre = sp.special.genlaguerre(n-l-1, 2*l+1)
    fac = np.math.factorial
    norm = np.sqrt( (2/n)**3 * fac(n-l-1) / ( 2*n*fac(n+l)) )
    
    return norm * np.exp(-x/2) * x**l * laguerre(x)


# aaltofunktio pallokoordinaateissa
def wavefunction_spherical(n,l,m,distance,angle_around_z,angle_from_z):

    return radial_function(n,l,distance)*angular_function(l,m,angle_around_z,angle_from_z)


# aaltofunktio karteesisissa koordinaateissa
def wavefunction_cartesian(n,l,m,x,y,z):

    try:
        n_x = len(x)
    except:
        n_x = 1
    xyz = np.zeros([n_x, 3])
    xyz[:,0] = x 
    xyz[:,1] = y 
    xyz[:,2] = z
    
    spherical = to_spherical(xyz)
    
    return wavefunction_spherical(n,l,m,spherical[:,0],spherical[:,1],spherical[:,2])



def plot_radial_function(n,l,r_max=30, psi_max=0.5, show=1):

    plt.close('all')
    r_n = 1001
    r = np.linspace(0,r_max,r_n)
    psi = radial_function(n,l,r)

    plt.plot(r, np.zeros(len(r)), 'k')
    plt.plot(r, psi, lw=2, c='k')
    plt.xlim(0, r_max)
    plt.ylim(-psi_max, psi_max)
    plt.xlabel("$r$ ($a$)")
    plt.ylabel("$\\psi$ ($a^{-1/2}$)")

    plt.savefig("radiaalifunktio."+KUVAFORMAATTI)
    plt.close('all')
        
    
def plot_radial_probability(n,l,r_max=30, rho_max=0.6, show=1):

    plt.close('all')
    r_n = 1001
    r = np.linspace(0,r_max,r_n)
    rho = np.abs( radial_function(n,l,r) )**2 * r**2

    plt.clf()
    plt.plot(r, rho, 'k')
    plt.xlim(0, r_max)
    plt.ylim(0, rho_max)
    plt.xlabel("$r$ ($a$)")
    plt.ylabel("$\\rho$ ($a^{-1}$)")
    
    ax = plt.gca()
    zero = np.zeros(r_n)
    ax.fill_between(r, zero, rho, color=(0.5, 0.70, 0.70, 0.5) )

    plt.savefig("radiaalifunktio_tiheys_"+ket(n,l,m)+"."+KUVAFORMAATTI)
    plt.close('all')


def plot_revolving_function(l,m, theta = np.pi/2, R=4, show=1):

    plt.close('all')
    phi_n = 1001
    phi = np.linspace(0, 2*np.pi, phi_n)
    psi = angular_function(l,m,theta,phi)
    dtheta = 0.0
    
    # in case psi is 0 at this theta
    while np.max( np.abs(psi) ) < 0.001:
        dtheta += 0.01*np.pi
        psi = angular_function(l,m,theta+dtheta,phi)

    scale = np.max( np.abs(psi) )

    psi = psi / scale * 0.75 #* 0.28
    
    zero = np.ones(len(phi))*R
    
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    plt.plot( phi, zero, 'k', lw=1)
    plt.plot( phi, zero+psi.real, lw=2, c='k')
    plt.plot( phi, zero+psi.imag, lw=2, c='r')
    
    plt.axis('off')
    
    ax = plt.gca()
    ax.fill_between(phi, zero, zero+np.abs(psi), color=(0.78, 0.70, 0.60, 1) )

    plt.savefig("kulmafunktio_xy_"+ket(n,l,m)+"."+KUVAFORMAATTI)
    plt.close('all')


def plot_standing_wave_function(l,m, phi = 0, R=4, show=1):

    plt.close('all')
    theta_n = 1001
    halfpoint = 501
    theta = np.linspace(0, 2*np.pi, theta_n)
    first_half = theta[:halfpoint]
    second_half = 2*np.pi-theta[halfpoint:]
    
    first_psi = angular_function(l,m,first_half,phi) * 1.5
    second_psi = angular_function(l,m,second_half,phi+np.pi) * 1.5
    
    psi = np.concatenate([first_psi, second_psi])
    zero = np.ones(theta_n)*R
    
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    plt.plot( theta+np.pi/2, zero, 'k', lw=1)
    plt.plot( theta+np.pi/2, zero+psi.real, lw=2, c='k')
    plt.plot( theta+np.pi/2, zero+psi.imag, lw=2, c='r')
    
    plt.axis('off')
    
    ax = plt.gca()
    ax.fill_between(theta+np.pi/2, zero, zero+np.abs(psi), color=(0.78, 0.70, 0.60, 1) )

    plt.savefig("kulmafunktio_xz_"+ket(n,l,m)+"."+KUVAFORMAATTI)
    plt.close('all')


def draw_random_distances(n, l, amount=10000, r_max=30, gf=0):

    r_n = 10001
    r = np.linspace(0,r_max,r_n)
    if gf == 0:
        rho = np.abs( radial_function(n,l,r) )**2 
    elif gf == 1:
        rho = np.abs( radial_function(n,l,r) )**2 * r
    else:
        rho = np.abs( radial_function(n,l,r) )**2 * r**2
    rho /= np.sum(rho)
    
    rng = np.random.default_rng()
    R = rng.choice(r, p=rho, size=amount)
    
    return R


def draw_random_theta_angles(l, m, phi=0, amount=10000, gf=0):

    theta_n = 10001
    theta = np.linspace(0,np.pi,theta_n)
    if gf == 0:
        rho = np.abs( angular_function(l,m,theta,phi) )**2 
    else:
        rho = np.abs( angular_function(l,m,theta,phi) )**2 * np.sin(theta)
    rho /= np.sum(rho)

    rng = np.random.default_rng()
    THETA = rng.choice(theta, p=rho, size=amount)
    #THETA += rng.normal(size=amount) * np.pi/theta_n * 0.5
    
    return THETA


def draw_random_phi_angles(l, m, theta=np.pi/2, amount=10000):

    phi_n = 10001
    phi = np.linspace(0,2*np.pi,phi_n)
    rho = np.abs( angular_function(l,m,theta,phi) )**2
    rho /= np.sum(rho)
    
    rng = np.random.default_rng()
    PHI = rng.choice(phi, p=rho, size=amount)
    #PHI += rng.normal(size=amount) * 2*np.pi/phi_n * 0.5
    
    return PHI


def draw_random_spherical_coordinates(n,l,m, amount=10000, r_max=30):

    R = draw_random_distances(n, l, amount, r_max, 2)
    THETA = draw_random_theta_angles(l, m, 0, amount, 1)
    PHI = draw_random_phi_angles(l, m, 0.12345, amount=amount)
    
    return R, THETA, PHI


def draw_random_cartesian_coordinates(n,l,m, amount=10000, r_max=30):
    
    R, T, P = draw_random_spherical_coordinates(n,l,m,amount,r_max)    

    spherical = np.zeros([amount,3])
    spherical[:,0] = R
    spherical[:,1] = T
    spherical[:,2] = P
    
    return to_cartesian(spherical)
    

def plot_random_points_in_xz_plane(n,l,m, amount=10000, r_max=30, phi=0):

    plt.close('all')
    rng = np.random.default_rng()
    R = draw_random_distances(n,l,amount,r_max,gf=1)
    
    THETA = draw_random_theta_angles(l,m,phi,amount,gf=0)
    THETA *= rng.choice([-1,1],size=amount)
        
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    plt.plot( np.pi/2+THETA, R, 'o', ms=0.5, c=(0,0,0,0.2))
    
    plt.axis('off')
        
    plt.savefig("pistetiheys_xz_"+ket(n,l,m)+".png",dpi=400)
    plt.close('all')
    
    
def plot_density_in_xz_plane(n,l,m, r_max=30, rho_max=1e-6):

    plt.close('all')
    nx = 101
    x = np.linspace(-r_max, r_max, nx)    
    z = np.linspace(-r_max, r_max, nx)
    
    X,Z = np.meshgrid(x,z)
    P = np.zeros( [nx, nx] )

    for i in range(nx):
        for j in range(nx):
            P[i,j] = np.abs( wavefunction_cartesian(n,l,m,X[i,j],0,Z[i,j]) )**2
            if P[i,j] > rho_max:
                P[i,j] = rho_max

    #plt.contourf( X, Z, P, 100 )
    plt.pcolormesh( X, Z, P, cmap='magma' )
    
    plt.axis('off')
    plt.gca().set_aspect('equal')
        
    plt.savefig("todennakoisyys_xz_"+ket(n,l,m)+".png",dpi=400)
    plt.close('all')
    
    
def plot_random_points(n,l,m, amount=10000, r_max=30):

    plt.close('all')
    xyz = draw_random_cartesian_coordinates(n,l,m,amount,r_max)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect( (1,1,1) )
    #plt.axes('off')
    ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],color=(0,0,0,0.2),s=0.2)

    plt.savefig("pistetiheys_3d_"+ket(n,l,m)+".png",dpi=400)
    plt.close('all')


def wavefunction_grid(n,l,m,x_max=30,delta=0.1):

    n_x = int(2*x_max/delta+1)
    x = np.array( np.linspace(-x_max, x_max, n_x).tolist()*n_x**2 )
    y = np.array( np.repeat(np.linspace(-x_max, x_max, n_x), n_x).tolist()*n_x )
    z = np.repeat( np.linspace(-x_max, x_max, n_x), n_x**2)
    
    psi = wavefunction_cartesian(n,l,m,x,y,z)
    
    return x, y, z, psi

     

def plot_angular_real_part(l,m, R=2):

    plot_limit = 2
    plt.clf()

    theta_n = 101
    phi_n = 201
    theta = np.linspace(0, np.pi, theta_n)
    phi = np.linspace(0, 2*np.pi, phi_n)
    THETA, PHI = np.meshgrid(theta, phi)
    
    # angular_function(l, m, angle_from_z, angle_around_z)
    A = R + np.real( angular_function(l,m,THETA,PHI) )
    
    X = A * np.sin(THETA) * np.cos(PHI)
    Y = A * np.sin(THETA) * np.sin(PHI)
    Z = A * np.cos(THETA)
    
    X0 = R * np.sin(THETA) * np.cos(PHI)
    Y0 = R * np.sin(THETA) * np.sin(PHI)
    Z0 = R * np.cos(THETA)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlim(-plot_limit, plot_limit)
    ax.set_ylim(-plot_limit, plot_limit)
    ax.set_zlim(-plot_limit, plot_limit)
    ax.set_box_aspect( (1,1,1) )

    plot1 = ax.plot_surface(
    X, Y, Z, rstride=1, cstride=1, facecolors = matplotlib.cm.coolwarm((A-R+1)/2),
    linewidth=0, antialiased=True, alpha=1)


    plt.axis('off')
    
    plt.savefig("kulmafunktio_reaaliosa_"+ket(n,l,m)+".png",dpi=400)
    plt.close('all')


def plot_angular_density(l,m, R=2):

    plot_limit = R
    plt.clf()

    theta_n = 101
    phi_n = 201
    theta = np.linspace(0, np.pi, theta_n)
    phi = np.linspace(0, 2*np.pi, phi_n)
    THETA, PHI = np.meshgrid(theta, phi)
    
    A = R + np.abs( angular_function(l,m,THETA,PHI) )**2
    
    X = A * np.sin(THETA) * np.cos(PHI)
    Y = A * np.sin(THETA) * np.sin(PHI)
    Z = A * np.cos(THETA)
    
    X0 = R * np.sin(THETA) * np.cos(PHI)
    Y0 = R * np.sin(THETA) * np.sin(PHI)
    Z0 = R * np.cos(THETA)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlim(-plot_limit, plot_limit)
    ax.set_ylim(-plot_limit, plot_limit)
    ax.set_zlim(-plot_limit, plot_limit)
    ax.set_box_aspect( (1,1,1) )
    

    plot1 = ax.plot_surface(
    X, Y, Z, rstride=1, cstride=1, facecolors = matplotlib.cm.coolwarm((A-R+1)/2),
    linewidth=0, antialiased=True, alpha=1)


    plt.axis('off')
    
    plt.savefig("kulmafunktio_tiheys_"+ket(n,l,m)+".png",dpi=400)
    plt.close('all')


def draw_angular_real_part_frame(frame, THETA, PHI, PSI, dphi):
    
    plot_limit = 1.8
    plt.clf()
    ax = plt.axes(projection ='3d') 
    ax.set_xlim(-plot_limit, plot_limit)
    ax.set_ylim(-plot_limit, plot_limit)
    ax.set_zlim(-plot_limit, plot_limit)
    ax.set_box_aspect( (1,1,1) )
    
    R = 2
    A = R + np.real( PSI * np.exp(-2j*np.pi*frame*dphi) )
    X = A * np.sin(THETA) * np.cos(PHI)
    Y = A * np.sin(THETA) * np.sin(PHI)
    Z = A * np.cos(THETA)
    
    X0 = R * np.sin(THETA) * np.cos(PHI)
    Y0 = R * np.sin(THETA) * np.sin(PHI)
    Z0 = R * np.cos(THETA)
    
    plot2 = ax.plot_surface(
    X, Y, Z, rstride=1, cstride=1, facecolors = matplotlib.cm.coolwarm((A-R+1)/2),
    linewidth=0, antialiased=True)

    plt.axis('off')


def calculate_E_and_L(n,l,m):
    E = -13.60 / n**2
    L = np.sqrt( l * (l+1) )
    Lz = m
    
    print("""
    tila """+ket(n,l,m)+f"""
    energia                        {E:8.2f} eV
    kulmaliikemäärä                {L:8.2f} h-viiva
    kulmaliikemäärän z-komponentti {Lz:8d} h-viiva
""")


def main():
    plot_random_points(n,l,m,r_max=R_MAX)
    plot_revolving_function(l,m)
    plot_standing_wave_function(l,m)
    plot_angular_real_part(l,m)
    plot_radial_probability(n,l,r_max=R_MAX)
    plot_random_points_in_xz_plane(n,l,m,N_PTS,r_max=R_MAX)
    plot_density_in_xz_plane(n,l,m,r_max=R_MAX,rho_max=RHO_MAX)
    
    calculate_E_and_L(n,l,m)

if __name__ == "__main__":
    main()
