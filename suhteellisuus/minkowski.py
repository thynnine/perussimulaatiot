#
# Visualisoidaan relativistinen liike.
#
# Ohjelma piirtää tiedostossa parametrit.py
# määriteltyjen ohjeiden mukaisesti ennalta
# simuloidut liikeradat kahdessa inertiaalikoordinaatistossa.
#
# Teemu Hynninen 2025
#
import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.animation import FuncAnimation
import warnings
warnings.filterwarnings("ignore")

# Lorentz-muuntaa annetut koordinaatit uuteen koordinaatistoon,
# jonka nopeus alkuperäisen suhteen on v.
def lorentz_muunnos(x, t, v, c=300.0):
    gamma = 1 / np.sqrt(1 - v**2 / c**2)
    xB = gamma * (x - v * t)
    tB = gamma * (t - (v * x) / c**2)
    return xB, tB

# funktio lorentz_muunnos loppuu


# Muuntaa joukon tapahtumia koordinatistosta A koordinaatistoon B,
# missä v on A:n nopeus B:n suhteen.
def muunna_joukko_AB(tapahtumat, v, c=300.0):

    B = []
    for e in tapahtumat:
        x, t = e[0], e[1]
        B += [ lorentz_muunnos(x, t, v, c) ]

    return B

# funktio muunna_joukko_AB loppuu


# Muuntaa joukon tapahtumia koordinatistosta B koordinaatistoon A,
# missä v on A:n nopeus B:n suhteen.
def muunna_joukko_BA(tapahtumat, v, c=300.0):

    return muunna_joukko_AB(tapahtumat, -v, c=300.0)

# funktio muunna_joukko_BA loppuu


# Etsii annetulta liikeradalta eli x,t-koordinaattivektorilta,
# missä kappale on annetulla ajan hetkellä t.
# Tekee lineaarisen interpolaation kahden lähimmän ajan hetken väliltä.
# Jos aikaväli on annetun datan ulkopuolella, palautetaan ääretön.
def paikka_radalta( rata, t ):

    for i in range(len(rata)-1):
        x1, t1 = rata[i,0], rata[i,1] 
        x2, t2 = rata[i+1,0], rata[i+1,1]
        
        if t1 <= t and t2 > t:
            return x1 + (x2-x1) * (t-t1)/(t2-t1)
            
    return np.inf

# funktio paikka_radalta loppuu


# Etsii suoran viivan leikkauspisteen suoran x=vakio kanssa.
def t_leikkaus( viiva, x ):
    x1, t1, x2, t2 = viiva
    k = (t2-t1) / (x2-x1)
    t = k * (x - x1) + t1
    return t

# funktio t_leikkaus loppuu


# Etsii suoran viivan leikkauspisteen suoran t=vakio kanssa.
def x_leikkaus( viiva, t ):
    x1, t1, x2, t2 = viiva
    k = (x2-x1) / (t2-t1)
    x = k * (t - t1) + x1
    return x

# funktio x_leikkaus loppuu


# Laskee annettujen ratakäyrien nopeuden numeerisesti ja piirtää
# nopeuden ajan funktiona.
def piirra_kuvaajat(kayrat, v, c=300.0, koko=4, n_points=501, koordinaatisto_A=True):

    plt.close('all')

    tmin = -koko
    tmax = koko
    dt = (tmax-tmin) / (n_points-1)
    
    if koordinaatisto_A:
        kirjain = 'A'
    else:
        kirjain = 'B'
        
    rata_data = []
    if koordinaatisto_A:
        rata_data = copy.copy(radat)
    else:
        for rata in radat:
            rata_data += [ np.array( muunna_joukko_AB( rata, v, c ) ) ]    
    
    xs = []
    ts = []
    vs = []
    
    for i in range(n_points):
        t = tmin + i*dt
        ts.append( t )
        
        uudet_paikat = []
        for rata in rata_data:
            uudet_paikat.append( paikka_radalta( rata, t ) )
            
        xs.append( uudet_paikat )

    xs = np.array(xs)
    
    for i in range(1,len(xs[:,0])-1):
        vs += [ (xs[i+1,:] - xs[i-1,:]) / (2*dt) ]
        
    vs = np.array(vs)    
    
    # nopeus
    for i in range(len(vs[0,:])):
        valku = vs[0,i]
        vloppu = vs[-1,i]
        info = ""
        if np.abs(valku) <= 1.1*c:
            info += "$v_{x,alku} = $"+f"{valku:5.1f}"+" m/μs"
        if np.abs(vloppu) <= 1.1*c:
            if len(info) > 0:
                info += " | "
            info += "$v_{x,loppu} = $"+f"{vloppu:5.1f}"+" m/μs "
        plt.plot( ts[1:-1], vs[:,i], '-', c=VARIT[i%len(VARIT)], label=info )
        
    plt.plot( [tmin, tmax], [c, c], ':', c='gray' )
    plt.plot( [tmin, tmax], [-c, -c], ':', c='gray' )
    plt.legend()
    plt.xlabel("$t_{("+kirjain+")}$ (μs)")
    plt.ylabel("$v_{x,("+kirjain+")}$ (m/μs)")
    plt.xlim(tmin, tmax)
    plt.ylim(-1.1*c, 1.1*c)
    plt.grid()
    plt.savefig("nopeus_"+kirjain+"."+KUVAFORMAATTI)
    plt.close('all')
   
# funktio piirra_kuvaajat loppuu     
    
    
# Piirtää nopeuden ajan funktiona sekä koordinaatistossa A että B.
def piirra_kuvaajat_AB(kayrat, v, c=300.0, koko=4, n_points=501):
    piirra_kuvaajat(kayrat, v, c, koko, n_points, koordinaatisto_A=True)
    piirra_kuvaajat(kayrat, v, c, koko, n_points, koordinaatisto_A=False)

# funktio piirra_kuvaajat_AB loppuu


# Piirtää annetut tapahtumat ja ratakäyrät Minkowski-kuvaajaan
# sekä koordinaatiston A että B näkökulmasta.
# Tässä v on B:n nopeus A:n suhteen ja koko on piirtoalueen koko
# aika-akselilla mikrosekunneissa.
def piirra_minkowski_koordinaatisto(tapahtumat, kayrat, v, c=300.0, koko=4):
    
    tmin = -koko
    tmax = koko
    xmin = tmin*c
    xmax = tmax*c
    labelcount = 5
    
    plt.close('all')
    fig, ax = plt.subplots()
    
    # koordinaatisto A
    viivat = []
    tyylit = []
    x_leikkuut = []
    x_arvot = []
    t_leikkuut = []
    t_arvot = []
    for i in range(-koko,koko+1):
        ( x1, t1, x2, t2 ) = ( i*c, tmin, i*c, tmax )
        ( x3, t3, x4, t4 ) = ( xmin, i, xmax, i )
        viivat += [ ( x1, t1, x2, t2 ) ]
        viivat += [ ( x3, t3, x4, t4 ) ]
        
        xx = x_leikkaus( ( x1, t1, x2, t2 ), tmin )
        if xx >= xmin and xx <= xmax:           
            x_leikkuut += [ xx ]
            x_arvot += [ i*c ]
            
        tt = t_leikkaus( ( x3, t3, x4, t4 ), xmin )
        if tt >= tmin and tt <= tmax:           
            t_leikkuut += [ tt ]
            t_arvot += [ i ]
            
        if i != 0:
            tyylit += [ ':', ':' ]
        else:
            tyylit += [ '-', '-']

    for i in range(len(viivat)):
        plt.plot( [ viivat[i][0], viivat[i][2] ], [ viivat[i][1], viivat[i][3] ], lw = 1, ls = tyylit[i], c='black' )

    tapahtumat = np.array( tapahtumat )

    if len(tapahtumat) > 0:
        ts = tapahtumat[:,1]
        xs = tapahtumat[:,0]
        plt.plot( xs, ts, 'ro' )
    
    for i in range(len(kayrat)):
        kk = np.array( kayrat[i] )
        ts = kk[:,1]
        xs = kk[:,0]
        plt.plot( xs, ts, '-', color=VARIT[i%len(VARIT)] )
        
    
    m = len(x_leikkuut) // labelcount + 1
    x_arvot = np.array( x_arvot , dtype=int )
    t_arvot = np.array( t_arvot , dtype=int )
    x0 = np.where(x_arvot == 0)[0][0]
    t0 = np.where(x_arvot == 0)[0][0]
    xoff = x0%m
    toff = t0%m
    plt.xticks( x_leikkuut[xoff::m] )
    plt.yticks( t_leikkuut[toff::m] )
    ax.set_xticklabels( x_arvot[xoff::m] )
    ax.set_yticklabels( t_arvot[toff::m] )    
    
    plt.ylabel("$t_{(A)}$ (μs)")
    plt.xlabel("$x_{(A)}$ (m)")
    plt.ylim(tmin, tmax)
    plt.xlim(xmin, xmax)
    ax.set_aspect(c)
    plt.savefig("koordinaatisto_A."+KUVAFORMAATTI)
    
    plt.close('all')
    fig, ax = plt.subplots()
    
    # koordinaatisto B
    viivat = []
    tyylit = []
    x_leikkuut = []
    x_arvot = []
    t_leikkuut = []
    t_arvot = []
    
    
    n = 20  
    for i in range(-n*koko,n*koko+1):
        x1, t1 = lorentz_muunnos( i*c, tmin*n, -v, c )
        x2, t2 = lorentz_muunnos( i*c, tmax*n, -v, c )
        x3, t3 = lorentz_muunnos( xmin*n, i, -v, c )
        x4, t4 = lorentz_muunnos( xmax*n, i, -v, c )
        viivat += [ ( x1, t1, x2, t2 ) ]
        viivat += [ ( x3, t3, x4, t4 ) ]
        
        xx = x_leikkaus( ( x1, t1, x2, t2 ), tmin )
        if xx >= xmin and xx <= xmax:           
            x_leikkuut += [ xx ]
            x_arvot += [ i*c ]
            
        tt = t_leikkaus( ( x3, t3, x4, t4 ), xmin )
        if tt >= tmin and tt <= tmax:           
            t_leikkuut += [ tt ]
            t_arvot += [ i ]
         
        if i != 0:
            tyylit += [ ':', ':' ]
        else:
            tyylit += [ '-', '-']     

    for i in range(len(viivat)):
        plt.plot( [ viivat[i][0], viivat[i][2] ], [ viivat[i][1], viivat[i][3] ], lw = 1, ls = tyylit[i], c='black' )

    if len(tapahtumat) > 0:
        ts = tapahtumat[:,1]
        xs = tapahtumat[:,0]
        plt.plot( xs, ts, 'ro' )
    
    for i in range(len(kayrat)):
        kk = np.array( kayrat[i] )
        ts = kk[:,1]
        xs = kk[:,0]
        plt.plot( xs, ts, '-', color=VARIT[i%len(VARIT)] )
    
    m = len(x_leikkuut) // labelcount + 1
    x_arvot = np.array( x_arvot , dtype=int )
    t_arvot = np.array( t_arvot , dtype=int )
    x0 = np.where(x_arvot == 0)[0][0]
    t0 = np.where(x_arvot == 0)[0][0]
    xoff = x0%m
    toff = t0%m
    plt.xticks( x_leikkuut[xoff::m] )
    plt.yticks( t_leikkuut[toff::m] )
    ax.set_xticklabels( x_arvot[xoff::m] )
    ax.set_yticklabels( t_arvot[toff::m] )    
    
    plt.ylabel("$t_{(B)}$ (μs)")
    plt.xlabel("$x_{(B)}$ (m)")
    plt.ylim(tmin, tmax)
    plt.xlim(xmin, xmax)
    ax.set_aspect(c)
    plt.savefig("koordinaatisto_B."+KUVAFORMAATTI)
    
    plt.close('all')
 
# funktio piirra_minkowski_koordinaatisto loppuu   


# Piirtää animaation ratakäyrien kuvaamien hiukkasten liikkeestä.
def piirra_animaatio(v,c=300.0,koko=4,n_frames=161,koordinaatisto_A=True):    

    tmin = -koko
    tmax = koko
    xmin = tmin*c
    xmax = tmax*c
    dt = (tmax-tmin) / (n_frames-1)
    
    if koordinaatisto_A:
        kirjain = 'A'
    else:
        kirjain = 'B'
            
    plt.close('all')
    fig, ax = plt.subplots()
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)
    ax.set_aspect('equal')
    ax.set_xlabel('$x_{('+kirjain+')}$ (m)')
    ax.set_ylabel('')
    ax.set_yticks([0])
    
    R = 100

    # Tallennetaan radat listaan rata_data.
    # Jos ei vaihdeta koordinaatistoa, tämä on vain sama data
    # kuin listassa radat.
    # Jos vaihdetaan koordinaatistoa, lasketaan ratapisteet
    # Lorentz-muunnoksella.
    rata_data = []
    if koordinaatisto_A:
        rata_data = copy.copy(radat)
    else:
        for rata in radat:
            rata_data += [ np.array( muunna_joukko_AB( rata, v, c ) ) ]    
    
    # Sitten etsitään, missä radat ovat milläkin ajan hetkellä
    # valitussa koordinaatistossa.
    xs = []
    for i in range(n_frames):
        t = tmin + i*dt
        
        uudet_paikat = []
        for rata in rata_data:
            uudet_paikat.append( paikka_radalta( rata, t ) )
            
        xs.append( uudet_paikat )

    xs = np.array(xs)


    # luodaan listat, johon ympyrät tallennetaan
    circle_patches = []
    time_text = ax.text(xmin + 0.05*(xmax - xmin), xmax - 0.1*(xmax - xmin), '', 
                    fontsize=12, color='black')

    
    for i in range(len(radat)):
        paikka = Circle((0, 0), R, fill=True, color=VARIT[i%len(VARIT)])
        ax.add_patch(paikka)
        circle_patches.append(paikka)
    

    def alusta():
        for circle in circle_patches:
            circle.set_center((0, 0))
        time_text.set_text('')
        return circle_patches

    # siirretään ympyrät oikeille paikoille
    def paivita(frame):

        t = tmin + frame*dt        
        time_text.set_text('$t_{('+kirjain+')} = $'+f'{t:6.2f} μs')

        for i in range(len(xs[0])):
            paikka = ( xs[frame,i], 0 )
            circle_patches[i].center = paikka
            circle_patches[i].set_center(paikka)

        return circle_patches + [time_text]
    
    # luodaan animaatio
    ani = FuncAnimation(fig, paivita, frames=n_frames, init_func=alusta, blit=True, interval=20)
    if TALLENNA_ANIMAATIO:
        print("tallennetaan animaatio "+kirjain)
        ani.save("animaatio_"+kirjain+".mp4")
    else:
        plt.show()
    plt.close(fig)

# funktio piirra_animaatio loppuu


# Piirtää animaation sekä koordinaatiston A etät B näkökulmasta.
def piirra_animaatio_AB(v,c=300.0,koko=4,n_frames=161):
    piirra_animaatio(v, c, koko, n_frames, koordinaatisto_A = True)
    piirra_animaatio(v, c, koko, n_frames, koordinaatisto_A = False)

# funktio piirra_animaatio_AB loppuu


# Lukee liikeradan koordinaatit annetuista tiedostoista.
def lue_radat( tiedostot ):

    radat = []

    for tiedosto in tiedostot:
        try:
            rata = np.loadtxt(tiedosto, delimiter=',', skiprows=0)
            radat += [rata]
        except:
            print("Tiedoston "+tiedosto+" lukeminen ei onnistunut.")

    return radat

# funktio lue_radat loppuu


# Listaa tapahtumien koordinaatit sekä koordinaatistoissa A että B.
def listaa_tapahtumat( tapahtumat ):

    print()
    for i in range(len(tapahtumat)):
        xA, tA = tapahtumat[i]
        xB, tB = lorentz_muunnos( xA, tA, v, c )
        print(f"tapahtuma {i+1:3d} | A: ({xA:10.2f} m, {tA:10.4f} μs ) | B: ({xB:10.2f} m, {tB:10.4f} μs )")

    print()

# funktio listaa_tapahtumat loppuu


# Pääohjelma
def main():
    if len(tapahtumat) > 0:
        listaa_tapahtumat(tapahtumat)

    piirra_minkowski_koordinaatisto(tapahtumat, radat, v, koko=8)
    if len(radat) > 0:
        piirra_kuvaajat_AB(radat, v, koko=8)
        piirra_animaatio_AB(v, koko=8)

# funktio main loppuu


# luetaan parametrit toisesta tiedostosta
import parametrit as par

KUVAFORMAATTI = par.KUVAFORMAATTI
TALLENNA_ANIMAATIO = par.TALLENNA_ANIMAATIO
VARIT = ['tab:blue', 'black', 'tab:orange', 'tab:brown', 'tab:olive',
'tab:gray', 'tab:purple', 'tab:red', 'tab:cyan', 'tab:green' ]

c = 300.0
v = par.VXBA

tapahtumat = par.TAPAHTUMAT_A + muunna_joukko_BA(par.TAPAHTUMAT_B, v, c)
radat = lue_radat(par.RADAT)

if __name__ == "__main__":
    main()