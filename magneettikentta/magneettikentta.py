#
# Laskee monikulmion muotoisen virtasilmukan magneettikentän.
#
# Teemu Hynninen 2025
#
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# haetaan parametrit
import parametrit as par

EPSILON = 8.8541878*10**-12; # sähkövakio
MU = 1.256637061*10**-6; # magneettivakio
KUVAFORMAATTI = par.KUVAFORMAATTI
KULMAT = par.KULMAT
VIRTA = par.VIRTA
for i in range(len(KULMAT)):
    KULMAT[i] = np.array(KULMAT[i])
DL = par.OSIEN_KOKO
ALKU = np.array( par.ALKUPISTE )
LOPPU = np.array( par.LOPPUPISTE )
XMAX = par.X_MAKSIMI
B_SKAALA = par.B_SKAALA


# laskee vektorin pituuden
def pituus(vektori):
    return np.sqrt(vektori@vektori)

# funktio pituus loppuu


# laskee pisteessä rl olevan, nopeudella v liikkuvan varauksen q
# magneettikentän pisteessä r
def pistevarauksen_magneettikentta(q, v, rq, r):

    # vektori varauksesta tarkastelupisteeseen ja sen pituus
    #
    # Lisätään pituuteen ihan pieni numero, jotta jos tarkastelupiste
    # on sama kuin varauksen paikka, ei päädytä jakamaan nollalla.
    r_vektori = r - rq
    r = pituus(r_vektori)+0.000001
    
    # Biot'n ja Savartin laki
    return MU * q * np.cross( v, r_vektori ) / ( 4.0 * np.pi * r**3 ) 

# funktio pistevarauksen_magneettikentta loppuu


# laskee pisteestä r_alku pisteeseen r_loppu kulkevan suoran johtimen,
# jossa kulkee virta i, magneettikentän pisteessä r
def suoran_johtimen_magneettikentta(r_alku, r_loppu, virta, r):

    # jaetaan johdin pieniin osiin, lasketaan kunkin osan kenttä
    # ja summataan nämä yhteen    
    dr = r_loppu - r_alku
    L = pituus( dr ) # johtimen pituus
    n = int( L / DL + 0.5 ) # osien lukumäärä
    dL = dr / n # yhtä osaa kuvaava vektori
    
    B = np.zeros(3)
    
    # käydään läpi kaikki osat
    for i in range(n):
    
        # osan paikkavektori
        r_johdin = r_alku + ( 0.5 + i ) * dL
        
        # osan magneettikenttä
        B += pistevarauksen_magneettikentta( virta, dL, r_johdin, r )

    return B

# funktio suoran_johtimen_magneettikentt loppuu


# laskee virtasilmukan magneettikentän pisteessä r.
# Silmukka on monikulmio, ja se määritellään antamalla
# lista sen kulmapisteiden koordinaateista, järjestyksessä.
def silmukan_magneettikentta(kulmat, virta, r):

    B = np.zeros(3)
    
    # käydään läpi kaikki monikulmion sivut
    for i in range(len(kulmat)):
        alkupiste = kulmat[i-1]
        loppupiste = kulmat[i]

        B += suoran_johtimen_magneettikentta( alkupiste, loppupiste, virta, r )

    return B

# funktio silmukan_magneettikentta loppuu


# laskee magneettikentän suoralla, joka kulkee pisteestä
# r_alku pisteeseen r_loppu
def laske_magneettikentta_suoralla(r_alku, r_loppu):

    # jaetaan suora pisteisiin ja lasketaan kenttä kussakin pisteessä
    dr = r_loppu - r_alku
    L = pituus( dr )
    n = 201
    dL = dr / n
    
    B_data = []
    r_data = []
    
    # käydään läpi kaikki pisteet
    for i in range(n):
    
        # tarkastelupisteen paikka
        r = r_alku + ( 0.5 + i ) * dL

        # kenttä tarkastelupisteessä
        B = silmukan_magneettikentta(KULMAT, VIRTA, r)

        r_data += [r]
        B_data += [B]
    
    return np.array(r_data), np.array(B_data)

# funktio laske_magneettikentta_suoralla loppuu


# Piirretään joukko kuvaajia ja tallennetaan ne tiedostoihin.
def piirra_kuvaajat(r ,B):
    
    L = pituus( r[-1] - r[0] )
    l = np.linspace( 0, L, len(r[:,0]) )
    Bx = B[:,0]
    By = B[:,1]
    Bz = B[:,2]

    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlabel("$r$ (m)")
    plt.ylabel("$B_x$ (V)")
    plt.plot(l, Bx, c='black', lw=2)
    plt.grid()   
    plt.savefig("Bx."+KUVAFORMAATTI)
    
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlabel("$r$ (m)")
    plt.ylabel("$B_y$ (V)")
    plt.plot(l, By, c='black', lw=2)
    plt.grid()     
    plt.savefig("By."+KUVAFORMAATTI)
    
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlabel("$r$ (m)")
    plt.ylabel("$B_z$ (V)")
    plt.plot(l, Bz, c='black', lw=2)
    plt.grid()     
    plt.savefig("Bz."+KUVAFORMAATTI)

# funktio piirra_kuvaajat loppuu


# piirretään kuva virtasilmukasta, jonka kenttää lasketaan
def piirra_silmukka():

    plt.close('all')
    fig, ax = plt.subplots()
    for i in range(len(KULMAT)):
        x0 = KULMAT[i-1][0]
        y0 = KULMAT[i-1][1]
        x1 = KULMAT[i][0]
        y1 = KULMAT[i][1]
        plt.plot([x0,x1], [y0,y1], 'k-', lw=3)
    
    ax.set_aspect('equal')    
    plt.savefig("silmukka."+KUVAFORMAATTI)

# funktio piirra_silmukka loppuu



# Piirretään kuva systeemistä.
#
def piirra_systeemi(kulmat, virta, xmax=XMAX, bmax=XMAX/10):

    plt.close('all')
    fig, ax = plt.subplots()
    
    # kuinka monessa pisteessä nuolet piirretään
    ngrid = int( (2*xmax) // bmax + 0.01 )
    # kuinka kaukana toisistaan nämä pisteet ovat
    dx = (2*xmax) / ngrid
    
    xs = np.linspace(-xmax+0.5*dx, xmax-0.5*dx, ngrid)
    ys = np.linspace(-xmax+0.5*dx, xmax-0.5*dx, ngrid)
    B = [ ]
    
    # käydään läpi kaikki pisteet, joihin kenttä piirretään,
    # ja lasketaan kenttä niissä
    for i in range(ngrid):
        for j in range(ngrid):
        
                xP = xs[ j ]
                zP = ys[ i ]
                rP = np.array( [ xP, 0, zP ] )
                        
                bP = silmukan_magneettikentta(kulmat, virta, rP) * B_SKAALA
                    
                if pituus(bP) > dx:
                    bP = bP / pituus(bP) * dx
                B += [ ( np.array( (xP, zP) ), np.array( (bP[0], bP[2]) ) ) ]
        
    # piirretään sähkökenttä nuolina
    for rP, bP in B:
        
        kentta = FancyArrowPatch(rP-0.5*bP, rP+0.5*bP, 
                    arrowstyle='->', color=(0.5,0,0.5,1), mutation_scale=2.5,
                    shrinkA=0.2,shrinkB=0.2)
        ax.add_patch(kentta)

    for i in range(len(KULMAT)):
        x0 = KULMAT[i-1][0]
        z0 = KULMAT[i-1][2]
        x1 = KULMAT[i][0]
        z1 = KULMAT[i][2]
        plt.plot([x0,x1], [z0,z1], 'k-', lw=3)
        
    # tallennetaan kuva
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-xmax, xmax)
    ax.set_aspect('equal')
    plt.grid(False)
    plt.xlabel('x (m)')
    plt.ylabel('z (m)')
    plt.savefig(f"kentta."+KUVAFORMAATTI)
    plt.close('all')

# funktio piirra_systeemi loppuu


# Pääohjelma.
def main():
    r, B = laske_magneettikentta_suoralla(ALKU, LOPPU)
    piirra_kuvaajat(r, B)
    piirra_silmukka()
    piirra_systeemi(KULMAT, VIRTA)
    
# funktio main loppuu


# Ajetaan pääohjelma.
if __name__ == '__main__':
    main()

