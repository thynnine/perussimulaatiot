#
# Laskee ja piirtää sähkökenttiä ja potentiaaleja.
#
# Teemu Hynninen 2025
#
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.colors import ListedColormap
from scipy.interpolate import interpn

# haetaan parametrit
import parametrit as par

EPSILON = 8.8541878*10**-12; # sähkövakio
Q = par.VARAUKSET
KUVAFORMAATTI = par.KUVAFORMAATTI    
PISTEET = par.TARKASTELUPISTEET
RP = par.PISTEIDEN_KOKO
RESO = par.POTENTIAALIN_RESOLUUTIO
SAHKOKENTTA_POTENTIAALISTA = par.SAHKOKENTTA_POTENTIAALISTA
Y_MIN = par.Y_KESKUS-par.KUVAN_LEVEYS/2
Y_MAX = par.Y_KESKUS+par.KUVAN_LEVEYS/2
X_MIN = par.X_KESKUS-par.KUVAN_LEVEYS/2
X_MAX = par.X_KESKUS+par.KUVAN_LEVEYS/2
E_SKAALA = par.NUOLTEN_MITTAKAAVA
E_MAX = par.NUOLTEN_ETAISYYS
V_MIN = par.PIENIN_POTENTIAALI
V_MAX = par.SUURIN_POTENTIAALI
JANA = np.array( par.KUVAAJAN_PAATEPISTEET )
PIIRRA_E = par.PIIRRA_SAHKOKENTTA
PIIRRA_V = par.PIIRRA_POTENTIAALI


# laskee vektorin pisteestä (x0,y0) pisteeseen (x1,y1)
def siirtymavektori(x0, y0, x1, y1):
    return np.array( [x1-x0, y1-y0] )

# funktio siirtymavektori loppuu


# laskee vektorin pituuden
def pituus(vektori):
    return np.sqrt(vektori@vektori)

# funktio pituus loppuu


# laskee pisteessä (xq,yq) olevan varauksen q sähkökentän pisteessä (x,y)
def pistevarauksen_sahkokentta(q, xq, yq, x, y):

    # vektori varauksesta tarkastelupisteeseen ja sen pituus
    #
    # Lisätään pituuteen ihan pieni numero, jotta jos tarkastelupiste
    # on sama kuin varauksen paikka, ei päädytä jakamaan nollalla.
    r_vektori = siirtymavektori(xq,yq,x,y)
    r = pituus(r_vektori)+0.000001
    
    return q / (4.0 * np.pi * EPSILON * r**3 ) * r_vektori 

# funktio pistevarauksen_sahkokentta loppuu


# laskee pisteessä (xq,yq) olevan varauksen q potentiaalin pisteessä (x,y)
def pistevarauksen_potentiaali(q, xq, yq, x, y):

    r_vektori = siirtymavektori(xq,yq,x,y)
    r = pituus(r_vektori)+0.000001
    
    return q / (4.0 * np.pi * EPSILON * r )
    
# funktio pistevarauksen_potentiaali loppuu


# laskee kaikkien varausten potentiaalien summan pisteessä (x,y)
def laske_potentiaali_pisteessa( varaukset, x, y ):

    # alustetaan muuttuja, johon potentiaali tallennetaan
    V = 0.0
     
    # käydään läpi kaikki varaukset
    for q, (xq, yq) in varaukset:
        
        # lasketaan yhteen potentiaalien summa
        V += pistevarauksen_potentiaali(q, xq, yq, x, y)

    # palautetaan lopputulos
    return V

# funktio laske_potentiaali_pisteessa loppuu


# Lasketaan potentiaalin arvot säännöllisin välein ja tallennetaan ne.
# Näitä taulukoita käytetään sekä potentiaalin tasa-arvokäyrien piirtämiseen
# että potentiaalin laskemiseen.
V = None
X = None
Y = None
XY = None
def alusta_potentiaali( Q ):
    global X, Y, XY, V

    # Jos potentiaali on jo laskettu, ei tehdä mitään
    if V is None:
    
        xs = np.linspace(X_MIN, X_MAX, RESO)
        ys = np.linspace(Y_MIN, Y_MAX, RESO)
        X, Y = np.meshgrid(xs, ys)
        XY = (ys, xs)
        V = np.zeros([RESO, RESO])
    
        for i in range(RESO):
            for j in range(RESO):
        
                xP = X[ i, j ]
                yP = Y[ i, j ]
            
                V[ i, j ] = laske_potentiaali_pisteessa(Q, xP, yP)

# funktio alusta_potentiaali loppuu
    
     
# Laskee potentiaalin pisteessä (x,y).
#
# Funktio ei varsinaisesti laske potentiaalia itse vaan se käyttää
# valmista taulukkoa V, johon on tallennettu valmiiksi potentiaalin arvoja
# säännöllisin välein xy-tasossa. Potentiaalin arvo vain interpoloidaan
# tästä valmiista taulukosta scipyn funktiolla interpn.
#
# Tämä säästää aikaa jos pistevarauksia on monta, koska potentiaalia ei 
# tarvitse laskea aina uudestaan niiden avulla. Lisäksi tällä tavalla
# potentiaali voidaan laskea ilman tietoa sen luoneista varauksista.
# Menetelmän huono puoli on se, että interpolointi ei ole tarkkaa,
# joten lopputulokseen tulee pieniä numeerisia epätarkkuuksia.
#
# Jos potentiaali on välin [V_MIN, V_MAX] ulkopuolella, palautetaan
# V_MIN tai V_MAX. Näin siksi, että potentiaali kasvaa rajatta
# pistevarausten lähellä, ja ilman tällaista katkaisua potentiaalin
# kuvaajasta ei näkisi juuri mitään.
def potentiaali( x, y ):
    global XY, V
    
    # Jos taulukkoa V ei ole vielä luotu, luodaan se.
    if V is None:
        alusta_potentiaali( Q )
    
    try:
        delta = 0.01
        
        # interpoloidaan potentiaali valmiista taulukosta
        potentiaali = interpn( XY, V, np.array([y,x]) )[0]
        
        # katkaistaan potentiaali, jos se on liian suuri tai pieni
        if potentiaali > V_MAX+delta:
            return V_MAX
        elif potentiaali < V_MIN-delta:
            return V_MIN
        else:
            return potentiaali
    except:
        return 0     

# funktio potentiaali loppuu


# Piirretään kuva systeemistä.
#
# varaukset: systeemin pistevaraukset
# pistet: tarkastelupisteet, jotka merkitään kuvaan
# emax: tätä voimakkaammat kentät piirretään kaikki yhtä pitkinä nuolina
# tasot: montako tasa-arvokäyrää potentiaalille piirretään
# piirra_E: piirretäänkö sähkökenttä nuolina?
# piirra_V: piirretäänkö potentiaali tasa-arvokäyrinä
def piirra_systeemi(varaukset, pisteet, emax=E_MAX, tasot = 11, 
                    piirra_E=PIIRRA_E, piirra_V=PIIRRA_V):

    plt.close('all')
    fig, ax = plt.subplots()
    
    # kuinka monessa pisteessä nuolet piirretään
    ngrid = int( (X_MAX-X_MIN+0.01) // emax )
    # kuinka kaukana toisistaan nämä pisteet ovat
    dx = (X_MAX-X_MIN) / ngrid
    
    xs = np.linspace(X_MIN+0.5*dx, X_MAX-0.5*dx, ngrid)
    ys = np.linspace(Y_MIN+0.5*dx, Y_MAX-0.5*dx, ngrid)
    C = np.linspace(V_MIN, V_MAX, tasot)
    E = [ ]
    
    # käydään läpi kaikki pisteet, joihin sähkökenttä piirretään,
    # ja lasketaan kenttä niissä
    if piirra_E:
        for i in range(ngrid):
            for j in range(ngrid):
        
                xP = xs[ j ]
                yP = ys[ i ]
                rP = ( xP, yP )
                        
                if SAHKOKENTTA_POTENTIAALISTA:
                    eP = par.sahkokentta_potentiaalista(xP, yP) * E_SKAALA
                else:
                    eP = par.sahkokentta(Q, xP, yP) * E_SKAALA
                    
                if pituus(eP) > dx:
                    eP = eP / pituus(eP) * dx
                E += [ (rP, eP) ]
        
    # piirretään potentiaali
    if piirra_V:
        N = tasot
        vals = np.ones((3*N, 4))
        vals[0*N:1*N, 0] = np.linspace(75/256, 1, N) # red
        vals[0*N:1*N, 1] = np.linspace(56/256, 200/256, N) # green
        vals[0*N:1*N, 2] = np.linspace(0, 0, N) # blue
        vals[1*N:2*N, 0] = np.linspace(1, 1, N)
        vals[1*N:2*N, 1] = np.linspace(200/256, 1, N)
        vals[1*N:2*N, 2] = np.linspace(0, 150/256, N)
        vals[2*N:, 0] = np.linspace(1, 1, N)
        vals[2*N:, 1] = np.linspace(1, 1, N)
        vals[2*N:, 2] = np.linspace(150/256, 1, N)
        cmap = ListedColormap(vals) 
 
        cont = ax.contourf(X, Y, V, C, cmap=cmap)
        cont2 = ax.contour(X, Y, V, C, colors=[(0,0,0,0.3)]*tasot,linewidths=0.5,
            levels=C)
        cbar = fig.colorbar(cont)
        cbar.add_lines(cont2)
        cbar.ax.set_ylabel('V (V)')

    # piirretään varaukset värillisinä ympyröinä
    for q, xy in varaukset:        
        if q > 0:
            vara = plt.Circle(xy, RP, fill=True, color='red')
        elif q < 0:
            vara = plt.Circle(xy, RP, fill=True, color='blue')
        ax.add_patch(vara)
        
    # piirretään tarkastelupisteet ympyröinä
    for xy in pisteet:
        pt = plt.Circle(xy, RP, fill=False, color='black')
        ax.add_patch(pt)
        
    # piirretään sähkökenttä nuolina
    for rP, eP in E:
        
        ekentta = FancyArrowPatch(rP-0.5*eP, rP+0.5*eP, 
                    arrowstyle='->', color=(0,0,0.5,1), mutation_scale=2.5,
                    shrinkA=0.2,shrinkB=0.2)
        ax.add_patch(ekentta)

    # tallennetaan kuva
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect('equal')
    plt.grid(False)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.savefig(f"kentta."+KUVAFORMAATTI)
    plt.close('all')

# funktio piirra_systeemi loppuu


# Piirretään joukko kuvaajia ja tallennetaan ne tiedostoihin.
def piirra_kuvaaja(varaukset,x0,y0,x1,y1,n_points=501,
                   piirra_E=PIIRRA_E,piirra_V=PIIRRA_V):

    plt.close('all')
    fig, ax = plt.subplots()
    
    xs = np.linspace(x0,x1,n_points)
    ys = np.linspace(y0,y1,n_points)    
    dr = np.sqrt( (x1-x0)**2 + (y1-y0)**2 )
    rs = np.linspace(0, dr, n_points)
    
    exs = []
    eys = []
    es = []
    vs = []
    
    for x,y in zip(xs, ys):
        if SAHKOKENTTA_POTENTIAALISTA:
            E = par.sahkokentta_potentiaalista( x, y )    
        else:
            E = par.sahkokentta(varaukset, x, y)
        exs += [ E[0] ]
        eys += [ E[1] ]
        es += [ pituus(E) ]
        vs += [ potentiaali(x, y) ]
    es = np.array(es)
    vs = np.array(vs)

    emax = np.max(es)
    vmax = np.max(vs)
    if emax == 0:
        emax = 1
    if vmax == 0:
        vmax = 1

    plt.xlabel("r (m)")
    plt.ylabel("Ex (V/m)")
    plt.plot(rs, exs, c='blue', lw=2)
    plt.grid()
    
    plt.xlim(0, dr)
    plt.ylim(-1.1*emax, 1.1*emax)
    
    ax.text(0.02*dr, 1.0*emax, f'({x0:4.2f} m, {y0:4.2f} m)', 
                    fontsize=12, color='black')
    ax.text(0.70*dr, 1.0*emax, f'({x1:4.2f} m, {y1:4.2f} m)', 
                    fontsize=12, color='black')
    
    if piirra_E:
        plt.savefig("Ex_xy."+KUVAFORMAATTI)
    plt.close('all')
    fig, ax = plt.subplots()   
    
    plt.xlabel("r (m)")
    plt.ylabel("Ey (V/m)")
    plt.plot(rs, eys, c='blue', lw=2)
    plt.grid()
    
    plt.xlim(0, dr)
    plt.ylim(-1.1*emax, 1.1*emax)
    
    ax.text(0.02*dr, 1.0*emax, f'({x0:4.2f} m, {y0:4.2f} m)', 
                    fontsize=12, color='black')
    ax.text(0.70*dr, 1.0*emax, f'({x1:4.2f} m, {y1:4.2f} m)', 
                    fontsize=12, color='black')
    
    if piirra_E:
        plt.savefig("Ey_xy."+KUVAFORMAATTI)
    plt.close('all')
    fig, ax = plt.subplots()   
    
    plt.xlabel("r (m)")
    plt.ylabel("|E| (V/m)")
    plt.plot(rs, es, c='blue', lw=2)
    plt.grid()
    
    plt.xlim(0, dr)
    plt.ylim(0, 1.1*emax)
    
    ax.text(0.02*dr, 1.02*emax, f'({x0:4.2f} m, {y0:4.2f} m)', 
                    fontsize=12, color='black')
    ax.text(0.70*dr, 1.02*emax, f'({x1:4.2f} m, {y1:4.2f} m)', 
                    fontsize=12, color='black')
    
    if piirra_E:
        plt.savefig("E_xy."+KUVAFORMAATTI)
    plt.close('all')
    fig, ax = plt.subplots()    
    
    plt.xlabel("r (m)")
    plt.ylabel("V (V)")
    plt.plot(rs, vs, c='orange', lw=2)
    plt.grid()
    
    plt.xlim(0, dr)
    plt.ylim(0, 1.1*vmax)
    
    ax.text(0.02*dr, 1.02*vmax, f'({x0:4.2f} m, {y0:4.2f} m)', 
                    fontsize=12, color='black')
    ax.text(0.70*dr, 1.02*vmax, f'({x1:4.2f} m, {y1:4.2f} m)', 
                    fontsize=12, color='black')
    
    if piirra_V:
        plt.savefig("V_xy."+KUVAFORMAATTI)
    plt.close('all')

# funktio piirra_kuvaajat loppuu


# laskee sähkökentän tarkastelupisteissä ja kertoo tulokset
def laske_sahkokentta_pisteissa(varaukset, pisteet):
    print()
    print("Sähkökenttä tarkastelupisteissä:")
    for xy in pisteet:
        if SAHKOKENTTA_POTENTIAALISTA:
            E = par.sahkokentta_potentiaalista( xy[0], xy[1] )    
        else:
            E = par.sahkokentta( varaukset, xy[0], xy[1] )
        print( f"(x,y) = ({xy[0]:7.3f} m, {xy[1]:7.3f} m): E = ({E[0]:9.3f} V/m, {E[1]:9.3f} V/m)" )
    print()

# funktio laske_sahkokentta_pisteissa loppuu


# laskee potentiaalin tarkastelupisteissä ja kertoo tulokset
def laske_potentiaali_pisteissa(pisteet):
    print()
    print("Potentiaali tarkastelupisteissä:")
    for xy in pisteet:
        V = potentiaali( xy[0], xy[1] )
        print( f"(x,y) = ({xy[0]:7.3f} m, {xy[1]:7.3f} m): V = {V:9.3f} V" )
    print()

# funktio laske_potentiaali_pisteissa loppuu


# Pääohjelma.
def main():
    if PIIRRA_E:
        laske_sahkokentta_pisteissa(Q, PISTEET)
    if PIIRRA_V:
        laske_potentiaali_pisteissa(PISTEET)
    piirra_kuvaaja(Q, JANA[0,0], JANA[0,1], JANA[1,0], JANA[1,1])
    piirra_systeemi(Q, PISTEET)
    
# funktio main loppuu


# Ajetaan pääohjelma.
if __name__ == '__main__':
    main()

