#
# Simuloi vaimennetun värähtelijän liikettä.
#
# Teemu Hynninen 2025
#
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.animation import FuncAnimation

#
# Haetaan tiedot tiedostosta parametrit.py
#
import parametrit as par
M = par.M
K = par.K
B = par.B
X0 = par.DX_ALKU
VX0 = par.VX_ALKU
F = np.abs(par.VOIMA)
OMEGA = 2.0 * np.pi * par.TAAJUUS

T_LOPPU = par.T_LOPPU         
ANIMAATIO_DT = par.ANIMAATIO_DT    
ANIMAATIO_LOPPU = par.ANIMAATIO_LOPPU
ANIMAATIONOPEUS = par.ANIMAATIONOPEUS
KUVAFORMAATTI = par.KUVAFORMAATTI    
TALLENNA_ANIMAATIO = par.TALLENNA_ANIMAATIO

DT = 0.001 # aika-askeleen pituus (pieni = tarkka mutta hidas simulaatio)
# lähellä resonanssia tarvitaan lisää tarkkuutta
if F > 0.0 and np.abs( np.sqrt( K/M )/( OMEGA ) - 1 ) < 0.01:
    DT = 0.0002
if ANIMAATIO_DT < DT:
    DT = ANIMAATIO_DT       
Y_MIN = -1.0     # koordinaatiston y-akselin minimi
Y_MAX =  1.0     # koordinaatiston y-akslein maksimi
X_MIN = -1.0     # koordinaatiston y-akselin minimi
X_MAX =  1.0     # koordinaatiston y-akslein maksimi
V_SKAALA = 1.0   # nopeusvektoreiden skaalaus
A_SKAALA = 1.0   # kiihtyvyysvektoreiden skaalaus
R = 0.1

# johdetut parametrit
N_ASKEL = int(T_LOPPU / DT)
ANIMAATIOASKEL = int(ANIMAATIO_DT / DT)
N_ANIMAATIO = int(ANIMAATIO_LOPPU / ANIMAATIO_DT )

# tallennetaan data kuvien piirtämistä varten
t_historia = [  ]
x_historia = [  ]
v_historia = [  ]
a_historia = [  ]




# Piirretään animaatio systeemin liikkeestä.
def piirra_animaatio():

    plt.close('all')
    fig, ax = plt.subplots()
    
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    
    n_frames = np.min( [ len(x_historia), N_ANIMAATIO ] )

    # luodaan listat, johon ympyrät ja nuolet tallennetaan
    circle_patches = []
    arrow_patches = []
    time_text = ax.text(X_MIN + 0.05*(X_MAX - X_MIN), Y_MAX - 0.1*(Y_MAX - Y_MIN), '', 
                    fontsize=12, color='black')

    
    paikka = Circle((0, 0), R, fill=False, color='black')
    ax.add_patch(paikka)
    circle_patches.append(paikka)

    nopeus = FancyArrowPatch((0, 0), (0, 0), 
                arrowstyle='->', color='red', mutation_scale=10)
    ax.add_patch(nopeus)
    arrow_patches.append(nopeus)
    
    kiihtyvyys = FancyArrowPatch((0, 0), (0, 0), 
                arrowstyle='->', color='orange', mutation_scale=10)
    ax.add_patch(kiihtyvyys)
    arrow_patches.append(kiihtyvyys)
    

    def alusta():
        #for circle in circle_patches:
        #    circle.set_center((0, 0))
        #for arrow in arrow_patches:
        #    arrow.set_positions((0, 0), (0, 0))
        time_text.set_text('')
        return circle_patches + arrow_patches

    # siirretään ympyrät ja nuolet oikeille paikoille
    def paivita(frame):

        t = t_historia[frame]
        x = x_historia[frame]
        v = v_historia[frame]
        a = a_historia[frame]
        
        time_text.set_text(f't = {t:.2f} s')
        
        paikka = np.array( x )
        nopeus = np.array( V_SKAALA*v )
        kiihtyvyys = np.array( A_SKAALA*a )
        
        circle_patches[0].center = paikka
        circle_patches[0].set_center(paikka)

        arrow_patches[0].set_positions(paikka, paikka+nopeus)
        arrow_patches[1].set_positions(paikka, paikka+kiihtyvyys)


        return circle_patches + arrow_patches + [time_text]

    # luodaan animaatio
    ani = FuncAnimation(fig, paivita, frames=n_frames, init_func=alusta, blit=True, interval=int(1000*ANIMAATIO_DT/ANIMAATIONOPEUS))
    if TALLENNA_ANIMAATIO:
        print("tallennetaan animaatio")
        ani.save("animaatio.mp4")
    else:
        plt.show()
    plt.close(fig)

# funktio piirra_animaatio loppuu



# Piirretään joukko kuvaajia ja tallennetaan ne tiedostoihin.
def piirra_kuvaajat():

    xs = np.array(x_historia)
    vs = np.array(v_historia)
    aa = np.array(a_historia)
    ts = np.array(t_historia)

    plt.close('all')

    plt.xlabel("t (s)")
    plt.ylabel("x (m)")
    plt.plot(ts,xs[:,1],c='black', lw=1)
    plt.grid()
    plt.savefig("x_t."+KUVAFORMAATTI)
    plt.clf()

    plt.xlabel("t (s)")
    plt.ylabel("$v_x$ (m/s)")
    plt.plot(ts,vs[:,1],c='red', lw=1)
    plt.grid()
    plt.savefig("vx_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("$a_x$ (m/s$^2$)")
    plt.plot(ts,aa[:,1],c='orange', lw=1)
    plt.grid()
    plt.savefig("ax_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.close('all')
    
# funktio piirra_kuvaajat loppuu


# Näytetään väliaikatietoja simulaation etenemisestä.
def edistyminen(n, n_max, kerrat=10):
    
    if n % (n_max // kerrat) == 0 :  
        p = int(np.round((100*n)/n_max,0))   
        print(f"simuloitu {p:3d} %")

# funktio edistyminen loppuu


# Laskee kiihtyvyyden.
def kiihtyvyys(x, v, t):

    return ( -K*x-B*v + F*np.cos(OMEGA*t)*np.array([0.0, 1.0]) ) / M

# funktio kiihtyvyys loppuu


# Simuloidaan ajan kulumista eli annetaan kappaleiden liikkua.
def simuloi(x, v, t = 0, askeleita = N_ASKEL):

    a = kiihtyvyys(x, v, t)

    # alustetaan nopeus simulaatiota varten
    v += (0.5 * DT * a)
    
    # otetaan haluttu määrä pieniä askeleita ajassa eteenpäin
    for askel in range(askeleita+1):
                
        # tallennetaan tiedot tasaisin välein
        if askel % ANIMAATIOASKEL == 0:
            t_historia.append( t )
            x_historia.append( copy.copy(x) )
            v_historia.append( v - (0.5 * a * DT) )
            a_historia.append( copy.copy(a) )
            
            
        # päivitetään aika sekä kappaleiden paikat ja nopeudet
        t = np.round(t+DT,5)
        x += v * DT    
        a = kiihtyvyys(x, v, t)
        v = ( ( 1 - B/(2*M)*DT )*v + kiihtyvyys(x, np.zeros(2), t)*DT ) / ( 1 + B/(2*M)*DT )
        
        # kerrotaan käyttäjälle, miten menee
        edistyminen( askel, askeleita+1 )
        
    # viimeistellään nopeus simulaation lopuksi
    v -= (0.5 * DT * a)

    return t
    
# funktio simuloi loppuu
    

# laskee, kuinka suurena animaatio piirretään
def skaalaa(x,v,a):
    global R, X_MAX, X_MIN, Y_MAX, Y_MIN, V_SKAALA, A_SKAALA
    maksimi = np.max(np.abs(x))
    vmaksimi = np.max(np.abs(v))
    amaksimi = np.max(np.abs(a))
    
    X_MAX =  1.2*maksimi
    X_MIN = -1.2*maksimi
    Y_MAX =  1.2*maksimi
    Y_MIN = -1.2*maksimi
    R     =  0.1*maksimi

    V_SKAALA = 0.5 * maksimi / vmaksimi
    A_SKAALA = 0.5 * maksimi / amaksimi

# funktio skaalaa loppuu


# raportoi amplitudin simulaation alussa ja lopussa
def kerro_amplitudi(x, dt=2.0):

    xa = np.array(x)
    N = np.min( [ int( dt/ANIMAATIO_DT+0.01 ), len(x) ] )
    
    A_alku = 0.5 * ( np.max( xa[:N,1] ) - np.min( xa[:N,1] ) )
    A_loppu = 0.5 * ( np.max( xa[-N:,1] ) - np.min( xa[-N:,1] ) )
    print( f"Amplitudi ensimmäisen {dt:5.1f} s aikana: {A_alku:8.3} m" )
    print( f"Amplitudi viimeisen   {dt:5.1f} s aikana: {A_loppu:8.3} m" )

# funktio kerro_amplitudi loppuu


# Pääohjelma.
def main():

    simuloi( np.array([0.0, X0]) , np.array([0.0, VX0]) )
    skaalaa( x_historia,v_historia,a_historia )
    piirra_kuvaajat()
    kerro_amplitudi(x_historia)
    piirra_animaatio()

# funktio main loppuu


# Ajetaan pääohjelma.
if __name__ == '__main__':
    main()

