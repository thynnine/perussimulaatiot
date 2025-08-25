#
# Simuloi liikettä xy-tasossa.
#
# Käyttäjä määrittelee parameterit.py-tiedostossa
# kappaleen kiihtyvyyden, ja tämä ohjelma laskee
# annettua kiihtyvyyttä vastaavan liikeradan.
# Tasossa on esteitä ja maaleja, ja ohjelma pitää
# kirjaa, osutaanko niihin.
#
# Teemu Hynninen 2025
#

import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.animation import FuncAnimation

#
# Haetaan lentotiedot tiedostosta parametrit.py
#
import parametrit as par
A = par.A
T_LOPPU = par.T_LOPPU         
ANIMAATIO_DT = par.ANIMAATIO_DT    
KUVA_DT = par.KUVA_DT
KUVAFORMAATTI = par.KUVAFORMAATTI    
TALLENNA_ANIMAATIO = par.TALLENNA_ANIMAATIO


R_ALUS   = 10.0   # aluksen säde
R_ESTE   = 40.0   # esteen säde
R_MAALI  = 10.0   # maalin säde

# esteiden keskipisteiden koordinaatit
ESTEET = [
(  -70.0, -340.0 ),
( -220.0, -180.0 ),
( -350.0,  -40.0 ),
(   50.0, -230.0 ),
(  270.0, -110.0 ),
(  180.0, -360.0 ),
(  -40.0,  360.0 ),
( -210.0,   40.0 ),
( -250.0,  260.0 ),
(   20.0,  100.0 ),
(  120.0,   30.0 ),
(  320.0,  260.0 )
]

# maalien keskipisteiden koordinaatit
MAALIT = [
( -210.0, -240.0 ),
( -300.0,  -60.0 ),
(   50.0, -330.0 ),
(  250.0, -180.0 ),
( -330.0,  350.0 ),
(  -50.0,  200.0 ),
(  200.0,   50.0 ),
(  350.0,  350.0 )
]

DT = 0.01         # aika-askeleen pituus (pieni = tarkka mutta hidas simulaatio)
Y_MIN = -500.0    # koordinaatiston y-akselin minimi
Y_MAX = 500.0     # koordinaatiston y-akslein maksimi
X_MIN = -500.0    # koordinaatiston y-akselin minimi
X_MAX = 500.0     # koordinaatiston y-akslein maksimi
V_SKAALA = 4.0    # nopeusvektoreiden skaalaus
A_SKAALA = 20.0   # kiihtyvyysvektoreiden skaalaus
A_MAX = 10.0001   # kiihtyvyyden sallittu maksimi

# johdetut parametrit
N_ASKEL = int(T_LOPPU / DT)
ANIMAATIOASKEL = int(ANIMAATIO_DT / DT)
KUVAASKEL = int(KUVA_DT / DT)

# aluksen ominaisuudet
x = np.zeros( 2 )
v = np.zeros( 2 )
a = np.zeros( 2 )
keratty = np.zeros( len(MAALIT) , dtype='int')

# tallennetaan data kuvien piirtämistä varten
t_historia = [  ]
x_historia = [  ]
v_historia = [  ]
a_historia = [  ]
k_historia = [  ]



# Piirretään kuva systeemistä yhdellä ajan hetkellä.
def piirra_systeemi(x, v, a, t, rata=True):

    plt.close('all')
    fig, ax = plt.subplots()

    paikka = plt.Circle(x, R_ALUS, fill=False, color='black')
    ax.add_patch(paikka)
    nopeus = FancyArrowPatch(x, x+V_SKAALA*v, 
            arrowstyle='->', color='red', mutation_scale=10)
    ax.add_patch(nopeus)
    kiihtyvyys = FancyArrowPatch(x, x+A_SKAALA*a, 
            arrowstyle='->', color='orange', mutation_scale=10)
    ax.add_patch(kiihtyvyys)
    k_vari = ['blue', 'gray']

    for es in ESTEET:
        este = plt.Circle(es, R_ESTE, fill=False, color='red')
        ax.add_patch(este)
        
    for i in range(len(MAALIT)):
        ma = MAALIT[i]
        maali = plt.Circle(ma, R_MAALI, fill=False, color=k_vari[ keratty[i] ])
        ax.add_patch(maali)
    
    if rata:
        ratadata = np.array(x_historia)
        xs = ratadata[:,0]
        ys = ratadata[:,1]
        plt.plot(xs,ys,c='black',lw=1,ls=':')
    
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect('equal')
    plt.grid(False)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title(f't = {t:.2f} s')
    plt.savefig(f"systeemi_t{t:.2f}."+KUVAFORMAATTI)
    plt.close('all')

# funktio piirra_systeemi loppuu


# Piirretään animaatio systeemin liikkeestä.
def piirra_animaatio():

    plt.close('all')
    fig, ax = plt.subplots()
    
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    
    n_frames = len(x_historia)

    # luodaan listat, johon ympyrät ja nuolet tallennetaan
    circle_patches = []
    arrow_patches = []
    time_text = ax.text(X_MIN + 0.05*(X_MAX - X_MIN), Y_MAX - 0.1*(Y_MAX - Y_MIN), '', 
                    fontsize=12, color='black')

    
    paikka = Circle((0, 0), R_ALUS, fill=False, color='black')
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
    
    for es in ESTEET:
        este = plt.Circle(es, R_ESTE, fill=False, color='red')
        ax.add_patch(este)
        circle_patches.append(este)
        
    for ma in MAALIT:
        maali = plt.Circle(ma, R_MAALI, fill=False, color='blue')
        ax.add_patch(maali)
        circle_patches.append(maali)
    

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
        k = k_historia[frame]
        
        k_vari = ['blue', 'gray']
        
        time_text.set_text(f't = {t:.2f} s')

        
        paikka = np.array( x )
        nopeus = np.array( V_SKAALA*v )
        kiihtyvyys = np.array( A_SKAALA*a )
        
        circle_patches[0].center = paikka
        circle_patches[0].set_center(paikka)
        
        n_este = len(ESTEET)
        n_maali = len(MAALIT)
        for i in range(n_maali):
            circle_patches[i+n_este+1].set_color( k_vari[ k[i] ] )

        arrow_patches[0].set_positions(paikka, paikka+nopeus)
        arrow_patches[1].set_positions(paikka, paikka+kiihtyvyys)


        return circle_patches + arrow_patches + [time_text]

    # luodaan animaatio
    ani = FuncAnimation(fig, paivita, frames=n_frames, init_func=alusta, blit=True, interval=int(200*ANIMAATIO_DT))
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
    plt.plot(ts,xs[:,0],c='black', lw=2)
    plt.savefig("x_t."+KUVAFORMAATTI)
    plt.clf()

    plt.xlabel("t (s)")
    plt.ylabel("y (m)")
    plt.plot(ts,xs[:,1],c='black', lw=2)
    plt.savefig("y_t."+KUVAFORMAATTI)
    plt.clf()
    
    xmax = X_MAX
    xmin = X_MIN
    ymin = Y_MIN
    ymax = Y_MAX
    fig, ax = plt.subplots()
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    plt.plot(xs[:,0],xs[:,1],c='black', lw=2)
    plt.savefig("x_y."+KUVAFORMAATTI)
    plt.clf()

    plt.xlabel("t (s)")
    plt.ylabel("$v_x$ (m/s)")
    plt.plot(ts,vs[:,0],c='red', lw=2)
    plt.savefig("vx_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("$v_y$ (m/s)")
    plt.plot(ts,vs[:,1],c='red', lw=2)
    plt.savefig("vy_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("$v$ (m/s)")
    plt.plot(ts,np.sqrt(vs[:,0]**2+vs[:,1]**2),c='red', lw=2)
    plt.savefig("v_t."+KUVAFORMAATTI)
    plt.clf()

    plt.xlabel("t (s)")
    plt.ylabel("$a_x$ (m/s$^2$)")
    plt.plot(ts,aa[:,0],c='orange', lw=2)
    plt.savefig("ax_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("$a_y$ (m/s$^2$)")
    plt.plot(ts,aa[:,1],c='orange', lw=2)
    plt.savefig("ay_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("$a$ (m/s$^2$)")
    plt.plot(ts,np.sqrt(aa[:,0]**2+aa[:,1]**2),c='orange', lw=2)
    plt.ylim(0,11)
    plt.savefig("a_t."+KUVAFORMAATTI)
    plt.clf()

# funktio piirra_kuvaajat loppuu


# Näytetään väliaikatietoja simulaation etenemisestä.
def edistyminen(n, n_max, kerrat=10):
    
    if n % (n_max // kerrat) == 0 :  
        p = int(np.round((100*n)/n_max,0))   
        print(f"simuloitu {p:3d} %")

# funktio edistyminen loppuu



# Lasketaan kiihtyvyys.
def kiihtyvyys(t, A):

    a = np.zeros( 2 )

    for i in range(len(A)):
        if t < A[i][1]:
            a = np.array( A[i][0] )
            break
    
    return a
   
# funktio kiihtyvyys loppuu


 
# Tarkastetaan, ollaanko törmätty esteeseen.
def tormasi(x):
    for es in ESTEET:
        rvec = x - es
        r = np.sqrt( rvec@rvec )
        
        if r < R_ALUS+R_ESTE:
            return True
    
    return False
    
# funktio tormasi loppuu


# Tarkastetaan, ollaanko päästy maaliin.
def maalissa(x):
    for i in range(len(MAALIT)):
        rvec = x - MAALIT[i]
        r = np.sqrt( rvec@rvec )
        
        if r < R_ALUS+R_MAALI:
            keratty[i] = 1
            break
    
    if (keratty == 1).all():
        return True
    
    return False
  
# funktio maalissa loppuu  


# Lasketaan kuljetun matkan pituus
def matka():

    l = 0.0
    for i in range(len(x_historia)-1):
        dr = x_historia[i]-x_historia[i+1]
        l += np.sqrt( dr @ dr )
        
    return l

# funktio matka loppuu


# Kerrotaan simulaation lopuksi, miten meni.
def lopeta(viesti, t):

    print()
    print(viesti+f", kun t = {t:0.2f} s.")
    l = matka()
    print(f"Kuljettu matka oli {l:4.1f} m")
    print(f"Keskivauhti oli {(l/t):3.1f} m/s")
    print("Osuit "+str(np.sum(keratty))+" / "+str(len(keratty))+" maaliin.")
    print()
 
# funktio lopeta loppuu   


# Simuloidaan ajan kulumista eli annetaan kappaleiden liikkua.
def simuloi(x, v, t = 0, askeleita = N_ASKEL):

    a = kiihtyvyys(t, A)

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
            k_historia.append( copy.copy(keratty) )
            
        # piirretään kuva tasaisin välein
        if askel % KUVAASKEL == 0:
            piirra_systeemi(x, v - (0.5 * a * DT), a, t)
            
        # päivitetään aika sekä kappaleiden paikat ja nopeudet
        t = np.round(t+DT,5)
        x += v * DT    
        a = kiihtyvyys(t, A)
        v += (a * DT)
        
        # kerrotaan käyttäjälle, miten menee
        edistyminen( askel, askeleita+1 )
        
        if tormasi(x):
            lopeta("Törmäsit esteeseen", t)
            piirra_systeemi(x, v - (0.5 * a * DT), a, t)
            return t
            
        if maalissa(x):
            lopeta("Pääsit maaliin", t)
            piirra_systeemi(x, v - (0.5 * a * DT), a, t)
            return t
        
        # testataan ettei kiihtyvyys ole liian suuri
        a_abs = np.sqrt(a@a)
        if a_abs > A_MAX:
            lopeta("Kiihtyvyys kasvoi liian suureksi (a = {a_abs:.2f} m/s2)", t)
            piirra_systeemi(x, v - (0.5 * a * DT), a, t)
            return t
        
    # viimeistellään nopeus simulaation lopuksi
    v -= (0.5 * DT * a)

    lopeta("Aika loppui", t)
    return t
    
# funktio simuloi loppuu
    


# Pääohjelma.
def main():

    simuloi(x, v)
    piirra_kuvaajat()
    piirra_animaatio()

# funktio main loppuu


# Ajetaan pääohjelma.
if __name__ == '__main__':
    main()

