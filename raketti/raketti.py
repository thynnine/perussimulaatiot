#
# Pienen raketin simulaatio
#
# Teemu Hynninen 2025
#
# Tässä ykinkertaisessa simulaatiossa on yksi suuri ja joukko pienempiä
# kappaleita, jotka ovat aluksi toisissaan kiinni.
# Suuren kappaleen voi asettaa ampumaan pieniä kappaleita yksi
# kerrallaan alaspäin, mikä matkii rakettimoottorin toimintaa.
# Kaikille kappaleille voi myös antaa samanaikaisesti satunnaiset
# impulssit, mikä matkii raketin räjähtämistä.
#

#
# SIMULAATIOPARAMETRIT
#
import parametrit

M_ISO    = parametrit.M_ISO  # suuren kappaleen massa, kg
M_PIENI  = parametrit.M_PIENI  # pienien kappaleiden massa, kg
N_PIENET = parametrit.N_PIENET  # pienien kappaleiden lukumäärä

YKSIULOTTEINEN = parametrit.YKSIULOTTEINEN # onko simulaatio yksiulotteinen? (jos ei, se on kaksiulotteinen)
VY_ALKU = parametrit.VY_ALKU # alkunopeuden x-komponentti, m/s (vain 2D)
VX_ALKU = parametrit.VX_ALKU # alkunopeuden y-komponentti, m/s
G = parametrit.G         # putoamiskiihtyvyys, m/s2 (0 = ei painovoimaa, 9.8 = painovoima alaspäin)

T_LOPPU = parametrit.T_LOPPU  # simulaation kesto, s

T_RAJAHDYS  = parametrit.T_RAJAHDYS  # räjähdyksen tapahtumishetki, s (negatiivinen luku poistaa räjähdyksen)
DV_RAJAHDYS = parametrit.DV_RAJAHDYS  # räjähdyksen voimakkuus: pienten kappaleiden nopeuden muutosten keskihajonta, m/s

T_KIIHDYTYS   =  parametrit.T_KIIHDYTYS  # suihkutuksen kesto, s (negatiivinen luku poistaa suihkutuksen)
DVX_KIIHDYTYS =  parametrit.DVX_KIIHDYTYS  # suihkutettavien hiukkasten nopeuden muutos x-suunnassa, m/s (positiivinen = vasemmalle, vain 2D)
DVY_KIIHDYTYS =  parametrit.DVY_KIIHDYTYS  # suihkutettavien hiukkasten nopeuden muutos y-suunnassa, m/s (positiivinen = alaspäin)

#
# Piirtämistä ohjaavat parametrit - näitäkin voi muuttaa.
#

ANIMAATIO_DT = parametrit.ANIMAATIO_DT        # tulosten tallennusväli, s (pieni = tarkat kuvaajat ja animaatio)
KUVA_DT = parametrit.KUVA_DT              # kuvien piirtoväli, s (pieni = paljon kuvia lyhyin aikavälein)
KUVAFORMAATTI = parametrit.KUVAFORMAATTI      # tallennettavien kuvien formaatti (png tai pdf)
TALLENNA_ANIMAATIO = parametrit.TALLENNA_ANIMAATIO # tallennetaanko animaatio tiedostoon? (jos ei, näytetään erillisessä ikkunassa)

KAIKKI_NOPEUDET = parametrit.KAIKKI_NOPEUDET # piirretäänkö myös pienten kappaleiden nopeusvektorit?
KOKONAISSUUREET = parametrit.KOKONAISSUUREET   # piirretäänkö massakeskipiste sekä kokonaisliikemäärä ja -energia?

X_CENTER = parametrit.X_CENTER # koordinaatiston keskipiste x-akselilla
Y_MIN    = parametrit.Y_MIN     # koordinaatiston y-akselin minimi
Y_MAX    = parametrit.Y_MAX    # koordinaatiston y-akslein maksimi
V_SKAALA = parametrit.V_SKAALA  # nopeusvektoreiden skaalaus

#
# SIMULAATIOPARAMETRIT LOPPUVAT TÄHÄN 
#




# aika-askeleen pituus (pieni = tarkka mutta hidas simulaatio)
DT = 0.01    
# suihkutuksen nopeus: aikaväli suihkutettavien kappaleiden irtoamisten välillä, s
DT_KIIHDYTYS  = 0.05  
  

# kappaleiden koko kuvissa
R_ISO = 1.0   # ison pallon säde kuvissa
R_PIENI = 0.4 # pienen pallon säde kuvissa

# johdetut parametrit
N_OSAT = N_PIENET + 1
N_ASKEL = int(T_LOPPU / DT)
RAJAHDYSASKEL = int(T_RAJAHDYS / DT)
KIIHDYTYSASKEL = int(DT_KIIHDYTYS / DT)
KIIHDYTYSLOPPUASKEL = int(T_KIIHDYTYS / DT)
ANIMAATIOASKEL = int(ANIMAATIO_DT / DT)
KUVAASKEL = int(KUVA_DT / DT)


# ladataan kirjastoja
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.animation import FuncAnimation


# kappaleiden ominaisuudet taulukoituna
m = np.array( [M_ISO] + [M_PIENI] * N_PIENET )
x = np.zeros( [N_OSAT, 2] )
v = np.array( [ [VX_ALKU, VY_ALKU] ] * N_OSAT )
n_jaljella = N_PIENET

# tallennetaan data kuvien piirtämistä varten
t_historia = [  ]
x_historia = [  ]
v_historia = [  ]



# Lasketaan massakeskipisteen koordinaatit annetun vektorin koordinaateista.
def laske_massakeskipiste(x, massa = m):

    x_mkp = np.zeros(2)
    
    for i in range(N_OSAT):
        x_mkp += massa[i] * x[i]
        
    x_mkp /= np.sum(massa)
    
    return x_mkp

# funktio laske_massakeskipiste loppuu   



# Piirretään kuva systeemistä yhdellä ajan hetkellä.
def piirra_systeemi(x, v, t):

    plt.close('all')
    cs = []
    vs = []
    for i in range(len(x[:,0])):
        cs += [ x[i,:] ]
        if KAIKKI_NOPEUDET or i == 0:
            vs += [ v[i,:] ]
        else:
            vs += [ [ 0, 0 ] ]
            
    cs = np.array(cs)
    vs = np.array(vs)
    rs = np.array( [R_ISO] + [R_PIENI] * N_PIENET )

    plt.clf()
    fig, ax = plt.subplots()

    for (keskus, sade, nuoli) in zip(cs, rs, vs):
        ympyra = plt.Circle(keskus, sade, fill=False, color='black')
        ax.add_patch(ympyra)
        #ax.arrow(keskus[0], keskus[1], V_SKAALA*nuoli[0], V_SKAALA*nuoli[1],
        #     head_width=0.1, head_length=0.15, fc='red', ec='red')

        arrow = FancyArrowPatch(keskus, keskus+V_SKAALA*nuoli, 
                            arrowstyle='->', color='red', mutation_scale=10)
        ax.add_patch(arrow)
        
    if KOKONAISSUUREET:
        mkp = plt.Circle( laske_massakeskipiste(x), R_PIENI, fill=True, color='brown' )
        ax.add_patch( mkp )

                    
    xmax = 0.5*(Y_MAX-Y_MIN)+X_CENTER
    xmin = -0.5*(Y_MAX-Y_MIN)+X_CENTER
    ymin = Y_MIN
    ymax = Y_MAX
    
    time_text = ax.text(xmin + 0.05*(xmax - xmin), ymax - 0.1*(ymax - ymin), f't = {t:.2f} s', 
                    fontsize=12, color='black')
                    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    plt.grid(False)
    plt.xlabel('')
    plt.ylabel('x (m)')
    #plt.title(f't = {t:.2f} s')
    plt.savefig(f"raketti_t{t:.2f}."+KUVAFORMAATTI)
    plt.close('all')

# funktio piirra_systeemi loppuu



# Piirretään animaatio systeemin liikkeestä.
def piirra_animaatio():

    plt.close('all')
    fig, ax = plt.subplots()
    
    xmax = 0.5*(Y_MAX-Y_MIN)+X_CENTER
    xmin = -0.5*(Y_MAX-Y_MIN)+X_CENTER
    ymin = Y_MIN
    ymax = Y_MAX
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    
    n_frames = len(x_historia)
    rs = np.array( [R_ISO] + [R_PIENI] * N_PIENET )

    # luodaan animaatiossa tarvittavat ympyrät ja nuolet
    circle_patches = []
    arrow_patches = []
    time_text = ax.text(xmin + 0.05*(xmax - xmin), ymax - 0.1*(ymax - ymin), '', 
                    fontsize=12, color='black')

    for i in range(N_OSAT):
        circle = Circle((0, 0), rs[i], fill=False, color='black')
        ax.add_patch(circle)
        circle_patches.append(circle)

        arrow = FancyArrowPatch((0, 0), (0, 0), 
                            arrowstyle='->', color='red', mutation_scale=10)
        ax.add_patch(arrow)
        arrow_patches.append(arrow)
    
    if KOKONAISSUUREET:
        mkp = plt.Circle( (0,0), R_PIENI, fill=True, color='brown' )
        ax.add_patch( mkp )
        circle_patches.append(mkp)

    # asetataan graafisten olioiden ominaisuudet nollaan
    def alusta():
        for circle in circle_patches:
            circle.set_center((0, 0))
        for arrow in arrow_patches:
            arrow.set_positions((0, 0), (0, 0))
        time_text.set_text('')
        return circle_patches + arrow_patches

    # siirretään ympyrät ja nuolet oikeille paikoille
    def paivita(frame):

        t = t_historia[frame]
        x = x_historia[frame]
        v = v_historia[frame]
        
        time_text.set_text(f't = {t:.2f} s')

        for i in range(N_OSAT):
            keskus = np.array( x[i,:] )
            if KAIKKI_NOPEUDET or i == 0:
                nuoli = np.array( V_SKAALA*v[i,:] )
            else:
                nuoli = np.array( [ 0, 0 ] )
            circle_patches[i].center = keskus

            # Update circle
            circle_patches[i].set_center(keskus)

            # Update arrow
            arrow_patches[i].set_positions(keskus, keskus+nuoli)

        if KOKONAISSUUREET:
            circle_patches[N_OSAT].set_center( laske_massakeskipiste(x) )

        return circle_patches + arrow_patches + [time_text]

    # luodaan animaatio
    ani = FuncAnimation(fig, paivita, frames=n_frames, init_func=alusta, blit=True, interval=int(1000*ANIMAATIO_DT))
    if TALLENNA_ANIMAATIO:
        print("tallennetaan animaatio")
        ani.save("animaatio.mp4")
    else:
        plt.show()
    plt.close(fig)

# funktio piirra_animaatio loppuu



# Piirretään joukko kuvaajia ja tallennetaan ne tiedostoihin.
def piirra_kuvaajat():

    plt.close('all')
    xs = np.array(x_historia)
    vs = np.array(v_historia)
    ts = np.array(t_historia)

    varit = [ 'red' ] + (N_PIENET*[ 'black' ])
    paksuudet = [ 3 ] + (N_PIENET*[ 1 ])
    nimet = ['suuri', 'pieni'] + N_PIENET*['']

    plt.close('all')

    plt.xlabel("t (s)")
    plt.ylabel("x (m)")
    for i in range(N_OSAT-1,-1,-1):            
        plt.plot(ts,xs[:,i,0],c=varit[i], lw=paksuudet[i], label=nimet[i])
    if KOKONAISSUUREET:
        mkp = np.zeros(len(ts))
        for i in range(len(ts)):
            mkp[i] = laske_massakeskipiste( xs[i,:,:] )[0]
        plt.plot(ts, mkp, c='brown', lw=1, ls=':', label='massakeskipiste')
    plt.legend()
    plt.savefig("x_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("y (m)")
    for i in range(N_OSAT-1,-1,-1):            
        plt.plot(ts,xs[:,i,1],c=varit[i], lw=paksuudet[i], label=nimet[i])
    if KOKONAISSUUREET:
        mkp = np.zeros(len(ts))
        for i in range(len(ts)):
            mkp[i] = laske_massakeskipiste( xs[i,:,:] )[1]
        plt.plot(ts, mkp, c='brown', lw=1, ls=':', label='massakeskipiste')
    plt.legend()
    plt.savefig("y_t."+KUVAFORMAATTI)
    plt.clf()
    
    xmax = 0.5*(Y_MAX-Y_MIN)+X_CENTER
    xmin = -0.5*(Y_MAX-Y_MIN)+X_CENTER
    ymin = Y_MIN
    ymax = Y_MAX
    fig, ax = plt.subplots()
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    for i in range(N_OSAT-1,-1,-1):            
        plt.plot(xs[:,i,0],xs[:,i,1],c=varit[i], lw=paksuudet[i], label=nimet[i])
    if KOKONAISSUUREET:
        mkp = np.zeros([len(ts),2])
        for i in range(len(ts)):
            mkp[i] = laske_massakeskipiste( xs[i,:,:] )
        plt.plot(mkp[:,0], mkp[:,1], c='brown', lw=1, ls=':', label='massakeskipiste')
    plt.legend()
    plt.savefig("x_y."+KUVAFORMAATTI)
    plt.clf()

    plt.xlabel("t (s)")
    plt.ylabel("$v_x$ (m/s)")
    for i in range(N_OSAT-1,-1,-1):  
        plt.plot(ts,vs[:,i,0],c=varit[i], lw=paksuudet[i], label=nimet[i])
    if KOKONAISSUUREET:
        mkp = np.zeros(len(ts))
        for i in range(len(ts)):
            mkp[i] = laske_massakeskipiste( vs[i,:,:] )[0]
        plt.plot(ts, mkp, c='brown', lw=1, ls=':', label='massakeskipiste')
    plt.legend()
    plt.savefig("vx_t."+KUVAFORMAATTI)
    plt.clf()

    plt.xlabel("t (s)")
    plt.ylabel("$v_y$ (m/s)")
    for i in range(N_OSAT-1,-1,-1):  
        plt.plot(ts,vs[:,i,1],c=varit[i], lw=paksuudet[i], label=nimet[i])
    if KOKONAISSUUREET:
        mkp = np.zeros(len(ts))
        for i in range(len(ts)):
            mkp[i] = laske_massakeskipiste( vs[i,:,:] )[1]
        plt.plot(ts, mkp, c='brown', lw=1, ls=':', label='massakeskipiste')
    plt.legend()
    plt.savefig("vy_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("$p_x$ (kgm/s)")
    total = np.zeros(len(ts))
    for i in range(N_OSAT-1,-1,-1):  
        plt.plot(ts,m[i]*vs[:,i,0],c=varit[i], lw=paksuudet[i], label=nimet[i])
        total += m[i]*vs[:,i,0]
    if KOKONAISSUUREET:
        plt.plot(ts,total,c='brown', lw=1, ls=':', label='kokonais')
    plt.legend()
    plt.savefig("px_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("$p_y$ (kgm/s)")
    total = np.zeros(len(ts))
    for i in range(N_OSAT-1,-1,-1):  
        plt.plot(ts,m[i]*vs[:,i,1],c=varit[i], lw=paksuudet[i], label=nimet[i])
        total += m[i]*vs[:,i,1]
    if KOKONAISSUUREET:
        plt.plot(ts,total,c='brown', lw=1, ls=':', label='kokonais')
    plt.legend()
    plt.savefig("py_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("K (J)")
    total = np.zeros(len(ts))
    for i in range(N_OSAT-1,-1,-1):  
        plt.plot(ts,0.5*m[i]*(vs[:,i,0]**2+vs[:,i,1]**2),c=varit[i], lw=paksuudet[i], label=nimet[i])
        total += 0.5*m[i]*(vs[:,i,0]**2+vs[:,i,1]**2)
    if KOKONAISSUUREET:
        plt.plot(ts,total,c='brown', lw=1, ls=':', label='kokonais')
    plt.legend()
    plt.savefig("K_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("U (J)")
    total = np.zeros(len(ts))
    for i in range(N_OSAT-1,-1,-1):  
        plt.plot(ts,G*m[i]*xs[:,i,1],c=varit[i], lw=paksuudet[i], label=nimet[i])
        total += G*m[i]*xs[:,i,1]
    if KOKONAISSUUREET:
        plt.plot(ts,total,c='brown', lw=1, ls=':', label='kokonais')
    plt.legend()
    plt.savefig("U_t."+KUVAFORMAATTI)
    plt.clf()

    plt.xlabel("t (s)")
    plt.ylabel("K+U (J)")
    total = np.zeros(len(ts))
    for i in range(N_OSAT-1,-1,-1):  
        plt.plot(ts, G*m[i]*xs[:,i,1] + 0.5*m[i]*(vs[:,i,0]**2+vs[:,i,1]**2),c=varit[i], lw=paksuudet[i], label=nimet[i])
        total += G*m[i]*xs[:,i,1] + 0.5*m[i]*(vs[:,i,0]**2+vs[:,i,1]**2)
    if KOKONAISSUUREET:
        plt.plot(ts,total,c='brown', lw=1, ls=':', label='kokonais')
    plt.legend()
    plt.savefig("E_t."+KUVAFORMAATTI)
    plt.clf()
    plt.close('all')

# funktio piirra_kuvaajat loppuu


# Tallennetaan rakettia (vain isoa kappaletta) kuvaava data tiedostoon.
def tallenna_data():

    xs = np.array(x_historia)[:,0,:]
    vs = np.array(v_historia)[:,0,:]
    ts = np.array(t_historia)
    ps = m[0]*vs
    ks = 0.5 * m[0] * ( vs[:,0]**2 + vs[:,1]**2 )
    us = G * m[0] * xs[:,1]
    
    all_data = np.zeros( [len(ts), 9] )
    all_data[:,0] = ts
    all_data[:,1:3] = xs
    all_data[:,3:5] = vs
    all_data[:,5:7] = ps
    all_data[:,7] = ks
    all_data[:,8] = us

    np.savetxt('raketti_data.csv', all_data, fmt='%.3f', delimiter=',',
        header='t, x, y, vx, vy, px, py, K, U')

# funktio tallenna_data loppuu


# Näytetään väliaikatietoja simulaation etenemisestä.
def edistyminen(n, n_max, kerrat=10):
    
    if n % (n_max // kerrat) == 0 :  
        p = int(np.round((100*n)/n_max,0))   
        print(f"simuloitu {p:3d} %")

# funktio edistyminen loppuu


# Simuloidaan ajan kulumista eli annetaan kappaleiden liikkua.
def simuloi(x, v, t = 0, askeleita = N_ASKEL, n_jaljella = n_jaljella):

    if YKSIULOTTEINEN:
        v[:,0] = 0.0

    g = G * np.array( [0,1] )

    # alustetaan nopeus simulaatiota varten
    v += (0.5 * DT * -g)
    
    # otetaan haluttu määrä pieniä askeleita ajassa eteenpäin
    for askel in range(askeleita+1):
                
        # tallennetaan tiedot tasaisin välein
        if askel % ANIMAATIOASKEL == 0:
            t_historia.append( t )
            x_historia.append( copy.copy(x) )
            v_historia.append( v + (0.5 * g * DT) )
            
        # piirretään kuva tasaisin välein
        if askel % KUVAASKEL == 0:
            piirra_systeemi(x, v + (0.5 * g * DT), t)
            
        # työnnetään rakettia ylöspäin suihkuttamalla pakokaasua alaspäin
        if askel % KIIHDYTYSASKEL == 0 and n_jaljella > 0 and askel > 0 and askel <= KIIHDYTYSLOPPUASKEL:
            kiihdytys(v, n_jaljella, np.array( [DVX_KIIHDYTYS, DVY_KIIHDYTYS] ) )
            n_jaljella -= 1
            
        # räjäytetään raketti, kun aika on sopiva
        if askel == RAJAHDYSASKEL:
            rajahdys(v, n_jaljella, DV_RAJAHDYS)
            n_jaljella = 0

            
        # päivitetään aika sekä kappaleiden paikat ja nopeudet
        t = np.round(t+DT,5)
        if YKSIULOTTEINEN:
            v[:,0] = 0.0
        x += v * DT
        v += (-g * DT)
        
        # kerrotaan käyttäjälle, miten menee
        edistyminen( askel, askeleita+1 )
        
    # viimeistellään nopeus simulaation lopuksi
    v += (0.5 * DT * g)

    return x, v, t

# funktio simuloi loppuu
    
    
    
# Räjäytetään raketti.
# Tämä tapahtuu niin, että raketin mukana kulkeville pienille 
# kappaleille annetaan satunnaisia impulsseja.
def rajahdys(v, n_osia, sigma):
    
    rng = np.random.default_rng()
    
    summa = np.zeros(2)
    
    for i in range(1, n_osia+1):
        dv = rng.normal(size=2)*sigma
        
        # Annetaan impulssi aina niin päin, että
        # suuren kappaleen saama kokonaisimpulssi
        # on pieni. Ei näin tarvitsisi tehdä, mutta
        # tällä tavalla ei käy niin, että kaikki
        # pienet kappaleet sattumalta lentävät samaan
        # suuntaan.
        for j in range(2):
            if dv[j] * summa[j] > 0:
                dv[j] = -dv[j]
            
        summa += dv
        
        v[i] += dv
        v[0] -= M_PIENI / M_ISO * dv

# funktio rajahdys loppuu
  
    

# Kiihdytetään rakettia.
# Tämä tapahtuu niin, että yksi raketin mukana kulkevista
# pienistä kappaleista ammutaan suurella vauhdilla
# alaspäin, jolloin raketti saa impulssin ylöspäin.
def kiihdytys(v, n_osia, dv):

    v[n_osia] -= dv
    dv_raketti = dv * M_PIENI / (M_ISO + M_PIENI * (n_osia-1))

    for i in range(n_osia):
        v[i] += dv_raketti
  
# funktio kiihdytys loppuu
  


# Pääohjelma.
def main():

    simuloi(x, v)
    tallenna_data()
    piirra_kuvaajat()
    piirra_animaatio()

# funktio main loppuu



# Ajetaan pääohjelma.
if __name__ == '__main__':
    main()

