#
# Simuloi kahden kaasusäiliön hakautumista tasapainoon.
#
# Teemu Hynninen 2025
#

import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.animation import FuncAnimation
import parametrit as par

TA = par.TA
TB = par.TB
VA = par.VA
VB = par.VB
nA = par.nA
nB = par.nB
cA = par.cA
cB = par.cB

LUKITTU = par.VAKIO_V
TERMOSTAATTI = par.VAKIO_T

T_LOPPU = par.T_LOPPU
ANIMAATIO_DT = par.ANIMAATIO_DT
KUVA_DT = par.KUVA_DT
KUVAFORMAATTI = par.KUVAFORMAATTI
TALLENNA_ANIMAATIO = par.TALLENNA_ANIMAATIO


DT = 0.01     # aika-askeleen pituus (pieni = tarkka mutta hidas simulaatio)
K  = 2.0      # tämä vaikuttaa lämmönsiirron nopeuteen
M  = 2.0e-8   # tämä vaikuttaa männän liikkeen nopeuteen
R = 8.314     # kaasuvakio

# johdetut parametrit
N_ASKEL = int(T_LOPPU / DT)
ANIMAATIOASKEL = int(ANIMAATIO_DT / DT)
KUVAASKEL = int(KUVA_DT / DT)

# tallennetaan data kuvien piirtämistä varten
t_historia = [  ]
A_historia = [  ]
B_historia = [  ]



# Piirretään kuva systeemistä yhdellä ajan hetkellä.
def piirra_systeemi(A, B, t):

    plt.close('all')
    fig, ax = plt.subplots()

    width = 100
    sailio = plt.Rectangle( (0,0), width, width//2, fill=False, color='black')
    ax.add_patch(sailio)
    splitter = ( A["V"] / (A["V"] + B["V"]) ) * width - 2
    manta = plt.Rectangle( (splitter,0), 4, width//2, fill=True, color='black')
    ax.add_patch(manta)    
    
    temp_max = 500
    ta = np.min( [temp_max, A['T']] )
    lampoA = plt.Rectangle( (0.05*width,5), 2, width//2 * ta/temp_max , fill=True, color='red')
    ax.add_patch(lampoA)  
    tb = np.min( [temp_max, B['T']] )
    lampoB = plt.Rectangle( (0.95*width-2,5), 2, width//2 * tb/temp_max , fill=True, color='red')
    ax.add_patch(lampoB)  
    
    A_text = ax.text(0.05*width, 0.6*width, 
                    f"$T_A$ = {A['T']:.1f} K \n"+\
                    f"$n_A$ = {A['n']:.2f} mol \n"+\
                    f"$V_A$ = {A['V']:.3f} m$^3$ \n"+\
                    f"$p_A$ = {A['P']/1000:.1f} kPa \n"+\
                    f"$\Delta S_A$ = {A['S']-A['S0']:.2f} J/K \n"+\
                    f"$E_A$ = {A['E']/1000:.2f} kJ", 
                    fontsize=12, color='black')
    
    B_text = ax.text(0.65*width, 0.6*width, 
                    f"$T_B$ = {B['T']:.1f} K \n"+\
                    f"$n_B$ = {B['n']:.2f} mol \n"+\
                    f"$V_B$ = {B['V']:.3f} m$^3$ \n"+\
                    f"$p_B$ = {B['P']/1000:.1f} kPa \n"+\
                    f"$\Delta S_B$ = {B['S']-B['S0']:3.2f} J/K \n"+\
                    f"$E_B$ = {B['E']/1000:.2f} kJ", 
                    fontsize=12, color='black')
    
    
    time_text = ax.text(0.37*width, 0.52*width, 
                    f't = {t:5.2f} s',
                    fontsize=12, color='black')
    ax.set_xlim(-5, width+5)
    ax.set_ylim(-5, width+5)
    ax.set_aspect('equal')
    plt.grid(False)
    plt.axis('off')
    #plt.title(f't = {t:.2f} s')
    plt.savefig(f"systeemi_t{t:.2f}."+KUVAFORMAATTI)
    plt.close('all')
    
# funktio piirra_systeemi loppuu


# Piirretään animaatio systeemin liikkeestä.
def piirra_animaatio():

    plt.close('all')
    fig, ax = plt.subplots()
    
    n_frames = len(t_historia)

    width = 100
    A = A_historia[0]
    B = B_historia[0]
    
    sailio = plt.Rectangle( (0,0), width, width//2, fill=False, color='black')
    ax.add_patch(sailio)
    splitter = ( A["V"] / (A["V"] + B["V"]) ) * width - 2
    manta = plt.Rectangle( (splitter,0), 4, width//2, fill=True, color='black')
    ax.add_patch(manta)    
    
    temp_max = 500
    ta = np.min( [temp_max, A['T']] )
    lampoA = plt.Rectangle( (0.05*width,5), 2, width//2 * ta/temp_max , fill=True, color='red')
    ax.add_patch(lampoA)  
    tb = np.min( [temp_max, B['T']] )
    lampoB = plt.Rectangle( (0.95*width-2,5), 2, width//2 * tb/temp_max , fill=True, color='red')
    ax.add_patch(lampoB)  
    
    A_text = ax.text(0.05*width, 0.6*width, 
                    f"$T_A$ = {A['T']:.1f} K \n"+\
                    f"$n_A$ = {A['n']:.2f} mol \n"+\
                    f"$V_A$ = {A['V']:.3f} m$^3$ \n"+\
                    f"$p_A$ = {A['P']/1000:.1f} kPa \n"+\
                    f"$\Delta S_A$ = {A['S']-A['S0']:.2f} J/K \n"+\
                    f"$E_A$ = {A['E']/1000:.2f} kJ", 
                    fontsize=12, color='black')
    
    B_text = ax.text(0.65*width, 0.6*width, 
                    f"$T_B$ = {B['T']:.1f} K \n"+\
                    f"$n_B$ = {B['n']:.2f} mol \n"+\
                    f"$V_B$ = {B['V']:.3f} m$^3$ \n"+\
                    f"$p_B$ = {B['P']/1000:.1f} kPa \n"+\
                    f"$\Delta S_B$ = {B['S']-B['S0']:3.2f} J/K \n"+\
                    f"$E_B$ = {B['E']/1000:.2f} kJ", 
                    fontsize=12, color='black')
                    
                    
    time_text = ax.text(0.37*width, 0.52*width, 
                    f't = 0 min',
                    fontsize=12, color='black')
    
    
    patches = [lampoA, lampoB, manta]
    texts = [A_text, B_text, time_text]
    
    ax.set_xlim(-5, width+5)
    ax.set_ylim(-5, width+5)
    ax.set_aspect('equal')
    plt.grid(False)
    plt.axis('off')


    def alusta():
        return patches

    # siirretään ympyrät ja nuolet oikeille paikoille
    def paivita(frame):

        t = t_historia[frame]
        A = A_historia[frame]
        B = B_historia[frame]
        
        texts[0].set_text(f"$T_A$ = {A['T']:.1f} K \n"+\
                    f"$n_A$ = {A['n']:.2f} mol \n"+\
                    f"$V_A$ = {A['V']:.3f} m$^3$ \n"+\
                    f"$p_A$ = {A['P']/1000:.1f} kPa \n"+\
                    f"$\Delta S_A$ = {A['S']-A['S0']:.2f} J/K \n"+\
                    f"$E_A$ = {A['E']/1000:.2f} kJ"
        )
        
        texts[1].set_text(f"$T_B$ = {B['T']:.1f} K \n"+\
                    f"$n_B$ = {B['n']:.2f} mol \n"+\
                    f"$V_B$ = {B['V']:.3f} m$^3$ \n"+\
                    f"$p_B$ = {B['P']/1000:.1f} kPa \n"+\
                    f"$\Delta S_B$ = {B['S']-B['S0']:3.2f} J/K \n"+\
                    f"$E_B$ = {B['E']/1000:.2f} kJ"
        )
        
        texts[2].set_text(f't = {t:5.2f} s')
        
        ta = np.min( [temp_max, A['T']] )
        h = width//2 * ta/temp_max
        patches[0].set_height( h )
        
        tb = np.min( [temp_max, B['T']] )
        h = width//2 * tb/temp_max
        patches[1].set_height( h )
        
        splitter = ( A["V"] / (A["V"] + B["V"]) ) * width - 2
        patches[2].set_xy( (splitter, 0) )
        
        return patches + texts

    # luodaan animaatio
    ani = FuncAnimation(fig, paivita, frames=n_frames, init_func=alusta, blit=True, interval=20)
    if TALLENNA_ANIMAATIO:
        print("tallennetaan animaatio")
        ani.save("animaatio.mp4")
    else:
        plt.show()
    plt.close(fig)

# funktio piirra_animaatio loppuu


# Piirretään joukko kuvaajia ja tallennetaan ne tiedostoihin.
def piirra_kuvaajat():

    def tee_vektori(historia, suure):
        vektori = []
        for hetki in historia:
            vektori += [ hetki[suure] ]            
        return np.array( vektori )

    AT = tee_vektori(A_historia, "T")
    BT = tee_vektori(B_historia, "T")
    AV = tee_vektori(A_historia, "V")
    BV = tee_vektori(B_historia, "V")
    AP = tee_vektori(A_historia, "P")
    BP = tee_vektori(B_historia, "P")
    AS = tee_vektori(A_historia, "S")
    BS = tee_vektori(B_historia, "S")
    AE = tee_vektori(A_historia, "E")
    BE = tee_vektori(B_historia, "E")

    plt.close('all')
    
    plt.xlabel("t (s)")
    plt.ylabel("T (K)")
    plt.grid()
    plt.plot(t_historia, AT, c='black', lw=2, label="kaasu A")
    plt.plot(t_historia, BT, c='red', lw=2, label="kaasu B")
    plt.legend()
    plt.savefig("T_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("V (m$^3$)")
    plt.grid()
    plt.plot(t_historia, AV, c='black', lw=2, label="kaasu A")
    plt.plot(t_historia, BV, c='red', lw=2, label="kaasu B")
    plt.plot(t_historia, AV+BV, c='brown', lw=2, label="summa")
    plt.legend()
    plt.savefig("V_t."+KUVAFORMAATTI)
    plt.clf()

    plt.xlabel("t (s)")
    plt.ylabel("p (kPa)")
    plt.grid()
    plt.plot(t_historia, AP/1000, c='black', lw=2, label="kaasu A")
    plt.plot(t_historia, BP/1000, c='red', lw=2, label="kaasu B")
    plt.legend()
    plt.savefig("P_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("$\Delta$S (J/K)")
    plt.grid()
    plt.plot(t_historia, AS - AS[0], c='black', lw=2, label="kaasu A")
    plt.plot(t_historia, BS - BS[0], c='red', lw=2, label="kaasu B")
    plt.plot(t_historia, ( AS+BS ) - ( AS[0]+BS[0] ), c='brown', lw=2, label="summa")
    plt.legend()
    plt.savefig("S_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("t (s)")
    plt.ylabel("E (J)")
    plt.grid()
    plt.plot(t_historia, AE, c='black', lw=2, label="kaasu A")
    plt.plot(t_historia, BE, c='red', lw=2, label="kaasu B")
    plt.plot(t_historia, AE+BE, c='brown', lw=2, label="summa")
    plt.legend()
    plt.savefig("E_t."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("$T_A$ (K)")
    plt.ylabel("$\Delta$S (J/K)")
    plt.grid()
    plt.scatter(AT, AS - AS[0], s=1, c='b', label="kaasu A")
    plt.scatter(AT, BS - BS[0], s=1, c='r', label="kaasu B")
    plt.scatter(AT, (AS+BS) - (AS[0]+BS[0]), s=1, c='brown', label="summa")
    plt.legend()
    plt.savefig("S_TA."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("$T_B$ (K)")
    plt.ylabel("$\Delta$S (J/K)")
    plt.grid()
    plt.scatter(BT, AS - AS[0], s=1, c='b', label="kaasu A")
    plt.scatter(BT, BS - BS[0], s=1, c='r', label="kaasu B")
    plt.scatter(BT, (AS+BS) - (AS[0]+BS[0]), s=1, c='brown', label="summa")
    plt.legend()
    plt.savefig("S_TB."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("$p_A$ (kPa)")
    plt.ylabel("$\Delta$S (J/K)")
    plt.grid()
    plt.scatter(AP/1000, AS - AS[0], s=1, c='b', label="kaasu A")
    plt.scatter(AP/1000, BS - BS[0], s=1, c='r', label="kaasu B")
    plt.scatter(AP/1000, (AS+BS) - (AS[0]+BS[0]), s=1, c='brown', label="summa")
    plt.legend()
    plt.savefig("S_PA."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("$p_B$ (kPa)")
    plt.ylabel("$\Delta$S (J/K)")
    plt.grid()
    plt.scatter(BP/1000, AS - AS[0], s=1, c='b', label="kaasu A")
    plt.scatter(BP/1000, BS - BS[0], s=1, c='r', label="kaasu B")
    plt.scatter(BP/1000, (AS+BS) - (AS[0]+BS[0]), s=1, c='brown', label="summa")
    plt.legend()
    plt.savefig("S_PB."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("$V_A$ (m$^3$)")
    plt.ylabel("$\Delta$S (J/K)")
    plt.grid()
    plt.scatter(AV, AS - AS[0], s=1, c='b', label="kaasu A")
    plt.scatter(AV, BS - BS[0], s=1, c='r', label="kaasu B")
    plt.scatter(AV, (AS+BS) - (AS[0]+BS[0]), s=1, c='brown', label="summa")
    plt.legend()
    plt.savefig("S_VA."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("$V_B$ (m$^3$)")
    plt.ylabel("$\Delta$S (J/K)")
    plt.grid()
    plt.scatter(BV, AS - AS[0], s=1, c='b', label="kaasu A")
    plt.scatter(BV, BS - BS[0], s=1, c='r', label="kaasu B")
    plt.scatter(BV, (AS+BS) - (AS[0]+BS[0]), s=1, c='brown', label="summa")
    plt.legend()
    plt.savefig("S_VB."+KUVAFORMAATTI)
    plt.clf()
    
    plt.xlabel("V (m$^3$)")
    plt.ylabel("p (kPa)")
    plt.grid()
    plt.scatter(AV, AP/1000, s=1, c='b', label="kaasu A")
    plt.scatter(BV, BP/1000, s=1, c='r', label="kaasu B")
    plt.legend()
    plt.savefig("P_V."+KUVAFORMAATTI)
    plt.clf()
    
    plt.close('all')
    
# funktio piirra_kuvaajat loppuu


# laskee paineen
def paine(T, V, n):

    return n * R * T / V

# laskee entropian
def entropia(T, V, n, c):

    return n * R * np.log( V ) + n * c * np.log( T )
    
# laskee energian
def energia(T, n, c):

    return n * c * T

# laskee lämpötilan
def lampotila(E, n, c):

    return E / (n * c)


# kaasujen ominaisuudet
GAS_A = {
"T" : TA,
"V" : VA,
"n" : nA,
"c" : cA,
"P" : paine(TA, VA, nA),
"S" : entropia(TA, VA, nA, cA),
"E" : energia(TA, nA, cA),
"S0": entropia(TA, VA, nA, cA)
}

GAS_B = {
"T" : TB,
"V" : VB,
"n" : nB,
"c" : cB,
"P" : paine(TB, VB, nB),
"S" : entropia(TB, VB, nB, cB),
"E" : energia(TB, nB, cB),
"S0": entropia(TB, VB, nB, cB)
}


# Näytetään väliaikatietoja simulaation etenemisestä.
def edistyminen(n, n_max, kerrat=10):
    
    if n % (n_max // kerrat) == 0 :  
        p = int(np.round((100*n)/n_max,0))   
        print(f"simuloitu {p:3d} %")
   
# funktio edistyminen loppuu


# Simuloidaan ajan kulumista.
def simuloi(A, B, t = 0, askeleita = N_ASKEL):
    
    # otetaan haluttu määrä pieniä askeleita ajassa eteenpäin
    for askel in range(askeleita+1):                
                
        # tallennetaan tiedot tasaisin välein
        if askel % ANIMAATIOASKEL == 0:
            t_historia.append( t )
            A_historia.append( copy.copy(A) )
            B_historia.append( copy.copy(B) )
            
        # piirretään kuva tasaisin välein
        if askel % KUVAASKEL == 0:
            piirra_systeemi(A, B, t)
            
        # päivitetään aika sekä kaasujen ominaisuudet
        t = np.round(t+DT,5)      
        deltaT = A["T"] - B["T"]
        deltaP = A["P"] - B["P"]
        avrP = 0.5*(A["P"]+B["P"])
        
            
        if LUKITTU:
            dV = 0
        else:
            dV =  M * deltaP * DT
                    
        dW = avrP * dV
        
        A["V"] += dV
        B["V"] -= dV
        
        if TERMOSTAATTI:
            A["P"] = A["n"] * R * A["T"] / A["V"]
            B["P"] = B["n"] * R * B["T"] / B["V"]
            
        else:
            dQ = -K * deltaT * DT
            A["E"] += dQ - dW
            B["E"] -= dQ - dW
            A["T"] = lampotila( A["E"], A["n"], A["c"] )
            B["T"] = lampotila( B["E"], B["n"], B["c"] )        
            A["P"] = paine( A["T"], A["V"], A["n"] )
            B["P"] = paine( B["T"], B["V"], B["n"] )
        
        A["S"] = entropia( A["T"], A["V"], A["n"], A["c"] )
        B["S"] = entropia( B["T"], B["V"], B["n"], B["c"] )

# funktio simuloi loppuu
     

    
    


# Pääohjelma.
def main():

    simuloi(GAS_A, GAS_B)
    piirra_kuvaajat()
    piirra_animaatio()
    
# funktio main loppuu




# Ajetaan pääohjelma.
if __name__ == "__main__":
    main()

