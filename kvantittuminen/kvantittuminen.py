#
# Simuloi kvanttimekaanisia seisovia aaltoja (ominaistiloja)
# ja niiden superpositioita erilaisissa potentiaalikuopissa.
#
# Teemu Hynninen 2025
#

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e
from scipy.linalg import eigh_tridiagonal
from matplotlib.animation import FuncAnimation

import parametrit as par

eV = 1.60218e-19  # elektronivoltti jouleissa
nm = 1e-9         # nanometri metreissä

L = par.KUOPAN_LEVEYS*nm
N = np.array( par.N ) - 1
KERTOIMET = np.array( par.KERTOIMET )
FPS = 30

if len(N) < len(KERTOIMET):
    KERTOIMET = KERTOIMET[:len(N)]
elif len(N) > len(KERTOIMET):
    print("kertoimia on vähemmän kuin kvanttilukuja")
    quit()

U_max_plot = par.KUVAN_KORKEUS*eV
T_SKAALA = par.SIMULAATION_NOPEUS
T_MAX = par.T_LOPPU

POTENTIAALIN_TYYPPI = par.POTENTIAALIN_TYYPPI
U0 = par.POTENTIAALIN_KORKEUS*eV

# ratkaisun resoluutio eli pisteiden lukumäärä x-akselilla:
# jos ratkaisuista tulee kummallisia, syy voi olla huono numeerinen
# tarkkuus, jolloin tätä voi kokeilla kasvattaa
Nx = 4001

KUVAFORMAATTI = par.KUVAFORMAATTI
TALLENNA_ANIMAATIO = par.TALLENNA_ANIMAATIO

# laskee potentiaalienergian U(x)
def potentiaalienergia(x):
    
    def porras(x, askeleen_korkeus):
    
        x0 = L / 2
        return np.piecewise(x, [x < x0, x >= x0], [0, askeleen_korkeus])
        
    def kaksi_kuoppaa(x, vallin_leveys, vallin_korkeus):
    
        x0 = L / 2
        dx = vallin_leveys / 2
        return np.piecewise(x, [x < x0-dx, ((x >= x0-dx) & (x < x0+dx)), x >= x0+dx], [0, vallin_korkeus, 0])

    def kuoppa_kuopassa(x, vallin_leveys, vallin_korkeus):
    
        x0 = L / 2
        dx = vallin_leveys / 2
        return np.piecewise(x, [x < x0-dx, ((x >= x0-dx) & (x < x0+dx)), x >= x0+dx], [vallin_korkeus, 0, vallin_korkeus])
        
    def paraabeli(x, k):
    
        x0 = L / 2
        return k*(x-x0)**2
            
    def maki(x, k):
        
        return k*x
        
    
    if POTENTIAALIN_TYYPPI == 1:
        return porras(x, U0)
    elif POTENTIAALIN_TYYPPI == 2:
        return kuoppa_kuopassa(x, 0.2*L, U0)
    elif POTENTIAALIN_TYYPPI == 3:
        return kaksi_kuoppaa(x, 0.05*L, U0)
    elif POTENTIAALIN_TYYPPI == 4:
        return paraabeli(x, U0 / (L/2)**2 )
    elif POTENTIAALIN_TYYPPI == 5:
        return maki(x, U0 / L )
    else:
        print("""
Potentiaalienergian tyypin pitää olla jokin näistä:
1: porras
2: pieni kuoppa suuremman sisällä
3: kaksi kuoppaa, joiden välissä on valli
4: paraabeli
5: suoraan nouseva mäki
""")
        quit()
    
# funktio potentiaalienergia loppuu
    

# ratkaisee ominaisfunktiot ja -energiat numeerisesti
def ominaistilat(L):
 
    dx = L / (Nx - 1)
    x = np.linspace(0, L, Nx)

    U = potentiaalienergia(x)

    # lasketaan Hamiltonin matriisi
    diagonaali = (hbar**2 / (m_e * dx**2)) + U  # diagonaalielementit
    off_diagonaali = -hbar**2 / (2 * m_e * dx**2) * np.ones(Nx - 1) 

    # ominaisratkaisut
    energiat, aaltofunktiot = eigh_tridiagonal( diagonaali, off_diagonaali )

    return energiat, aaltofunktiot
   
# funktio ominaistilat loppuu 
    
    
# normittaa annetut superpositiokertoimet
def normita_kertoimet(kertoimet):

    summa = np.sum( kertoimet @ kertoimet )
    return kertoimet / np.sqrt(summa)
    
# funktio normita_kertoimet loppuu
    
    
# piirtää kuvan potentiaalienergiasta sekä aaltofunktiosta
def piirra_systeemi(L, energiat, aaltofunktiot, kvanttiluvut, kertoimet, tiedosto="aaltofunktio"):

    plt.close('all')
    
    a = normita_kertoimet( kertoimet )
    dx = L / (Nx - 1)
    x = np.linspace(0, L, Nx)

    psi = np.zeros( Nx )
    E_odotus = 0.0
    
    U = potentiaalienergia(x)

    for i in range(len(kvanttiluvut)):
        
        n = kvanttiluvut[i]
        a_n = a[i]
        p_n = a_n**2
        psi_n = aaltofunktiot[:, n]
        E_n = energiat[n]
    
        #norm = np.sqrt(np.sum(psi_n**2) * dx)    
        psi += a_n * psi_n
        E_odotus += p_n * E_n        

    norm = np.sqrt(np.sum(psi**2) * dx) 
    psi = psi / norm

    psi_max = np.max( np.abs(psi) )

    fig, ax1 = plt.subplots()
    
    ax1.set_xlabel("$x$ ($10^{-9}$ m)")
    
    psi_plot_min = -psi_max - E_odotus / U_max_plot * 4 * psi_max
    psi_plot_max = psi_max*4 - E_odotus / U_max_plot * 4 * psi_max
    try:
        psi_suuruus = int( np.log10(psi_plot_max) )
    except:
        psi_suuruus = 0
    
    
    ax2 = ax1.twinx()
    ax1.grid(True, axis='x', lw=0.5, ls=':', zorder=-1)
    ax2.grid(True, axis='y', lw=0.5, ls=':', zorder=-1)
    
    ax2.set_ylabel("$U$ (eV)", color='orange')
    ax2.plot(x / nm, U / eV, color='orange', linestyle='-', lw=2, label="U(x)", zorder=10)
    ax1.set_xlim(0,L/nm)
    ax2.set_ylim(-0.25*U_max_plot / eV, U_max_plot / eV)
    ax2.tick_params(axis='y', labelcolor='orange')

    if len(kvanttiluvut) > 1:
        for n in kvanttiluvut:
    
            ax2.axhline( y = energiat[n] / eV, color='gray', linestyle='--', lw=1, zorder=4)
        
    ax2.axhline( y = E_odotus / eV, color='gray', linestyle='--', zorder=5)

    ax1.set_ylabel('$\psi$ ($10^'+f'{psi_suuruus:d}'+'$ m$^{-1/2}$)', color='k')
    ax1.plot(x / nm, psi / 10**psi_suuruus , color='k', lw=2, zorder=15 )    
    ax1.set_ylim(psi_plot_min / 10**psi_suuruus, psi_plot_max / 10**psi_suuruus )
    ax1.tick_params(axis='y', labelcolor='k')
    
    fig.tight_layout()
    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)
    
    ax1.set_zorder(ax2.get_zorder()+1)
    ax1.set_frame_on(False)
    plt.savefig(tiedosto+"."+KUVAFORMAATTI)
    
    plt.close('all')
   
# funktio piirrä systeemi loppuu 


# piirtää animaation aaltofunktiosta
def piirra_animaatio(L, energiat, aaltofunktiot, kvanttiluvut, kertoimet):

    print("piirretään animaatio...")
    
    
    plt.close('all')
    a = normita_kertoimet( kertoimet )
    dx = L / (Nx - 1)
    x = np.linspace(0, L, Nx)
    
    U = potentiaalienergia(x)
    
    psi = np.zeros( Nx )
    E_odotus = 0.0
        
    for i in range(len(kvanttiluvut)):
        
        n = kvanttiluvut[i]
        a_n = a[i]
        p_n = a_n**2
        psi_n = aaltofunktiot[:, n]
        E_n = energiat[n]
        
        norm = np.sqrt(np.sum(psi_n**2) * dx)    
        psi += a_n * psi_n / norm
        E_odotus += p_n * E_n

    psi_max = np.max( np.abs(psi) )
    
    kesto = T_MAX
    num_frames = int(kesto * FPS / T_SKAALA)
    t = np.linspace(0, kesto * 10**-15, num_frames)
    
    fig, ax1 = plt.subplots()
    
    psi_plot_min = -psi_max - E_odotus / U_max_plot * 4 * psi_max
    psi_plot_max = psi_max*4 - E_odotus / U_max_plot * 4 * psi_max
    try:
        psi_suuruus = int( np.log10(psi_plot_max) )
    except:
        psi_suuruus = 0

    ax1.set_ylabel('$\psi$ ($10^'+f'{psi_suuruus:d}'+'$ m$^{-1/2}$)', color='k')
    ax1.set_ylim(psi_plot_min / 10**psi_suuruus, psi_plot_max / 10**psi_suuruus )
    ax1.tick_params(axis='y', labelcolor='k')

    ax2 = ax1.twinx()
    ax1.set_xlabel("$x$ ($10^{-9}$ m)")
    ax2.set_ylabel("$U$ (eV)", color='orange')
    ax2.plot(x/nm, U / eV, color='orange', linestyle='-', lw=2, label="U(x)", zorder=10)
    ax2.set_xlim(0,L/nm)
    ax2.set_ylim(-0.25*U_max_plot / eV, U_max_plot / eV)
    ax2.tick_params(axis='y', labelcolor='orange')
    ax2.tick_params(axis='x', labelcolor='k')

    ax2.axhline(y=E_odotus / eV, color='gray', linestyle='--', zorder=5)

    ax1.grid(True, axis='x', lw=0.5, ls=':', zorder=-1)
    ax2.grid(True, axis='y', lw=0.5, ls=':', zorder=-1)
    fig.tight_layout()
    
    line_real, = ax1.plot([], [], color='black', label='Re[$\\psi(x)$]',zorder=15)
    line_imag, = ax1.plot([], [], color='red', label='Im[$\\psi(x)$]',zorder=14)
    
    ax1.set_zorder(ax2.get_zorder()+1)
    ax1.set_frame_on(False)
    
    time_text = ax2.text(0.01*L/nm, -0.23*U_max_plot / eV, f't ={0:6.0f} fs', fontsize=10, color='black')
    
    def alusta():
        line_real.set_data(x, np.zeros_like(x))
        line_imag.set_data(x, np.zeros_like(x))
        return line_real, line_imag
    
    def animoi(frame):
        time = t[frame]
        
        psi_t = np.zeros( Nx, dtype=complex )

        for i in range(len(kvanttiluvut)):
        
            n = kvanttiluvut[i]
            a_n = a[i]
            p_n = a_n**2
            psi_n = aaltofunktiot[:, n]
            E_n = energiat[n]
        
            norm = np.sqrt(np.sum(psi_n**2) * dx)
            psi_t += a_n * psi_n / norm * np.exp(-1j * E_n * time / hbar) / 10**psi_suuruus
        
        line_real.set_data(x/nm, np.real(psi_t) )
        line_imag.set_data(x/nm, np.imag(psi_t) )
        
        time_text.set_text(f't ={time*10**15:8.2f} fs')
        
        for coll in ax1.collections:
            coll.remove()
        ax1.fill_between(x/nm, np.abs(psi_t) , color='peru', alpha=0.3, step='mid')
        
        kuva_askel = int(par.KUVA_DT * FPS / T_SKAALA + 0.1)
    
        if frame%kuva_askel == 0:
            prosentti = np.min( [ int(100 * frame / (FPS * T_MAX / T_SKAALA) ), 100 ] )
            print(f"simuloitu: {prosentti:3d} %")
            plt.savefig(f"vangittu_elektroni_{int(time*10**15+0.1):d}fs."+par.KUVAFORMAATTI)
        
        return line_real, line_imag, ax1.collections[0], time_text
    
    ani = FuncAnimation(fig, animoi, frames=num_frames, init_func=alusta, blit=True, interval=1000 / FPS)
    ax1.legend()
    
    if TALLENNA_ANIMAATIO:
        ani.save('vangittu_elektroni.mp4', fps=FPS)
    else:
        plt.show()
    
    plt.close('all')

# funktio piirrä animaatio loppuu


# tulostaa annettuja kvanttilukuja vastaavat energiat
# ja superposition energian odotusarvon
def kerro_energiat(energiat, aaltofunktiot, kvanttiluvut, kertoimet):

    E_odotus = 0.0
    a = normita_kertoimet( kertoimet )
    
    for i in range(len(kvanttiluvut)):
        
        n = kvanttiluvut[i]
        a_n = a[i]
        p_n = a_n**2
        psi_n = aaltofunktiot[:, n]
        E_n = energiat[n]
        
        E_odotus += p_n * E_n
    
        print(f"tilan |{n+1:3d}) energia: {E_n/eV:8.2f} eV")
        
    if len(kvanttiluvut) > 1:
        print(f"superposition energian odotusarvo: {E_odotus/eV:8.3} eV")

# funktio kerro_energiat loppuu


# pääohjelma
def main():
    e, psi = ominaistilat(L)
    kerro_energiat(e, psi, N, KERTOIMET)    
    #piirra_systeemi(L, e, psi, N, KERTOIMET)
    piirra_animaatio(L, e, psi, N, KERTOIMET)

# funktio main loppuu    
    
    
if __name__ == "__main__":
    main()