#
# Simuloi kvanttimekaanisen elektronin liikettä erilaisissa systeemeissä,
# kun elektroni ei muodosta seisovia aaltoja.
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import parametrit as par


# Asetetaan luonnonvakioille arvoksi 1, ja käsitellään pituudet nanometreissä,
# sillä laskenta on numeerisesti tarkempaa, kun lukuarvot eivät ole todella suuria tai pieniä.
# Tämä muuttaa mittakaavaa mutta ei vaikuta muuten tuloksiin. 
# Skaalataan lopuksi tulokset SI-yksiköihin seuraavilla kertoimilla.
hbar = 1.0
m = 1.0
x_skaala, x_exp = 1, -9          # nm -> m
psi_skaala, psi_exp = 3.162, 4   # nm^-1/2 -> m^-1/2
p_skaala, p_exp = 1.055, -25     # hbar / nm -> kgm/s
phi_skaala, phi_exp = 3.079, 12  # (hbar / nm)^-1/2 -> (kgm/s)^-1/2
t_skaala, t_exp = 8.634, -15     # nm^2 / ( hbar / me ) -> s
u_skaala = 0.0761996             # hbar^2 / (me nm^2) -> eV

DT = par.SIMULAATION_NOPEUS/t_skaala # aika-askel
L = 300 # simulaatioalueen leveys, nm

# Laskee aaltoa vaimentavan potentiaali simulaatioalueen reunoille.
# Tämä vähentää aaltojen heijastumista reunoilta takaisin simulaatioalueelle.
def absorptiopotentiaali(x, width=30.0, strength=1.0):
    Vabs = np.zeros_like(x, dtype=complex)
    edge = np.max(np.abs(x)) - width
    mask = np.abs(x) > edge
    Vabs[mask] = -1j * strength * ((np.abs(x[mask]) - edge) / width)**2
    return Vabs

# Laskee potentiaalienergian.
def potentiaalienergia(x):
    
    U_type = par.POTENTIAALIN_TYYPPI
    U0 = par.POTENTIAALIN_KORKEUS

    if U_type == 1:
        # vapaa hiukkanen, U = 0
        U = 0*x 

    elif U_type == 2:
        # mäki, U = k x
        U = x * U0 / (L/2) / u_skaala
        
    elif U_type == 3:
        # paraabeli
        U = x**2 * U0 / (L/2)**2 / u_skaala
        
    elif U_type == 4:
        # askel
        U = 0.5*( U0 * np.tanh(0.2*x) + U0 ) / u_skaala

    elif U_type == 5:
        # valli
        sigma = 2
        U = U0 * np.exp(-(x)**2 / (2 * sigma**2)) / u_skaala
        
        
    else:
        print("""
Potentiaalienergian tyypin pitää olla jokin näistä:
1: nolla (vapaa hiukkanen)
2: mäki 
3: paraabeli 
3: askel
4: valli
""")
        quit()
        
    return U, U + absorptiopotentiaali(x)


# simuloidaan aaltofunktion kehitys ajan funktiona ja piirretään samalla animaatio
def simuloi():
    global psi

    # paikka- ja liikemääräavaruudet
    x_min, x_max = -300, 300
    x_plotmin, x_plotmax = x_min/2, x_max/2
    N = 2**14 # simulaatiopisteiden määrä
    x = np.linspace(x_min, x_max, N)
    dx = x[1] - x[0]
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    p = hbar * k
    dp = p[1] - p[0]
    psimax = 0.5 * psi_skaala
    phimax = 3.0 * phi_skaala

    # aika
    dt = DT
    tmax = par.T_LOPPU/t_skaala
    tsteps = int(tmax/dt+0.1)+1

    # potentiaalienergia
    V, V_total = potentiaalienergia(x)

    # aaltofunktio alussa: Gaussin pulssi
    x0 = par.X_ALKU
    sigma = par.DX_ALKU
    p0 = par.P_ALKU / p_skaala
    p_max = 5.0
    psi0 = (1 / (2*np.pi * sigma**2))**0.25 * np.exp(-(x - x0)**2 / (4 * sigma**2)) * np.exp(1j * p0 * x / hbar)
    psi = psi0.copy()

    # aaltolukuesitys
    psi_p = np.fft.fftshift(np.fft.fft(psi)) * dx / np.sqrt(2 * np.pi * hbar)
    p_shifted = np.fft.fftshift(p)
        
    # e odotusarvo - simulaatio on liian epätarkka energian keskihajonnan laskemiseksi
    mu_K = np.real( np.sum( np.abs(psi_p)**2 * (p_shifted ** 2) / (2 * m) ) * dp )
    mu_U = np.real( np.sum( np.abs(psi)**2 * V_total ) * dx )    
    print()
    print("elektronin energian odotusarvo alussa:")
    print(f"K = {mu_K*u_skaala:8.4f} eV")
    print(f"U = {mu_U*u_skaala:8.4f} eV")
    print(f"E = {(mu_K+mu_U)*u_skaala:8.4f} eV")
    print()

    # simulaation aikakehitysoperaattorit
    expV_half = np.exp(-1j * V_total * dt / (2 * hbar))
    expT = np.exp(-1j * (p**2) * dt / (2 * m * hbar))

    # tallennetaan x ja p
    mu_x_list = []
    sigma_x_list = []
    mu_p_list = []
    sigma_p_list = []

    # luodaan kuvaajat
    fig, (ax_x, ax_p) = plt.subplots(1, 2, figsize=(14, 5), tight_layout=True)

    # aaltofunktion kuvaaja
    line_abs_x, = ax_x.plot([], [], color="brown", linewidth=0.1, label="|ψ(x)|")
    line_re_x, = ax_x.plot([], [], color="black", linewidth=0.4, label="Re[ψ(x)]")
    line_im_x, = ax_x.plot([], [], color="red", linewidth=0.4, label="Im[ψ(x)]")
    
    ax2 = ax_x.twinx()
    ax2.set_ylabel('$U$ (eV)', color='orange')
    ax2.set_ylim( -psimax, psimax )
    ax2.tick_params(axis='y', labelcolor='orange')
    ax2.plot(x, V*u_skaala , color="orange", linestyle="-", label="U(x)")
    
    fill_x = None
    ax_x.set_xlim(x_plotmin, x_plotmax)
    ax_x.set_ylim(-psimax, psimax)
    ax_x.set_xlabel("$x$ ($10^{"+f"{x_exp}"+"}$ m)")
    ax_x.set_ylabel("$\psi$ ($10^{"+f"{psi_exp}"+"}$ m$^{-1/2}$)")
    ax_x.set_title("")
    ax_x.grid(True)
    ax_x.legend()

    # aaltolukuesityksen kuvaaja
    line_abs_p, = ax_p.plot([], [], color="saddlebrown", linewidth=0.1, label="|a(p)|")
    #line_re_p, = ax_p.plot([], [], color="black", linewidth=0.2, label="Re[a(p)]")
    #line_im_p, = ax_p.plot([], [], color="red", linewidth=0.2, label="Im[a(p)]")
    fill_p = None
    ax_p.set_xlim(-p_max, p_max)
    ax_p.set_ylim(-phimax, phimax)
    ax_p.set_xlabel("$p_x$ ($10^{"+f"{p_exp}"+"}$ kgm/s)")
    ax_p.set_ylabel("$a$ ($10^{"+f"{phi_exp}"+"}$ (kgm/s)$^{-1/2}$)")
    ax_p.set_title("")
    ax_p.grid(True)
    ax_p.legend()


    time_text = ax2.text(-0.49*L, -0.9*psimax, f't ={0:6.0f} fs', fontsize=10, color='black')
    

    # animaatio
    def init():
        global fill_x, fill_p
        line_abs_x.set_data([], [])
        line_re_x.set_data([], [])
        line_im_x.set_data([], [])
        fill_x = ax_x.fill_between([], [], [], color="peru", alpha=0.3)

        line_abs_p.set_data([], [])
        #line_re_p.set_data([], [])
        #line_im_p.set_data([], [])
        fill_p = ax_p.fill_between([], [], [], color="sandybrown", alpha=0.3)
    
        return line_abs_x, line_re_x, line_im_x, fill_x, line_abs_p, fill_p
        #return line_abs_x, line_re_x, line_im_x, fill_x, line_abs_p, line_re_p, line_im_p, fill_p

    def animate(i):
        global psi, fill_x, fill_p

        # poimitaan aaltofunktion itseisarvo sekä reaali- ja imaginääriosat
        abs_psi = np.abs(psi)
        re_psi = np.real(psi)
        im_psi = np.imag(psi)
    
        # samoin aaltolukuesitykselle
        psi_p = np.fft.fftshift(np.fft.fft(psi)) * dx / np.sqrt(2 * np.pi * hbar)
        p_shifted = np.fft.fftshift(p)
        abs_psi_p = np.abs(psi_p)
        #re_psi_p = np.real(psi_p)
        #im_psi_p = np.imag(psi_p)
        
        # päivitetään kuvaajat
        line_abs_x.set_data(x * x_skaala, abs_psi * psi_skaala)
        line_re_x.set_data(x * x_skaala, re_psi * psi_skaala)
        line_im_x.set_data(x * x_skaala, im_psi * psi_skaala)

        for coll in ax_x.collections:
            coll.remove()
        fill_x = ax_x.fill_between(x * x_skaala, 0, abs_psi * psi_skaala, color="peru", alpha=0.3)

        line_abs_p.set_data(p_shifted * p_skaala, abs_psi_p * phi_skaala)
        #line_re_p.set_data(p_shifted, re_psi_p)
        #line_im_p.set_data(p_shifted, im_psi_p)

        for coll in ax_p.collections:
            coll.remove()
        fill_p = ax_p.fill_between(p_shifted * p_skaala, 0, abs_psi_p * phi_skaala, color="sandybrown", alpha=0.3)

        ax_x.set_title("aaltofunktio")
        ax_p.set_title("liikemääräesitys")
    
        # tulostetaan väliaikatietoja
        t = i * DT * t_skaala
        time_text.set_text(f't ={int(t):6d} fs')
        kuva_askel = int( par.KUVA_DT/t_skaala / DT + 0.1 )
        if i%kuva_askel == 0:
            prosentti = np.min( [ int(100*i*dt/par.T_LOPPU*t_skaala+0.1), 100 ] )
            print(f"simuloitu: {prosentti:3d} %")
            plt.savefig(f"elektroniaalto_{int(t+0.1):d}fs."+par.KUVAFORMAATTI)

        # x odotusarvo ja keskihajonta
        mu_x = np.sum(np.abs(psi)**2 * x) * dx
        x2 = np.sum(np.abs(psi)**2 * x**2 ) * dx
        sigma_x = np.sqrt(x2 - mu_x**2)

        mu_x_list.append(mu_x)
        sigma_x_list.append(sigma_x)

        # p odotusarvo ja keskihajonta
        psi_p = np.fft.fftshift(np.fft.fft(psi)) * dx / np.sqrt(2 * np.pi * hbar)
        prob_p = np.abs(psi_p)**2
        mu_p = np.sum(p_shifted * prob_p) * dp
        p2 = np.sum(p_shifted**2 * prob_p) * dp
        sigma_p = np.sqrt(p2 - mu_p**2)

        mu_p_list.append(mu_p)
        sigma_p_list.append(sigma_p)

        # annetaan ajan kulua hiukan eteenpäin
        psi = expV_half * psi
        psi_k = np.fft.fft(psi)
        psi_k = expT * psi_k
        psi = np.fft.ifft(psi_k)
        psi = expV_half * psi
        
        # poimitaan itseisarvo, reaaliosa, imaginääriosa
        abs_psi = np.abs(psi)
        re_psi = np.real(psi)
        im_psi = np.imag(psi)

        # liikemääräesitys
        psi_p = np.fft.fftshift(np.fft.fft(psi)) * dx / np.sqrt(2 * np.pi * hbar)
        p_shifted = np.fft.fftshift(p)
        abs_psi_p = np.abs(psi_p)
        #re_psi_p = np.real(psi_p)
        #im_psi_p = np.imag(psi_p)

        return line_abs_x, line_re_x, line_im_x, fill_x, line_abs_p, fill_p, time_text
        #return line_abs_x, line_re_x, line_im_x, fill_x, line_abs_p, line_re_p, line_im_p, fill_p


    print("simuloidaan ja piirretään samalla animaatio")
    ani = animation.FuncAnimation(
        fig, animate, init_func=init, frames=tsteps, interval=20, blit=True
    )

    if par.TALLENNA_ANIMAATIO:
        ani.save("elektroniaalto.mp4")
    else:
        plt.show()
    plt.close("all")
    
    return np.array(mu_x_list), np.array(sigma_x_list), np.array(mu_p_list), np.array(sigma_p_list)

def piirra_keskiarvo_ja_hajonta(times, mu_list, sigma_list, label, filename):
    mu_arr = np.array(mu_list)
    sigma_arr = np.array(sigma_list)

    lowers = mu_arr - sigma_arr
    uppers = mu_arr + sigma_arr

    plt.figure(figsize=(8, 5))
    plt.plot(times, mu_arr, label="$\\mu$", color="black")
    plt.fill_between(times, lowers, uppers, color="sandybrown", alpha=0.5, label="$\\pm \\sigma$")
    plt.xlabel("$t$ ($10^{"+f"{t_exp}"+"}$ s)")
    plt.ylabel(label)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def main():
    mu_x_list, sigma_x_list, mu_p_list, sigma_p_list = simuloi()
    ajat = np.arange(len(mu_x_list)) * DT * t_skaala
    piirra_keskiarvo_ja_hajonta(ajat, mu_x_list*x_skaala, sigma_x_list*x_skaala, "$x$ ($10^{"+f"{x_exp}"+"}$ m)", "x_t."+par.KUVAFORMAATTI)
    piirra_keskiarvo_ja_hajonta(ajat, mu_p_list*p_skaala, sigma_p_list*p_skaala, "$p_x$ ($10^{"+f"{p_exp}"+"}$ kgm/s)", "p_t."+par.KUVAFORMAATTI)


if __name__ == "__main__":
    main()

