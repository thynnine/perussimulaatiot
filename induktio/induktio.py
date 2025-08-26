#
# Simuloi sähkömagneettisen induktion magneetin liikkuessa käämin läpi.
#
# Teemu Hynninen 2025
#
import copy
import numpy as np
import matplotlib.pyplot as plt

# haetaan parametrit
import parametrit as par

EPSILON = 8.8541878*10**-12; # sähkövakio
MU = 1.256637061*10**-6; # magneettivakio
G = 9.8; # putoamiskiihtyvyys
Lm = par.MAGNEETIN_PITUUS
Lk = par.KAAMIN_PITUUS
R = par.KAAMIN_SADE
N = par.KAAMIN_KIERROKSET
MUm = par.DIPOLIMOMENTTI
Z0 = par.ALKUKORKEUS
DT = par.DT_KUVA
T_LOPPU = par.T_LOPPU
KUVAFORMAATTI = par.KUVAFORMAATTI

# Laskee magneettivuon renkaan läpi, kun magneettinen dipoli mu
# on etäisyydellä z renkaan keskipisteen yläpuolella.
def magneettivuo(z, mu):

    return MU*mu/2 * R**2 / (z**2 + R**2)**(3/2)

# funktio magneettivuo loppuu


# Laskee magneettivuon käämin läpi, kun magneetin keskipiste
# on etäisyydellä z käämin keskipisteen yläpuolella.
def kokonaisvuo(z):

    # jaetaan magneetti ja käämi näin moneen osaan
    nm = int(1000*Lm+1)
    nk = int(1000*Lk+1)
    
    vuo = 0.0
    
    # jaetaan magneetti ja käämi osiin ja käydään läpi kaikki osat
    for im in range(nm):
        for ik in range(nk):

            # osan etäisyys magneetin tai käämin keskipisteestä
            dzm = ( 0.5 - (0.5 + im) / nm ) * Lm
            dzk = ( 0.5 - (0.5 + ik) / nk ) * Lk 
                       
            # magneetin osan vuo käämin osan läpi
            vuo += N / nk * magneettivuo( z + dzm + dzk, MUm / nm )

    return vuo

# funktio kokonaisvuo loppuu


# Laskee induktiojännitteen käämissä, kun magneetin keskipiste
# on etäisyydellä z käämin keskipisteen yläpuolella ja liikkuu
# vauhdilla v alaspäin.
def jannite(z, v, dt = DT):
    
    phi0 = kokonaisvuo(z + v*dt) # vuo hetki sitten
    phi1 = kokonaisvuo(z - v*dt) # vuo hetken päästä
    return -( phi1 - phi0 ) / (2*dt) # Faradayn induktiolaki

# funktio jannite loppuu


# Simuloidaan koe, jossa magneetti putoaa käämin läpi.
def simuloi():

    t = 0.0
    z = Z0
    v = 0.0
    phi = kokonaisvuo(z)
    dv = jannite(z, v)
    
    aikahistoria = [t]
    paikkahistoria = [z]
    nopeushistoria = [v]
    vuohistoria = [phi]
    jannitehistoria = [dv]
    
    n_askel = int( T_LOPPU / DT + 0.5 )
    
    for i in range(n_askel):
    
        t += DT
        z -= v*DT + 0.5*G*DT**2
        v += G*DT
        phi = kokonaisvuo(z)
        dv = jannite(z, v)
        
        aikahistoria += [t]
        paikkahistoria += [z]
        nopeushistoria += [v]
        vuohistoria += [phi]
        jannitehistoria += [dv]
        
        
    return aikahistoria, paikkahistoria, nopeushistoria, vuohistoria, jannitehistoria

# funktio simuloi loppuu


# Piirretään joukko kuvaajia ja tallennetaan ne tiedostoihin.
def piirra_kuvaajat(t, z, v, phi, dv):

    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlabel("$t$ (s)")
    plt.ylabel("$\\Delta V$ (V)")
    plt.plot(t, dv, c='blue', lw=2)
    plt.grid()     
    plt.ylim([-2,2])
    plt.savefig("dV_t."+KUVAFORMAATTI)
    
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlabel("$z$ (m)")
    plt.ylabel("$\\Delta V$ (V)")
    plt.plot(z, dv, c='blue', lw=2)
    plt.grid()     
    plt.ylim([-2,2])
    plt.savefig("dV_z."+KUVAFORMAATTI)
    
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlabel("$t$ (s)")
    plt.ylabel("$\\Phi_B$ (Wb)")
    plt.plot(t, phi, c='red', lw=2)
    plt.grid()     
    plt.savefig("Phi_t."+KUVAFORMAATTI)
    
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlabel("$z$ (m)")
    plt.ylabel("$\\Phi_B$ (Wb)")
    plt.plot(z, phi, c='red', lw=2)
    plt.grid()     
    plt.savefig("Phi_z."+KUVAFORMAATTI)
    
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlabel("$t$ (s)")
    plt.ylabel("$z$ (m)")
    plt.plot(t, z, c='black', lw=2)
    plt.grid()     
    plt.savefig("z_t."+KUVAFORMAATTI)
    
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlabel("$t$ (s)")
    plt.ylabel("$v$ (m/s)")
    plt.plot(t, v, c='black', lw=2)
    plt.grid()     
    plt.savefig("v_t."+KUVAFORMAATTI)

# funktio piirra_kuvaajat loppuu



# Pääohjelma.
def main():
    t, z, v, phi, dv = simuloi()
    piirra_kuvaajat(t, z, v, phi, dv)
    
# funktio main loppuu


# Ajetaan pääohjelma.
if __name__ == '__main__':
    main()

