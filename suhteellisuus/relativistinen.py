#
# Suhteellisuusteoreettinen simulaatio hiukkasen liikkeestä.
#
# Simulaatio tapahtuu tiedoston parametrit.py ohjeiden mukaisesti.
# Hiukkaseen voi kohdistaa vakiovoiman ja sen voi hajottaa kahdeksi
# pienemmäksi hiukkaseksi.
#
# Teemu Hynninen 2025
#
import parametrit as par
import numpy as np
import copy
import warnings
warnings.filterwarnings("ignore")

XALKU_A = par.XALKU_A
VALKU_A = par.VALKU_A
F_A = par.F_A
T_KIIHDYTYS_ALKU_A = par.T_KIIHDYTYS_ALKU_A
T_KIIHDYTYS_LOPPU_A = par.T_KIIHDYTYS_LOPPU_A
T_HAJOAMINEN_A = par.T_HAJOAMINEN_A
M1 = par.M1
M2 = par.M2
M3 = par.M3
TMIN_A = par.TMIN_A
TMAX_A = par.TMAX_A

m = np.array( [M1,M2,M3] ) # massat
t = TMIN_A # aika

# tallennetaan paikka ja nopeus 3-vektoriin
# 1. alkio: alkuperäinen hiukkanen
# 2. ja 3. alkio: mahdollisessa hajoamisessa syntyneet hiukkaset
x = np.array( [ XALKU_A, np.inf, np.inf ] ) # paikka
vx = np.array( [ VALKU_A, 0.0, 0.0 ] ) # nopeus
c = 300 # valonnopeus, m/μs

x_historia = []
t_historia = []
dt = 0.001
historia_dt = 0.1
askeleita = int( (TMAX_A-TMIN_A) // dt + 1 )
historia_askel = int( historia_dt // dt )


# Lorentzin tekijä
def gamma(v, c=300):
    return 1.0/np.sqrt(1.0-v**2/c**2)


# Nopeuden Lorentz-muunnos
def v_muunnos(vx, vxAB, c=300):

    return ( vx - vxAB ) / ( 1 - (vx * vxAB)/c**2 )


# Hajottaa hiukkasen kahdeksi pienemmäksi.
def hajoa(x, vx):
    
    x[1] = x[0]
    x[2] = x[0]
    x[0] = np.inf
    
    # loppunopeudet massakeskipistekoordinaatistossa
    # hiukkanen 2 lähtee negatiiviseen suuntaan,
    # hiukkanen 3 positiiviseen suuntaan
    vx2B = -c * np.sqrt( (M1**4 + (M2**2 - M3**2)**2 - 2*M1**2 * (M2**2 + M3**2)) / (M1**2 + M2**2 - M3**2)**2 )
    vx3B = c * np.sqrt( (M1**4 + (M2**2 - M3**2)**2 - 2*M1**2 * (M2**2 + M3**2)) / (M1**2 - M2**2 + M3**2)**2 )
    
    # loppunopeudet laboratorion koordinaatistossa
    vx[1] = v_muunnos( vx2B, -vx[0] )
    vx[2] = v_muunnos( vx3B, -vx[0] )
    vx[0] = 0.0
    
    
# Laskee hiukkasiin kohdistuvan kiihtyvyyden.
def kiihtyvyys(t, v):
    if t >= T_KIIHDYTYS_ALKU_A and t <= T_KIIHDYTYS_LOPPU_A:
        a = np.zeros(3)
        for i in range(len(a)):
        
            # Pitää tarkistaa, onko mukana massattomia hiukkasia.
            if m[i] > 0:
                a[i] = F_A / (gamma(v[i])**3 * m[i])
            else:
                a[i] = 0.0
                
        return a
    
    return 0.0


# Simuloi hiukkasten liike.
def simuloi(t, x, vx):

    a = kiihtyvyys(t, vx)   
     
    # alustetaan nopeus simulaatiota varten
    vx += (0.5 * dt * a)
    
    # otetaan haluttu määrä pieniä askeleita ajassa eteenpäin
    for askel in range(askeleita + 1):
                
        # tallennetaan tiedot tasaisin välein
        if askel % historia_askel == 0:
            t_historia.append( t )
            x_historia.append( copy.copy(x) )            
            
            
        # päivitetään aika sekä paikka ja nopeus
        t = np.round(t+dt,5)
        x += vx * dt   
        a = kiihtyvyys(t, vx)
        vx += (a * dt)
        
        
        # hajotetaan hiukkanen, jos aika on sopiva
        if t <= T_HAJOAMINEN_A and t+dt > T_HAJOAMINEN_A:
            hajoa(x, vx)
        
    # viimeistellään nopeus simulaation lopuksi
    vx -= (0.5 * dt * a)


# Tallenna simulaatiodata tiedostoon.
def tallenna_data():

    xs = np.array(x_historia)
    ts = np.array(t_historia)
    
    for i in range(3):
        all_data = []
        
        for j in range( len(ts) ):

            test = xs[j,i]
            if not ( np.isnan( test ) or np.isinf( test ) ):
                all_data += [ ( xs[j,i], ts[j] ) ]
        

        if len(all_data) > 0:
            np.savetxt(f'liikerata_{i+1:d}.csv', np.array(all_data), fmt='%.3f', delimiter=',',
                header='x, t')


# pääohjelma
def main():
    simuloi(t, x, vx)
    tallenna_data()


# ajetaan pääohjelma
if __name__ == "__main__":
    main()