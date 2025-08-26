#
# Tekee Fourier-synteesin eli
# syntetisoi äänen annetusta taajuusdatasta
# laskemalla summan
#
#  u(t) = sum_i A(i) cos( 2 pi f(i) t + phi(i) )
#
# Teemu Hynninen 2025
#
import numpy as np
import scipy.io.wavfile
import matplotlib.pyplot as plt
import syntetisaattori_parametrit as par

kesto = 4  # äänen kesto sekunteina
bps = 12800  # datapisteiden määrä yhden sekunnin aikana
taajuudet  = par.TAAJUUDET
amplitudit = par.AMPLITUDIT
vaiheet = par.VAIHEET

# amplitudeja ja vaiheita tarvitaan yhtä monta kuin taajuuksia
# lisää nollia, jos arvoja puuttuu
while len(amplitudit) < len(taajuudet):
    amplitudit += [0]
while len(vaiheet) < len(taajuudet):
    vaiheet += [0]

# Varmistetaan, että amplitudit ovat tarpeeksi
# pienet, jottei äänestä tule kovin voimakas.
#
# Älä muuta tätä!
amplitudit = np.abs( np.array( amplitudit, dtype=float ) ) 
amplitudit /= 5*np.max( amplitudit )

# aika
t = np.linspace(0,kesto,bps*kesto)

# alustetaan aaltofunktio nollaksi
aani = np.zeros( len(t) )

# Fourier-synteesi: summataan harmoniset värähtelyt
for i in range( len(taajuudet) ):
    a = amplitudit[i]
    f = taajuudet[i]
    phi = vaiheet[i]
    
    aani += a * np.cos( 2*np.pi*f*t + phi )

# skaalataan ääni niin, että se voimistuu ja heikkenee vähitellen
aani *= np.exp( -( t - kesto/2 )**2 / (0.2*kesto)**2 )

# tallennetaan ääni
scipy.io.wavfile.write("synteettinen_aani.wav", bps, aani)
