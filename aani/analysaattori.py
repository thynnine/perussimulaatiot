#
# Analysoi wav-muodossa olevan äänitiedoston ja piirtää
# kuvaajan äänen värähtelystä sekä spektristä.
#
# Teemu Hynninen 2025
#
import numpy as np
import scipy.io.wavfile
import matplotlib.pyplot as plt

import analysaattori_parametrit as par
NIMI = par.NIMI
MAKSIMITAAJUUS = par.MAKSIMITAAJUUS
KUVAFORMAATTI = par.KUVAFORMAATTI

# hakee äänitiedoston ja palauttaa siinä esiintyvät taajuudet ja näiden amplitudit
def laske_spektri():

    input_data = scipy.io.wavfile.read(NIMI+".wav")
    audio = input_data[1]
    aika = np.linspace(0,len(audio)/input_data[0],len(audio))

    try:
        # muutetaan stereoääni monoksi
        audio = audio.mean(axis=1)
    except:
        pass

    # kuinka paljon dataa on
    n = len(audio)
    n_erilliset = int(np.ceil((n+1)/2.0))
    
    # lasketaan taajuudet
    f = np.arange(0, n_erilliset, 1.0) * (input_data[0]/n);
    
    # lasketaan Fourier-muunnos
    p = np.fft.fft(audio) 

    # otetaan vain informatiivinen osuus
    p = p[0:n_erilliset]
    
    # spektri on fourier-muunnoksen itseisarvo tai neliö
    # otetaan nyt vain itseisarvo
    p = np.abs(p)

    # skaalataan 
    p = p / n

    # nollapiste tulee tuplana, joten skaalataan kaikki muut pisteet
    if n % 2 > 0:
        p[1:len(p)] = p[1:len(p)] * 2
    else:
        p[1:len(p) -1] = p[1:len(p) - 1] * 2 
        
    # palautetaan taajuudet ja niitä vastaavat amplitudit sekä alkuperäinen värähtely
    return f, p, aika, audio


# piirtää värähtelystä ja sen spektristä kuvaajat
def piirra_kuvaajat(f, p, t, a):

    if MAKSIMITAAJUUS > 4000:
        delta_f = 200
    elif MAKSIMITAAJUUS > 2000:
        delta_f = 100
    else:
        delta_f = 50
    n_ticks = int( MAKSIMITAAJUUS // delta_f + 1.1 )
    f_ticks = np.linspace(0, (n_ticks-1)*delta_f, n_ticks)

    plt.plot(t, a/np.max(np.abs(a)), color='k', lw=1)
    plt.xlim(0,t[-1])
    plt.ylim(-1.1, 1.1)
    plt.xlabel('aika (t)')
    plt.ylabel('värähtely')
    plt.grid()
    plt.savefig(NIMI+"_oskillaatio."+KUVAFORMAATTI)
    plt.clf()

    plt.plot(t, a/np.max(np.abs(a)), color='k', lw=1)
    plt.xlim(t[-1]/2-0.01,t[-1]/2+0.01)
    plt.ylim(-1.1, 1.1)
    plt.xlabel('aika (t)')
    plt.ylabel('värähtely')
    plt.grid()
    plt.savefig(NIMI+"_oskillaatio_lyhyt."+KUVAFORMAATTI)
    plt.clf()
    

    fig, ax = plt.subplots()
    ax.set_xticks(f_ticks,minor=True)
    plt.plot(f, p/np.max(p), color='k', lw=1)
    plt.xlim(0,MAKSIMITAAJUUS)
    plt.ylim(0,np.max(1.1 * p/np.max(p)))
    plt.xlabel('taajuus (Hz)')
    plt.ylabel('amplitudi')
    plt.grid()
    plt.savefig(NIMI+"_amplitudi."+KUVAFORMAATTI)
    plt.clf()

    power = 10 * np.log10(p**2)
    power = 80 + power - np.max(power)
    
    fig, ax = plt.subplots()
    ax.set_xticks(f_ticks,minor=True)
    plt.plot(f, power, color='k', lw=1)
    plt.xlim(0,MAKSIMITAAJUUS)
    plt.ylim(0,100)
    plt.xlabel('taajuus (Hz)')
    plt.ylabel('voimakkuus (dB)')
    plt.grid()
    plt.savefig(NIMI+"_voimakkuus."+KUVAFORMAATTI)
    plt.clf()
    

# pääohjelma
def main():
    f, p, t, a = laske_spektri()
    piirra_kuvaajat(f, p, t, a)
    
    
if __name__ == "__main__":
    main()