# Interferenssi
#
# Simuloidaan valon kulku erilaisten rakojen ja hilojen läpi.
#
# Kun monokromaattinen valoaalto kulkee kapeista raoista, se
# diffraktoituu ja eri rakojen kautta kulkeneet aallot
# interferoivat keskenään. Lopputulos on se, että valo ei kulje raoista
# suoraan kuten valon sädemalli ennustaa vaan hajaantuu eri kulmiin.
# Erityisesti valo ei hajaannu miten tahansa vaan eri reittejä kulkeneet
# aallot vahvistavat toisiaan vain niissä suunnissa, joissa aallot
# ovat samassa vaiheessa.
#
# Tämä ohjelma laskee annetuista raoista kulkeneen valon intensiteetin
# (kirkkauden) eri havaintokulmissa ja piirtää valon rakojen taakse 
# asetetulle varjostimelle muodostaman kuvion.
#
# Lasku tehdään vain yhdessä ulottuvuudessa.
# Käyttäjä voi siis laittaa valon kulkemaan kuinka leveistä raoista
# tahansa, mutta raot ovat kaikki pitkiä suorakulmioita.
# Niiden muotoa ei voi valita.
#
# Teemu Hynninen 2024

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
 


#
# SIMULAATIOPARAMETREJÄ, JOITA MUUTTAMALLA VOI OHJATA SIMULAATIOTA
#

# Kokeessa käytetyn valon ominaisuuksia.
#
# Näitä voi ja kannattaakin muuttaa.
#
# kirkkaus
# - piirrettävän kuvion kirkkaus
# - oltava väliltä (0,1)
# - tämä on vain piirtämisessä käytetty parametri, jolla ei
#   ole täsmällistä fysikaalista merkitystä
# - suuri luku tarkoittaa kirkasta valoa
# - jos kuviossa on himmeitä, heikosti näkyviä osia, kokeile
#   suurempaa arvoa
# - jos kuvio näyttää olevan melkein pelkkää valkoista, kokeile
#   pienempää arvoa
# - 0.5 ... 0.7 on usein hyvä arvo
#
# aallonpituus
# - valon aallonpituus nanometreissä
# - näkyvän valon aallonpituusalue on noin 350 nm ... 750 nm
# - interferenssikuvion leveys riippuu aallonpituuden ja rakojen
#   geometrian (leveyksien ja etäisyyksien) suhteesta
# - mitä lyhyempi aallonpituus, sitä pienempiin kulmiin
#   valo diffraktoituu
#
kirkkaus = 0.9 # yleensä 0.7 on hyvä arvo, mutta voi kokeilla 0.01 ... 0.99
aallonpituus = 532 # nm, valon aallonpituus, 400 ... 700


# Rakojen geometrian ominaisuuksia.
#
# Näitä voi ja kannattaakin muuttaa.
#
# Ohjelmaan on määritelty valmiiksi joukko kokeita, ja näillä
# parametreillä voit valita niistä haluamasi.
# 
# Toki voit laittaa valon kulkemaan myös ihan millaisesta raosta haluat.
# Se onnistuu kirjoittamalla tiedoston lopussa sopivat arvot muuttujiin
# rakojen_keskipisteet sekä rakojen_leveydet.
# 
# tapaus
# - valmiiksi määriteltyjä mielenkiintoisia kokeita.
#   1: kaksi kapeaa rakoa (kaksoisrakokoe eli Youngin koe)
#   2: kaksi kapeaa rakoa kaksinkertaisella etäisyydellä
#   3: neljä kapeaa rakoa
#   4: monta kapeaa rakoa (diffraktiohila)
#   5: yksi kapea rako (voimakas diffraktio)
#   6: yksi melko leveä rako (yhden raon diffraktio)
#   7: yksi rako kaksinkertaisella leveydellä
#   8: yksi hyvin leveä rako (valo kulkee suoraan, sädemalli toimii)
#   9: kaksi melko leveää rakoa (yhdistetty diffraktio ja interferenssi)
#  10: kaksi melko leveää rakoa kaksinkertaisella etäisyydellä
#  11: monta leveää rakoa
#  12: epäsäännöllinen joukko rakoja
#  13: epäsäännöllinen joukko rakoja
#
# w
# - kapeiden rakojen leveys
# - oletusarvoisesti valon aallonpituutta pienempi
# - Jos w on pienempi kuin aallonpituus, täsmällinen arvo ei
#   juurikaan vaikuta interferenssikuvion muotoon.
# - mitä kapeampi rako, sitä vähemmän valoa pääsee läpi
#   ja sitä himmeämpi valo varjostimella nähdään
#
# d
# - kapeiden rakojen keskipisteiden välinen etäisyys
# - leveiden rakojen leveys on tästä puolet
# - oletusarvoisesti valon aallonpituutta suurempi
# - Interferenssikuvion leveys on verrannollinen aallonpituuden ja
#   rakojen etäisyyden suhteeseen. Ts. mitä pienempi d, sitä
#   leveämpiä interferenssikuvioita nädhdään.
#
tapaus = 7
w = 300 # nm rakojen leveys
d = 50000 # nm rakojen välinen etäisyys

#
# MUITA PARAMETREJÄ, JOITA EI PITÄISI OLLA TARVETTA MUUTTAA
#

# Näkyvän valon spektrin rajat.
# Näitä ei pitäisi olla syytä muuttaa.
#
# - spektri (ei täsmälleen oikein, mutta lähellä todellista)
#   400 nm: violetti
#   450 nm: sininen
#   500 nm: turkoosi (syaani)
#   550 nm: vihreä
#   580 nm: keltainen
#   600 nm: oranssi
#   650 nm: punainen
MIN_AALLONPITUUS = 350.0 # nm, ultraviolettiraja - lyhyempiä aallonpituuksia ei piirretä
MAX_AALLONPITUUS = 750.0 # nm, infrapunaraja - pidempiä aallonpituuksia ei piirretä
PUNA_AALLONPITUUS = 630.0 # nm, punaisen valon aallonpituus
LILA_AALLONPITUUS = 380.0 # nm, violetin valon aallonpituus

# Rakojen eli aaltojen lähteiden geometriaan liittyviä vakioita.
# Näitä ei pitäisi olla syytä muuttaa.
MAX_ALKUPISTE = 50000.0 # nm, rakojen suurin mahdollinen paikkakoordinaatti
DELTA_ALKUPISTE = 10.0 # nm, laskennassa käytetty resoluutio, oltava tarpeeksi pieni
N_ALKUPISTE = int((2*MAX_ALKUPISTE)/DELTA_ALKUPISTE+1) # rakopisteiden määrä

# Havaintopisteiden geometriaan liittyviä vakoita.
# Näitä ei pitäisi olla syytä muuttaa.
# Poikkeuksena on MAX_KULMA, jota voi muuttaa välillä (0, 90) astetta.
# Se määrää piirrettyjen kuvien reunojen paikat.
MAX_KULMA = 5 # astetta, suurin mahdollinen havaintokulma
DELTA_KULMA = 0.005 # astetta, tulosten resoluutio
N_KULMA = int(2*MAX_KULMA/DELTA_KULMA+1) # havaintopisteiden lukumäärä
N_Y = 1000 # interferenssikuvion pikselien määrä pystysuunnassa
# Y_HAJONTA
# - intensiteettikuviot pehmennetään pystysuunnassa skaalaamalla
#   intensiteetti Gaussin funktiolla, joka keskipiste on kohdassa
#   0.5 (kuvan keskipiste) ja hajonta tämä arvo
# - pystysuunnassa piirtorajat ovat [0,1], joten hajonnan pitäisi
#   olla tätä suuruusluokkaa, jotta saadaan hyvä kuva
# - mitä pienempi hajonta, sitä vähemmän kuva leviää pystysuunnassa
# - jos None, intensiteetti piirretään pystysuunnassa vakiona
# - matkii pystysuuntaista diffraktiota, mutta ei oikeasti perustu
#   fysiikkaan vaan on pelkästään graafinen parametri
Y_HAJONTA = 0.01



#
# Lasketaan aaltojen vaiheet eri havaintokulmissa.
# 
# Funktio laskee vaiheet käyttäen kaukana olevan havaintopisteen approksimaatiota.
# Tällöin x-akselin pisteestä x kulmaan theta lähtevän aallon matkaero
# origosta lähtevään aaltoon verrattuna on
#
# Delta r = x sin( theta ), 
#
# jolloin jos origosta lähteneen aallon vaiheeksi asetetaan nolla, pisteestä x
# lähteneen aallon vaihe on
#
# phi = 2 pi Delta r / lambda = 2 pi ( x / lambda ) sin( theta ),
#
# missä lambda on aallonpituus.
#
# Funktio laskee nämä vaihetekijät kaikille annetuille alkupisteille x ja kulmille theta.
# Jos pisteitä on Nx ja kulmia Nk kpl, laskettavia vaihetekijöitä on yhteensä
# Nx kertaa Nk kpl. Tulokset palautetaan (Nx,Nk)-matriisina.
#
# kulmat: numpy-vektori, johon on listattu havaintokulmat theta.
# alkupisteet: numpy-vektori, johon on listattu alkupisteet x.
# aallonpituus: aallonpituus reaalilukuna.
# kulma_asteissa: jos True, kulmat oletetaan annetun asteissa, muuten radiaaneissa.
#
def vaiheet(kulmat, alkupisteet, aallonpituus, kulma_asteissa=True):

    if kulma_asteissa:
        muuntokerroin = np.pi/180
    else:
        muuntokerroin = 1
        
    vaiheet = 2*np.pi/aallonpituus * np.outer( alkupisteet, np.sin( kulmat*muuntokerroin ) )
    
    return vaiheet


#
# Lasketaan intensiteettijakauma.
#
# Havaintopisteeseen saapuvien samalla taajuudella värähtelevien aaltojen
# summa kannattaa kuvata vaiheenosoittimella (tai kompleksiluvulla).
# Aaltojen summaa kuvaa nimittäin vaiheenosoitin, joka on saapuvien
# aaltojen vaiheenosoittimien vektorisumma (tai kompleksilukujen summa).
#
# Jos aallon amplitudi on A ja vaihe phi, sitä kuvaavan vaiheenosoittimen
# (tai kompleksiluvun) komponentit ovat
# 
# [ A cos phi, A sin phi ].
#
# Laskemalla näitä komponentteja yhteen, saamme summa-aallon komponentit [Sx, Sy].
# Kun nämä tiedetään, aallon amplitudi saadaan näistä Pythagoraan lauseella.
# Aallon intensiteetti on puolestaan verrannollinen amplitudin neliöön, joka
# on siis komponenttien neliösumma,
#
# I ~ Sx^2 + Sy^2.
#
# kulmat: numpy-vektori, johon on listattu havaintokulmat theta.
# alkupisteet: numpy-vektori, johon on listattu alkupisteet x.
# amplitudit: numpy-vektori, johon on listattu kunkin aallon amplitudi
#             samassa järjestyksessä kuin alkupisteet
# aallonpituus: aallonpituus reaalilukuna.
#
def intensiteettijakauma(kulmat, alkupisteet, amplitudit, aallonpituus):

    # lasketaan vaihematriisi
    #
    # Tämän alkio (i,j) kertoo alkupisteestä x(i) kulmaan theta(j)
    # kulkeneen aallon vaiheen havaintopisteessä.
    vaihekulmat = vaiheet(kulmat, alkupisteet, aallonpituus)

    # lasketaan vaiheenosoittimia kuvaava matriisi
    #
    # Matriisissa on kulmien lukumäärän verran 2-alkioisia vektoreita.
    # Kunkin vektorin komponentit kertovat havaintokulmaan muodostuneen
    # summa-aallon vaiheenosoittimien x- ja y-komponentit
    # (tai kompleksilukuesityksessä reaali- ja imaginäärikomponentit).
    vaiheenosoittimet = np.zeros( [2, N_KULMA] )    
    vaiheenosoittimet[0,:] = amplitudit @ np.cos( vaihekulmat )
    vaiheenosoittimet[1,:] = amplitudit @ np.sin( vaihekulmat )
    
    # lasketaan intensiteettijakauma
    #
    # Jakauma lasketaan vektorina, jonka komponentit kuvaavat
    # intensiteettia kussakin havaintokulmassa.
    #
    # Intensiteetti saadaan vaiheenosoittimen pituuden neliönä.
    amplitudin_nelio = vaiheenosoittimet[0,:]**2 + vaiheenosoittimet[1,:]**2
    amplitudin_nelio /= (N_ALKUPISTE)**2
    
    return amplitudin_nelio

    
#
# Piirretään intensiteetin kuvaaja I(theta).
#
# kulmat: numpy-vektori, johon on listattu havaintokulmat theta.
# intensiteetti: numpy-vektori, johon on listattu intensiteetti kussakin kulmassa I.
#
def piirra_intensiteettijakauman_kuvaaja(kulmat, intensiteetti):

    #
    # TEHTÄVÄ A:
    #
    # Piirrä kuvaaja, jonka vaaka-akselilla kulma ja pystyakselilla intensiteetti.
    # Piirrettävä data on tallennettu vektoreihin 
    # kulmat (vaaka-akseli, x) sekä
    # intensiteetti (pysty-akseli, y).
    #
    # Moduuli matplotlib.pyplot on jo haettu nimellä plt, joten voit kutsua
    # piirtokomentoja kirjoittamalla plt.plot(x,y) tms.
    #
    # Muista näyttää kuvaaja komennolla plt.show() tai kirjoita se tiedostoon
    # komennolla plt.savefig().
    #
    
    #
    # TEHTÄVÄ B:
    #
    # Muokkaa kuvaajaa.
    # - rajaa x-akseli välille [-MAX_KULMA, MAX_KULMA] (plt.xlim)
    # - rajaa y-akseli välille [0, y_maksimi] (plt.ylim)
    # - anna x-akselille nimeksi "havaintokulma, $\\theta$ ($^\circ$)" (plt.xlabel)
    # - anna y-akselille nimeksi "suhteellinen intensiteetti, $I$" (plt.ylabel)
    # - piirrä kuvaaja mustana viivana (plt.plot(..., c="black"))
    # - piirrä kuvaaja paksulla viivalla (plt.plot(..., lw=3))
    # - aseta kuvaajaan taustaviivoitus (plt.grid)
    #
    plt.clf() # tyhjennetään kuva, jos siellä oli jo jotakin muuta
    y_maksimi = 1.1*np.max(intensiteetti) # lasketaan sopiva yläraja y-akselille

    plt.xlim(-MAX_KULMA,MAX_KULMA)
    plt.ylim(0, y_maksimi)
    plt.plot(kulmat, intensiteetti, c="black", lw=2)
    plt.xlabel("havaintokulma, $\\theta$ ($^\circ$)")
    plt.ylabel("suhteellinen intensiteetti, $I$")
    plt.grid()
    plt.savefig("intensiteettijakauma.pdf")
    plt.show()
    
    
#
# Sekoitetaan kaksi väriä keskenään.
#
# Funktio ottaa kaksi RGBA-vektorina määriteltyä väriä ja laskee niiden
# perusteella uuden RGBA-vektorin annetussa suhteessa.
#
# Värit pitää antaa nelikomponenttisena listana tai vektorina, mittä kukin
# komponentti on väliltä 0,...,1. Komponentit ovat järjestyksessä
# red, green, blue, alpha eli 
# punainen, vihreä, sininen, läpinäkymättömyys.
#
# eka: ensimmäinen väri
# toka: toinen väri
# suhde: sekoitussuhde - jos 0, palautetaan eka, jos 1, palautetaan toka
#
def sekoita_vareja(eka, toka, suhde):
    
    sekoitus = [0,0,0,0]
    for i in range(4):
        eka_osa = eka[i]
        toka_osa = toka[i]
        sekoitus[i] = eka_osa*(1-suhde)+toka_osa*suhde
    return sekoitus

    
#
# Valitaan väri aallonpituuden mukaan.
#
# Näkyvän valon väri riippuu sen taajuudesta, joka puolestaan riippuu aallonpituudesta.
# Lyhyin aallonpituus, jonka ihminen näkee, on noin 400 nm, ja se nähdään violettina.
# Pisin aallonpituus on puolestaan noin 800 nm punaisena.
#
# Tämä funktio etsii interferenssikuvion piirtämiseksi värin, joka vastaa simulaatiossa
# käytettyä aallonpituutta. Funktio on varsin karkea eikä tulos ole fysikaalisesti tarkka,
# mutta lopputulos on likimain oikein.
#
# Funktio palauttaa neljän värin listan liukuvärien määrittelyä varten.
# Listan ensimmäinen väri on musta. Kaksi seuraavaa ovat aallonpituutta vastaava väri.
# Viimeinen väri on valkoinen, aallonpituutta vastaava väri tai näiden jokin sekoitus
# aallonpituudesta ja kirkkaudesta riippuen.
#
def varijakauma(aallonpituus):

    musta = (0,0,0,1)
    valkoinen = (1,1,1,1)
        
    if aallonpituus <= PUNA_AALLONPITUUS and aallonpituus >= LILA_AALLONPITUUS:    
        suhde = ( aallonpituus - LILA_AALLONPITUUS ) / (PUNA_AALLONPITUUS - LILA_AALLONPITUUS)
        cmap = matplotlib.colormaps['gist_rainbow_r']
        vari = cmap( suhde )
        huippuvari = valkoinen
        jakauma = [musta, vari, vari, huippuvari]

    elif aallonpituus > PUNA_AALLONPITUUS and aallonpituus < MAX_AALLONPITUUS:
        suhde = 1-( aallonpituus - PUNA_AALLONPITUUS ) / (MAX_AALLONPITUUS - PUNA_AALLONPITUUS)
        vari = (suhde,0,suhde*0.16,1)
        huippuvari = (suhde,suhde,suhde,1)
        huippuvari = sekoita_vareja(vari, huippuvari, kirkkaus)
        huippuvari = sekoita_vareja(huippuvari, vari, suhde)
        jakauma = [musta, vari, vari, huippuvari]
        
    elif aallonpituus < LILA_AALLONPITUUS and aallonpituus > MIN_AALLONPITUUS:
        suhde = ( aallonpituus - MIN_AALLONPITUUS ) / (LILA_AALLONPITUUS - MIN_AALLONPITUUS)
        vari = (suhde,0,0.75*suhde,1)
        huippuvari = (suhde, suhde, suhde, 1)
        huippuvari = sekoita_vareja(vari, huippuvari, kirkkaus)
        huippuvari = sekoita_vareja(huippuvari, vari, suhde)
        jakauma = [musta, vari, vari, huippuvari]

    else:
        jakauma = [musta, musta, musta, musta]

    return jakauma
    
    
#
# Funktio interferenssikuvion pystysuuntaista pehmentämistä varten.
#
# Tämä on normittamaton Gaussin funktio annetulla hajonnalla.
# Funktion keskipiste on kohdassa y = 0.5, jossa funktio saa arvon 1.
#
def pystyjakauma(y, hajonta=Y_HAJONTA):
    if hajonta is None:
        return 1.0
    else:
        return np.exp(-(0.5-y)**2/(2*hajonta**2))
    
    
#
# Funktio interferenssikuvion kirkkauden skaalaamiseksi.
#
# Piirtofunktiot normeeraavat käyttämänsä väriskaalan automaattisesti
# suurimman intensiteetin mukaan. Tämän takia kuvan kirkkain alue
# piirtyy normaalisti aina valkoisena.
#
# Kuitenkin todellisuudessa käy niin, että jos valon intensiteetti on
# pieni, se näyttää kaikkialla himmeätä. Jotta piirretyt kuviot
# toimisivat suunnilleen tähän tapaan, pienen intensiteetin kuvia
# pitää skaalata himmeämmiksi. Tämä funktio laskee tätä varten
# maksimi-intensiteetistä riippuvan kertoimen.
#
def skaalaa_kirkkaus(intensiteetti):
    
    maksimi = np.max(intensiteetti)
    skaalaus = 1.0
        
    if maksimi < 1e-4:
        x = -np.log10(maksimi)
        skaalaus = 1+10*(4-x)**2
        
    return skaalaus
    
    
#
# Piirretään kuva varjostimelle piirtyvästä kuviosta.
#
# Valoton alue piirretään mustana ja valaistu alue likimain sellaisena kuin
# miten valo sen valaisee.
#
# kulmat: numpy-vektori, johon on listattu havaintokulmat theta.
# intensiteetti: numpy-vektori, johon on listattu intensiteetti havaintokulmissa.
# aallonpituus: aallonpituus reaalilukuna.
#
def piirra_interferenssikuvio(kulmat, intensiteetti, aallonpituus):

    # Kuvan piirtämistä varten valmistellaan ensin vaaka- ja pystyakselien
    # arvot sisältävät vektorit.
    xs = kulmat - 0.5*DELTA_KULMA
    xs = np.append(xs, [xs[-1]+DELTA_KULMA])    
    ys = np.linspace(0,1,N_Y)
    eys = pystyjakauma(ys)
    
    # Kootaan sitten intensiteettiä kuvaavat arvot matriisiin.
    # Matriisin alkio (i,j) kertoo intensiteetin pisteessä [x(j),y(i)].
    maksimi = np.max(intensiteetti) * skaalaa_kirkkaus(intensiteetti)    
    z = np.outer( eys[:-1], intensiteetti )
    
    # Valitaan väriskaala aallonpituuden mukaan.
    varit = varijakauma(aallonpituus)

    # Luodaan liukuväriesitys kuvaamaan intensiteettiä nollasta (ei valoa, siis musta)
    # maksimiin (kirkas valo, siis valkoinen tai kirkas väri).
    askeleet = [0, 0.35-0.3*kirkkaus, (1-0.5*kirkkaus), 1] 
    liukuvari = LinearSegmentedColormap.from_list("valo", list(zip(askeleet, varit)))
    
    # Jos intensiteettijakaumassa on hyvi korkea piikki ja väritys skaalataan sen
    # mukaan, himmeämpiä yksityiskohtia kuten diffraktiosivumaksimeja ei näe.
    # Nämä saadaan näkyviin asettamalla väriskaalan yläraja maksimia pienempään
    # intensiteettiin.
    #
    # Vastaavasti jos kirkkaus on pieni, yläraja halutaan maksimi-intensiteettiä ylemmäs,
    # jotta mikään alue ei piirry kovin kirkkaana.
    #
    # Tämä skaalaustekijä huolehtii, että näin tehdään.
    # Katkaisuraja asetetaan maksimiin, kun kirkkaus = 0.5.
    # Kun kirkkaus -> 1, katkaisuraja -> 0 (yhä suurempi osuus piirretään valkoisena).
    # Kun kirkkaus -> 0, katkaisuraja -> ääretön (yhä suurempi osuus mustana). 
    katkaisuraja = (1/kirkkaus-1)*maksimi
    
    # Piirretään kuva.
    plt.clf()
    p=plt.pcolormesh(xs,ys,z, cmap=liukuvari, vmin=0, vmax=katkaisuraja)
    plt.xlabel("havaintokulma, $\\theta$ ($^\circ$)")
    plt.savefig("interferenssikuvio.png")
    plt.show()


#
# Piirretään kuva aukoista, joiden läpi valo kulkee.
#
# Kuvaan piirretään aukot valkoisena ja valon pysäyttävä alue mustana.
#
# paikat: numpy-vektori, johon on listattu kaikki alkupisteet.
# aukot: numpy-vektori, johon on listattu kutakin paikkaa kohti,
#        onko pisteessä aukko (1) vai este (0)
#
def piirra_raot(paikat, aukot):
    
    # Kuvan piirtämistä varten valmistellaan ensin vaaka- ja pystyakselien
    # arvot sisältävät vektorit.
    xs = paikat - 0.5*DELTA_ALKUPISTE
    xs = np.append(xs, xs[-1]+DELTA_ALKUPISTE)
    xs = xs/1000
    
    ny = 10
    ys = np.linspace(0,10,ny)
    
    # Kootaan sitten aukkoja kuvaavat arvot matriisiin.
    # Matriisin alkio (i,j) kertoo onko pisteessä [x(j),y(i)] aukko vai este.
    z = np.zeros([ny-1,N_ALKUPISTE])
    for j in range(ny-1):
        if j > 0 and j < ny-2:
            z[j,:] = aukot[:]

    # piirretään kuva
    plt.clf()
    varit = ["black",  "white"]
    askeleet = [0, 1]
    varikartta = LinearSegmentedColormap.from_list("aukko", list(zip(askeleet, varit)))
    p=plt.pcolormesh(xs,ys,z, cmap=varikartta)
    plt.xlabel("raot, ($\mu$m)")
    plt.savefig("raot.png",dpi=600)
    plt.show()
    
    
# 
# Kertoo, onko kysytyssä paikassa aukko.
#
# Funktiolle annetaan paikkakoordinaatti ja tiedot rakojen paikoista.
# Funktio palauttaa 1, jos kysytty koordinaatti osuu rakoon, ja 0 muuten.
#
# paikka: x-koordinaatti
# keskipisteet: lista tai vektori rakojen keskipisteiden koordinaateista.
# leveydet: rakojen leveydet lukuna tai listana, jos leveydet eivät ole samat.
#
def on_aukko(paikka, keskipisteet, leveydet):

    # jos leveys oli yksi luku, tehdään siitä vektori, jossa luku toistuu
    # yhtä monta kertaa kuin rakoja on
    if isinstance(leveydet,int) or isinstance(leveydet,float):
        leveydet = np.ones( len(keskipisteet) )*leveydet
            
    # käydään läpi kaikki raot
    for k, l in zip(keskipisteet, leveydet):

        # jos paikka on lähellä raon keskipistettä, palautetaan 1
        if abs(paikka-k) <= l/2:
            return 1

    # jos paikka oli tarpeeksi kaukana kaikista raoista, palautetaan 0
    return 0
    

#
# Luo vektorin, joka kertoo kustakin pisteestä lähtevän aallon amplitudin.
#
# Funktiolle annetaan kaikkien mahdollisten alkupisteiden vektori sekä tiedot
# rakojen paikoista. Se luo uuden vektorin, jossa on yhtä monta alkiota kuin
# mahdollisia alkupisteitä on. Tämän vektorin kukin alkio kertoo vastaavasta
# pisteestä lähteneen aallon amplitudin.
#
# Jos ko. pisteessä on aukko, josta valo pääsee kulkemaan, amplitudiksi tulee 1.
# Jos pisteessä ei ole aukkoa, sieltä ei voi tulla valoa, ja amplitudi on 0.
#
# alkupisteet: numpy-vektori, johon on listattu alkupisteet x.
# keskipisteet: lista tai vektori rakojen keskipisteiden koordinaateista.
# leveydet: rakojen leveydet lukuna tai listana, jos leveydet eivät ole samat.
#
def generoi_amplitudijakauma( alkupisteet, keskipisteet, leveydet ):

    # luodaan sopivan pituinen vektori täynnä nollia
    amplitudit = np.zeros(N_ALKUPISTE)

    # käydään läpi kaikki mahdolliset alkupisteet
    for i in range(N_ALKUPISTE):
    
        # tarkastetaan, onko tässä pisteessä aukko, ja
        # aukkokohtiin asetetaan arvoksi 1
        amplitudit[i] = on_aukko( alkupisteet[i], keskipisteet, leveydet )
        
    return amplitudit
    

#
# Muuntaa rakojen reunapisteiden vektorin niiden keskipisteiden ja leveyksien vektoreiksi.
#
# Ohjelma käsittelee rakoja kahtena vektorina:
# - yhdessä on rakojen keskipisteiden koordinaatit
# - toisessa on rakojen leveydet
#
# Raot voisi kuitenkin määritellä myös antamalla vektorina rakojen
# reunapisteiden koordinaatit. Tämä funktio laskee annettuja reunapisteitä
# vastaavat keskipisteet ja leveydet.
#
# reunat: vektori, johon on listattu rakojen reunapisteiden koordinaatit.
#
def laske_aukot_reunoista(reunat):

    n_aukot = len(reunat) // 2
    
    keskikohdat = np.zeros(n_aukot)
    leveydet = np.zeros(n_aukot)
    
    for i in range(n_aukot):
        vasen = reunat[2*i]
        oikea = reunat[2*i+1]
        
        keskikohdat[i] = 0.5*(vasen+oikea)
        leveydet[i] = oikea-vasen
        
    return keskikohdat, leveydet
    

# 
# Arpoo fotonin osumakohdan annetun intensiteettijakauman mukaisesti.
#
# Kvanttimekaanisessa mallissa valo kulkee raoista sähkömagneettisena aaltona,
# mutta aallon vuorovaikuttaessa aineen kanssa varjostimella se paikallistuu 
# johonkin tiettyyn paikkaan. Tämä on satunnaisprosessi, jossa paikan
# todennäköisyysjakauma on aallon intensiteettijakauman.
#
# Tämä funktio valitsee satunnaisesti paikkoja mahdollisten havaintopisteiden
# joukosta. Pisteet valitaan kaksiulotteiselta tasolta.
# Pisteiden x-koordinaattien jakauma noudattaa annettua intensiteettiä.
# Pisteiden y-koordinaattien jakauma noudattaa normaalijakaumaa.
# Funktio palauttaa x- ja y-koordinaatit kahtena erillisenä vektorin.
#
# intensiteetti: numpy-vektori, johon on listattu intensiteetti havaintokulmissa.
# lukumäärä: näin monta pistettä arvotaan.
# 
def arvo_paikka_intensiteettijakaumasta(intensiteetti, lukumaara):

    # käytetään numpyn satunnaislukueneraattoria
    rng = np.random.default_rng()
    
    # mahdollisten kulmien arvot
    mahdolliset_x = np.linspace( -MAX_KULMA, MAX_KULMA, N_KULMA )

    # x-suuntainen todennäköisyystiheys saadaan normittamalla intensiteettijakauma
    #
    # Normitus tarkoittaa jakauman skaalaamista niin, että kaikkien todennäköisyyksien
    # summaksi tulee yksi.
    #
    fx = intensiteetti / np.sum(intensiteetti)

    # arvotaan pyydetty lukumäärä satunnaisia x-arvoja (havaintokulmia)
    X = rng.choice(mahdolliset_x, p=fx, size=lukumaara)
    
    # Lisätään arvoihin vielä pientä satunnaiskohinaa.
    #
    # Jos tätä ei tee, vektoriin X voi tulla vain vektoriin mahdolliset_x
    # listatut lukuarvot. Erityisesti jos arvotaan paljon pisteitä, moni piste
    # alkaa osua täsmälleen samaan x-koordinaattiin, mitä ei pitäisi käydä.
    #
    # Oikeasti satunnaisluvut pitäisi arpoa jatkuvasta jakaumasta, mutta
    # sen toteuttaminen vaatisi kertymäfunktion laskemisen ja siis hiukan enemmän
    # työtä. Tässä käytetty menetelmä on hyvin helppo ja varsin hyvä approksimaatio.
    X += rng.normal(scale=DELTA_KULMA, size=lukumaara)
    
    # y-suunnassa mahdolliset arvot ovat välillä 0,1
    mahdolliset_y = np.linspace(0,1,N_Y)
    
    # y-suunnassa pisteet asettuvat funktion pystyjakauma määräämällä tavalla
    fy = pystyjakauma(mahdolliset_y)

    # todennäköisyystiheyden normitus
    fy /= np.sum(fy)
    
    # arvotaan y-arvot ja lisätään niihinkin vähän satunnaispoikkeamaa.
    Y = rng.choice(mahdolliset_y, p=fy, size=lukumaara)
    Y += rng.normal(scale=1/N_Y, size=lukumaara)
    
    return X, Y
    
    
#
# Piirtää kuvan intensiteettijakaumaa vastaavista fotoneista.
#
# Valon vuorovaikutus aineen kanssa noudattaa kvanttimekaanisia sääntöjä.
# Valo on klassisesti sähkömagneettinen aalto, ja tämä malli kuvaa esimerkiksi
# valon kulun rakojen läpi täysin oikein. Kuitenkin kun valo sitten osuu johonkin pintaan,
# se ei itse asiassa valaise sitä tasaisesti kuten meistä näyttää.
# Sen sijaan valoaallon vuorovaikutus aineen kanssa saa valon luovuttamaan energiaansa
# pinnalle tietyissä pisteissä. Koska vuorovaikutus on paikallinen, valo näyttää siinä
# käyttäytyvän kuin hiukkanen. Tällöin puhutaan fotoneista.
#
# Jos valon intensiteetti on suuri, näitä paikallisia vuorovaikutuksia on hyvin paljon
# eikä niitä erota. Silloin valo näyttää valaisevan pintoja tasaisesti ja klassinen malli
# toimii oikein hyvin.
#
# Jos valon intensiteetti on pieni, yksittäiset vuorovaikutuspisteet alkavat erottua.
# Silloin valon muodostamasta kuviosta tulee rakeinen.
#
# Tämä funktio piirtää kuvan valon piirtämästä kuviosta tässä jälkimmäisessä tapauksessa.
#
# intensiteetti: numpy-vektori, johon on listattu intensiteetti havaintokulmissa.
#
def piirra_fotonikuvio(intensiteetti):

    # Arvotaan fotonien osumapisteitä annetun intensiteettijakauman mukaisesti.
    #
    # Valitaan piirrettävien pisteiden määrä valon intensiteetin mukaaan.
    # Jos intensiteetti on suuri, piirretään enemmän pisteitä.
    # Tämä skaalaus ei ole fysikaalisesti oikein vaan se tehdään vain
    # jotta kuvista tulisi hyvän näköisiä.
    skaalaus = skaalaa_kirkkaus(intensiteetti)
    N = int(150000 * (1/(1-kirkkaus) - 1) / skaalaus )
    X, Y = arvo_paikka_intensiteettijakaumasta(intensiteetti, N)
    
    # piirrettävien pisteiden ominaisuuksia
    pisteen_koko = 800 / N_KULMA
    pisteen_vari = (1,1,1,0.2)


    # tyhjennetään kuva varmuuden vuoksi
    plt.clf()

    #
    # TEHTÄVÄ C:
    #
    # Piirrä kuva fotoneista.
    # Fotonien koordinaatit on tallennettu vektoreihin X ja Y,
    # joten kuvion piirtäminen onnistuu tavallisella plot-komennolla,
    # kunhan kuva piirretään nimenomaan pisteinä eikä murtoviivana.
    # Tämä onnistuu antamalla plot-funktiolle kolmantena argumenttina 'o'.
    #
    # Lisäksi jotta kuvasta saataisiin järkevän näköinen, pisteiden koko ja väri
    # tulisi asettaa muuttujien pisteen_koko ja pisteen_vari mukaisesti ja
    # kuvan taustaväri mustaksi.
    # Pisteiden muuttaminen onnistuu antamalla plot-funktiolle avainsana-argumentit
    # c = haluttu väri (c, color) sekä
    # ms = haluttu koko (ms, markersize)
    # eli tässä tapauksessa siis
    # plt.plot(..., c = pisteen_vari, ms = pisteen_koko).
    #
    # Taustavärin muuttaminen on vähän hankalampaa. Se onnistuu niin, että
    # ensin pyydetään matplotlibiltä kuvaajan akselien tiedot funktiolla 
    #
    # akselit = plt.gca() (gca tarkoittaa 'get current axes')
    #
    # ja sitten muutetaan sieltä taustaväri komennolla
    #
    # akselit.set_facecolor("black").
    #
    # Lopuksi muokaa kuvaa vielä seuraavasti:
    # - aseta vaakasuunnassa piirtorajoiksi [-MAX_KULMA, MAX_KULMA] (plt.xlim)
    # - aseta pystysuunnassa piirtorajoiksi [0, 1] (plt.ylim)
    # - aseta vaaka-akselin nimeksi "havaintokulma, $\\theta$ ($^\circ$)" (plt.xlabel)

    plt.plot(X,Y,'o', c=pisteen_vari, ms=pisteen_koko)
    plt.xlim(-MAX_KULMA,MAX_KULMA)
    plt.ylim(0,1)
    plt.xlabel("havaintokulma, $\\theta$ ($^\circ$)")  
    ax = plt.gca()
    ax.set_facecolor('black') 
    plt.savefig("fotonikuvio.png",dpi=300) 
    plt.show()


#
# Valitsee kokeen valmiista listasta.
#
# Käyttäjä voi simuloida interferenssikuvion ihan millaisilla raoilla haluaa, 
# mutta tähän funktioon on määritelty valmiiksi joukko erityisen mielenkiintoisia
# tapauksia.
#
# tapaus: kokonaisluku 1, ..., 13 valitsee kokeen tyypin
# w: rakojen leveys tietyissä kokeissa
# d: rakojen välinen etäisyys tietyissä kokeissa
#
def valmistele_koeasetelma(tapaus, w, d):

    if tapaus == 1: # kaksoisrako
        leveydet = w 
        keskipisteet = [-d/2, d/2] 
    
    elif tapaus == 2: # kaksoisrako tuplaetäisyydellä
        leveydet = w 
        keskipisteet = [-d, d] 
    
    elif tapaus == 3: # neloisrako
        leveydet = w 
        keskipisteet = np.array([-3,-1,1,3])*d/2 

    elif tapaus == 4: # hila
        leveydet = w # nm
        nn = 500
        keskipisteet = np.linspace(-nn,nn,nn+1)*d/2 

    elif tapaus == 5: # yksi kapea rako
        leveydet = w
        keskipisteet = [0]
    
    elif tapaus == 6: # yksi melko leveä rako
        leveydet = d/2
        keskipisteet = [0]

    elif tapaus == 7: # yksi tuplaleveä rako
        leveydet = d
        keskipisteet = [0]
    
    elif tapaus == 8: # yksi hyvin leveä rako
        leveydet = 2*MAX_ALKUPISTE
        keskipisteet = [0]
    
    elif tapaus == 9: # kaksi melko leveää rakoa
        leveydet = d/2
        keskipisteet = [-d/2,d/2]
    
    elif tapaus == 10: # kaksi melko leveää rakoa tuplaetäisyydellä
        leveydet = d/2
        keskipisteet = [-d, d]
    
    elif tapaus == 11: # leveden rakojen hila
        leveydet = d/2
        keskipisteet = np.linspace(-100,100,51)*d/2
    
    elif tapaus == 12: # sekalaisia rakoja
        leveydet = [800, 1200, 1600]
        keskipisteet = [-4000, 0, 6000]
    
    elif tapaus == 13: # voi antaa myös rakojen reunapisteet
        reunat = [-6000,-3000, -2000,1000, 4000,5500]
        keskipisteet, leveydet = laske_aukot_reunoista(reunat)
    
    return keskipisteet, leveydet
    
    
#
# Pääohjelma alkaa tästä
#
   
# valitaan koeasetelma tiedoston alussa annettujen parametrien mukaisesti 
rakojen_keskipisteet, rakojen_leveydet = valmistele_koeasetelma(tapaus, w, d)    

# luodaan vektori, johon on koottu valon lähtöpisteiden x-koordinaatit
alkupisteet = np.linspace( -MAX_ALKUPISTE, MAX_ALKUPISTE, N_ALKUPISTE )

# luodaan vektori, johon on koottu valon amplitudi lähtöpisteissä
amplitudit = generoi_amplitudijakauma(alkupisteet, rakojen_keskipisteet, rakojen_leveydet)

# luodaan vektori, johon on koottu valon interferenssikuvioon laskettavat kulmat
havaintokulmat = np.linspace( -MAX_KULMA, MAX_KULMA, N_KULMA )

# piirretään kuva raoista, joiden läpi valo kulkee
piirra_raot(alkupisteet, amplitudit)

# lasketaan intensiteetti havaintokulmissa
intensiteetti = intensiteettijakauma(havaintokulmat, alkupisteet, amplitudit, aallonpituus)

# piirretään kuvaaja intensiteettijakaumasta
piirra_intensiteettijakauman_kuvaaja(havaintokulmat, intensiteetti)

# piirretään kuva valon piirtämästä interferenssikuviosta
piirra_interferenssikuvio(havaintokulmat, intensiteetti, aallonpituus)

# piirretään kuva valon piirtämästä kuviosta, jos intensiteetti on niin pieni, 
# että yksittäiset fotonit erottuvat
piirra_fotonikuvio(intensiteetti)
