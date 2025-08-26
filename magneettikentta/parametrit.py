#
# Kokeen fysiikkaan liittyviä parametrejä.
# Kaikki parametrit on annettu SI-perusyksiköissä.
# Esim. paikkakoordinaatit ovat siis metreissä.
#
VIRTA  = 10.0    # sähkövirta

# virtasilmukan kulmien koordinaatit
# järjestyksessä virran kulkusuunnassa
#
# tämä esimerkki luo neliön muotoisen silmukan
KULMAT = [
(  0.001,       0, 0 ),
(      0,   0.001, 0 ),
( -0.001,       0, 0 ),
(      0,  -0.001, 0 )
]

#
# Kuvaajien piirtoon liittyviä parametrejä
#

# kentän komponenttien kuvaajat piirretään suoraa pitkin 3D-avaruudessa
ALKUPISTE     = ( 0.0, 0.0, 0.0 )  # piirtosuoran alkupiste
LOPPUPISTE    = ( 0.0, 0.0, 0.01 ) # piirtosuoran loppupiste

# kenttä piirretään vektorikuvaajana xz-tasossa
B_SKAALA      = 50     # nuolten mittakaava vektorikuvaajassa
X_MAKSIMI     = 0.01   # piirtoalueen koko vektorikuvaajassa (maksimietäisyys origosta)
KUVAFORMAATTI = "pdf"  # kuvatiedoston tyyppi (png tai pdf)

# OSIEN_KOKO määrittelee, kuinka pieniin osiin silmukka jaetaan laskun aikana.
# Pieni arvo antaa tarkan tuloksen mutta silloin lasku kestää pitkään.
# Suuri arvo antaa epätarkan tuloksen nopeasti.
# Hyvä kompromissi on esim. noin kymmenesosa silmukan sivun pituudesta.
OSIEN_KOKO = 10**-4
