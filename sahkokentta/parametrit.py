import numpy as np

nC = 1.0 * 10**-9 ; # tallennetaa "nano" muuttujaan

#
# lista systeemissä olevista pistevarauksista
#
# Listan pitää olla muotoa
#
# [
# [ q1, (x1, y1) ],
# ...
# [ qn, (xn, yn) ]
# ]
#
# missä 
# q1, ..., qn ovat varausten suuruudet coulombeina
# x1, ..., xn ovat varausten x-koordinaatit metreinä ja
# y1, ..., yn ovat vastaavat y-kordinaatit.
# 
VARAUKSET = [ 
[  1.0*nC, (  0.10,  0.20 )],
[ -0.5*nC, (  0.30, -0.40 )]
]

#
# parametrejä kuvien piirtämistä varten
#
# Ohjelma piirtää systeemistä vektorikuvaajan
# neliön muotoiselta alueelta xy-tasosta.
# Kaikki mitat ovat metreissä.
#
X_KESKUS = 0.0                     # kuvan keskipisteen x-koordinaatti
Y_KESKUS = 0.0                     # kuvan keskipisteen y-koordinaatti
KUVAN_LEVEYS = 1.00                # kuvan leveys
NUOLTEN_ETAISYYS = 0.02            # sähkökenttänuolten välinen etäisyys kuvassa
NUOLTEN_MITTAKAAVA = 4.0*10**-5    # nuolten pituusskaala (pieni: lyhyitä nuolia)
SUURIN_POTENTIAALI =  500          # potentiaalin tasa-arvokäyrien maksimi
PIENIN_POTENTIAALI = -500          # potentiaalin tasa-arvokäyrien minimi
POTENTIAALIN_RESOLUUTIO = 201      # potentiaalin tarkkuus (suuri: tarkempi mutta lasku kestää kauemmin)
SAHKOKENTTA_POTENTIAALISTA = False # Lasketaanko sähkökenttä potentiaalista (jos ei, lasketaan suoraan pistevarauksista)
PISTEIDEN_KOKO = 0.005             # kuinka suurina varaukset ja tarkastelupisteet merkitään kuvaan
KUVAFORMAATTI = "pdf"              # kuvatiedoston formaatti (pdf tai png)
PIIRRA_SAHKOKENTTA = True          # Piirretäänkö sähkökenttä?
PIIRRA_POTENTIAALI = False         # Piirretäänkö potentiaali?

#
# Ohjelma piirtää myös sähkökentän ja potentiaalin kuvaajat suoralla.
# Aseta tässä suoran alku- ja loppupisteen koordinaatit.
#
KUVAAJAN_PAATEPISTEET = [ (0.0, 0.0), (0.0, 0.25) ]

#
# Ohjelma laskee sähkökentän ja potentiaalin näissä pisteissä.
# Pisteet merkitään myös kuvaan.
#
TARKASTELUPISTEET = [
(  0.00,  0.10 ),
( -0.20, -0.30 )
]


import sahkokentta as sahko

#
# Tehtävä 1: toteuta funktio, joka laskee systeemin sähkökentän
# annetussa pisteessä pistevarausten kenttien summana.
# Argumentti 'varaukset' on lista pistevarauksista (kuten VARAUKSET).
# Argumentit x ja y ovat tarkastelupisteen koordinaatit.
#
# Käytössä on jo valmiiksi toteutettu funktio
# sahko.pistevarauksen_sahkokentta( q, xq, yq, x, y ),
# joka laskee pisteessä (xq,yq) olevan varauksen q
# sähkökentän pisteessä (x,y).
# Funktion palauttaa tuloksenaan vektorin [Ex,Ey].
#
def sahkokentta( varaukset, x, y ):

    # Luodaan muuttuja, johon sähkökenttä tallennetaan
    # ja alustetaan se nollavektoriksi [0.0, 0.0].
    E_kokonais = np.zeros(2)
     
     # Käydään silmukassa läpi kaikki varaukset.
     # Tässä muuttuja q saa vuorollaan arvokseen
     # kunkin varauksen suuruuden ja muuttujat
     # (xq, yq) varauksen koordinaatit
    for q, (xq, yq) in varaukset:
        
        # Tämä laskee yhden pistevarauksen kentän pisteessä (x,y),
        # mutta ei tee sillä mitään. Korjaa asia!
        sahko.pistevarauksen_sahkokentta(q, xq, yq, x, y)

    # Palautetaan lopuksi muuttujaan E_kokonais tallennettu tulos.
    # Aluksi se on nolla, koska muuttujalle ei tehty mitään.
    return E_kokonais


#
# Tehtävä 2: toteuta funktio, joka laskee systeemin sähkökentän
# annetussa pisteessä potentiaalista.
# Argumentit x ja y ovat tarkastelupisteen koordinaatit.
#
# Käytössä on jo valmiiksi toteutettu funktio
# sahko.potentiaali( x, y )
# joka laskee potentiaalin pisteessä (x,y).
# Funktio palauttaa tuloksenaan skalaarin V.
#
def sahkokentta_potentiaalista( x, y ):

    # Luodaan muuttuja, johon sähkökenttä tallennetaan
    # ja alustetaan se nollavektoriksi [0.0, 0.0].
    E_kokonais = np.zeros(2)
     
    # Tämä laskee potentiaalin pisteessä (x,y),
    # mutta ei tee sillä mitään. Korjaa asia!
    sahko.potentiaali( x, y )
     
    # Palautetaan lopuksi muuttujaan E_kokonais tallennettu tulos.
    # Aluksi se on nolla, koska muuttujalle ei tehty mitään.
    return E_kokonais
        

