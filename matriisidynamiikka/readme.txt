
Jatkuvien funktioiden simulointi

Teemu Hynninen 2024
CC-SA-BY-4.0

Tämä ohjelma simuloi jatkuvien suureiden dynamiikkaa.

Jatkuvia suureita ovat fysiikassa esimerkiksi lämpötila, konsentraatio ja
aaltofunktio. Näitä kuvataan fysiikan teoriassa osittaisdifferentiaaliyhtälöillä.
Tietokoneella niitä ei voi kuitenkaan käsitellä aivan sellaisenaan vaan ne on 
muutettava muotoon, jota voi käsitellä numeerisesti. Tähän on olemassa monia
menetelmiä, mutta käytännössä funktio pitää muuttaa jonkinlaiseksi lukujen
kokoelmaksi - siis vektoriksi v - ja yhtälöt pitää muuttaa säännöksi, joka
kertoo miten nämä luvut muuttuvat - matriisiksi M. Käytännössä systeemin tila
tulevaisuudessa lasketaan säännöllä 
v(t+dt) = Mv(t) - kv(t-dt), 
missä k on vakio: 1 aalloille ja 0 lämpötilalle.

Ohjelmassa on graafinen ikkuna, joka kuvaa simuloitavan suureen
tasossa. Positiiviset arvot piirretään punaisena, negatiiviset sinisenä ja
nolla valkeana. Aluksi funktio on nolla, mutta hiirellä klikkaamalla käyttäjä
voi luoda suureeseen pulsseja (ykkösnapilla positiivisia, kakkosnapilla 
negatiivisia).

Ohjelma lukee erillisestä tiedostosta suureen dynamiikkaa kuvaavan matriisin.
Tarjolla on seitsemän mahdollista valmista tiedostoa.
Itse ohjelma on aina täsmälleen sama, mutta vaihtamalla dynamiikkamatriisia
suureelle voidaan valita esimerkiksi lämpötilan tai aaltojen dynamiikka.
Matriisi vastaa myös siitä, mitä systeemin reunoilla tapahtuu.

Simulaatio ajetaan suorittamalla tiedostossa dynamiikka.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Matriisin määrittelevän tiedoston nimi on tallennettu muuttujaan
DYNAMIIKKA, ja systeemin fysiikkaa muutetaan vaihtamalla tiedostoa.

