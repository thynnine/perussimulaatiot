
Kovien pallojen simulaatio

Teemu Hynninen 2024
CC-SA-BY-4.0

Tämä ohjelma simuloi kovien pallojen törmäyksiä systeemissä, joka näyttää
hiukan snooker-pöydältä.

Simulaatio on graafinen ja sitä ohjataan hiirellä (tai näppäimistöllä).
Pelaaja voi lyödä valkoista palloa. Lyönti tapahtuu valkoisen ympyrän suuntaan.
Ympyrä seuraa hiirtä tai sitä voi siirtää nuolinäppäimillä.
Yksi klikkaus hiirellä tai välilyönnin painaminen aloittaa lyönnin.
Se näkyy siten, että valkoinen ympyrä alkaa pienentyä.
Toinen klikkaus lähettää valkoisen pallon matkaan, ja se saa sitä suuremman
nopeuden, mitä pienempi valkoinen ympyrä oli. Jos ympyrä kutistuu pisteeksi,
lyönti lähtee automaattisesti maksiminopeudella.

Simulaatio perustuu kovien pallojen malliin. Simulaattori laskee, mitkä palloista
törmäävät seuraavaksi ja kauanko törmäykseen on aikaa. Sitten se antaa ajan kulua
(pienissä askeleissa animaation piirtämistä varten) kunnes tämä törmäys tapahtuu,
ja laskee pallojen nopeudet törmäyksen jälkeen olettaen törmäykset täysin
elastisiksi.

Oletusarvoisesti simulaatiossa siis pelataan biljardia.
Käyttäjä voi kuitenkin poistaa simulaatiosta dissipatiiviset vuorovaikutukset, jolloin
pallot eivät koskaan pysähdy. Lisäksi pallojen kokoa ja määrää voi muuttaa.
Kun systeemiin laittaa paljon pieniä ikuisesti liikkuvia palloja, siitä tuleekin 
kaksiulotteinen kaasu. 

Tällaiselle simulaatiolla käyttäjä voi laskea ja piirtää vauhtijakauman painamalla 
näppäintä 'p' (niin kuin plot). Kolmiulotteisen kaasun hiukkasten vauhdit noudattavat 
Maxwell-Boltzmann-jakaumaa, ja kaksiulotteisen kaasun hiukkaset samankaltaista jakaumaa
f(v) = a v exp( -a/2 v^2 ).
Ohjelma piirtää sekä tämän teoreettisen jakauman että histogrammin simulaatiossa olevien
pallojen hetkellisistä vauhdeista.

Simulaatio kirjoittaa valkoisen pallon xy-koordinaatit ja kokonaisenergian ajan 
funktiona tiedostoihin.

Simulaatio ajetaan suorittamalla tiedostossa pallot.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Simulaation ominaisuuksia voi kuitenkin muuttaa tiedoston alussa olevien
parametrien avulla. Nämä ovat:

SKAALA: montako pikseliä 1 cm etäisyys on kuvassa
R: pallojen säde (cm)
LEVEYS: systeemin leveys (cm)
KORKEUS: systeemin korkeus (cm)
DT: animaation kuvien aikaväli (ms)
MAX_VOIMA: suurin mahdollinen nopeus valkoiselle pallolle
DELTA_VOIMA: vaikuttaa siihen, kuinka nopeasti lyönnin voimakkuus latautuu
KITKA: dissipatiivisten vuorovaikutusten (vierimisvastus, ilmanvastus) voimakkuus
KERROKSIA: montako kerrosta punaisia palloja alkukolmiossa on
ERO: pallojen väliin jäävä tila alkukolmiossa
KARKI: alkukolmion kärjen paikka suhteessa pöydän leveyteen (esim. 0.5 tarkoittaa keskellä pöytää)
VARI: pallojen väri rgb-heksaformaatissa (jos puna-vihreä ei erotu, muuta tätä)

Kokeita

Jos haluat pelata biljardia, aseta:
R = 3 * SKAALA 
KITKA = 5.0 
KERROKSIA = 5 
ERO = 0 * SKAALA 
KARKI = 0.75 

Jos haluat simuloida kaksiulotteista ideaalikaasua, aseta:
R = 1 * SKAALA
KITKA = 0 
KERROKSIA = 15
ERO = 2 * SKAALA
KARKI = 0.5