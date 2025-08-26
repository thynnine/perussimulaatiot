
Sähkökenttä

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi varausten luoman sähkökentän ja potentiaalin.
Se koostuu tiedostoista
- sahkokentta.py
- parametrit.py

Simulaatio laskee annettujen pistevarausten sähkökentän ja potentiaalin
xy-tasossa.

Simulaatio ajetaan suorittamalla tiedostossa sahkokentta.py oleva ohjelma.
Tässä tiedostossa oleva ohjelma toimii eikä sitä tarvitse muuttaa.
Ohjelmaan on kuitenkin määritelty funktiot vain yksittäisten
pisteiden sähkökentän laskemiseen.
Koko systeemin sähkökentän laskemista varten on määritelty
tiedostoon parametrit.py kaksi funktiota:

Funktion sahkokentta( varaukset, x, y ) on tarkoitus laskea sähkökenttä
summaamalla kaikkien pistevarausten kentät yhteen.
Funktioon on kirjoitettu valmiiksi silmukka, joka käy läpi kaikki
varaukset, sekä funktiokutsu pistevarauksen_sahkokentta(), joka
laskee kunkin yksittäisen varauksen kentän.
Näitä kenttiä ei kuitenkaan tallenneta mihinkään, ja käyttäjän tehtävänä
on täydentää ohjelma niin, että se laskee kaikki nämä sähkökenttävektorit
yhteen.

Funktion sahkokentta_potentiaalista( x, y ) on tarkoitus laskea sähkökenttä
systeemin potentiaalifunktion avulla. Tämä funktio ei saa argumentiksi
varausten listaa, joten yksittäisten varausten sähkökenttien summaaminen
ei nyt ole mahdollista. Sen sijaan on tarkoitus käyttää funktiota
potentiaali(), joka laskee systeemin potentiaalin annetussa pisteessä.
Käyttäjän tehtävä on selvittää, miten sähkökenttä ja potentiaali liittyvät
toisiinsa ja toteuttaa numeerisesti sähkökentän x- ja y-komponenttien 
likimääräinen lasku potentiaalin kautta.

Tiedostossa parametrit.py määritellään myös joukko muita muuttujia:

VARAUKSET: pistevaraukset (C) ja niiden koordinaatit (m)
 (Varaukset pitää antaa listojen listana, missä kukin alkio kertoo
  varauksen suuduuren sekä sen koordinaatit xy-tasossa.
  Esim. VARAUKSET = [ [1, (2,3)], [4, (5,6)] ]
  tarkoittaisi sitä, että 1 C varaus on pisteessä (x,y) = (2 m, 3 m)
  ja 4 C varaus on pisteessä (5 m, 6 m).)
X_KESKUS: systeemistä piirretyn kuvan keskipiste x-akselilla (m)
Y_KESKUS: systeemistä piirretyn kuvan keskipiste y-akselilla (m)
KUVAN_LEVEYS: systeemistä piirretyn kuvan leveys (m)
 (Kuva on neliön muotoinen, joten tämä on myös korkeus.)
NUOLTEN_ETAISYYS: sähkökenttänuolten keskipisteiden etäisyys kuvassa (m)
 (Sähkökenttä lasketaan neliöhilasta valituissa pisteissä, ja kuhunkin
  pisteeseen piirretään siinä olevaa kenttää kuvaava nuoli. Tämä arvo
  määrää, kuinka tiheässä nämä pisteet ovat. Huom. jos tähän laittaa
  hyvin pienen arvon, laskettavaa on paljon ja lasku kestää pitkään.)
NUOLTEN_MITTAKAAVA: sähkökenttänuolten mittakaava
 (Sähkökentän yksikkö on V/m, mutta systeemin kuva piirretään koordinaatistoon,
  jossa pituudet ovat metrejä. Sähkökenttää kuvaavat nuolet voidaan siis
  piirtää millaiseen mittakaavaan tahansa. Tämä arvo säätää nuolten kokoa.
  Jos arvo on esim. 0.001, 1 V/m sähkökenttä piirretään 0.001 m pituisena nuolena.
  Huom. nuolten maksimipituus on muuttujan NUOLTEN_ETAISYYS arvo, vaikka 
  sähkökenttä olisi voimakkaampi ja nuoli pitäisi piirtää tätä pidempänä.)
SUURIN_POTENTIAALI: potentiaalin tasa-arvokäyrien maksimi (V)
PIENIN_POTENTIAALI: potentiaalin tasa-arvokäyrien minimi (V)
POTENTIAALIN_RESOLUUTIO: potentiaalin laskentapisteiden määrä
 (Kuten sähkökenttä, myös potentiaali lasketaan neliöhilasta poimituissa pisteissä.
  Ohjelma laskee potentiaalin ensin näissä pisteissä ja tallentaa tulokset.
  Kun siltä sen jälkeen kysyy potentiaalin arvoa missä tahansa pisteessä,
  se interpoloi tuloksen. Tämä muuttuja määrää, kuinka monta pistettä x- ja y-
  akseleiden suunnassa on. Mitä suurempi arvo, sitä tarkemmin potentiaalin
  interpolaatio toimii, mutta sitä pidempään taulukon luominen kestää.)
SAHKOKENTTA_POTENTIAALISTA: määrää, käytetäänkö sähkökentän laskemiseen potentiaalia
 (Kuten yllä kerrottiin, koko systeemin sähkökentän laskemiseen on tiedostossa
  parametrit.py kaksi funktiota, joita kumpaakin pitää hiukan täydentää.
  Jos tämä muuttuja on tosi (True), ohjelma käyttää funktiota sahkokentta(),
  kun sen täytyy laskea sähkökenttä.
  Jos muuttuja on epätosi (False), käytetään funktiota sahkokentta_potentiaalista().
PISTEIDEN_KOKO: varausten ja tarkastelupisteiden säde (m)
 (Pistevaraukset ja mahdolliset muut pisteet piirretään ympyröinä.)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
PIIRRA_SAHKOKENTTA: valitsee, piirretäänkö sähkökenttä
 (Jos tosi (True), sähkökenttä piirretään. Jos epätosi (False), ei piirretä.)
PIIRRA_POTENTIAALI: valitsee, piirretäänkö potentiaali
 (Jos tosi (True), potentiaali piirretään. Jos epätosi (False), ei piirretä.
  Näiden valintojen avulla käyttäjä saa halutessaan kuvan sähkökentästä,
  potentiaalista tai molemmista.)
KUVAAJAN_PAATEPISTEET: xy-tasosta valitun suoran päätepisteet
 (Ohjelma laskee sähkökentän ja potentiaalin xy-tasossa ja piirtää niistä
  tasokuvan. Tämän lisäksi ohjelma voi piirtää potentiaalin sekä kentän komponenttien
  kuvaajat suoralla. Käyttäjä voi asettaa tämän suoran miten haluaa.
  Esim. jos KUVAAJAN_PAATEPISTEET = [ (0,0), (1,0) ]
  kuvaajat piirretään positiivisella x-akselilla välillä x = [0,1], y = 0.)
TARKASTELUPISTEET: pisteet, joiden kenttä ja potentiaali listataan erikseen
 (Tämän pitää olla lista xy-koordinaateista.
  Ohjelma piirtää kuvia sähkökentästä ja potentiaalista laajalla alueella.
  Käyttäjä voi myös pyytää ohjelmalta kentän ja potentiaalin lukuarvot
  haluamissaan pisteissä. Pisteet piirretään myös kuviin.)
  