
Raketti

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi osiin hajoavan kappaleen liikettä.
Se koostuu tiedostoista
- raketti.py
- parametrit.py

Simulaation alussa systeemissä on yksi kappale, joka koostuu
yhdestä suuresta ja mahdollisesti useasta pienestä osasta.
Aluksi nämä osat ovat kaikki kiinni toisissaan muodostaen yhden
kappaleen. Kappaleen voi kuitenkin määrätä räjähtämään, jolloin
kaikki pienet osat lentävät siitä pois satunnaisilla nopeuksilla,
tai kiihdyttämään raketin tavoin, jolloin pienet osat ammutaan
kappaleesta pois yksi kerrallaan määrätyllä nopeudella.

Simulaatio ajetaan suorittamalla tiedostossa raketti.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

M_ISO: suuren osan massa (kg)
M_PIENI: kunkin pienen osan massa (kg)
N_PIENET: pienten osien lukumäärä
 (Esim. jos N_PIENET = 2, kappale koostuu kolmesta osasta, 1 suuresta ja 2 pienestä.)
YKSIULOTTEINEN: Jos totta (True), simulaatio on yksiulotteinen y-suunnassa.
 Jos epätotta (False), simulaatio tapahtuu xy-tasossa.
VX_ALKU: kappaleen alkunopeuden x-komponentti (m/s)
 (Jos simulaation on yksiulotteinen, tämä ei vaikuta mihinkään.)
VY_ALKU: kappaleen alkunopeuden y-komponentti (m/s)
G: putoamiskiihtyvyys (m/s)
 (Painovoima vaikuttaa alaspäin eli negatiiviseen y-suuntaan, jos G > 0.
  Jos G = 0, simulaatiossa ei ole painovoimaa.)
T_LOPPU: simulaation lopetusaika (s)
T_RAJAHDYS: kappaleen räjähdyksen hetki (s)
 (Jos kappaleessa on suuria ja pieniä osia tällä hetkellä, kaikille pienille
  osille annetaan satunnaisesti valittu nopeus. Jos tämä on negatiivinen
  tai suurempi kuin T_LOPPU, räjähdystä ei koskaan tapahdu.)
DV_RAJAHDYS: räjähdyksessä arvottujen nopeuksien keskihajonta (m/s)
 (Jos kappale räjäytetään, pienille osille arvotaan nopeudet normaalijakaumasta,
  jonka keskihajonta on tämän parametrin arvo. Suuri arvo tarkoittaa siis
  voimakasta räjähdystä.)
T_KIIHDYTYS: raketin kiihdytyksen kesto (s)
 (Raketti toimii suihkuttamalla kaasuja taaksepäin, jolloin kaasut työntävät
  rakettia eteenpäin. Tätä voi simuloida antamalla kappaleen ampua pieniä osiaan
  alaspäin yksi kerrassaan. Tämä arvo määrää, kuinka kauan näin tehdään.
  Kiihdytys loppuu myös siinä tapauksessa, että pienet osat loppuvat, mikä vastaa
  polttoaineen loppumista. Jos tämä arvo on negatiivinen, kiihdytystä ei tehdä.)
DVX_KIIHDYTYS: suihkutettavien hiukkasten nopeuden muutos x-suunnassa (m/s)
 (Kiihdytyksessä kappaleesta pois ammuttujen pienten osien nopeuden muutos
  x-suunnassa. Positiivinen arvo tarkoittaa osien saavan nopeutta
  vasemmalle, mikä työntää rakettia oikealle. 
  Tämä ei tee mitään, jos simulaatio on yksiulotteinen.)
DVXY_KIIHDYTYS: suihkutettavien hiukkasten nopeuden muutos y-suunnassa (m/s)
 (Kiihdytyksessä kappaleesta pois ammuttujen pienten osien nopeuden muutos
  y-suunnassa. Positiivinen arvo tarkoittaa kappaleiden saavan nopeutta
  alaspäin, mikä työntää rakettia ylöspäin.)
ANIMAATIO_DT: kuinka usein suureet tallennetaan analyysiä varten (s)
KUVA_DT: kuinka usein systeemistä piirretään erillinen kuva (s)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
TALLENNA_ANIMAATIO: valitsee, näytetäänkö animaatio ruudulla vai tallennetaanko se
 (Animaation tallentaminen vaatii, että tietokoneelle on asennettuna jokin
  ohjelma, joka osaa kirjoittaa videotiedoston. Jos sellaista ei ole, anna
  ohjelman vain näyttää animaatio tietokoneen ruudulta.)
KAIKKI_NOPEUDET: valitsee, piirretäänkö kuvaajiin kaikkien pienten osien nopeudet
KOKONAISSUUREET: valitsee, piirretäänkö kuvaajiin kokonaisliikemäärä ja -energia
X_CENTER: piirretyn kuvan keskipisteen x-koordinaatti (m)
Y_MIN: piirretyn kuvan y-akselin alin piste (m)
Y_MAX: piirertyn kuvan y-akselin korkein piste (m)
V_SKAALA: nopeusvektorien skaalaustekijä
 (Koordinaatiston yksiköt ovat metrejä, mutta nopeuden yksikkö on m/s.
  Niinpä nopeusvektorin voi piirtää kuvaan kuinka pitkänä haluaa.
  Jos tämä arvo on esim. 0.5, 1 m/s nopeus piirretään kuvaan 0.5 m -pituisena.)
  
Kokeita
  
Jos haluat simuloida kitkattomalla pinnalla liikkuvan kappaleen räjähtävän 
erotuksen, aseta 
N_PIENET = 1
G = 0 
0 < T_RAJAHDYS < T_LOPPU
T_KIIHDYTYS < 0

Jos haluat simuloida saman heitetylle kappaleelle, aseta 
G = 9.8
VY_ALKU > 0

Jos haluat simuloida raketin ampumisen taivaalle, aseta
G = 9.8
T_KIIHDYTYS > 0
DVY_KIIHDYTYS > 0
