
Elektronin kvanttimekaaninen aaltofunktio dynaamisena pulssina

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi elektronin liikettä kvanttimekaanisena aaltopulssina.
Se koostuu tiedostoista
- elektroniaalto.py
- parametrit.py

Simulaation alussa elektroniaalto muodostaa Gaussin pulssin,
jolle epätarkkuustulo saavuttaa minimin. Käyttäjä saa valita 
pulssin keskipisteen paikan, leveyden sekä nopeuden.
Elektronin potentiaalienergian voi myös valita
viidestä erilaisesta vaihtoehdosta. 

Simulaatio laskee aallon liikkeen tässä ympäristössä ja
piirtää animaation sekä elektronin aaltofunktiosta että sen
liikemääräesityksestä eli aaltofunktion Fourier-muunnoksesta. 
Ohjelma piirtää myös kuvaajat elektronin paikan ja liikemäärän 
jakaumien aikakehityksestä.

Simulaatio ajetaan suorittamalla tiedostossa elektroniaalto.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

X_ALKU: elektronin paikan odotusarvo alussa (nm)
DX_ALKU: elektronin paikan epätarkkuus alussa (nm)
P_ALKU: elektronin liikemäärän odotusarvo alussa (10^-25 kgm/s)
POTENTIAALIN_TYYPPI: potentiaalin valinta kokonaislukuna
 (Simulaatioon on määritelty valmiiksi viisi erilaista potentiaaliprofiilia,
  ja tämä arvo, jonka pitäisi olla jokin luvuista 1...5, valitsee niistä yhden.
  1: U = 0, vapaa hiukkanen
  2: U = kx, lineaarisesti muuttuva potentiaalienergia, "mäki"
  3: U = kx^2, paraabelin muotoinen potentiaalienergia, harmoninen värähtelijä
  4: U = 0, x < 0, U = U0, x > 0, askelfunktio
  5: U > 0, |x| < x0, U = 0 muuten, kapea valli.)
POTENTIAALIN_KORKEUS: potentiaalienergia systeemistä piirretyn kuvan reunalla (eV)
 (Simulaatiosta piirretään 300 nm levyinen alue, jonka reunat ovat pisteissä
  x = -150 nm ja x = 150 nm. Potentiaalien 2, 3 ja 4 tapauksessa potentiaalienergia
  saa tämän muuttujan arvon pisteessä x = 150 nm eli kuvan oikeassa laidassa. 
  Potentiaalin 5 tapuksesssa tämä on vallin korkeus eli potentiaalin arvo origossa.)
SIMULAATION_NOPEUS: aika-askeleen pituus (fs)
 (Tämä vaikuttaa sekä piirrettyjen animaatioiden että itse laskennan nopeuteen.
  Mitä pienempi tämä arvo on, sitä tarkemmin lasketaan ja sitä hitaammin animaatiot
  kulkevat.)
T_LOPPU: simulaation lopetusaika (fs)
KUVA_DT: kuinka usein systeemistä piirretään erillinen kuva (s)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
TALLENNA_ANIMAATIO: valitsee, näytetäänkö animaatio ruudulla vai tallennetaanko se
 (Animaation tallentaminen vaatii, että tietokoneelle on asennettuna jokin
  ohjelma, joka osaa kirjoittaa videotiedoston. Jos sellaista ei ole, anna
  ohjelman vain näyttää animaatio tietokoneen ruudulta.)
  
