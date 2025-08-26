
Suhteellisuusteoria

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi liikettä ja visualisoi sen
kahdessa eri inertiaalikoordinaatistossa huomioiden
suppean suhteellisuusteorian.
Se koostuu tiedostoista
- relativistinen.py
- minkowski.py
- parametrit.py

Simulaatiossa on aluksi yksi hiukkanen.
Käyttäjä saa asettaa hiukkasen paikan ja nopeuden.
Hiukkaseen voi myls kohdistaa vakiovoiman ja hiukkasen
voi antaa hajota kahdeksi pienemmäksi hiukkaseksi.

Simulaatio ajetaan suorittamalla tiedostossa relativistinen.py oleva ohjelma.
Ohjelma kirjoittaa simulaation tulokset tiedostoihin nimeltä
liikerata_1.csv jne. Tuloksista piirretään kuvaajia ja animaatioita
suorittamalla samassa hakemistossa tiedoston minkowski.py ohjelma.
Kumpikin ohjelma toimii sellaisenaan eikä niitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

XALKU_A: hiukkasen 1 paikka aluksi koordinaatistossa A (m)
 (Huom. alkuhetki ei ole t = 0 s vaan muuttujan TMIN_A arvo.)
VALKU_A: hiukkasen 1 alkunopeus koordinaatistossa A (m/μs)
 (Huom. ei saa ylittää valonnopeutta 300 m/μs.)
TMIN_A: ajanoton alkuhetki koordinaatistossa A (μs)
TMAX_A: simulaation loppuhetki koordinaatistossa A (μs)
T_HAJOAMINEN_A: hiukkasen 1 hajoamisen hetki koordinaatistossa A (μs)
 (Jos tämä ei ole arvojen TMIN_A ja TMAX_A välissä, hiukkanen ei hajoa.)
T_KIIHDYTYS_ALKU_A: hetki, jolloin hiukkaseen 1 kohdistetaan voima koordinaatistossa A (μs)
T_KIIHDYTYS_LOPPU_A: hetki, jolloin voiman vaikutus loppuu koordinaatistossa A (μs)
F_A: hiukkaseen 1 kohdistuva voima koordinaatistossa A (N)
M1: hiukkasen 1 massa (ng)
M2: hajoamisessa syntyneen hiukkasen 2 massa (ng)
M3: hajoamisessa syntyneen hiukkasen 3 massa (ng)
 (Syntyvien hiukkasten massat eivät voi olla yhdessä suuremmat kuin hiukkasen 1 massa oli.)
VXBA: koordinaatiston B nopeus A:n suhteen (m/μs)
TAPAHTUMAT_A: lista tapahtumien koordinaateista A:ssa (m, μs)
 (Tapahtumat eivät vaikuta hiukkassimulaatioon mitenkään.
  Nämä koordinaatit piirretään Minkowski-kuvaajiin.
  Niitä voi siis käyttää visualisoinnissa apuvälineenä.)
TAPAHTUMAT_B: lista tapahtumien koordinaateista B:ssä (m, μs)
RADAT: lista tiedostoista, joihin simulaatiodata on tallennettu
 (Listassa voi olla kuinka monta tiedostoa tahansa.
  Tiedostojen pitää olla samassa hakemistossa kuin missä
  ohjelma minkowski.py ajetaan. Ohjelma lukee nämä tiedostot
  ja piirtää niiden kuvaamat ratakäyrät samaan kuvaajaan.)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
TALLENNA_ANIMAATIO: valitsee, näytetäänkö animaatio ruudulla vai tallennetaanko se
 (Animaation tallentaminen vaatii, että tietokoneelle on asennettuna jokin
  ohjelma, joka osaa kirjoittaa videotiedoston. Jos sellaista ei ole, anna
  ohjelman vain näyttää animaatio tietokoneen ruudulta.)
  
