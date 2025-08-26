
Äänen spektri

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämän paketin ohjelmat analysoivat oikeaa ääntä ja
syntetisoivat keinotekoista ääntä.
Se koostuu tiedostoista
- analysaattori.py
- analysaattori_parametrit.py
- syntetisaattori.py
- syntetisaattori_parametrit.py

Ääni on jaksollisista värähtelyistä koostuva aalto, ja tällaisten
värähtelyjen kuvaamisessa Fourier-muunnos on tärkeä työkalu.
Mikä tahansa jaksollinen funktio voidaan esittää harmonisten funktioiden
A cos( 2 pi f t + phi ) summana ja mikä tahansa funktio niiden integraalina.

Jos meillä on todellinen funktio kuten paineen muutokset ääniaallossa,
Fourier-analyysi tarkoittaa tässä summassa tai integraalissa esiintyvien
harmonisten komponenttien amplitudien A, taajuuksien f ja vaiheiden phi
etsimistä. Jos meillä puolestaan on tiedossa nämä suureet, alkuperäisen funktion
laskeminen niiden avulla on Fourier-synteesiä.
Tässä paketissa onkin kaksi ohjelmaa: analysaattori ja syntetisaattori.

Analysaattorille annetaan wav-formaattiin tallennettu ääni, jolloin ohjelma
laskee äänen Fourier-muunnoksen eli etsii amplitudit A taajuuden f funktiona.
Ohjelma piirtää näiden perusteella äänen spektrin |A(f)| sekä lineaarisella
että logaritmisella (desibeli) asteikolla. Tietoa komponenttien vaiheista phi ei
tallenneta. Ohjelma piirtää myös aallon aiheuttaman värähtelyn kuvaajan 
ajan funktiona.

Syntetisaattorille annetaan lista harmonisten komponenttien amplitudeista,
taajuuksista ja vaiheista, jolloin ohjelma laskee näitä vastaavien
harmonisten värähtelyjen summan ja tallentaa sitä vastaavan äänen wav-tiedostona.
Tämän tiedoston voi syöttää takaisin analysaattorille.

Analyysi ajetaan suorittamalla tiedostossa analysaattori.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Analyysiä ohjataan muuttamalla tiedostossa analysaattori_parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

NIMI: tiedoston nimi
 (Anna nimi ilman päätettä. Jos nimi on esim. "huilu", ääni luetaan
  tiedostosta huilu.wav ja kuvaajat piirretään tiedostoihin nimeltä
  huilu_*.*)
MAKSIMITAAJUUS: kuvaajien taajuusakselin maksimi (Hz)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf

Synteesi ajetaan suorittamalla tiedostossa syntetisaattori.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Synteesiä ohjataan muuttamalla tiedostossa syntetisaattori_parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

TAAJUUDET: lista taajuuksista (Hz)
AMPLITUDIT: lista amplitudeista
VAIHEET: lista vaiheista (rad)
 (Listoissa tulisi olla yhtä monta alkiota. Jos ei ole, TAAJUUDET määrää
  summaan tulevien harmonisten komponenttien lukumäärän.
  Jos muissa listoissa on vähemmän alkioita, puuttuvien arvojen tilalle
  laitetaan nollia. Amplitudin tapauksessahan tämä tarkoittaa sitä, että
  komponentti ei vaikuta summmaan mitenkään.)

