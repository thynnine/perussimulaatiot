
Kaasujen tasapaino

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi kahden kaasusäiliön tilanmuuttujia
ajan kuluessa.
Se koostuu tiedostoista
- kaasu.py
- parametrit.py

Simulaatiossa on kaksiosainen kaasusäiliö.
Säiliö on eristetty ympäristöstään, mutta säiliön osia, A ja B, 
erottaa seinämä, joka voi käyttäjän valinnoista riippuen liikkua tai
johtaa lämpöä. Simulaatio laskee, mitä kaasuille tapahtuu ajan kuluessa.
Kaasujen oletetaan noudattavan ideaalikaasumallia.

Simulaatio ajetaan suorittamalla tiedostossa kaasu.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä niitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

TA: kaasun A alkulämpötila (K)
TB: kaasun B alkulämpötila (K)
VA: kaasun A alkutilavuus (m3)
VB: kaasun B alkutilavuus (m3)
nA: kaasun A ainemäärä (mol)
nB: kaasun B ainemäärä (mol)
cA: kaasun A ominaislämpökapasiteetti (J/molK)
cB: kaasun B ominaislämpökapasiteetti (J/molK)
VAKIO_V: valitsee, pidetäänkö kaasujen tilavuudet vakiona
 (Jos tosi (True), kaasujen välinen seinä pysyy paikoillaan ja tilavuus on vakio.
  Jos epätosi (False), seinä liikkuu, jos kaasujen paine työntää sitä.
  Huom. seinän liike on dissipatiivista, ja sen liike-energia muuttuu
  lopulta lämpöenergiaksi. Jos näin ei olisi, seinä päätyisi värähtelemään 
  edestakaisin, mikä ei olisi realistista.)
VAKIO_T: valitsee, pidetäänkö kaasujen lämpötila vakiona
 (Jos tosi (True), ulkoinen termostaatti pitää kaasut aina vakiolämpötilassa.
  Jos epätosi (False), kaasujen lämpötila muuttuu, jos niiden energia
  muuttuu työn tai lämmön kautta.)
T_LOPPU: simulaation lopetusaika (s)
ANIMAATIO_DT: kuinka usein suureet tallennetaan analyysiä varten (s)
KUVA_DT: kuinka usein systeemistä piirretään erillinen kuva (s)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
TALLENNA_ANIMAATIO: valitsee, näytetäänkö animaatio ruudulla vai tallennetaanko se
 (Animaation tallentaminen vaatii, että tietokoneelle on asennettuna jokin
  ohjelma, joka osaa kirjoittaa videotiedoston. Jos sellaista ei ole, anna
  ohjelman vain näyttää animaatio tietokoneen ruudulta.)
  
