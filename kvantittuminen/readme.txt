
Elektronin kvanttimekaaninen aaltofunktio ominaistiloina

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi elektronin liikettä kvanttimekaanisissa ominaistiloissa
sekä ominaistilojen superpositioissa.
Se koostuu tiedostoista
- kvantittuminen.py
- parametrit.py

Simulaatio ratkaisee elektronin ominaistilat ja -energiat systeemissä
ja piirtää kuvan sekä animaation pyydetyn tilan aaltofunktiosta.
Pyydetty tila voi olla ominaistila |n) tai ominaistilojen superpositio
c(1)|1) + c(2)|2) + c(3)|3) + ...

Simulaatio ajetaan suorittamalla tiedostossa kvantittuminen.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:


N: tilojen kvanttiluvut n
 (Tämän pitää olla lista kokonaisluvuista. Jos simuloidaan ominaistila,
  listassa on vain yksi luku, tilan kvanttiluku. Jos simuloidaan ominaistilojen
  superpositio, listassa on oltava kaikkien näiden tilojen kvanttiluvut.)
KERTOIMET: tilojen superpositiokertoimet c(n)
 (Tämän pitää olla lista reaaliluvuista. Jos simuloidaan ominaistila,
  tämä voi olla mikä tahansa yksi luku. Jos simuloidaan ominaistilojen 
  superpositio, listassa on oltava kaikkien tilojen kertoimet.
  Ohjelma normittaa kertoimet automaattisesti, joten vain niiden
  suhteellisilla suuruuksilla on merkitystä.)
KUOPAN_LEVEYS: systeemin leveys (nm)
 (Simuloitavan systeemin ulkopuolella potentiaalienergia on ääretön
  eli aaltofunktio on nolla.)
KUVAN_KORKEUS: systeemistä piirrettyjen kuvien energian maksimi (eV)
 (Kuviin piirretään potentiaalienergian kuvaaja. Lisäksi aaltofunktio
  piirretään niin, että sen nollakohta asetetaan energian tai energian
  odotusarvon korkeudelle energia-akselilla.)
POTENTIAALIN_TYYPPI: potentiaalin valinta kokonaislukuna
 (Simulaatioon on määritelty valmiiksi viisi erilaista potentiaaliprofiilia,
  ja tämä arvo, jonka pitäisi olla jokin luvuista 1...5, valitsee niistä yhden.
  1: askelfunktio
  2: pienempi kuoppa suuren kuopan keskellä
  3: kapea valli kuopan keskellä
  4: paraabeli, harmoninen potentiaali
  5: lineaarinen potentiaali, mäki
POTENTIAALIN_KORKEUS: potentiaalikuopan korkeus (eV)
 (Potentiaalien 1, 2 ja 3 kohdalla tämä on potentiaalienergiassa nähtävän askelman korkeus.
  Potentiaalien 4 ja 5 kohdalla tämä on potentiaalienergian arvo kuopan oikealla reunalla.)
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
  
