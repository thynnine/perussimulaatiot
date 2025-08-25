
Lentokilpailu

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi tasossa liikkuvan kappaleen dynamiikkaa.
Se koostuu tiedostoista
- lento.py
- parametrit.py

Simulaatiossa on kappale, drooni, joka pääsee liikkumaan xy-tasossa.
Droonia ohjataan määräämällä sen kiihtyvyyden
x- ja y-komponentit (paloittain vakioina) ajan funktioina.
Simulaatiossa on myös väistettäviä esteitä sekä kerättäviä maaleja,
ja käyttäjän tehtävänä on ohjelmoida drooni kulkemaan kaikkien
maalien kautta mahdollisimman nopeasti. Jos drooni törmää matkalla
esteeseen, simulaatio loppuu siihen. Simulaatiossa on myös
maksimiaika, jolloin se lopetetaan joka tapauksessa jollei maaleja ole
saatu siihen mennessä kerättyä.

Simulaatio ajetaan suorittamalla tiedostossa lento.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

A: kiihtyvyyden komponentit (m/s)
 (Komponentit pitää antaa listojen listana, jossa kukin alkio määrittelee 
  kiihtyvyyden komponentit tiettyyn aikaan asti.
  Esim. A = [ [ (1,2), 3 ], [ (4,5), 6 ] ]
  tarkoittaisi sitä, että kiihtyvyys on 
  ax = 1 m/s2, ay = 2 m/s2 aikavälillä 0...3 s,
  ax = 4 m/s2, ay = 5 m/s2 aikavälillä 3...6 s
  ja nolla 6 s jälkeen.)

T_LOPPU: simulaation lopetusaika, jos maaliin ei päästä (s)
ANIMAATIO_DT: kuinka usein suureet tallennetaan analyysiä varten (s)
KUVA_DT: kuinka usein systeemistä piirretään erillinen kuva (s)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
TALLENNA_ANIMAATIO: valitsee, näytetäänkö animaatio ruudulla vai tallennetaanko se
 (Animaation tallentaminen vaatii, että tietokoneelle on asennettuna jokin
  ohjelma, joka osaa kirjoittaa videotiedoston. Jos sellaista ei ole, anna
  ohjelman vain näyttää animaatio tietokoneen ruudulta.)
  
