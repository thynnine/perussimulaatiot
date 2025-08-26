
Oskillaattori

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi vaimennetun harmonisen värähtelijän liikettä.
Se koostuu tiedostoista
- oskillaattori.py
- parametrit.py

Simulaatiossa on yksi harmonisessa potentiaalissa liikkuva kappale.
Kappale voisi olla esimerkiksi kiinni jousessa.
Jos kappaleeseen ei kohdistu muita ulkoisia vuorovaikutuksia,
se jää värähtelemään harmonisesti.
Jos kappaleeseen kohdistetaan dissipatiivinen vuorovaikutus kuten
ilmanvastus, se menettää mekanisen energiansa ja pysähtyy.
Tämä on vaimennettua värähtelyä.
Jos kappaleeseen kohdistetaan ulkoinen harmonisesti värähtelevä voima,
kappale päätyy loppujen lopuksi värähtelemään tämän ulkoisen voiman taajuudella.
Tämä on pakotettua värähtelyä.
Värähtelyn amplitudi riippuu siitä, kuinka hyvin ulkoinen voima sopii yhteen
värähtelijän ominaisuuksien kanssa. Jos vastaavuus on hyvä, amplitudi voi
kasvaa suureksi. Tämä on resonanssi.

Simulaatio ajetaan suorittamalla tiedostossa oskillaattori.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

M: kappaleen massa (kg)
K: jousivakio (N/m)
B: ilmanvastuskerroin (kg/s)
DX_ALKU: kappaleen poikkeama tasapainosta aluksi (m)
VX_ALKU: kappaleen nopeus aluksi (m/s)
VOIMA: ulkoisen voiman amplitudi (N)
TAAJUUS: ulkoisen voiman taajuus (Hz)
T_LOPPU: simulaation lopetusaika (s)
ANIMAATIO_LOPPU: animaation lopetusaika (s)
 (Joissakin tapauksissa voi olla syytä simuloida hyvin pitkä aika.
  Pitkän animaation tuottaminen vie kuitenkin paljon laskentaresursseja,
  joten silloin kannattaa lopettaa animaatio ennen simulaation loppua.)
ANIMAATIO_DT: kuinka usein suureet tallennetaan analyysiä varten (s)
ANIMAATIONOPEUS: animaatiossa kuluneen ajan suhde todelliseen aikaan
 (Jos tämä arvo on esimerkiksi 2 ja animaation kesto on 10 s simulaatioaikaa,
  siitä luodaan todellisuudessa 5 sekuntia kestävä nopeutettu animaatio.)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
TALLENNA_ANIMAATIO: valitsee, näytetäänkö animaatio ruudulla vai tallennetaanko se
 (Animaation tallentaminen vaatii, että tietokoneelle on asennettuna jokin
  ohjelma, joka osaa kirjoittaa videotiedoston. Jos sellaista ei ole, anna
  ohjelman vain näyttää animaatio tietokoneen ruudulta.)
  
Kokeita
  
Jos haluat simuloida harmonisen värähtelyn, aseta
B = 0
VOIMA = 0

Jos haluat simuloida vaimennetun värähtelyn, aseta
B > 0
VOIMA = 0

Jos haluat simuloida pakotetyn värähtelyn, aseta
B > 0
VOIMA > 0
Huom. jos B on pieni, värähtelijällä voi kestää hyvin pitkään ennen kuin se
saavuttaa lopullisen tasaisen värähtelyn tilan. Aja siis ensin pitkä simulaatio
ja katso kuvaajasta, onko värähtelyn amplitudi asettunut jo johonkin vakioarvoon.
Jos ei ole, aja vielä pidempi simulaatio.
