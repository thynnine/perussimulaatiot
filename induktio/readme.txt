
Sähkömagneettinen induktio

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi käämin läpi putoavan magneetin aiheuttaman induktiojännitteen.
Se koostuu tiedostoista
- induktio.py
- parametrit.py

Simulaatiossa on pystyasennossa pidetty sauvamagneetti, joka päästetään
vapaaseen pudotukseen suoran, poikkileikkaukseltaan ympyrän muotoisen
käämin läpi. Simulaatio laskee käämin läpi kulkevan magneettivuon
sekä käämiiin indusoituneen jännitteen ajan ja magneetin paikan funktiona.

Simulaatio ajetaan suorittamalla tiedostossa induktio.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

MAGNEETIN_PITUUS: sauvamagneetin pituus (m)
KAAMIN_PITUUS: käämin pituus (m)
KAAMIN_SADE: käämin poikkileikkauken säde (m)
KAAMIN_KIERROKSET: montako kertaa johdin on kierretty silmukalle käämin ympäri
DIPOLIMOMENTTI: sauvamagneetin magneettisen dipolimomentin suuruus (Am)
ALKUKORKEUS: magneetin ja käämin keskipisteiden välinen etäisyys alussa (m)
T_LOPPU: simulaation lopetusaika (s)
KUVA_DT: kuinka usein suureet tallennetaan analyysiä varten (s)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
