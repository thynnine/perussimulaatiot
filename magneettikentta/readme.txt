
Magneettikenttä

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi virtasilmukan luoman magneettikentän.
Se koostuu tiedostoista
- magneettikentta.py
- parametrit.py

Simulaatio laskee siihen määritellyn monikulmion muotoisen virtasilmukan
magneettikentän 3D-avaruudessa.

Simulaatio ajetaan suorittamalla tiedostossa magneettikentta.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

VIRTA: silmukassa kulkeva sähkovirta (A)
KULMAT: silmukan kulmapisteiden koordinaatit (m)
 (Silmukka määritellään antamalla sen kulmapisteiden xyz-koordinaatit
  virran kulkusuunnan mukaisessa järjestyksessä listana.
  Esim. KULMAT = [ (1,0,0), (0,1,0), (-1,0,0), (0,-1,0)]
  määrittelee origokeskisen neliön muotoisen silmukan.
  Huom. pisteiden ei tarvitse olla tasossa vaan silmukalla
  voi olla millainen kolmiulotteinen rakenne tahansa.)
ALKUPISTE: avaruuden suoran alkupisteen xyz-koordinaatit (m)
LOPPUPISTE: avaruuden suoran loppupisteen xyz-koordinaatit (m)
 (Ohjelma piirtää silmukan luoman magneettikentän komponenttien
  kuvaajat suoraa pitkin. Käyttäjää saa asettaa suoran 3D-avaruuteen
  miten haluaa määrittelemällä sen päätepisteet.)
B_SKAALA: kenttänuolten piirtämisen mittakaava
 (Magneettikentän yksikkö on T, mutta systeemin kuva piirretään koordinaatistoon,
  jossa pituudet ovat metrejä. Kenttää kuvaavat nuolet voidaan siis
  piirtää millaiseen mittakaavaan tahansa. Tämä arvo säätää nuolten kokoa.
  Jos arvo on esim. 0.001, 1 T sähkökenttä piirretään 0.001 m pituisena nuolena.)
X_MAKSIMI: systeemistä piirrettävän kuvan koko (m)
 (Systeemistä piirretään kuva xz-tasossa niin, että origo on kuvan keskellä.
  Tämä arvo määrää, kuinka kauas origosta kuva ulottuu.)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
OSIEN_KOKO: numeerisen integroinnin jakovälin pituus (m)
 (Ohjelma määrittää magneettikentän laskemalla numeerisen integraalin
  Biot'n ja Savartin laista pitkin virtasilmukkaa. Silmukka koostuu
  useista lyhyistä suorista johtimista, ja ohjelma jakaa laskua varten
  nämä johtimet pieniin osiin. Tämä arvo on näiden osien pituus.
  Mitä pienempi arvo on, sitä tarkempi lasku on, mutta sitä pidempään
  se myös kestää. Arvon pitäisi olla pienempi kuin lähimmän tarkastelupisteen
  etäisyys mistään silmukkaan kuuluvasta johtimesta.)