
Valon kulku raoista

Teemu Hynninen 2024
CC-SA-BY-4.0

Tämä ohjelma simuloi valon kulkua kapeissa raoissa.

Pienistä raoista kulkeva valo kokee aaltoilmiöitä kuten diffraktiota ja interferenssiä.
Tämä simulaatio laskee valon intensiteettijakauman tällaisten rakojen takana ja piirtää
sekä intensiteetin kuvaajan että kuvan valon piirtämästä kuviosta.
Käytännössä valon piirtämä kuvio on rakojen Fourier-muunnos, ja sen ohjelma laskee.

Simulaatio ajetaan suorittamalla tiedostossa interferenssi.py oleva ohjelma.
Valon ja rakojen ominaisuuksia säädetään tiedoston alussa olevien parametrien avulla:

kirkkaus: piirrettävien kuvien kirkkaus välillä [0,1]
aallonpituus: valon aallonpituus (nm)
tapaus: kokonaisluku, joka valitsee millaisista raoista valo kulkee
 (1: kaksi kapeaa rakoa
  2: kaksi kapeaa rakoa kaksinkertaisella etäisyydellä
  3: neljä kapeaa rakoa
  4: monta kapeaa rakoa (diffraktiohila)
  5: yksi kapea rako (voimakas diffraktio)
  6: yksi melko leveä rako
  7: yksi rako kaksinkertaisella leveydellä
  8: yksi hyvin leveä rako
  9: kaksi melko leveää rakoa
  10: kaksi melko leveää rakoa kaksinkertaisella etäisyydellä
  11: monta leveää rakoa
  12: joukko rakoja: määrittele rakojen leveydet ja keskipisteet funktiossa valmistele_koeasetelma()
  13: joukko rakoja: määrittele rakojen reunapisteet funktiossa valmistele_koeasetelma().)
w: rakojen leveys tapauksissa 1, 2, 3, 4 (nm)
d: rakojen välinen etäisyys tapauksissa 1, 2, 3, 4 ja raon leveys tapauksessa 7 (nm)
