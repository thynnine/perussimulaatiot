#
# värähtelijän ominaisuudet
#
# Simuloidaan jouseen kiinnitetyn kappaleen värähtelyä.
# Kaikki arvot on annettu SI-perusyksiköissä.
#
M = 1.00  # massa
K = 10.0  # jousivakio
B = 0.01  # ilmanvastuskerroin

#
# alkutilan ominaisuudet
#
DX_ALKU = 0.10 # poikkeama tasapainosta alussa
VX_ALKU = 0.00 # nopeus alussa

#
# ulkoisen voiman ominaisuudet
#
# Jos VOIMA = 0, värähtelijä liikkuu itsekseen.
# Jos VOIMA > 0, värähtelijää työntää liikkeelle
# harmonisesti värähtelevä voima.
#
VOIMA   = 0.0 # voiman maksimi (amplitudi)
TAAJUUS = 1.0 # voiman taajuus


#
# kuvien piirtoon liittyviä parametrejä
#
T_LOPPU         = 10.0       # simulaation loppuaika
ANIMAATIO_LOPPU =  2.0       # animaation loppuaika
ANIMAATIO_DT    = 0.02       # tulosten tallennusväli (pieni = tarkat kuvaajat ja animaatio)
ANIMAATIONOPEUS =  1.0       # montako sekuntia simulaatioaikaa animaatiossa näytetään yhden todellisen sekunnin aikana. 1: reaaliaika, <1: hidastettu, >1: nopeutettu
KUVAFORMAATTI   = "pdf"      # tallennettavien kuvien formaatti (png tai pdf)
TALLENNA_ANIMAATIO = True    # tallennetaanko animaatio tiedostoon?
