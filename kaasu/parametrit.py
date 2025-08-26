#
# Alkutilan ominaisuudet
#
TA = 300.0      # A:n lämpötila (K)
TB = 300.0      # B:n lämpötila (K)
VA = 0.050      # A:n tilavuus (m3)
VB = 0.050      # B:n tilavuus (m3)
nA = 2.00       # A:n ainemäärä (mol)
nB = 2.00       # B:n ainemäärä (mol)
cA = 10.000     # A:n ominaislämpökapasiteetti (J/molK)
cB = 10.000     # B:n ominaislämpökapasiteetti (J/molK)

VAKIO_V = False   # pidetäänkö tilavuus vakiona?
VAKIO_T = False   # pidetäänkö lämpötila vakiona?


#
# Näitäkin parametrejä voi tarvittaessa muuttaa:
#

T_LOPPU      = 100.0       # simulaation loppuaika
ANIMAATIO_DT =   0.1       # tulosten tallennusväli (pieni = tarkat kuvaajat ja animaatio)
KUVA_DT      =   5.0       # kuvien piirtoväli (pieni = paljon kuvia lyhyin aikavälein)
KUVAFORMAATTI = "pdf"      # tallennettavien kuvien formaatti (png tai pdf)
TALLENNA_ANIMAATIO = False # tallennetaanko animaatio tiedostoon?
