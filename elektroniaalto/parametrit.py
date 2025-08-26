#
# elektronin tila alussa
#
# Ohjelma luo elektronin tilaan, jossa sen aaltofunktio muodostaa
# Gaussin pulssin eli paikka on normaalijakautunut.
#
X_ALKU  = -50.0  # elektronin paikan keskiarvo aluksi, nm
DX_ALKU =   5.0  # elektronin paikan keskihajonta aluksi, nm
P_ALKU  =   1.0  # elektronin liikemäärän keskiarvo aluksi, 10^-25 kgm/s

#
# systeemin potentiaalienergia
#
# valitse potentiaalienergian tyyppi ja muoto
# 1: vapaa hiukkanen
# 2: mäki
# 3: paraabeli
# 4: askelfunktio
# 5: valli
POTENTIAALIN_TYYPPI   = 1
POTENTIAALIN_KORKEUS  = 0.1 # potentiaalin korkeus kuvan reunalla tai askelman kohdalla, eV

#
# simulaatio ja visualisointi
#
SIMULAATION_NOPEUS = 10    # simulaation aika-askelten pituus, fs (pieni: hidas ja tarkka simulaatio)
T_LOPPU            = 400   # simulaation kesto, fs
KUVA_DT            = 200   # kuvien tallentamisen väli, fs
TALLENNA_ANIMAATIO = False  # tallennetaanko animaatio tiedostoon?
KUVAFORMAATTI      = "pdf" # kuvien formaatti (pdf tai png)