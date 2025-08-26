#
# Parametrejä simulaatioihin.
# Nämä ovat kaikki koordinaatiston A arvoja.
#
# Suureiden yksiköt ovat 
# [x] = m
# [t] = μs
# [v] = m/μs (eli c = 300)
# [m] = ng (huom. nanogrammaa eli 10**-12 kg)
# [F] = N
#

XALKU_A =   0.0             # hiukkasen 1 alkupaikka koordinaatistossa A (hetkellä t = TMIN_A)
VALKU_A =   0.0             # hiukkasen 1 alkunopeus A:ssa (pitää olla |v| < 300)
TMIN_A  = -20.0             # simulaation alkuhetki A:ssa
TMAX_A  =  20.0             # simulaation loppuhetki A:ssa
T_HAJOAMINEN_A      = 100.0 # hiukkasen 1 hajoamisen hetki A:ssa (jos arvo ei ole välillä [TMIN_A, TMAX_A], hiukkanen ei hajoa)
T_KIIHDYTYS_ALKU_A  =  -2.0 # hiukkasen 1 kiihdytyksen alkuhetki A:ssa
T_KIIHDYTYS_LOPPU_A =   2.0 # hiukkasen 1 kiihdytyksen loppuhetki A:ssa
F_A =  0.0                  # hiukkasiin kohdistuva vakiovoima A:ssa
M1  =  1.0                  # hiukkasen 1 massa
M2  =  0.5                  # hajoamisessa syntyneen hiukkasen 2 massa
M3  =  0.5                  # hajoamisessa syntyneen hiukkasen 3 massa (pitää olla M1 >= M2 + M3)




#
# Parametrejä kuvien piirtämiseen
#

VXBA = 100.0  # koordinaatiston B nopeus koordinaatiston A suhteen (pitää olla |v| < 300)

# Tapahtumien koordinaatteja A:ssa, (xA, tA).
#
# Nämä eivät vaikuta simulaatioon mitenkään.
# Ohjelma minkowski.py vain laskee näiden
# Lorentz-muunnokset ja piirtää ne pisteinä
# Minkowski-diagrammeihin.
#
# Voit listata tähän mielestäsi tärkeitä
# tapahtumia, jolloin näet miten ne sijoittuvat
# eri koordinaatistoihin.
TAPAHTUMAT_A = [
    (  -300,  0.00),
    (     0,  1.00)
]

# tapahtumien koordinaatteja B:ssa, (xB, tB)
TAPAHTUMAT_B = [
    (  300,  0.00),
    (    0, -1.00)
]

# Jos et halua ylimääräisiä pisteitä kuviin, määrittele
# tapahtumat tyhjinä listoina:
# TAPAHTUMAT_A = [] 
# TAPAHTUMAT_B = []

# tiedostojen nimet, joista luetaan liikeratojen koordinaatit
#
# Tiedostojen pitää olla samassa hakemistossa kuin missä ohjelma ajetaan.
# Tiedostoja voi olla kuinka monta tahansa.
# Jos tiedostoja ei löydy, ohjelma ilmoittaa siitä ja jatkaa toimintaansa
# normaalisti.
RADAT = [ "liikerata_1.csv", "liikerata_2.csv", "liikerata_3.csv" ]


KUVAFORMAATTI = "pdf"      # kuvaajien formaatti (png, pdf)
TALLENNA_ANIMAATIO = False  # tallennetaanko animaatio tiedostoon?