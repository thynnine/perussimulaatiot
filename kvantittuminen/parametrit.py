# Ilmaistaan tila kvanttilukujen n ja kertoimien c avulla.
#
# Jos kyseessä on ominaistila |n), sitä kuvaa yksi kvanttiluku n = 1,2,3,...
# eikä kertoimella ole mitään merkitystä.
# Jos esim. haluat simuloida tilan |3), aseta
# N = [3]
# KERTOIMET = [1]
#
# Jos kyseessä on superpositiotila, se on ominaistilojen summa
# c(1)|1) + c(2)|2) + c(3)|3) + ...
# Tämä määritellään antamalla lista superpositiossa olevien tilojen
# kvanttiluvuista n ja niiden kertoimista c(n).
# Jos esim. haluat simuloida tilan 3/5 |1) + 4/5 |3), aseta
# N = [1, 3]
# KERTOIMET = [3/5, 4/5]
#
# Huom. sekä N että KERTOIMET pitää antaa listana, vaikka tiloja olisi vain yksi.
# Huom. ohjelma normittaa kertoimet automaattisesti. Edellisessä esimerkissä siis
# toimisi aivan yhtä hyvin myös
# KERTOIMET = [3, 4]
#
N = [1]
KERTOIMET = [1]

KUOPAN_LEVEYS = 1.0   # kuopan koko leveys, jonka ulkopuolella U = ääretön, nm
KUVAN_KORKEUS = 30    # kuvaajaan piirretty energian maksimiarvo, eV

SIMULAATION_NOPEUS =   2.0  # animaation nopeus: montako fs nähdään animaatiossa 1 s aikana
T_LOPPU            =  10.0  # simulaation lopetusaika, fs
KUVA_DT            =   5.0  # aika kuvien piirtämisen välillä, fs

KUVAFORMAATTI = "pdf"      # tallennettavan kuvan formaatti ('pdf' tai 'png')
TALLENNA_ANIMAATIO = True  # tallennetaanko animaatio tiedostoon?

# valitse potentiaalienergian tyyppi ja muoto
# 1: askelfunktio,    U(x) = 0 (x < L/2), U0 (x > L/2)
# 2: kuoppa kuopassa, U(x) = U0 (x < L/2-dL tai x > L/2 + dL), 0 (L/2-dL < x < L/2+dL), dL = L / 10
# 3: valli,           U(x) = 0 (x < L/2-dL tai x > L/2 + dL), U0 (L/2-dL < x < L/2+dL), dL = L / 40
# 4: paraabeli,       U(x) = k*(x - L/2)^2
# 5: mäki,            U(x) = k*x
# L on sama kuin KUOPAN_LEVEYS, eli L/2 on kuopan keskipiste
POTENTIAALIN_TYYPPI   =   1 
POTENTIAALIN_KORKEUS  =   5   # U0, eV

