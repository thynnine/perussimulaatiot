#
# Kiihtyvyys ajan funktiona.
#
# Kirjoita ratkaisu tähän seuraavassa formaatissa:
#
# A = [
# [ ( ax0, ay0 ), t0 ],
# [ ( ax1, ay1 ), t1 ],
# ...
# [ ( axn, ayn ), tn ]
# ]
#
# Tässä (ax0, ay0) ovat kiihtyvyyden x- ja y-skalaarikomponentit aikavälillä 0 < t < t0,
# (ax1, ay1) ovat komponentit välillä t0 < t < t1, jne.
# Viimeisen annetun ajan hetken jälkeen, t > tn, kiihtyvyys on nolla.
#
# Huom: Kaikki pilkut ja sulkeet ovat välttämättömiä!
# Jos kirjoitat ne väärin, ohjelma antaa virheilmoituksen.
# Lienee parasta ottaa alla annettu malli, kopioida sen rivejä
# tarpeeksi monta kertaa ja muuttaa niitä halutulla tavalla.
#
# Vinkki: voit tallentaa usein käytettyjä kiihtyvyysvektoreita muuttujiin.
# Esim. tämä toimii:
#
# N = ( 0.0,  1.0) # pohjoinen
# E = ( 1.0,  0.0) # itä
# S = ( 0.0, -1.0) # etelä
# W = (-1.0,  0.0) # länsi
#
# A = [
# [ N, 1.0 ],
# [ W, 3.0 ],
# [ S, 6.0 ]
# ]
#


A = [
[ (   0.0,   0.0 ),   0.5 ],  # kiihtyvyys on (0.0, 0.0) hetkeen 0.5 s asti
[ (  -5.0,  -8.0 ),   6.5 ],  # sen jälkeen kiihtyvyys on (-5.0, -8.0) hetkeen 6.5 s asti
[ (   0.0,  10.0 ),  20.0 ]   # sitten kiihtyvyys on (0.0, 10.0) hetkeen 20.0 s asti
] 


#
# Näitäkin parametrejä voi tarvittaessa muuttaa.
# Etenkin simulaation lopetusaika on liian pieni.
# Kannattaa kuitenkin pidentää simulaatiota vähitellen,
# jotta ei joudu turhaan odottelemaan pitkän simulaation 
# loppua.
#

T_LOPPU       = 20.0       # simulaation loppuaika - jos maaliin ei päästä, lopetetaan tässä kohtaa
ANIMAATIO_DT  = 0.5        # tulosten tallennusväli (pieni = tarkat kuvaajat ja animaatio)
KUVA_DT       = 10.0       # kuvien piirtoväli (pieni = paljon kuvia lyhyin aikavälein)
KUVAFORMAATTI = "pdf"      # tallennettavien kuvien formaatti (png tai pdf)
TALLENNA_ANIMAATIO = False # tallennetaanko animaatio tiedostoon?
