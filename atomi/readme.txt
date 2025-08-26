
Elektronin kvanttimekaaninen aaltofunktio vetyatomissa

Teemu Hynninen 2025
CC-SA-BY-4.0

Tämä ohjelma simuloi elektronin aaltofunktion vetyatomissa.
Se koostuu tiedostoista
- atomi.py
- parametrit.py

Simulaatio ratkaisee elektronin ominaistilan |n,l,m,s) aaltofunktion
vetyatomissa ja piirtää kuvia sekä aaltofunktiosta että paikan 
todennäköisyysjakaumasta (eli siis orbitaalista) eri tavoin.

Simulaatio ajetaan suorittamalla tiedostossa atomi.py oleva ohjelma.
Ohjelma toimii sellaisenaan eikä sitä tarvitse muuttaa mitenkään.
Simulaation kulkua ohjataan muuttamalla tiedostossa parametrit.py
määriteltyjen muuttujien arvoja. Nämä ovat:

N: tilan pääkvanttiluku n
L: tilan sivukvanttiluku l
M: tilan sivukvanttiluku m
R_MAKSIMI: maksimietäisyys ytimestä kuvia piirrettäessä (Bohrin sädettä a)
VOIMAKKUUS: orbitaalikuvien tummuus tai kirkkaus
 (Pieni arvo tarkoittaa sitä, että vain suuren todennäköisyyden alueet piirretään näkyviin.)
KUVAFORMAATTI: piirrettävien kuvien tiedostomuoto, png tai pdf
