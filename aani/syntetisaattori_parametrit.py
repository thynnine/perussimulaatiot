#
# Listaa syntetisoitavan äänen harmonisten komponenttien
#
#    A cos(2 pi f t + phi)
#
# o taajuudet f hertseissä
# o amplitudit A mielivaltaisissa yksiköissä 
#  - Amplitudi skaalataan, joten vain suhteellisilla arvoilla on väliä
#  - Jos annat amplitudille negatiivisen arvon, miinusmerkki jätetään huomioimatta.
# o vaihetekijät phi radiaaneissa
#
# Komponentteja voi olla kuinka monta tahansa.
#
# Jos amplitudeja tai vaiheita annetaan enemmän kuin taajuuksia, ylimääräiset jätetään huomioimatta.
# Jos niitä annetaan vähemmän kuin taajuuksia, puuttuvat oletetaan nolliksi.
#
TAAJUUDET  = [ 400, 800, 1200, 1600 ]
AMPLITUDIT = [ 1.0, 0.5, 0.3, 0.2 ]
VAIHEET    = [ 0.0, 0.0, 0.0, 0.0 ]
