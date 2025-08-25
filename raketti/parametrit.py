
M_ISO   = 1.00  # suuren kappaleen massa, kg
M_PIENI = 0.50  # pienien kappaleiden massa, kg
N_PIENET =   1  # pienien kappaleiden lukumäärä

YKSIULOTTEINEN = True # onko simulaatio yksiulotteinen? (jos ei, se on kaksiulotteinen)
VX_ALKU =  0.0  # alkunopeuden x-komponentti, m/s (vain 2D)
VY_ALKU = 30.0  # alkunopeuden y-komponentti, m/s
G = 9.8         # putoamiskiihtyvyys, m/s2 (0 = ei painovoimaa, 9.8 = painovoima alaspäin)

T_LOPPU =  4.0  # simulaation kesto, s

T_RAJAHDYS  = -2.0  # räjähdyksen tapahtumishetki, s (negatiivinen luku poistaa räjähdyksen)
DV_RAJAHDYS = 10.0  # räjähdyksen voimakkuus: pienten kappaleiden nopeuden muutosten keskihajonta, m/s

T_KIIHDYTYS   = -50.0  # suihkutuksen kesto, s (negatiivinen luku poistaa suihkutuksen)
DVX_KIIHDYTYS =   0.0  # suihkutettavien hiukkasten nopeuden muutos x-suunnassa, m/s (positiivinen = vasemmalle, vain 2D)
DVY_KIIHDYTYS = 100.0  # suihkutettavien hiukkasten nopeuden muutos y-suunnassa, m/s (positiivinen = alaspäin)

#
# Piirtämistä ohjaavat parametrit - näitäkin voi muuttaa.
#

ANIMAATIO_DT = 0.05        # tulosten tallennusväli, s (pieni = tarkat kuvaajat ja animaatio)
KUVA_DT     =     1        # kuvien piirtoväli, s (pieni = paljon kuvia lyhyin aikavälein)
KUVAFORMAATTI = "pdf"      # tallennettavien kuvien formaatti (png tai pdf)
TALLENNA_ANIMAATIO = False # tallennetaanko animaatio tiedostoon? (jos ei, näytetään erillisessä ikkunassa)

KAIKKI_NOPEUDET = False  # piirretäänkö myös pienten kappaleiden nopeusvektorit?
KOKONAISSUUREET = False  # piirretäänkö massakeskipiste sekä kokonaisliikemäärä ja -energia?

X_CENTER =  0.0  # koordinaatiston keskipiste x-akselilla
Y_MIN    = -5.0  # koordinaatiston y-akselin minimi
Y_MAX    = 45.0  # koordinaatiston y-akslein maksimi
V_SKAALA =  0.5  # nopeusvektoreiden skaalaus
    