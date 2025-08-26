#
# Kokeen fysiikkaan liittyviä parametrejä.
# Kaikki parametrit on annettu SI-perusyksiköissä.
# Esim. pituudet ovat siis metrejä.
#
MAGNEETIN_PITUUS  = 0.0    # magneetin pituus ( 0 = pistemäisen dipolin malli )
KAAMIN_PITUUS     = 0.0    # käämin korkeus ( 0 = tasomaisen silmukan malli )
KAAMIN_SADE       = 0.010  # käämin säde olettaen käämin poikkileikkaus ympyräksi
KAAMIN_KIERROKSET = 1      # montako kertaa johdin on kierretty käämin ympäri
DIPOLIMOMENTTI    = 0.010  # magneetin dipolimomentti
ALKUKORKEUS       = 0.10   # magneetin etäisyys käämistä aluksi

#
# Kuvaajien piirtoon liittyviä parametrejä
#
T_LOPPU       = 0.4    # simulaation ajallinen kesto
DT_KUVA       = 0.001  # mittauspisteiden välinen aikaero
KUVAFORMAATTI = "pdf"  # kuvatiedoston tyyppi (png tai pdf)
