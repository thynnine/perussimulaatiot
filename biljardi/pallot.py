#
# Biljardipeli
#
# Simuloidaan biljardipallojen liikettä käyttäen kovien pallojen mallia.
#
# Simulaatio käynnistyy, kun pelaaja lyö valkoista palloa.
# Tähdätä voi joko hiirellä tai nuolinäppäimillä. Lyönti tapahtuu klikkaamalla
# hiirellä tai painamalla välilyöntiä kaksi kertaa.
# Mitä pidempi näiden klikkausten väli on, sitä kovempaa palloa lyödään.
#
# Pallojen ja vallien väliset törmäykset simuloidaan täysin elastisina ja kitkattomina.
# Toisin sanoen kokonaisenergia ja -liikemäärä ovat törmäyksissä vakiot.
# Liikkuviin palloihin vaikuttaa vierimis- ja ilmanvastus, joten ne pysähtyvät
# aikanaan. (Vastusvoimat voi kuitenkin poistaa asettamalla muuttujan KITKA arvoksi
# nollan. Silloin systeemin kokonaisenergia pysyy vakiona ja simulaatio muuttuu
# biljardista ideaalikaasuksi.)
#
# Teemu Hynninen 2024
#

import tkinter as tk
import math
import random
import time

import matplotlib.pyplot as plt
import numpy as np

# vakioita
SKAALA = 5 # grafiikan mittakaava: 1 cm piirretään näin monena pikselinä.
R = 3 * SKAALA # pallon säde (cm)
LEVEYS = 224 * SKAALA # pelialueen leveys (cm)
KORKEUS = 112 * SKAALA # pelialueen korkeus (cm)
DT = 20 # kuvan päivityksen aikaväli (ms): DT 20 on 50 FPS päivitystaajuus
MAX_VOIMA = 400 * SKAALA # suurin sallittu lyöntivoima (eli suurin mahdollinen nopeus)
DELTA_VOIMA = 0.2 * SKAALA*DT # lyönnin latauksen nopeus
KITKA = 5.0 # vierimis- ja ilmanvastus: laita tähän nolla, niin pallot eivät pysähdy ikinä
KERROKSIA = 5 # kerrosten lukumäärä alkukolmiossa
ERO = 0 * SKAALA # pallojen välinen etäisyys alkukolmiossa, oltava >= 0
KARKI = 0.75 # alkukolmion kärjen suhteellinen paikka pöydällä

VARI = "#990000" # pallojen väri RGB-heksana: vaihda jos puna-vihreä ei erotu hyvin
TIEDOSTO_TK = "aika_energia.txt" # tiedoston nimi datan kirjoittamista varten
TIEDOSTO_XY = "liikerata.txt"    # tiedoston nimi datan kirjoittamista varten

# Ohjelman käyttämiä tunnisteita erilaisten törmäysten erottamiseksi.
# Älä muuta näitä arvoja.
VAAKA = 1
PYSTY = 2
OSUMA = 3
OFFSET = 3 # piirtämistä varten


# Kirjoittaa dataa tiedostoon.
def kirjoita_tiedostoon(x_lista, y_lista, tiedoston_nimi):

    # tallennetaan listojen pituus muuttujaan n
    n = len(x_lista)
    
    # avataan tiedosto lisäämistä varten (append: 'a')
    f = open(tiedoston_nimi, 'a')

    # toistetaan n kertaa
    for i in range(n):
        # poimitaan i:s arvo x- ja y-listoista
        data = str(x_lista[i])+" "+str(y_lista[i])+"\n"
        f.write(data)

    f.close()
    

# Luo tyhjän tiedoston datan kirjoittamista varten.
def alusta_tiedosto(tiedoston_nimi):

    f = open(tiedoston_nimi, 'w')
    f.close()
    

  






# Tämä osa ohjelmaa ohjaa tähtäintä
class Kursori:

    # Luodaan tähtäin.
    def __init__(self, x, y, canvas):
        self.kulma = 0.0 # tähtäyskulma positiivisesta x-suunasta
        self.keskus_x = x # keskuspisteen koordinaatit
        self.keskus_y = y
        self.r = 8*R # piirtoetäisyys keskuspisteestä
        self.s = R   # kursorin säde
        self.alue = canvas
        
        # piirretään kursori ympyränä
        self.kuva = self.alue.create_oval(
            self.keskus_x-self.r,
            self.keskus_y-self.r,
            self.keskus_x+self.r,
            self.keskus_y+self.r,
            outline="white"
        )
        
        
    # Piilotetaan kursori siirtämällä se ulos piirtoalueelta.
    def piilota(self):
        self.keskus_x = -LEVEYS
        self.keskus_y = -KORKEUS
        
    
    # Palautetaan kursori näkyviin asettamlla sille uusi keskuspiste.
    def palauta(self, x, y):
        self.keskus_x = x
        self.keskus_y = y
        self.s = R
        # tähdätään kokhti pöydän keskustaa
        self.tahtaa_kohti(LEVEYS/2, KORKEUS/2)
        
        
    # Aseteteaan tähtäyskulma niin, että se on kohti annettuja koordinaatteja x,y.
    def tahtaa_kohti(self,x,y):

        # keskuspiseen ja tähtäyspisteen koordinaattien erotus
        dx = self.keskus_x - x+OFFSET
        dy = self.keskus_y - y+OFFSET
        
        # kulma saadaan laskettua arkustangentilla, mutta eri neljännekset
        # täytyy huomioida eri tavoin
        if dx == 0:
            if dy > 0: # suoraan ylös
                self.kulma = 3*math.pi/2
            else: # suoraan alas
                self.kulma = math.pi/2
        elif dx > 0: # vasemmalle
            self.kulma = math.pi + math.atan(dy/dx)
        else: # oikealle
            if dy > 0: # yläviistoon, kun kulman pitää olla välillä [0, 2pi]
                self.kulma = 2*math.pi + math.atan(dy/dx)
            else: # alaviistoon
                self.kulma = math.atan(dy/dx)
                

    # asetetaan piirretty tähtäin kuvassa oikeaan paikkaan
    def siirry(self):
    
        # tähtäimen koordinaatit trogonometrialla:
        x = self.keskus_x + self.r * math.cos(self.kulma)
        y = self.keskus_y + self.r * math.sin(self.kulma)

        self.alue.coords(self.kuva, 
            x-self.s + OFFSET, y-self.s + OFFSET, 
            x+self.s + OFFSET, y+self.s + OFFSET
        )


    # käännetään tähtäyskulmaa hieman myötäpäivään
    def kaanny_myotapaivaan(self):    
        self.kulma += math.pi/180 + 0.001*(1-2*random.random())
        if self.kulma > 2*math.pi:
            self.kulma -= 2*math.pi
            
    # käännetään tähtäyskulmaa hieman vastapäivään
    def kaanny_vastapaivaan(self):
        self.kulma -= math.pi/180 + 0.001*(1-2*random.random())
        if self.kulma < 0.0:
            self.kulma += 2*math.pi





# Tämä osa ohjelmaa ohjaa pelissä liikkuvia palloja
class Pallo:

    # Tämä luo pallon ja tallentaa sen ominaisuudet
    def __init__(self, x, y, canvas, vari=VARI):
        # paikkakoordinaatit
        self.x = x
        self.y = y
        # nopeuden komponentit
        self.vx = 0
        self.vy = 0
        self.alue = canvas
        self.kuva = self.alue.create_oval(
            self.x-R, self.y-R, self.x+R, self.y+R, fill=vari, outline=vari, width=0
        )
        
        
    # Siirretään pallon kuva oikeaan paikkaan.
    def siirry(self):
        self.alue.coords(self.kuva, 
            self.x-R + OFFSET, self.y-R + OFFSET,
            self.x+R + OFFSET, self.y+R + OFFSET
        )
    
    
    # Annetaan pallon liikkua suoraan ajan dt verran.    
    def liiku_suoraan(self, dt):
        self.x += self.vx*dt
        self.y += self.vy*dt
        
        
    # Pienennetään pallon vauhtia. Määrä riippuu muuttujasta KITKA.
    def hidasta(self):
        vauhti = math.sqrt(self.vx**2 + self.vy**2)
        kitka = KITKA*0.1
        if vauhti > 20*SKAALA:
            self.vx -= kitka*self.vx*0.01
            self.vy -= kitka*self.vy*0.01
        elif vauhti > 0.2*SKAALA*kitka:
            self.vx -= 0.2*SKAALA*kitka/vauhti*self.vx
            self.vy -= 0.2*SKAALA*kitka/vauhti*self.vy
        else:
            self.vx = 0
            self.vy = 0
        

    # Lasketaan pallon liike-energia.
    def laske_energia(self):

        # K = 1/2 m v^2.
        # Massaksi on asetettu 0.2 kg, mikä on suunnilleen biljardipallon massa.
        # Nopeus on yksiköissä cm/s ja sisältää skaalauskertoimen, joka pitää poistaa.
        return 0.1*(self.vx**2+self.vy**2)/(100*SKAALA)**2

        
    # Lasketaan, milloin tämä ja toinen pallo osuvat.
    # Ks. esim. https://introcs.cs.princeton.edu/java/assignments/collisions.html
    def aika_osumaan(self, toinen):    
        dx = self.x - toinen.x
        dy = self.y - toinen.y
        dvx = self.vx - toinen.vx
        dvy = self.vy - toinen.vy
        dvdr = (dvx*dx + dvy*dy)
        d = dvdr**2 - (dvx*dvx + dvy*dvy)*(dx*dx + dy*dy - 4*R*R)
        
        if dvdr < 0 and d >= 0: # pallot törmäävät vain jos tämä ehto pätee
            return - ( dvdr + math.sqrt(d) ) / (dvx*dvx + dvy*dvy)
        else: # pallot eivät törmää koskaan
            return math.inf 
        
        
    # Lasketaan, milloin pallo osuu seinään.
    def aika_valliin(self):
        pysty = math.inf 
        vaaka = math.inf
        
        # liikutaan vasemmalle
        if self.vx < 0: 
            vaaka = (R - self.x)/self.vx # törmäysaika

        # liikutaan oikealle
        elif self.vx > 0:
            vaaka = (LEVEYS - R - self.x)/self.vx
        
        # liikutaan ylös
        if self.vy < 0: 
            pysty = (R - self.y)/self.vy

        # liikutaan alas
        elif self.vy > 0:
            pysty = (KORKEUS - R - self.y)/self.vy
            
        # palautetaan ensimmäisenä tapahtuva törmäys
        if vaaka < pysty: # törmätään sivuvalliin
            return VAAKA, vaaka
            
        elif pysty < vaaka: # törmätään ylä- tai alavalliin
            return PYSTY, pysty
            
        else: # pallo on paikoillaan, eikä törmää seiniin
            return 0, math.inf
        
        
        
    # Tuhotaan kuvaan piirretty pallo.        
    def tuhoudu(self):
        self.alue.delete(self.kuva)




# Tämä ohjaa varsinaista peliä.
class Peli:

    # Aloitetaan peli.
    def __init__(self, root):
        self.root = root
        self.root.title("Biljardi")
        self.kulma = 0.0
        self.lyodaan = False
        self.odottaa = True
        self.voima = 0.0
        self.aika = [0]
        self.energia = [0.0]
        self.aikaleima = 0
        
        # canvas kuvaa peli-ikkunaa, johon tapahtumat piirretään
        self.canvas = tk.Canvas(root, width=LEVEYS, height=KORKEUS, bg="green")
        self.canvas.pack()
        
        # valkoisen pallon alkukoordinaatit
        alku_x = LEVEYS*0.25
        alku_y = KORKEUS*0.5
        self.rata_x = [alku_x]
        self.rata_y = [alku_y]
        
        # luodaan valkoinen pallo
        self.kivi = Pallo(alku_x, alku_y, self.canvas, "white")
        # lisätään valkoinen ja punaiset pallot pallojen listaan
        self.pallot = [self.kivi] + self.aseta_pallot_kolmioon()
        # luodaan tähtäin  
        self.tahtain = Kursori(alku_x, alku_y, self.canvas)      
        
        # Nämä kertovat pelille, että sen pitää tehdä jotakin.
        # Jos pelaaja painaa vasenta nuolta, käännetään tähtäystä vasemmalle.
        self.canvas.bind("<Left>", self.vasemmalle)
        # Jos pelaaja painaa oikeaa nuolta, käännetään tähtäystä oikealle.
        self.canvas.bind("<Right>", self.oikealle)
        # Jos pelaaja klikkaa tai painaa välilyöntiä, lyödään palloa.
        self.canvas.bind("<space>", self.lyo)     
        self.canvas.bind("<Button 1>", self.lyo)     
        # Jos pelaaja liikuttaa hiirtä, käännetään tähtäin hiiren suuntaan.   
        self.canvas.bind("<Motion>", self.tahtaa)
        # Tämä ei vielä tee mitään...
        self.canvas.bind("<p>", self.piirra_nopeusjakauma)

        # Tämä aktivoi peli-ikkunan.
        self.canvas.focus_set()
        
        # Käynnistetään kello, joka liikuttaa animaatiota.
        self.kello()
       
       
    # Luo joukon kohdepalloja, asettelee ne kolmioon ja palauttaa listan palloista.
    def aseta_pallot_kolmioon(self, kerroksia=KERROKSIA):
    
        pallot = []
        
        huippu_x = LEVEYS*KARKI
        huippu_y = KORKEUS*0.5
        delta = 0.01
        ero = R+2*delta+ERO
        
        # luodaaan palloja kolmioon
        # lisätään pallojen koordinaatteihin hyvin pientä
        # satunnaisuutta, jotta avauslyönnistä tulee kaoottinen
        for i in range(kerroksia):
            for j in range(i+1):
                x = huippu_x + i*ero*math.sqrt(3) + delta*(1-2*random.random())
                y = huippu_y + (2*j-i)*ero + delta*(1-2*random.random())
                uusi_pallo = Pallo(x, y, self.canvas)
                pallot.append(uusi_pallo)
                
        return pallot
            
            
    # Testaa ovatko kaikki pallot liikkumatta.
    def kaikki_on_levossa(self):
    
        # käydään läpi kaikki pallot
        for p in self.pallot:
        
            # Jos yhdenkin pallon vauhti on selvästi nollasta 
            # poikkeava, pallot eivät ole levossa.
            if p.vx**2 + p.vy**2 > 0.0001:
                return False
                
        # Kaikki olivat paikoillaan.
        return True
    
    
    # Tämä osa hoitaa tietyin väliajoin tapahtuvat asiat.
    def kello(self):
        
        # simuloi, mitä tapahtuu seuraavaksi    
        self.simuloi(0.001*DT)
        
        # hidastetaan palloja vähän ja siirretään niiden kuvat oikeisiin paikkoihin
        for p in self.pallot:
            p.hidasta()
            p.siirry()
            
        # jos lyönti on aloitettu, kasvatetaan lyöntivoimaa
        if self.lyodaan:
        
            # voima ei ole vielä maksimissaan
            if self.voima < MAX_VOIMA:
            
                # kasvatetaan voimaa aluksi hitaasti
                if self.voima < 2*DELTA_VOIMA:
                    self.voima += 0.2*DELTA_VOIMA
                else: # ja sitten jo nopeammin
                    self.voima += DELTA_VOIMA
                    
                # näytetään voimistuminen pelaajalle pienentämällä tähtäintä
                osuus = (MAX_VOIMA - self.voima) / MAX_VOIMA
                self.tahtain.s = osuus*R
                
            else: # jos maksimivoima on saavutettu, lyödään
                self.lyo(None)
                
        # siirretään tähtäin oikeaan paikkaan
        self.tahtain.siirry()
                
        # tämä ehto on totta, jos lyönti on tehty ja pallot liikkuvat
        if not self.odottaa and not self.lyodaan:
        
            # tallennetaan simulaatiosta dataa tiedostoon kirjoittamista varten
            self.energia.append( round(self.laske_kokonaisenergia(),8) )
            self.aika.append( round(self.aika[-1]+0.001*DT,4) )
            self.rata_x.append( round(self.kivi.x, 3) )    
            self.rata_y.append( round(self.kivi.y, 3) )          
            
            # Tarkistetaan, ovatko pallot jo pysähtyneet lyönnin jäljiltä.
            # Jos ovat, tallennetaan simulaatiodata ja jäädään odottamaan
            # pelaajan seuraavaa lyöntiä.
            if self.kaikki_on_levossa():
                self.tallenna_data()
                self.odottaa = True
                self.tahtain.palauta(self.kivi.x, self.kivi.y)
                
            # Tallennetaan data jokaisen 100 aika-askeleen jälkeen,
            # jotta sitä ei tarvitse pitää muistissa.
            elif len(self.energia) > 100:
                self.tallenna_data()
        
        # tehdään sama uudestaan hetken päästä
        self.canvas.after( DT-1 , self.kello)


    # Tallennetaan simulaatiodataa tiedostoon.
    # Tyhjennetään samalla listat, joihin data oli kerätty,
    # jottei muisti täyty.
    def tallenna_data(self):
        kirjoita_tiedostoon(self.aika[:-1], self.energia[:-1], TIEDOSTO_TK)
        kirjoita_tiedostoon(self.rata_x[:-1], self.rata_y[:-1], TIEDOSTO_XY)
        self.energia = [self.energia[-1]]
        self.aika = [self.aika[-1]]
        self.rata_x = [self.rata_x[-1]]
        self.rata_y = [self.rata_y[-1]]


    #
    # Simuloidaan pallojen liikettä dt-pituisen aikajakson verran eteenpäin.
    #
    # Varsinainen fysiikka on tässä funktiossa.
    #
    def simuloi(self, dt):
    
        # Muuttujaan seuraava_aika lasketaan aika, joka vielä kuluu kunnes
        # seuraava törmäys tapahtuu. Koska simuloimme vain
        # aikajakson dt, kaikki sitä pidemmät odotusajat hylätään.
        seuraava_aika = dt
        
        # Muuttujaan tyyppi tallannetaan törmäyksen tyyppi.
        tyyppi = None
        
        # Muuttujaan pallot tallennetaan törmäävät pallot.
        pallot = []
        
        # Jos törmäys tapahtuu ajassa dt, tämä saa arvon True.
        jotakin_tapahtuu = False
        
        # Lasketaan, milloin kukin pallo osuu seuraavan kerran seinään,
        # jos jatkaa suoraviivaista liikettä.
        for p in self.pallot:
            suunta, aika = p.aika_valliin()
            
            # Tarkista, onko tämä törmäys se, joka tapahtuu seuraavaksi.
            # Ts. onko odotusaika siihen lyhyempi kuin toistaiseksi lyhyin
            # löydetty odotusaika.            
            if aika < seuraava_aika:
                seuraava_aika = aika
                jotakin_tapahtuu = True
                tyyppi = suunta
                pallot = [p]
                
        # Lasketaan, milloin kukin pallopari seuraavan kerran törmää.
        # Huom. Pallot eivät välttämättä törmää koskaan, jos eivät ole liikkeessä
        # kohti toisiaan. Jos näin on, funktio aika_osumaan() palauttaa äärettömän.
        npallot = len(self.pallot)

        # Käydään läpi kaikki pallot.
        for i in range(npallot):
            p1 = self.pallot[i]
            
            # Käydään läpi pallot, jotka ovat palloa p1 ennen pallojen listassa.
            # Näin jokainen pallopari tutkitaan täsmälleen kerran.
            for j in range(i):
                p2 = self.pallot[j]
                aika = p1.aika_osumaan(p2)
                
                # Tarkista, onko tämä törmäys se, joka tapahtuu seuraavaksi.
                if aika < seuraava_aika:
                    seuraava_aika = aika
                    jotakin_tapahtuu = True
                    tyyppi = OSUMA
                    pallot = [p1, p2]
        
                   
        # Liikutetaan palloja suoraan joko koko aikaväli dt 
        # tai vain seuraavaan tapahtumaan asti.
        # Tarvittava aikaväli on joka tapauksessa tallennettuna
        # muuttujaan seuraava_aika.
        for p in self.pallot:
            p.liiku_suoraan(seuraava_aika)
            
            
        # Jos jokin palloista törmää johonkin, simulaatio on nyt ajettu
        # täsmälleen törmäyshetkeen asti.
        # Ennen kuin simulaatiota voi jatkaa, törmäävän pallon
        # nopeutta täytyy muuttaa.
        # Törmäykset oletetaan täysin elastisiksi ja pallojen pyöriminen
        # jätetään huomioimatta (ts. ei ole kierrettä).
        if jotakin_tapahtuu:

            # elastinen osuma sivuvalliin kääntää liikkeen ympäri x-suunnassa
            if tyyppi == VAAKA:
                pallot[0].vx = -pallot[0].vx
            
            # elastinen osuma ylä- tai alavalliin kääntää liikkeen ympäri y-suunnassa
            elif tyyppi == PYSTY:
                pallot[0].vy = -pallot[0].vy
            
            # Osuma toiseen palloon antaa kummallekin pallolle impulssin,
            # jonka suuruus ja suunta riippuu osumakohdasta.
            # Ks. esim. https://introcs.cs.princeton.edu/java/assignments/collisions.html
            elif tyyppi == OSUMA:
                p1 = pallot[0]
                p2 = pallot[1]
                dx = p1.x - p2.x
                dy = p1.y - p2.y
                dvx = p1.vx - p2.vx
                dvy = p1.vy - p2.vy
                dvdr = (dvx*dx + dvy*dy)
                impulssix = dvdr*dx/(4*R*R)
                impulssiy = dvdr*dy/(4*R*R)
                                
                p1.vx -= impulssix
                p1.vy -= impulssiy
                p2.vx += impulssix
                p2.vy += impulssiy                
                
            # jos koko aika dt ei vielä kulunut loppuun, simuloidaan
            # jäljellä oleva aika loppuun
            aikaa_jaljella = dt - seuraava_aika 
            self.simuloi(aikaa_jaljella)
  
    #
    # Funktio simuloi() päättyy tähän.
    #
                
            

    # Käännetään tähtäystä vasemmalle.
    def vasemmalle(self, event):
        if self.odottaa and not self.lyodaan:
            self.tahtain.kaanny_vastapaivaan()


    # Käännetään tähtäystä oikealle.
    def oikealle(self, event):
        if self.odottaa and not self.lyodaan:
            self.tahtain.kaanny_myotapaivaan()


    # Käännetään tähtäystä kohti annettuja koordinaatteja.
    def tahtaa(self, event):
        if self.odottaa and not self.lyodaan:
            self.tahtain.tahtaa_kohti(event.x, event.y)


    # Lyödään palloa. Tämä toimii niin, että ensimmäinen kutsu
    # aloittaa lyönnin ja vasta toinen viimeistelee sen.
    # Kun lyönti on aloitettu, kello()-funktio kasvattaa lyönnin
    # voimaa kunnes lyönti lopetetaan. Niinpä
    # lyönti on sitä voimakkaampi, mitä pidempään näiden
    # kahden kutsun välillä odotetaan.
    def lyo(self, event):
    
        if self.odottaa: # aloita lyönti
            self.lyodaan = True
            self.odottaa = False
            self.voima = 0.0

        elif self.lyodaan: # lopeta lyönti
            self.lyodaan = False

            # piilotetaan tähtäin, jotta pallojen liike näkyy kivasti
            suunta = self.tahtain.kulma
            self.tahtain.piilota()
        
            # annetaan valkoiselle pallolle alkunopeus tähtäyssuuntaan
            self.kivi.vx = self.voima*math.cos(suunta)
            self.kivi.vy = self.voima*math.sin(suunta)
            self.voima = 0.0



    # Laskee yhteen kaikkien pallojen liike-energiat.
    def laske_kokonaisenergia(self):
    
        # tallennetaan energia tähän muuttujaan
        e = 0.0
    
        # käydään läpi kaikki pallot
        for p in self.pallot:
        
            # lisätään pallon liike-energia kokonaisenergiaan
            e += p.laske_energia()
            
        # palautetaan lopputulos
        return e
        
       
    # Piirretään pallojen nopeuksien jakauma.
    def piirra_nopeusjakauma(self,event):
            
        # Tämä ei ole oikein mielekäs funktio, jos kyseessä on biljardisimulaatio,
        # joten eipä tehdä mitään ellei kitkaa ole poistettu.
        if KITKA > 0:
            return
            
        # Tallennetaan kokonaisenergia ja vauhdit näihin muuttujiin.
        e = 0.0
        vauhdit = []
        
        # Käydään läpi kaikki pallot sekä lasketaan niiden vauhti ja liike-energia.
        for p in self.pallot:
            e += p.laske_energia()
            vauhti = math.sqrt( p.vx**2 + p.vy**2 )
            vauhdit.append(vauhti)
            
        # Jos kaikki energia olisi yhdellä pallolla, sen vauhti olisi tämä.
        max_vauhti = np.sqrt(10.0*e)
            
        # Skaalataan vauhdit suhteelliselle asteikolle jakamalla maksimivauhdilla.
        vauhdit = np.array( vauhdit ) / (max_vauhti*100*SKAALA+0.00001)
            
        # Luodaan vektori, jossa on tasavälisesti pisteitä vaaka-akselilta.
        # Kuvaaja piirretään murtoviivana, joka kulkee näiden pisteiden kautta. 
        #
        # Komento linspace() luo vektorin, jossa on lukuja tasavälisesti.
        # Ensimmäinen argumentti on vektorin ensimmäinen luku.
        # Toinen argumentti on vektorin viimeinen luku.
        # Kolmas argumentti määrää, montako lukua vektoriin tulee.
        x_arvot = np.linspace(0,0.5,40)

        # Statistisen fysiikan teoriasta voidaan johtaa tulos, jonka mukaan
        # kaasussa molekyylien vauhdit v jakautuvat satunnaisesti noudattaen
        # Maxwell-Boltzmann-jakaumaa. Tämän jakauman tiheysfunktio on
        #
        # f(v) = a v^(d-1) exp( -b v^2 ),
        # 
        # missä a ja b ovat systeemin lämpötilasta riippuvat vakiot ja
        # d on ulottuvuuksien lukumäärä.
        #
        # Simulaatiossa on d = 2 ulottuvuutta, ja tiheysfunktio on itse asiassa nyt
        #
        # f(v) = a v exp( -a/2 v^2 ).
        #
        # Koska skaalasimme jo vauhdit, vakio a riippuu ainoastaan
        # simulaatiossa olevien pallojen lukumäärästä.
        #
        a = 2*len(self.pallot)        
        teoreettinen_jakauma = a*x_arvot * np.exp( -a/2*x_arvot**2 )
        
        # piirretään
        plt.clf() # tyhjennetään kuva      
        plt.hist(vauhdit, bins=np.linspace(0,0.5,26), density=True, label="simulaatio")
        plt.plot(x_arvot, teoreettinen_jakauma, label="teoria")
        plt.xlabel("suhteellinen nopeus ($v/v_{max}$)")
        plt.ylabel("todennäköisyystiheys")
        plt.xlim(0,0.3)
        plt.legend()
        plt.grid()
        
        plt.savefig("vauhtijakauma.pdf")
        plt.show()  



# Tämä osa luo peli-ikkunan ja käynnistää pelin.
if __name__ == "__main__":
    alusta_tiedosto(TIEDOSTO_TK) # Luodaan datalle tyhjä tiedosto.
    alusta_tiedosto(TIEDOSTO_XY) # Luodaan datalle tyhjä tiedosto.
    moottori = tk.Tk()
    game = Peli(moottori)
    moottori.mainloop()