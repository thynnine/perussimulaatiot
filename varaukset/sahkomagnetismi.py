import tkinter as tk
import numpy as np
import copy
import time

# epsilon0 (F/m)
sahkovakio = 8.8541878128e-12
# mu0 (N/A^2)
magneettivakio = 1.25663706212e-6
# näistä johdetut vakiot
k_eps = 1/(4*np.pi*sahkovakio)
k_mu = magneettivakio/(4*np.pi)

# simulaatioparametrejä
DT = 30      # animaation jaksonaika (ms)
LEVEYS = 15  # havaintopisteiden määrä vaakasuunnassa
KORKEUS = 15 # havaintopisteiden määrä pystysuunnassa
RUUTU = 30   # ruudun koko pikseleissä
OHJE = 40    # ohjeen korkeus pikseleissä
R = 10       # varauksen säde
INERTIA = 10 # varauksen liikkeen hitaus - en suosittele muuttamaan

# Sähkö- ja magneettikentän voimakkuus divergoi lähellä pistevarausta.
# Jos kentät piirretään aidosti mittakaavassa, kuvasta ei näe mitään muuta kuin
# aivan pistevarauksen lähellä olevan kentän. Jotta kentän muoto näkyisi muuallakin,
# asetetaan kentille maksimivoimakkuus niin että piirtovaiheessa kaikki tätä
# voimakkaammat kentät piirretään samalla tavalla kuvaan.
E_MAX = 1e6   # sähkökentän maksimi
B_MAX = 1e-10 # magneettikentän maksimi
V_MAX = 1e8   # potentiaalin maksimi

# värejä
POSVARI = "#FF5555" # + varaus
NEGVARI = "#5555FF" # - varaus
EVARI = "#0000FF"   # E kenttä
BVARI = "#CC00CC"   # B kenttä
TAUSTAVARI = "#000000" # ohjeen tausta
OHJEVARI = "#FFFFFF"   # ohjeen teksti

# muita vakioita, älä muuta niitä
OFFSET = 3 
PLUSSA = 1
MIINUS = 2
VETO = 0


#
# Vektorilaskennan apufunktioita.
#

# lasketaan vektori alkupisteestä loppupisteeseen
def vektori(alku, loppu):
    return loppu-alku

# lasketaan vektorin pituus
def pituus(vektori):
    # Pythagoras: 
    # - ensin korotetaan komponentit neliöön **2
    # - sitten lasketaan neliöiden summa np.sum()
    # - lopuksi otetaan neliöjuuri np.sqrt()
    return np.sqrt(np.sum(vektori**2))
    
# lasketaan etäisyys alkupisteestä loppupisteeseen
def etaisyys(alku, loppu):
    return pituus(loppu-alku)


# tämä osa kuvaa varattuja hiukkasia
class Hiukkanen:

    # luodaan hiukkanen, jonka varaus on q ja paikkavektori r
    def __init__(self, r, q, canvas):
        self.r = r
        self.v = np.zeros(2)
        self.a = np.zeros(2)
        self.q = q
        self.kohde = copy.copy(r)
        
        # piirretään hiukkanen pallona, jonka väri valitaan varauksen mukaan
        self.alue = canvas
        if q > 0:
            self.vari = POSVARI
        else:
            self.vari = NEGVARI
        self.kuva = self.alue.create_oval(
            self.r[0]-R, self.r[1]-R, self.r[0]+R, self.r[1]+R, 
            fill=self.vari, outline=self.vari, width=0
        )
        self.seuraa = False


    # lasketaan tämän pistevarauksen luoma potentiaali tarkastelupisteessä
    def potentiaali(self, tarkastelupiste):
    
        # etäisyys hiukkasesta tarkastelupisteeseen
        # lisätään varmuuden vuoksi siihen jotakin pientä,
        # ettei se varmasti ole nolla
        r_pituus = etaisyys(self.r, tarkastelupiste)+0.0001 
    
        # Potentiaali on
        #
        # V = k q/r, missä 
        #
        # q on varaus,   
        # k = 1/(4 pi epsilon) ja 
        # r on etäisyys tarkastelupisteeseen
        #
        return k_eps * self.q / r_pituus
        
        
    # lasketaan tämän pistevarauksen luoma sähkökenttävektori tarkastelupisteessä
    def sahkokentta(self, tarkastelupiste):
    
        # vektori hiukkasesta tarkastelupisteeseen
        r_vektori = vektori(self.r, tarkastelupiste)
        
        # etäisyys hiukkasesta tarkastelupisteeseen
        # lisätään varmuuden vuoksi siihen jotakin pientä,
        # ettei se varmasti ole nolla
        r_pituus = etaisyys(self.r, tarkastelupiste)+0.0001        
        
        # hiukkasen varaus lyhyemmin kirjoitettuna
        q = self.q
        
        # Sähkökentän voimakkuus on Coulombin lain perusteella
        #
        # E = k q/r^2 , missä 
        #
        # q on varaus,        
        # k = 1/(4 pi epsilon) ja 
        # r on etäisyys tarkastelupisteeseen (tallennettu muuttujaan r_pituus).
        #
        # Kentän suunta saadaan laskemalla hiukkasesta tarkastelupisteeseen
        # osoittava yksikkövektori
        # r_hattu = r_vektori / r.
        # Huom. muuttujaan r_vektori on tallennettu valmiiksi oikea vektori,
        # mutta muuttujaa r_hattu ei ole määritelty vielä.
        #
        # Kenttävektori on siis kaikkiaan
        #
        # E_vektori = E rhattu = k q/r^2 r_hattu
        #                      = k q/r^2 r_vektori / r
        #                      = k q/r^3 r_vektori.
        #
              
        E_vektori = k_eps * self.q / r_pituus**3 * r_vektori
        
        return E_vektori
        
        
    # lasketaan tämän pistevarauksen luoma magneettikenttä tarkastelupisteessä
    def magneettikentta(self, tarkastelupiste):
    
        # vektori hiukkasesta tarkastelupisteeseen
        r_vektori = vektori(self.r, tarkastelupiste)
        
        # etäisyys hiukkasesta tarkastelupisteeseen
        # lisätään varmuuden vuoksi siihen jotakin pientä,
        # ettei se varmasti ole nolla
        r_pituus = etaisyys(self.r, tarkastelupiste)+0.0001
        
        # hiukkasen nopeusvektori
        v_vektori = self.v
        
        # hiukkasen varaus lyhyemmin kirjoitettuna
        q = self.q
        
        # Magneettikenttävektori on Biot-Savartin lain perusteella
        #
        # B_vektori = k q v_vektori x r_hattu / r^2
        #           = k q v_vektori x r_vektori / r^3,  missä
        #
        # q on varaus,   
        # k = mu / (4 pi),
        # v_vektori on hiukkasen nopeusvektori,
        # r_hattu on yksikkövektori hiukkasesta tarkastelupisteeseen
        # r on etäisyys tarkastelupisteeseen (r_pituus) ja
        # x on ristitulon kertomerkki ( Sen voi laskea funktiolla np.cross(v1,v2). )
        #
        b_vektori = k_mu * self.q * np.cross( v_vektori, r_vektori ) / r_pituus**3
        
        return b_vektori
        
    
    # Ohjaa hiukkasta kohti sille annettuja koordinaatteja.
    #
    # Käytännössä tämä toimii niin, että hiukkaseen kohdistetaan
    # jousivoimaa muistuttava kiihtyvyys, joka osoittaa kohti kohdepistettä
    # ja jonka voimakkuus on verrannollinen hiukkasen ja kohteen väliseen
    # etäisyyteen.
    def seuraa_kohdetta(self):
    
        # vektori ja etäisyys kohteeseen
        kohti_vektori = vektori(self.r, self.kohde)
        kohti_matka = etaisyys(self.r, self.kohde)
                
        # hyvin lähellä siirretään hiukkanen täsmälleen kohteeseen 
        if kohti_matka < 0.1:
            self.a = np.zeros(2)
            self.v = np.zeros(2)
            self.r =self.kohde

        # kauempana annetaan hiukkaselle kiihtyvyys kohti kohdetta
        else:
            self.a = kohti_vektori / INERTIA  
        
           
    # liikuttaa hiukkasta 
    def liiku(self):
    
        # muutetaan paikkaa nopeusvektorin suuntaan
        self.r += self.v
        
        # Muutetaan nopeutta kiihtyvyyden suuntaan.
        # Samalla hiukkaseen kohdistetaan nopeuteen verrannollinen jarruttava voima,
        # jotta se pysähtyy eikä jää pyörimään loputtomiin.
        self.v += self.a - 0.6*self.v
        
    
    # liikuta hiukkasta ja siirrä sen kuvaa näytöllä
    def paivita(self):
    
        # liikutaan
        self.seuraa_kohdetta()
        self.liiku()
    
        # asetetaan hiukkasen kuva oikeaan paikkaan
        self.alue.coords(self.kuva,
            self.r[0]-R, self.r[1]-R, self.r[0]+R, self.r[1]+R
        )
        
        # jos hiukkanen on valittu seuraamaan hiirtä, piirretään sille ääriviivat
        if self.seuraa:
            ulkovari = "black"
            w = 2
        else:
            ulkovari = self.vari
            w = 0
            
        self.alue.itemconfig(self.kuva,
            outline=ulkovari, width=w
        )
       

# tämä osa kuvaa havaintopisteitä
class Havaintopiste:

    # Luodaan havaintopiste.
    #
    # Pisteisiin liittyy potentiaalin, sähkökentän ja magneettikentän
    # graafinen esitys. Luodaan näihin tarvittavat graafiset elementit valmiiksi.
    def __init__(self, r, canvas):
        self.r = r
        self.alue = canvas
        
        # potentiaali piirretään neliönä, jonka väri muuttuu
        self.v = canvas.create_rectangle(self.r[0]-0.5*RUUTU+OFFSET, 
            self.r[1]-0.5*RUUTU+OFFSET, 
            self.r[0]+0.5*RUUTU+OFFSET, 
            self.r[1]+0.5*RUUTU+OFFSET, 
            fill="white", width=0)

        # sähkökenttävektori on kuvan tasossa ja piirretään nuolena            
        self.e_nuoli = canvas.create_line(self.r[0], self.r[1], self.r[0], self.r[1], fill=EVARI, arrow=tk.LAST, width=3) 

        # magneettikenttä on kuvaan kohtisuorassa ja piirretään ristinä tai ympyränä
        self.b_risti1 = canvas.create_line(self.r[0], self.r[1], self.r[0], self.r[1], fill=BVARI) 
        self.b_risti2 = canvas.create_line(self.r[0], self.r[1], self.r[0], self.r[1], fill=BVARI) 
        self.b_ympyra = canvas.create_oval(self.r[0], self.r[1], self.r[0], self.r[1], outline=BVARI) 


    # päivitetään havaintopisteen graafinen esitys vastaamaan
    # annettua sähkökenttää e, magneettikenttää b ja potentiaalia v
    def paivita(self, e, b, v):
    
        # havaintopisteen keskipiste
        x = self.r[0]+OFFSET
        y = self.r[1]+OFFSET
        
        # sähkökentän x- ja y-komponentit
        ex = e[0]
        ey = e[1]
        
        # magneettikentän z-komponentti
        bz = b[2]
        
        # suureiden itseisarvot
        eabs = pituus(e)
        babs = np.abs(bz)
        vabs = np.abs(v)
        
        # Jos kenttä on liian voimakas, skaalataan sen arvo
        # lähemmäs nollaa MAX-arvoon.
        #
        # Näin kannattaa tehdä, sillä muuten kentästä nähdään
        # vain pistevarauksien lähimmät pisteet.
        #
        if eabs > E_MAX:            
            ex *= E_MAX/eabs
            ey *= E_MAX/eabs
            eabs = E_MAX
            
        if babs > B_MAX:
            bz *= B_MAX/babs
            babs = B_MAX
            
        if vabs > V_MAX:
            v *= V_MAX/vabs
            vabs = V_MAX
        
        # sähkökenttävektorin x- ja y-komponentit kuvassa
        dx = ex*0.5*RUUTU/E_MAX
        dy = ey*0.5*RUUTU/E_MAX
        # vektorinuolen paksuus
        w = 2.0*eabs/E_MAX
        
        # päivitetään sähkökenttänuoli
        self.alue.coords(self.e_nuoli, 
            x-dx, y-dy, 
            x+dx, y+dy
        )
        self.alue.itemconfig(self.e_nuoli,
            width=w, arrowshape=(3*w,3*w,w)
        )
        
        # magneettikenttävektoria kuvaavan ristin tai ympyrän koko kuvassa
        bb = bz*0.4*RUUTU/B_MAX
        # viivojen paksuus
        w = 2.0*babs/B_MAX
        
        # jos Bz > 0 magneettikenttä osoittaa kuvan sisään,
        # ja kenttä piirretään ristinä
        if bb > 0:
        
            # päivitetään ristiä kuvaavat viivat
            self.alue.coords(self.b_risti1, 
                x-bb, y-bb, 
                x+bb, y+bb
            )
            self.alue.coords(self.b_risti2, 
                x-bb, y+bb, 
                x+bb, y-bb
            )
            # kutistetaan ympyrä pisteeksi
            self.alue.coords(self.b_ympyra, 
                x, y, 
                x, y
            )
            self.alue.itemconfig(self.b_risti1,
                width=w
            )
            self.alue.itemconfig(self.b_risti2,
                width=w
            )
            self.alue.itemconfig(self.b_ympyra,
                width=0
            )
            
        # jos Bz < 0, magneettikenttä osoittaa kuvasta ulos,
        # ja kenttä piirretään ympyränä
        else:
        
            # kutistetaan risti pisteeksi
            self.alue.coords(self.b_risti1, 
                x, y, 
                x, y
            )
            self.alue.coords(self.b_risti2, 
                x, y, 
                x, y
            )
            # päivitetään kenttää kuvaava ympyrä
            self.alue.coords(self.b_ympyra, 
                x+bb, y+bb, 
                x-bb, y-bb
            )
            self.alue.itemconfig(self.b_risti1,
                width=0
            )
            self.alue.itemconfig(self.b_risti2,
                width=0
            )
            self.alue.itemconfig(self.b_ympyra,
                width=w
            )
            
        # nollapotentiaalia kuvaavan värin RGB-arvot
        red   = 255
        green = 200
        blue  = 0
        
        # Muutetaan pisteen taustaväriä sen potentiaalin mukaan.
        #
        # Jos V < 0, muutetaan väriä tummemmaksi.
        if v < 0:
            red -= int(vabs*150/V_MAX)
            green -= int(vabs*120/V_MAX)
            
        # Jos V > 0, muutetaan väriä vaaleammaksi.
        else:
            green += int(vabs*50/V_MAX)
            blue += int(vabs*150/V_MAX)
        
        # päivitetään väri
        v_vari = '#%02x%02x%02x' % (red, green, blue)
        
        self.alue.itemconfig(self.v,
            fill=v_vari
        )


# Tämä ohjaa graafista ikkunaa.
class Avaruus:

    # Aloitetaan simulaatio.
    def __init__(self, root):
        self.root = root
        self.root.title("Sähkömagneettinen kenttä")
        
        # ohjausmoodi voi olla
        # - PLUSSA: luodaan positiivisia varauksia
        # - MIINUS: luodaan negatiivisia varauksia
        # - VETO: siirretään varauksia
        self.moodi = PLUSSA
        
        # canvas kuvaa ikkunaa, johon tapahtumat piirretään
        self.canvas = tk.Canvas(root, width=LEVEYS*RUUTU, height=KORKEUS*RUUTU+OHJE, bg=TAUSTAVARI)
        self.canvas.pack()
                
        # kirjoitetaan ikkunan alareunaan käyttöohje
        self.ohje = self.canvas.create_rectangle(-5, KORKEUS*RUUTU, LEVEYS*RUUTU+5, KORKEUS*RUUTU+OHJE+5, fill=TAUSTAVARI, tag="ohje")
        self.ohjeteksti = self.canvas.create_text(0.2*OHJE, KORKEUS*RUUTU+0.9*OHJE, anchor=tk.SW,
            text=self.anna_ohje(),font=("Futura", str(int(0.3*OHJE))), fill=OHJEVARI, tag="ohje")
        
        # luodaan havaintopisteiden taulukko     
        self.pisteet = [ ]
        self.luo_havaintoverkko()
        
        # lista varauksille - se on aluksi tyhjä
        self.hidut = [ ]
                
        # Nämä kertovat simulaatiolle, että sen pitää tehdä jotakin.  
        self.canvas.bind("<Button 1>", self.hiiren_nappi)       
        self.canvas.bind("<Motion>", self.veda)
        self.canvas.bind("<Button 2>", self.vaihda_moodi)
        self.canvas.bind("<space>", self.vaihda_moodi)

        # Tämä aktivoi ikkunan.
        self.canvas.focus_set()
        
        # Käynnistetään kello, joka liikuttaa animaatiota.
        self.kello()
       
    
    # vaihdetaan ohjausmoodia
    def vaihda_moodi(self, event):
    
        # Jos jokin hiukkanen seurasi hiirtä, lopetetaan seuraaminen
        # käskemällä jokaista hitua erikseen lopettamaan.
        for h in self.hidut:
            if h.seuraa:
                h.seuraa = False
    
        # vaihdetaan moodia
        self.moodi = (self.moodi+1)%3
        
        # vaihdetaan ohjetekstiä
        self.canvas.itemconfig(self.ohjeteksti, text=self.anna_ohje())
        
        
    # kirjoitetaan ohje kulloisenkin moodin mukaan
    def anna_ohje(self):
        if self.moodi == PLUSSA:
            teksti = "hiiri: luo ja tuhoa positiivisia varauksia"
        elif self.moodi == MIINUS:
            teksti = "hiiri: luo ja tuhoa negatiivisia varauksia"
        elif self.moodi == VETO:
            teksti = "hiiri: valitse varaus ja ohjaa se uuteen paikkaan"
            
        teksti += "\nvälilyönti: vaihda toimintoa"
        
        return teksti
    
        
    # Tehdään tämä kun hiiren nappia painetaan.
    # Tämä kutsuu toisia funktioita sen mukaan, mikä moodi on päällä.
    def hiiren_nappi(self, event):
        if self.moodi == VETO:
            self.poimi(event)
        else:
            self.luo(event)
            
            
    # Luodaan tai tuhotaan pistevaraus.
    #
    # Jos hiirellä klikattiin tyhjää kohtaa, sinne luodaan uusi varaus.
    #
    # Jos hiirellä klikattiin varausta, se tuhotaan.
    def luo(self, event):
    
        # poimitaan klikkauspisteen koordinaatit
        r = np.array( [event.x, event.y], dtype = 'float' )
        
        # Tarkistetaan, osuiko klikkaus varaukseen.
        # Jos osui, otetaan varauksen indeksi talteen,
        # sillä haluamme poistaa klikatut varaukset.
        tuhottavat = []
        for i in range(len(self.hidut)):
            h = self.hidut[i]
            
            # klikkaus osui, jos se oli tarpeeksi lähellä
            # varauksen koordinaatteja (piiretyn ympyrän sisällä)
            if etaisyys(r, h.r) < R:
                tuhottavat = [i] + tuhottavat
                
        # Jos tuhottavien listaan on kerätty varauksia,
        # käydään ne läpi ja poistetaan ne simulaatiosta.
        for i in tuhottavat:
            pois = self.hidut.pop(i) # poistetaan simulaatiosta
            self.canvas.delete( pois.kuva ) # poistetaan kuvasta
            
             
        # jos yhtään hiukkasta ei löytynyt tuhottavaksi,
        # luodaan uusi hiukkanen
        if len(tuhottavat) == 0:
            # uuden varauksen merkki riippuu moodista
            if self.moodi == PLUSSA:
                q = 1
            elif self.moodi == MIINUS:
                q = -1
            else:
                return

            # luodaan uusi hiukkanen
            uusi_hitu = Hiukkanen(r, q, self.canvas)
            self.hidut = [uusi_hitu] + self.hidut   
            
        self.canvas.tag_raise("ohje")
        
       
    # Valitaan simulaatiosta yksi hiukkanen, joka alkaa seurata hiirtä.
    def poimi(self,event):
    
        # poimitaan klikkauspisteen koordinaatit
        r = np.array( [event.x, event.y], dtype = 'float' )
             
        poimittu = False
        
        # Käydään läpi kaikki hiukkaset.
        for h in self.hidut:
        
            # Emme halua usean hiukkasen seuraavan hiirtä, koska silloin ne
            # päätyvät lopulta päällekäin.
            #
            # Niinpä jos hiukkanen seurasi jo hiirtä, käsketään sen lopettaa.
            if h.seuraa:
                h.seuraa = False
                
            # Jos kuitenkin klikkaus osui hiukkaseen, joka ei vielä seurannut,
            # käsketään sen aloittaa hiiren seuraaminen.
            #
            # Huom. koska haluamme vain yhden hiukkasen seuraavan, tämä
            # tehdään vain ensimmäiselle sopivalle hiukkaselle.
            # Muuttuja poimittu huolehtii tästä.
            elif etaisyys(r, h.r) < R and not poimittu:
                h.seuraa = True
                h.kohde = r
                poimittu = True
                

    # Jos jokin hiukkanen on määrätty seuraamaan hiirtä, kerrotaan
    # hiukkaselle, että sen kohde on siellä, missä hiiri on.
    def veda(self,event):
        for h in self.hidut:
            if h.seuraa:
                h.kohde = np.array( [event.x, event.y], dtype = 'float' )
                
                
                    
                   
    # Luodaan kuvaruudun täyttävä havaintopisteiden verkko.
    def luo_havaintoverkko(self):
    
        # vaakasuuntaan LEVEYS pistettä
        for i in range(LEVEYS):
            # pystysuuntaan KORKEUS pistettä
            for j in range(KORKEUS):
            
                # havaintopisteen koordinaatit pikseleissä
                x = (i+0.5)*RUUTU
                y = (j+0.5)*RUUTU
                
                # koordinaatit paikkavektorina
                r = np.array( [x,y] , dtype = 'float' )
                
                # luodaan havaintopiste graafisena oliona
                piste = Havaintopiste( r , self.canvas )
                self.pisteet.append(piste)
                
        
    # Lasketaan sähkömagneettinen kenttä kaikissa havaintopisteissä.
    def laske_kentta(self):

        # käydään läpi kaikki havaintopisteet
        for p in self.pisteet:
        
            # tallennetaan kentät näihin muuttujiin
            e = np.zeros(2) # sähkökenttä (xy-tason vektori)
            b = np.zeros(3) # magneettikenttä (3D-vektori)
            v = 0.0         # potentiaali (skalaari)

            # käydään läpi kaikki hiukkaset ja summataan
            # niiden tuottamat kentät yhteen
            #
            # havaintopisteen paikkavektori on p.r
            for h in self.hidut:
                e += h.sahkokentta(p.r)
                b += h.magneettikentta(p.r)
                v += h.potentiaali(p.r)
                
            # päivitetään kentän graafinen esitys
            p.paivita(e,b,v)
        
    # Siirretään varauksia sekä simulaatiossa että kuvassa.
    def aseta_varaukset(self):
        for h in self.hidut:
            h.paivita()
        

    # Tämä osa hoitaa tietyin väliajoin tapahtuvat asiat.
    def kello(self):
        
        # otetaan aikaa
        t0 = time.time()

        # simuloidaan
        self.aseta_varaukset()
        self.laske_kentta()

        # otetaan aikaa
        t1 = time.time()
        aikaa_kului = int((t1-t0)*1000)        
        
        # yritetään pitää tasainen päivitystaajuus
        if aikaa_kului < DT:
            odota = DT-aikaa_kului
        else:
            odota = 1
        
        # tehdään sama uudestaan hetken päästä
        self.canvas.after( odota , self.kello)



# Tämä osa luo peli-ikkunan ja käynnistää pelin.
if __name__ == "__main__":
    moottori = tk.Tk()
    simulaatio = Avaruus(moottori)
    moottori.mainloop()