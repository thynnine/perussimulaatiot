import tkinter as tk
import numpy as np

#
# Kokeile muuttaa systeemin dynamiikkaa.
#
# Systeemin dynamiikan määrittelevät parametrit luetaan tiedostosta. 
# Valitse sopiva tiedosto poistamalla rivin edestä kommenttimerkki.
#

#DYNAMIIKKA = "ei_dynamiikkaa.txt"
#DYNAMIIKKA = "virtaus_poistavilla_reunoilla.txt"
#DYNAMIIKKA = "virtaus_jaksollisilla_reunoilla.txt"
#DYNAMIIKKA = "diffuusio_poistavilla_reunoilla.txt"
#DYNAMIIKKA = "diffuusio_jaksollisilla_reunoilla.txt"
DYNAMIIKKA = "aalto_heijastavilla_reunoilla.txt"
#DYNAMIIKKA = "aalto_jaksollisilla_reunoilla.txt"

# muita simulaatioparametrejä
DT = 30      # animaation jaksonaika (ms)
RUUTU = 20   # ruudun koko pikseleissä
R = 3*RUUTU  # alkupulssin koko

# värejä 
#
# näitä voi halutessaan muuttaa
#
# värit on annettu muodossa [puna, viher, sini] kokonaislukuna 0, ..., 255
HUIPPU = np.array( [255, 0,   0  ], dtype='int' )
KESKI  = np.array( [255, 255, 255], dtype='int' )
POHJA  = np.array( [0,   0,   255], dtype='int' )

# ohjelman sisäisiä muuttujia, joihin ei ole syytä koskea
MAX = 10
KOKO = 40
OFFSET = 3 


#
# Muuttaa RGB-väriarvon heksadesimaalimuotoon, jota tkinter ymmärtää.
#
def rgb_heksaksi(rgb):    
    return '#%02x%02x%02x' % (rgb[0],rgb[1],rgb[2])
    

#
# Lukee tiedostosta simulaation dynamiikkaa kuvaavan matriisin.
#    
def lue_simulaation_tiedot(tiedoston_nimi):
    global KOKO

    # avataan ja luetaan tiedosto
    f = open(tiedoston_nimi)
    rivit = f.readlines()
    f.close()
    
    # ensimmäisellä rivillä on simulaatioalueen koko
    KOKO = int(rivit[0])
    
    # toisella rivillä on simulaatiossa tarvittava skalaarikerroin
    kerroin = float(rivit[1])
        
    # loput rivit määrittelevät dynamiikkamatriisin nollasta poikkeavat alkiot
    #
    # kullakin rivillä on ensin rivin ja sarakkeen indeksit sekä matriisialkion arvo
    matriisi = np.zeros([KOKO*KOKO, KOKO*KOKO])
    for rivi in rivit[2:]:
        osat = rivi.split() # jaetaan rivi osiin välilyöntien kohdalla
        i = int(osat[0])
        j = int(osat[1])
        m = float(osat[2])
        
        matriisi[i,j] = m

    # kerrotaan vähän, mitä tiedostosta luettiin     
    print()   
    print("Luettiin dynamiikkaparametrit tiedostosta "+str(tiedoston_nimi))
    print("  simulaation koko:    "+str(KOKO)+" x "+str(KOKO))
    print("  skalaarikerroin:     "+str(kerroin))
    print("  matriisin koko:      "+str(KOKO*KOKO)+" x "+str(KOKO*KOKO))
    print("  alkioita (ei nolla): "+str(len(rivit)-2))
    print("  matriisin alku:      ")
    kerrotaan = 5
    for i in range(kerrotaan):
        rivi = ""
        for j in range(kerrotaan):
            rivi += str( '{0:.2f}'.format(matriisi[i,j]) )+" "
        print("  "+rivi+" ... ")

    rivi = ""
    for j in range(kerrotaan):
        rivi += " ... "
    print("  "+rivi)
    print()
        
        
    return kerroin, matriisi

    

# tämä osa kuvaa havaintopisteitä
class Havaintopiste:

    # luodaan havaintopiste ja sen graafinen esitys    
    def __init__(self, x, y, canvas):
        self.x = x
        self.y = y
        self.arvo = 0
        self.alue = canvas
        
        # piste piirretään neliönä, jonka väri muuttuu
        self.v = canvas.create_rectangle(self.x*RUUTU+OFFSET, 
            self.y*RUUTU+OFFSET, 
            (self.x+1)*RUUTU+OFFSET, 
            (self.y+1)*RUUTU+OFFSET, 
            fill=rgb_heksaksi(KESKI), width=0)


    # päivitetään havaintopisteen graafinen esitys
    def paivita(self, arvo):
        
        self.arvo = arvo
        vari = np.zeros(3, dtype='int')
                
        # Muutetaan pisteen taustaväriä sen arvon mukaan.
        # 
        # Jos arvo on negatiivinen, se värjätään värien KESKI ja POHJA sekoituksella.
        if arvo < 0:
            suhde = min(MAX,-arvo)/MAX
            vari = np.int32((1-suhde)*KESKI + suhde*POHJA)
            
        # Jos arvo on positiivinen, se värjätään värien KESKI ja HUIPPU sekoituksella.
        else:
            suhde = min(MAX,arvo)/MAX
            vari += np.int32((1-suhde)*KESKI + suhde*HUIPPU)
        
        # muutetaan väri tkinterin käyttämään heksadesimaalimuotoon
        heksavari = rgb_heksaksi(vari)
        
        # värjätään hvaintopiste uudella värillä
        self.alue.itemconfig(self.v,
            fill=heksavari
        )


# Tämä ohjaa graafista ikkunaa.
class Avaruus:

    # Aloitetaan simulaatio.
    def __init__(self, root):
        self.root = root
        self.root.title("Dynamiikkaa matriiseilla")

        # luetaan simulaation ominaisuudet tiedostosta
        self.kerroin, self.matriisi = lue_simulaation_tiedot(DYNAMIIKKA)
        
        # canvas kuvaa ikkunaa, johon tapahtumat piirretään
        self.canvas = tk.Canvas(root, width=KOKO*RUUTU, height=KOKO*RUUTU)
        self.canvas.pack()

        # luodaan havaintopisteiden taulukko     
        self.pisteet = [ ]
        self.luo_havaintoverkko()
        
        # luodaan vektori, johon tallennetaan funktion arvot havaintopisteissä
        self.arvot = np.zeros([KOKO*KOKO])
        self.vanhat_arvot = np.zeros([KOKO*KOKO])
        
        # Nämä kertovat simulaatiolle, että sen pitää tehdä jotakin.  
        self.canvas.bind("<Button 1>", self.positiivinen_pulssi)    
        self.canvas.bind("<Button 2>", self.negatiivinen_pulssi)    
        self.canvas.bind("<space>", self.nollaa) 

        # Tämä aktivoi ikkunan.
        self.canvas.focus_set()
        
        # Käynnistetään kello, joka liikuttaa animaatiota.
        self.kello()
            

    # luodaan positiivinen pulssi
    def positiivinen_pulssi(self,event):
        self.pulssi(event,1)
        
        
    # luodaan negatiivinen pulssi
    def negatiivinen_pulssi(self,event):
        self.pulssi(event,-1)


    # luo väliaineeseen Gaussisen pulssin
    def pulssi(self, event, merkki):
    
        # hiiren klikkauksen koordinaatit
        x = event.x-OFFSET
        y = event.y-OFFSET
        
        # koordinaatteja vastaavan havaintopisteen indeksit
        i = x // RUUTU
        j = y // RUUTU
        
        # luodaan pulssi tyhjään vektoriin
        alku = np.zeros([KOKO*KOKO])
        
        # Gaussin funktion leveys
        r = (2*R) // RUUTU
        
        # Käydään läpi vektorin ne alkiot, jotka poikkeavat
        # klikkauspisteen indekseistä korkeintaan r verran.
        # Näihin lisätään Gaussin funktion arvo.
        #
        # Nyt siis klikkauspisteen x- ja y-indeksit ovat i ja j,
        # ja näistä poiketaan di ja dj.
        #
        for di in range(-r,r+1):
            for dj in range(-r, r+1):
            
                # tarkasteltavan alkion indeksi vektorissa
                k = i+di + (j+dj)*KOKO
                    
                # tarkastetaan, ettei mennä ulos simulaatiosta
                if i+di < KOKO and i+di >= 0 and j+dj < KOKO and j+dj >= 0:
                    alku[k] = merkki*5*MAX*np.exp(10/r**2*(-di**2-dj**2))  
                
        # aaltodynamiikkaa varten myös aikaisemmat arvot pitää päivittää
        sitten = 0.5*self.matriisi @ alku
        
        # lisätään edellä lasketun Gaussin funktion arvot simulaatioon
        self.vanhat_arvot += sitten
        self.arvot += alku
                
        
    # aseta väliaine kaikkialla nollaan
    def nollaa(self,event):
    
        self.arvot[:] = 0.0
        self.vanhat_arvot[:] = 0.0
                
                   
    # Luodaan kuvaruudun täyttävä havaintopisteiden matriisi listojen listana.
    def luo_havaintoverkko(self):
    
        self.pisteet = []
    
        # vaakasuuntaan KOKO pistettä
        for i in range(KOKO):
        
            pistesarake = []
            
            # pystysuuntaan KOKO pistettä
            for j in range(KOKO):
            
                # luodaan havaintopiste graafisena oliona
                piste = Havaintopiste( i, j , self.canvas )
                pistesarake.append(piste)
            
            self.pisteet.append(pistesarake)
                


    # Tämä osa hoitaa tietyin väliajoin tapahtuvat asiat.
    def kello(self):

        # simuloidaan
        self.simuloi()
        self.paivita_kuva()
        
        # tehdään sama uudestaan hetken päästä
        self.canvas.after( DT , self.kello)
        
                
        
    # käydään läpi kaikki havaintopisteet ja päivitetään
    # niiden graafinen esitys
    def paivita_kuva(self):
    
        for i in range(KOKO):
            for j in range(KOKO):
            
                k = i+j*KOKO
                arvo = self.arvot[k]
                self.pisteet[i][j].paivita(arvo)
        
        
        
    # Päivitä väliaine simulaattorin dynamiikkamatriisin mukaisesti
    def simuloi(self):
    
        M = self.matriisi
        V = self.arvot
        U = self.vanhat_arvot
        k = self.kerroin
    
        #
        # Edistä simulaatiota hiukan ajassa eteenpäin laskemalla
        # väliaineen uudet arvot hetken päästä tulevaisuudessa.
        # 
        # Periaate on hyvin yksinkertainen.
        # Väliaineen nykyiset ominaisuudet on vektorissa V.
        # Väliaineen ominaisuudet hetki sitten on vektorissa U.
        # Väliaineen dynamiikka on matriisissa M ja skalaarissa k.
        #
        # Väliaineen ominaisuudet hetken päästä saadaan laskemalla
        #
        # M V - k U
        #
        # Laske tämä lauseke ja tallenna se muuttujaan nimeltä uudet_arvot.
        #
        # Ole tarkkana siinä, että käytät oikeaa kertolaskua.
        # M V on matriisitulo kun taas k U on skalaarin ja vektorin tulo.
        #
        
        uudet_arvot = M @ V - k * U

        # tallennetaan vanhat arvot aaltodynamiikkaa varten
        self.vanhat_arvot = self.arvot
        self.arvot = uudet_arvot



# Tämä osa luo peli-ikkunan ja käynnistää pelin.
if __name__ == "__main__":
    moottori = tk.Tk()
    simulaatio = Avaruus(moottori)
    moottori.mainloop()