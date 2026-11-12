import tkinter as tk
import random
import time
import math

#Manipulation
zoom = 1

#Particules
N = 500                                                                 #Nombre de particules
rdt = 0.5                                                               #Ratio entre le nombr de particule de Deutérium et de Tritium

#Dimensions
V = (210 * N) / (2.4e23)#m^3                                            #Volume simulé
r = (V / (2 * math.pi**2 * (16 / 9)))**(1/3)                            #Petit rayon du tore
R = (16 / 9) * r                                                        #Grand rayon du tore
m_p = ((2 * math.pi * R) / 1600)                                        #Nombre de mètre par pixel
Largueur = 1600
Hauteur = 900
Taille = 2e-15#m                                                        #Taille d'un atome d'hydrogène
Taille_p = Taille / m_p                                                 #Taille d'un atome d'hydrogène en pixel

#Température
T = 1.5e8#K                                                             #Température à laquelle le plasma est maintenue

#Constantes
k_b = 1.380649e-23#J/K                                                  #Constante de Boltzmann en J/K
epsilon_0 = 8.854e-12#F/m                                               #Permittivité du vide en F/m
q = 1.602e-19#C                                                         #Charge élémentaire
M_d = 3.343583719e-27#kg                                                #Masse des atomes de Deutérium
M_t = 5.007356744e-27#kg                                                #Masse des atomes de Tritium

#Maxwell
V_MAX = 1.5e7#m/s                                                       #Vitesse maximum des particules (5% de c)
N_sub = 10000                                                           #Nombre de subdivision de l'axe des abcsisses de la courbe des vitesses
facteur_d = math.sqrt(2 / math.pi) * (M_d / (k_b * T))**(3/2)           #Facteur multiplicateur Deutérium pour la densité de probabilité maxwellienne
facteur_t = math.sqrt(2 / math.pi) * (M_t / (k_b * T))**(3/2)           #Facteur multiplicatif Tritium pour la densité de probabilité maxwellienne

#Traitement des données
with open(r"/home/thomas/Info", "a") as f:
    f.write(
        "Dimension :" + "\n"
        "\n"
        "Volume : " + str(V) + "\n"
        "Petit r : " + str(r) + "\n"
        "Grand R : " + str(R) + "\n"
        "Mètre par pixel : " + str(m_p) + "\n"
        "\n"
        "Paramètre :" + "\n"
        "\n"
        "Température : " + str(T) + "\n"
        "Nombre de particule : " + str(N) + "\n"
        "Taille d'un atome en pixel : " + str(Taille_p)
    )
start = time.perf_counter()

#Calcul du rayon de la section efficace en m
#Cette fonction a été obtenue à partir de l'approximation linéaire de la coubre DT du graphique des section efficaces page 3 du diaporama
def section_efficace(E):
    if E < 2.3:
        return 0.0
    elif 2.3 <= E < 4.3:
        sigma = (4.95e-4 * (E - 2.3) + 10**-5) * 1e-24
    elif 4.3 <= E < 10:
        sigma = ((23/1900) * (E - 4.3) + 10**-3) * 1e-24
    elif 10 <= E < 11.05:
        sigma = ((31/35) * (E - 10) + 0.07) * 1e-24
    elif 11.05 <= E < 13.7:
        sigma = ((100/53) * (E - 11.05) + 1) * 1e-24
    elif 13.7 <= E < 15:
        sigma = ((12/13) * (E - 13.7) + 6) * 1e-24
    elif 15 <= E < 18.2:
        sigma = ((-11/32) * (E - 15) + 7.2) * 1e-24
    elif 18.2 <= E < 290:
        sigma = ((-53/2720) * (E - 18.2) + 6.1) * 1e-24
    elif 290 <= E < 800:
        sigma = ((-61/51000) * (E - 290) + 0.8) * 1e-24
    elif 800 <= E < 1500:
        sigma = ((-9/70000) * (E - 800) + 0.19) * 1e-24
    elif 1500 <= E < 5000:
        sigma = ((-1/175000) * (E - 1500) + 0.1) * 1e-24
    return math.sqrt(sigma/math.pi)

#Calcul la probabiité que la particule aille à la vitesse de la particule soit comprise entre v et v+dv à l'aide de la formule parge 2 du diaporama
def maxwell(v, particule):
    if particule == 'D':
        masse = M_d
        facteur = facteur_d
    elif particule == 'T':
        masse = M_t
        facteur = facteur_t
    return facteur * v**2 * math.exp(-masse * v**2 / (2 * k_b * T))

#Calcul de la vitesse d'une particule à partir de la répartition maxwellienne des vitesses en m/s
def vitesse(particule):
    Probabilité = []
    Vitesse = []
    for i in range(N_sub + 1):
        Vitesse.append(i * V_MAX/N_sub)
        Probabilité.append(maxwell(i * V_MAX/N_sub, particule))
    return random.choices(Vitesse, weights=Probabilité, k=1)[0]

class Particule:
    #Définition de ce qu'est une particule pour la simulation
    def __init__(self, canvas, x, y, v, theta, type, masse):
        self.canvas = canvas
        self.masse = masse
        self.x = x
        self.y = y
        self.v = v
        self.vx = v * math.cos(theta)
        self.vy = v * math.sin(theta)
        self.type = type
        self.color = "red" if type == "D" else "blue"
        self.section = self.section_efficace()
        #Dessin de l'atome
        self.atome = self.canvas.create_oval(
            self.x - Taille_p, self.y - Taille_p,
            self.x + Taille_p, self.y + Taille_p,
            fill = self.color
        )
        #Dessin de la section efficace
        self.section_d = self.canvas.create_oval(
            self.x - self.section, self.y - self.section,
            self.x + self.section, self.y + self.section,
            outline = "green"
        )
    
    #Cette fonction permet de claculer la distance entre deux particules.
    #On ajuste la distance entre les particules sachant que nous somme dans un tore.
    def distance(self, other):
        dx = abs(self.x - other.x)
        dy = abs(self.y - other.y)
        if dx > Largueur / 2:
            dx -= Largueur
        if dy > Hauteur / 2:
            dy -= Hauteur
        dist = math.sqrt((dx * m_p)**2 + (dy * m_p)**2)
        return dist, dx, dy
    
    #Cette fonction permet de calculer la force exercer sur la particule self par toutes les autres particules de la simulation.
    def force_electrostatique(self, others):
        fx = 0.0
        fy = 0.0
        for other in others:
            if other is self:
                continue
            dist, dx, dy = self.distance(other)
            if dist > 10e-20:
                return 0.0, 0.0
            force = (q * q) / (4 * math.pi * epsilon_0 * dist**2)
            fx -= force * (dx * m_p) / dist
            fy -= force * (dy * m_p) / dist
        return fx, fy
    
    #Cette fonction permet de calculer l'energie cinéatique de la particule afin d'en déduire
    #sa section efficace.
    def section_efficace(self):
        vitesse = math.sqrt(self.vx**2 + self.vy**2)
        energie = 0.5 * self.masse * vitesse**2 * 6.242 * 10**15
        return section_efficace(energie)
    
    #Cette fonction régit le mouvement de la particules en fonction de toutes les autres
    def move(self, others):
        dt = 3e-15
        fx, fy = self.force_electrostatique(others)
        ax = fx / self.masse
        ay = fy / self.masse
        self.vx += ax * dt
        self.vy += ay * dt
        self.x = (self.x + (self.vx * dt) / m_p) % Largueur
        self.y = (self.y + (self.vy * dt) / m_p) % Hauteur
        self.section = self.section_efficace() / m_p
        #On met à jour la position des dessins
        self.canvas.coords(
            self.atome,
            self.x - Taille_p, self.y - Taille_p,
            self.x + Taille_p, self.y + Taille_p
        )
        self.canvas.coords(
            self.section_d,
            self.x - self.section, self.y - self.section,
            self.x + self.section, self.y + self.section,
        )

class Simulation:
    #On définit l'environement de la simulation
    def __init__(self, root):
        self.canvas = tk.Canvas(root, width=Largueur, height=Hauteur, bg="black")
        self.canvas.pack()
        self.fusion_count = 0
        self.fusion_label = tk.Label(root, text="Fusions: 0", fg="white", bg="black", font=("Arial", 12))
        self.fusion_label.place(relx=1.0, rely=1.0, x=-10, y=-10, anchor="se")
        self.particules = self.creation_particules()
        self.animation()

    #Cette fonction permet d'ajouter les particules à l'intérieur de la simulation
    def creation_particules(self):
        particules = []
        for _ in range(int(N)):
            x = random.uniform(0, Largueur)
            y = random.uniform(0, Hauteur)
            type = "D" if random.random() < rdt else "T"
            theta = random.uniform(0, 2 * math.pi)
            masse = M_d if type == 'D' else M_t
            v = vitesse(type)
            particules.append(Particule(self.canvas, x, y, v, theta, type, masse))
        return particules
    
    #Cette fonction permet la fusion des particules, elle compare la distance entre les deux particules avec les sections efficaces
    def fusion(self):
        for i in range(len(self.particules)):
            p1 = self.particules[i]
            for j in range(i+1, len(self.particules)):
                p2 = self.particules[i]
                if (p1.type == p2.type):
                    continue
                dist_px = p1.distance(p2)[0] / m_p
                if dist_px < (p1.section + p2.section):
                    self.fusion_count += 1
                    self.fusion_label.config(text=f"Fusions: {self.fusion_count}")
                    with open(r"/home/thomas/Fusion.txt", "a") as f:
                        f.write(str(self.fusion_count) + "\n")
                    with open(r"/home/thomas/Instant.txt", "a") as f:
                        f.write(str(time.perf_counter() - start) + "\n")
                    x = 0
                    y1 = random.uniform(0, Hauteur)
                    y2 = random.uniform(0, Hauteur)
                    v1 = vitesse("D")
                    v2 = vitesse("T")
                    theta1 = random.uniform(0, 2 * math.pi)
                    theta2 = random.uniform(0, 2 * math.pi)
                    p1.x = x
                    p1.y = y1
                    p1.theta = theta1
                    p1.v = v1
                    p2.x = x
                    p2.y = y2
                    p2.theta = theta2
                    p2.v = v2
    
    #Cette fonction permet la simulation, c'est elle qui met en mouvement la simulation
    def animation(self):
        self.fusion()
        for particule in self.particules:
            particule.move(self.particules)
        self.canvas.after(30, self.animation)

#Lancement de la simulation
root = tk.Tk()
root.title("Simulation de fusion nucléaire")
sim = Simulation(root)
root.mainloop()