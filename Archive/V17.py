import tkinter as tk
import random
import time
import math

#Manipulation
zoom = 1

#Particules
N = 100                                                                 #Nombre de particules
rdt = 0.5                                                               #Ratio entre le nombre de particule de Deutérium et de Tritium

#Dimensions
qVs = N / (8.4e22)
V = 840 * qVs                                                           #Volume simulé
r = (V / (2 * math.pi**2 * 3.1))**(1/3)                                 #Petit rayon du tore
R = 3.1 * r                                                             #Grand rayon du tore
m_p = ((2 * math.pi * R) / 1240)                                        #Nombre de mètre par pixel
Largueur = 1240
Hauteur = 400
Taille = 2e-15#m                                                        #Taille d'un atome d'hydrogène en mètre
Taille_p = Taille / m_p                                                 #Taille d'un atome d'hydrogène en pixel

#Température
T = 1.5e8#K                                                             #Température à laquelle le plasma est maintenue

#Constantes
k_b = 1.380649e-23#J/K                                                  #Constante de Boltzmann en J/K
epsilon_0 = 8.854e-12#F/m                                               #Permittivité du vide en F/m
e = 1.602e-19#C                                                         #Charge élémentaire en C
M_d = 3.343583719e-27#kg                                                #Masse des atomes de Deutérium en kg
M_t = 5.007356744e-27#kg                                                #Masse des atomes de Tritium en kg

#Conversion :

barn = 1e-24#cm2
cm = 0.01#m
jouletoeV = 6.242e18#eV

#Maxwell
V_MAX = 1.5e7#m/s                                                       #Vitesse maximum des particules (50% de c)
N_sub = 10000                                                           #Nombre de subdivision de l'axe des abcsisses de la courbe des vitesses
facteur_d = math.sqrt(2 / math.pi) * (M_d / (k_b * T))**(3/2)           #Facteur multiplicateur Deutérium pour la densité de probabilité maxwellienne
facteur_t = math.sqrt(2 / math.pi) * (M_t / (k_b * T))**(3/2)           #Facteur multiplicatif Tritium pour la densité de probabilité maxwellienne

#Traitement des données
with open(r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Info.txt", "a") as f:
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
    if E < 2.3 or E >= 5000:
        return 0.0
    elif 2.3 <= E < 4.3:
        sigma = 10**(7.36 * math.log(E) - 7.66) * barn
    elif 4.3 <= E < 10:
        sigma = 10**(5.03 * math.log(E) - 6.19) * barn
    elif 10 <= E < 11.05:
        sigma = 10**(19.03 * math.log(E) - 20.18) * barn
    elif 11.05 <= E < 13.7:
        sigma = 10**(10.24 * math.log(E) - 10.86) * barn
    elif 13.7 <= E < 15:
        sigma = 10**(2.01 * math.log(E) - 1.51) * barn
    elif 15 <= E < 18.2:
        sigma = 10**(-0.86 * math.log(E) + 1.87) * barn
    elif 18.2 <= E < 290:
        sigma = 10**(-0.73 * math.log(E) + 1.71) * barn
    elif 290 <= E < 800:
        sigma = 10**(-1.36 * math.log(E) + 3.26) * barn
    elif 800 <= E < 1500:
        sigma = 10**(-1.10 * math.log(E) + 2.50) * barn
    elif 1500 <= E < 5000:
        sigma = 10**(-0.19 * math.log(E) - 0.41) * barn
    return math.sqrt(sigma/math.pi) * cm

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
        self.section = self.section_efficace() / m_p
    
    #Cette fonction permet de claculer la distance entre deux particules.
    #On ajuste la distance entre les particules sachant que nous somme dans un tore.
    def distance(self, other):
        dx = abs(self.x - other.x)
        if dx > Largueur / 2:
            dx -= Largueur
        dy = abs(self.y - other.y)
        if dy > Hauteur / 2:
            dy -= Hauteur
        return dx * m_p, dy * m_p
    
    #Cette fonction permet de calculer la force exercée sur la particule self par toutes les autres particules de la simulation.
    def force_electrostatique(self, others):
        fx = 0.0
        fy = 0.0
        for other in others:
            if other is self:
                continue
            dx, dy = self.distance(other)
            if dx == 0 or dy == 0:
                fx += 0
                fy += 0
            else:
                fx -= e**2 / (4 * math.pi * epsilon_0 * dx**2)
                fy -= e**2 / (4 * math.pi * epsilon_0 * dy**2)
        return fx, fy
    
    #Cette fonction permet de calculer l'energie cinéatique de la particule afin d'en déduire sa section efficace.
    def section_efficace(self):
        energie = 0.5 * self.masse * (self.vx**2 + self.vy**2) * jouletoeV * 0.001
        return section_efficace(energie)
    
    #Cette fonction régit le mouvement de la particules en fonction de toutes les autres
    def move(self, others):
        dt = 3e-15
        fx, fy = self.force_electrostatique(others)
        ax = fx / self.masse
        ay = fy / self.masse
        self.vx += ax * dt
        self.vy += ay * dt
        self.v = math.sqrt(self.vx**2 + self.vy**2)
        self.x = (self.x + (self.vx * dt) / m_p) % Largueur
        self.y = (self.y + (self.vy * dt) / m_p) % Hauteur
        self.section = self.section_efficace() / m_p

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
    
    #Cette fonction permet la fusion des particules, elle compare la distance entre les deux particules avec la sommes de leurs sections efficaces
    def fusion(self):
        for i in range(len(self.particules)):
            p1 = self.particules[i]
            for j in range(i+1, len(self.particules)):
                p2 = self.particules[j]
                if (p1.type == p2.type):
                    continue
                dx, dy = p1.distance(p2)
                dist_px = (dx**2 + dy**2) / m_p
                if dist_px < (p1.section + p2.section):
                    self.fusion_count += 1
                    self.fusion_label.config(text=f"Fusions: {self.fusion_count}")
                    with open(r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Fusion.txt", "a") as f:
                        f.write(str(self.fusion_count) + "\n")
                    with open(r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Instant.txt", "a") as f:
                        f.write(str(time.perf_counter() - start) + "\n")
                    y = 200
                    if p1.type == "D":
                        x1 = 310
                        x2 = 930
                    else:
                        x1 = 930
                        x2 = 310
                    v1 = vitesse(p1.type)
                    v2 = vitesse(p2.type)
                    theta1 = random.uniform(0, 2 * math.pi)
                    theta2 = random.uniform(0, 2 * math.pi)
                    p1.x = x1
                    p1.y = y
                    p1.theta = theta1
                    p1.v = v1
                    p1.vx = p1.v * math.cos(theta1)
                    p1.vy = p1.v * math.sin(theta1)
                    p2.x = x2
                    p2.y = y
                    p2.theta = theta2
                    p2.v = v2
                    p2.vx = p2.v * math.cos(theta2)
                    p2.vy = p2.v - math.sin(theta2)
    
    #Cette fonction permet la simulation, c'est elle qui met en mouvement la simulation et permet la fusion des particules
    def animation(self):
        self.fusion()
        moy = 0
        for particule in self.particules:
            particule.move(self.particules)
            moy += abs(particule.v)
        moyenne = 1/N * moy
        with open(r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Vitesse.txt", "a") as f:
            f.write(str(moyenne) + "\n")
        with open(r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Instantv.txt", "a") as f:
            f.write(str(time.perf_counter() - start) + "\n")
        self.canvas.after(30, self.animation)

#Lancement de la simulation
root = tk.Tk()
root.title("Simulation de fusion nucléaire")
sim = Simulation(root)
root.mainloop()