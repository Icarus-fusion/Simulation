#Simulation d'une réaction de fusion nucléaire dans un espace en deux dimensions

import tkinter as tk                                                    #Module permettant la simulation
import random
import math

#Zoom
zoom = 1

#Dimensions:
Largueur = 1600                                                         #Largueur de la fenêtre en pixels
Hauteur = 900                                                           #Hauteur de la fenêtre en pixels
m_p = 5e-10 / 1                                                         #Mètre par pixel

#Particules
N = 1000                                                                #Nombre de particules
rdt = 0.5                                                               #Ratio entre Deutérium et Tritium
Taille = 1.4e-13                                                        #Taille des particules en pixels
#8,4e-16m

#Température
T = 1.5e8                                                               #Température initial du plasma

#Constantes
k_b = 1.380649e-23                                                      #Constante de Boltzmann en J/K
epsilon_0 = 8.854e-12                                                   #Permittivité du vide en F/m
q = 1.602e-19                                                           #Charge élémentaire en C
M_d = 3.343583719e-27                                                   #Masse des atomes de Deutérium en kg
M_t = 5.007356744e-27                                                   #Masse des atomes de Tritium en kg

#Maxwell
V_MAX = 1.5e7                                                           #Vitesse maximum des particules en m/s
N_sub = 10000                                                           #Nombre de subdivision de l'axe des abcsisses de la courbe des vitesses
facteur_d = math.sqrt(2 / math.pi) * (M_d / (k_b * T))**(3/2)           #Facteur multiplicateur Deutérium
facteur_t = math.sqrt(2 / math.pi) * (M_t / (k_b * T))**(3/2)           #Facteur multiplicatif Tritium

#Temps
temps = 0.00000000000001#%                                               #Pourcentage de la vitesse d'écoulement du temps
fr = 30#ms                                                               #Période de rafraichissement de la position des particules

def cross_section(E):                                                    #Calcul du rayon de la section efficace en m
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

def maxwell(v, particule):                                              #Calcul de la courbe de répartition maxwellienne des vitesses du plasma
    if particule == 'D':
        masse = M_d
        facteur = facteur_d
    elif particule == 'T':
        masse = M_t
        facteur = facteur_t
    return facteur * v**2 * math.exp(-masse * v**2 / (2 * k_b * T))

def vitesse(particule):                                                 #Calcul de la vitesse d'une particule à partir de la répartition maxwellienne des vitesses en m/s
    Probabilité = []
    Vitesse = []
    for i in range(N_sub + 1):
        Vitesse.append(i * V_MAX/N_sub)
        Probabilité.append(maxwell(i * V_MAX/N_sub, particule))
    return random.choices(Vitesse, weights=Probabilité, k=1)[0]

class Particle:
    def __init__(self, canvas, x, y, v, theta, p_type, masse):
        self.canvas = canvas
        self.masse = masse
        self.x = x
        self.y = y
        self.theta = theta
        self.v = v
        self.vx = v * math.cos(theta)
        self.vy = v * math.sin(theta)
        self.type = p_type
        self.color = "red" if p_type == "D" else "blue"
        #Section efficace initial
        self.section = self.cross_section() / m_p
        #Dessine l'atome
        self.id = self.canvas.create_oval(
            self.x - Taille, self.y - Taille,
            self.x + Taille, self.y + Taille,
            fill=self.color
        )
        #Dessine la section efficace
        self.id_cs = self.canvas.create_oval(
            self.x - self.section, self.y - self.section,
            self.x + self.section, self.y + self.section,
            outline="green"
        )
    
    def distance(self, other):
        dx = abs(self.x - other.x)
        dy = abs(self.y - other.y)
        #Ajustement des distances car nous sommes dans un tore
        if dx > Largueur / 2:
            dx = dx - Largueur
        if dy > Hauteur / 2:
            dy = dy - Hauteur
        dist_px = math.sqrt(dx**2 + dy**2)                              #Distance entre les deux particules en pixel
        dist = dist_px * m_p
        return dist, dx, dy

    def force_électrostatique(self, others):
        fx = 0.0
        fy = 0.0
        for other in others :
            if other is self:
                continue
            dist, dx, dy = self.distance(other)
            if dist > 10e-20:
                return 0.0, 0.0
            force = (q * q) / (4* math.pi * epsilon_0 * dist**2)
            fx = fx - force * (dx / (dist / m_p))
            fy = fy - force * (dy / (dist / m_p))
        return fx, fy
    
    def cross_section(self):
        vitesse = math.sqrt(self.vx**2 + self.vy**2)
        energy = 0.5 * self.masse * vitesse**2
        return cross_section(energy)

    def move(self, others):
        dt = fr * temps / 1000
        fx, fy = self.force_électrostatique(others)
        ax = fx / self.masse
        ay = fy / self.masse
        self.vx = self.vx + ax * dt
        self.vy = self.vy + ay * dt
        self.x = (self.x + ((self.vx * dt) / m_p)) % Largueur
        self.y = (self.y + ((self.vy * dt) / m_p)) % Hauteur
        self.section = self.cross_section() / m_p
        self.canvas.coords(
            self.id,
            self.x - Taille, self.y - Taille,
            self.x + Taille, self.y + Taille
        )
        self.canvas.coords(
            self.id_cs,
            self.x - self.section, self.y - self.section,
            self.x + self.section, self.y + self.section,
        )

class FusionSimulation:
    def __init__(self, root):
        self.canvas = tk.Canvas(root, width=Largueur, height=Hauteur, bg="black")
        self.canvas.pack()
        self.fusion_count = 0
        self.fusion_label = tk.Label(root, text="Fusions: 0", fg="white", bg="black", font=("Arial", 12))
        self.fusion_label.place(relx=1.0, rely=1.0, x=-10, y=-10, anchor="se")
        self.particules = self.create_particules()
        self.animate()

    def create_particules(self):
        particules = []
        for _ in range(int(N)):
            x = random.uniform(0, Largueur)
            y = random.uniform(0, Hauteur)
            p_type = "D" if random.random() < rdt else "T"
            theta = random.uniform(0, 2 * math.pi)
            masse = M_d if p_type == "D" else M_t
            v = vitesse("D" if random.random() < rdt else "T")
            particules.append(Particle(self.canvas, x, y, v, theta, p_type, masse))
        return particules

    def fusion(self):
        to_remove = set()
        for i in range(len(self.particules)):
            if i in to_remove:
                continue
            p1 = self.particules[i]
            for j in range(i+1, len(self.particules)):
                if j in to_remove:
                    continue
                p2 = self.particules[j]
                if (p1.type == p2.type):
                    continue
                dist = p1.distance(p2)[0]
                if dist < ((p1.section * m_p) + (p2.section * m_p)):
                    to_remove.add(i)
                    to_remove.add(j)
                    self.canvas.delete(p1)
                    self.canvas.delete(p1.id_atom)
                    self.canvas.delete(p1.id_ring)
                    self.canvas.delete(p2)
                    self.canvas.delete(p2.id_atom)
                    self.canvas.delete(p2.id_ring)
                    self.fusion_label.config(text=f"Fusions: {self.fusion_count}")
                    self.fusion_count += 1
        self.particules = [p for k, p in enumerate(self.particules) if k not in to_remove]

    def animate(self):
        self.fusion()
        for particule in self.particules:
            particule.move(self.particules)
        self.canvas.after(fr, self.animate)


root = tk.Tk()
root.title("Simulation fusion nucléaire")
sim = FusionSimulation(root)
root.mainloop()