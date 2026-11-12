import tkinter as tk
import random
import math
import time

# --- Constants ---
WIDTH = 1200
HEIGHT = 1000
NUM_PARTICLES = 100
DT_RATIO = 0.5
RADIUS = 4

k_b = 1.380649e-23
EPSILON_0 = 8.854e-12  
Q = 1.602e-19          
T = 1.5e8  # Kelvin
M_D = 3.343583719e-27
M_T = 5.007356744e-27
TIME_SCALE = 0.00000000000000000001

V_MAX = 1.5e7
N_STEPS = 10000
M_d = 3.343583719e-27#kg                                                #Masse des atomes de Deutérium
M_t = 5.007356744e-27#kg                                                #Masse des atomes de Tritium
facteur_d = math.sqrt(2 / math.pi) * (M_d / (k_b * T))**(3/2)           #Facteur multiplicateur Deutérium pour la densité de probabilité maxwellienne
facteur_t = math.sqrt(2 / math.pi) * (M_t / (k_b * T))**(3/2)           #Facteur multiplicatif Tritium pour la densité de probabilité maxwellienne

start = time.perf_counter()

CROSS_SECTION_TABLE = [
    (10e3 * 1.602e-19, 0.01),
    (50e3 * 1.602e-19, 0.5),
    (100e3 * 1.602e-19, 1.0),
    (200e3 * 1.602e-19, 0.6),
    (500e3 * 1.602e-19, 0.2),
]

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
    for i in range(N_STEPS + 1):
        Vitesse.append(i * V_MAX/N_STEPS)
        Probabilité.append(maxwell(i * V_MAX/N_STEPS, particule))
    return random.choices(Vitesse, weights=Probabilité, k=1)[0]

def compute_cross_section(E):
    """Linear interpolation of cross-section from energy (in J)"""
    for i in range(len(CROSS_SECTION_TABLE) - 1):
        (E1, σ1), (E2, σ2) = CROSS_SECTION_TABLE[i], CROSS_SECTION_TABLE[i + 1]
        if E1 <= E <= E2:
            return σ1 + (σ2 - σ1) * (E - E1) / (E2 - E1)
    return 0.0  # Outside defined range

def sigma_to_radius(sigma_barns):
    """Convert cross-section in barns to radius in pixels, physically scaled"""
    sigma_m2 = sigma_barns * 1e-28  # convert barns to m²
    radius_m = math.sqrt(sigma_m2 / math.pi)  # effective radius
    scale = 4e15  # 4 pixels = 1e-15 m
    return radius_m * scale

class Particle:
    def __init__(self, canvas, x, y, v, theta, particle_type, mass):
        self.canvas = canvas
        self.mass = mass
        self.v = v
        self.theta = theta
        self.vx = v * math.cos(theta) * 1e-4
        self.vy = v * math.sin(theta) * 1e-4
        self.type = particle_type
        self.color = "red" if particle_type == "D" else "blue"
        self.x = x
        self.y = y

        # Initial kinetic energy and cross section
        self.cross_section = self.compute_cross_section()
        self.radius_sigma = sigma_to_radius(self.cross_section)

        # Draw inner atom
        self.id_atom = self.canvas.create_oval(
            self.x - RADIUS, self.y - RADIUS,
            self.x + RADIUS, self.y + RADIUS,
            fill=self.color
        )

        # Draw cross section ring
        self.id_ring = self.canvas.create_oval(
            self.x - self.radius_sigma, self.y - self.radius_sigma,
            self.x + self.radius_sigma, self.y + self.radius_sigma,
            outline="green"
        )

    def compute_electrostatic_force(self, others):
        fx, fy = 0.0, 0.0
        for other in others:
            if other is self:
                continue
            dx = other.x - self.x
            dy = other.y - self.y

            # Apply toroidal boundary correction
            if dx > WIDTH / 2: dx -= WIDTH
            elif dx < -WIDTH / 2: dx += WIDTH
            if dy > HEIGHT / 2: dy -= HEIGHT
            elif dy < -HEIGHT / 2: dy += HEIGHT

            dist_px = math.sqrt(dx**2 + dy**2)
            meters_per_pixel = 1e-15 / 4
            r = dist_px * meters_per_pixel  # convert to meters

            if r < 1e-15:  # avoid singularity and unphysically strong forces
                continue

            force = (Q * Q) / (4 * math.pi * EPSILON_0 * r**2)

            # REVERSE DIRECTION: because same charges REPEL
            fx -= force * (dx / dist_px)
            fy -= force * (dy / dist_px)
        return fx, fy

    def compute_cross_section(self):
        speed = math.sqrt(self.vx**2 + self.vy**2)
        kinetic_energy = 0.5 * self.mass * speed**2
        return compute_cross_section(kinetic_energy)

    def move(self, others):
        dt = 30e-23
        fx, fy = self.compute_electrostatic_force(others)

        self.ax = fx / self.mass
        self.ay = fy / self.mass

        self.vx += self.ax * dt
        self.vy += self.ay * dt

        meters_per_pixel = 1e-15 / 4
        self.x = (self.x + (self.vx * dt) / meters_per_pixel) % WIDTH
        self.y = (self.y + (self.vy * dt) / meters_per_pixel) % HEIGHT

        self.cross_section = self.compute_cross_section()
        self.radius_sigma = sigma_to_radius(self.cross_section)

        self.canvas.coords(
            self.id_atom,
            self.x - RADIUS, self.y - RADIUS,
            self.x + RADIUS, self.y + RADIUS
        )
        self.canvas.coords(
            self.id_ring,
            self.x - self.radius_sigma, self.y - self.radius_sigma,
            self.x + self.radius_sigma, self.y + self.radius_sigma
        )


class FusionSimulation:
    def __init__(self, root):
        self.canvas = tk.Canvas(root, width=WIDTH, height=HEIGHT, bg="black")
        self.canvas.pack()
        self.fusion_count = 0
        self.fusion_label = tk.Label(root, text="Fusions: 0", fg="white", bg="black", font=("Arial", 12))
        self.fusion_label.place(relx=1.0, rely=1.0, x=-10, y=-10, anchor="se")
        self.particles = self.create_particles()
        self.animate()


    def create_particles(self):
        particles = []
        for _ in range(NUM_PARTICLES):
            x = random.uniform(0, WIDTH)
            y = random.uniform(0, HEIGHT)
            p_type = "D" if random.random() < DT_RATIO else "T"
            mass = M_D if p_type == "D" else M_T
            v = vitesse(p_type)
            angle = random.uniform(0, 2 * math.pi)
            particles.append(Particle(self.canvas, x, y, v, angle, p_type, mass))
        return particles
    
    def fusion(self):
        for i in range(len(self.particles)):
            p1 = self.particles[i]
            for j in range(i+1, len(self.particles)):
                p2 = self.particles[j]
                if (p1.type == p2.type):
                    continue
                dx = abs(p1.x - p2.x)
                dy = abs(p1.y - p2.y)
                if dx > WIDTH / 2:
                    dx = WIDTH - dx
                if dy > HEIGHT / 2:
                    dy = HEIGHT - dy
                dist = math.hypot(dx, dy)
                if dist < (p1.radius_sigma + p2.radius_sigma):
                    self.fusion_count += 1
                    self.fusion_label.config(text=f"Fusions: {self.fusion_count}")
                    with open(r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Fusion.txt", "a") as f:
                        f.write(str(self.fusion_count) + "\n")
                    with open(r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Instant.txt", "a") as f:
                        f.write(str(time.perf_counter() - start) + "\n")
                    p1.v = 0
                    p2.v = 0
                    x1 = 0
                    x2 = random.uniform(0, WIDTH)
                    y1 = random.uniform(0, HEIGHT)
                    y2 = 0
                    v1 = vitesse(p1.type)
                    v2 = vitesse(p2.type)
                    theta1 = random.uniform(0, 2 * math.pi)
                    theta2 = random.uniform(0, 2 * math.pi)
                    p1.x = x1
                    p1.y = y1
                    p1.ax = 0
                    p1.ay = 0
                    p1.theta = theta1
                    p2.x = x2
                    p2.y = y2
                    p2.ax = 0
                    p2.ay = 0
                    p2.theta = theta2
                    p1.v = v1
                    p2.v = v2

    def animate(self):
        self.fusion()
        moy = 0
        for particule in self.particles:
            particule.move(self.particles)
            moy += abs(particule.v)
        moyenne = 1/NUM_PARTICLES * moy
        with open(r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Vitesse.txt", "a") as f:
            f.write(str(moyenne) + "\n")
        with open(r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Instantv.txt", "a") as f:
            f.write(str(time.perf_counter() - start) + "\n")
        self.canvas.after(30, self.animate)


if __name__ == "__main__":
    root = tk.Tk()
    root.title("Fusion Simulation")
    sim = FusionSimulation(root)
    root.mainloop()