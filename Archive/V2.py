import tkinter as tk
import random
import math

# --- Constants ---
WIDTH = 600
HEIGHT = 600
NUM_PARTICLES = 50
DT_RATIO = 0.5  # Ratio of Deuterium to total
RADIUS = 4

# Physical constants
K_B = 1.380649e-23  # Boltzmann constant (J/K)
T = 1e8  # Temperature in Kelvin (100 million K)

# Masses (kg)
M_D = 3.343583719e-27  # Deuterium
M_T = 5.007356744e-27  # Tritium

# Time control (1 = real time, 0 = pause, 0.5 = half speed, etc.)
TIME_SCALE = 1  # Try changing this to 0, 1, or any float between

#0.0000000000000001 pour une visibilité de la vitesse

class Particle:
    def __init__(self, canvas, x, y, vx, vy, particle_type):
        self.canvas = canvas
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.type = particle_type
        self.color = "red" if particle_type == "D" else "blue"
        self.id = self.canvas.create_oval(
            self.x - RADIUS, self.y - RADIUS,
            self.x + RADIUS, self.y + RADIUS,
            fill=self.color
        )

    def move(self):
        self.x = (self.x + self.vx * TIME_SCALE) % WIDTH
        self.y = (self.y + self.vy * TIME_SCALE) % HEIGHT
        self.canvas.coords(
            self.id,
            self.x - RADIUS, self.y - RADIUS,
            self.x + RADIUS, self.y + RADIUS
        )

class FusionSimulation:
    def __init__(self, root):
        self.canvas = tk.Canvas(root, width=WIDTH, height=HEIGHT, bg="black")
        self.canvas.pack()
        self.particles = self.create_particles()
        self.animate()

    def maxwell_velocity(self, mass):
        std_dev = math.sqrt(K_B * T / mass)
        vx = random.gauss(0, std_dev)
        vy = random.gauss(0, std_dev)
        return vx * 1e10, vy * 1e10  # scale for visibility

    def create_particles(self):
        particles = []
        for _ in range(NUM_PARTICLES):
            x = random.uniform(0, WIDTH)
            y = random.uniform(0, HEIGHT)
            p_type = "D" if random.random() < DT_RATIO else "T"
            mass = M_D if p_type == "D" else M_T
            vx, vy = self.maxwell_velocity(mass)
            particles.append(Particle(self.canvas, x, y, vx, vy, p_type))
        return particles

    def animate(self):
        for particle in self.particles:
            particle.move()
        self.canvas.after(30, self.animate)

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Fusion Simulation - Time Control")
    sim = FusionSimulation(root)
    root.mainloop()
