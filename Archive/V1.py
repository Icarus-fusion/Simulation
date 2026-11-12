import tkinter as tk
import random
import math

# Simulation parameters
WIDTH = 600
HEIGHT = 600
NUM_PARTICLES = 50
DT_RATIO = 0.5  # Ratio of Deuterium to total
RADIUS = 4  # Radius for particle display
VELOCITY_SCALE = 1.5

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
        self.x = (self.x + self.vx) % WIDTH
        self.y = (self.y + self.vy) % HEIGHT
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

    def create_particles(self):
        particles = []
        for _ in range(NUM_PARTICLES):
            x = random.uniform(0, WIDTH)
            y = random.uniform(0, HEIGHT)
            angle = random.uniform(0, 2 * math.pi)
            speed = random.uniform(0.5, VELOCITY_SCALE)
            vx = speed * math.cos(angle)
            vy = speed * math.sin(angle)
            p_type = "D" if random.random() < DT_RATIO else "T"
            particles.append(Particle(self.canvas, x, y, vx, vy, p_type))
        return particles

    def animate(self):
        for particle in self.particles:
            particle.move()
        self.canvas.after(30, self.animate)

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Fusion Simulation - Deuterium & Tritium")
    sim = FusionSimulation(root)
    root.mainloop()
