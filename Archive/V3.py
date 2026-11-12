import tkinter as tk
import random
import math

# --- Constants ---
WIDTH = 600
HEIGHT = 600
NUM_PARTICLES = 50
DT_RATIO = 0.5
RADIUS = 4

# Physical constants
K_B = 1.380649e-23  # Boltzmann constant (J/K)
T = 1e8  # Temperature in Kelvin
M_D = 3.343583719e-27  # kg
M_T = 5.007356744e-27  # kg
TIME_SCALE = 0.1

# Maxwell speed sampling resolution
V_MAX = 2e6  # max speed (m/s) to consider in the sampling table
N_STEPS = 1000

def compute_maxwell_cdf(mass, temp):
    """Precompute the CDF of the Maxwellian distribution"""
    dv = V_MAX / N_STEPS
    v_vals = [i * dv for i in range(N_STEPS + 1)]
    f_vals = []
    factor = math.sqrt(2 / math.pi) * (mass / (K_B * temp)) ** 1.5
    for v in v_vals:
        f = factor * v ** 2 * math.exp(-mass * v ** 2 / (2 * K_B * temp))
        f_vals.append(f)
    # Normalize to form CDF
    total = sum(f_vals)
    cdf = []
    cumulative = 0
    for f in f_vals:
        cumulative += f
        cdf.append(cumulative / total)
    return v_vals, cdf

def sample_speed(v_vals, cdf):
    """Sample a speed using inverse transform sampling"""
    u = random.random()
    for i in range(len(cdf)):
        if u <= cdf[i]:
            return v_vals[i]
    return v_vals[-1]

class Particle:
    def __init__(self, canvas, x, y, v, theta, particle_type):
        self.canvas = canvas
        self.x = x
        self.y = y
        self.vx = v * math.cos(theta) * 1e-4  # scale to pixels/frame
        self.vy = v * math.sin(theta) * 1e-4
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
        # Precompute CDF tables
        self.vd_vals, self.cd_vals = compute_maxwell_cdf(M_D, T)
        self.vt_vals, self.ct_vals = compute_maxwell_cdf(M_T, T)
        self.particles = self.create_particles()
        self.animate()

    def create_particles(self):
        particles = []
        for _ in range(NUM_PARTICLES):
            x = random.uniform(0, WIDTH)
            y = random.uniform(0, HEIGHT)
            p_type = "D" if random.random() < DT_RATIO else "T"
            if p_type == "D":
                v = sample_speed(self.vd_vals, self.cd_vals)
            else:
                v = sample_speed(self.vt_vals, self.ct_vals)
            angle = random.uniform(0, 2 * math.pi)
            particles.append(Particle(self.canvas, x, y, v, angle, p_type))
        return particles

    def animate(self):
        for particle in self.particles:
            particle.move()
        self.canvas.after(30, self.animate)

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Fusion Simulation - Maxwellian Formula")
    sim = FusionSimulation(root)
    root.mainloop()
