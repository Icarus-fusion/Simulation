import matplotlib.pyplot as plt

def load_list(path):
    data = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line != "":
                data.append(float(line))
    return data

Instant_fichier = r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Instantv.txt"
Vitesse_fichier  = r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Vitesse.txt"

Instant = load_list(Instant_fichier)
Vitesse  = load_list(Vitesse_fichier)

plt.figure(figsize=(10, 5))
plt.plot(Instant, Vitesse)
plt.title("Évolution de la vitesse moyenne")
plt.xlabel("Temps (s)")
plt.ylabel("Vitesse (m/s)")
plt.grid(True)
plt.show()