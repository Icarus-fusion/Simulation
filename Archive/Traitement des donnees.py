import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import SGDRegressor

def load_list(path):
    data = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line != "":
                data.append(float(line))
    return data

Instant_fichier = r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Instant.txt"
Fusion_fichier  = r"C:\Users\Thomas\Documents\Cours\TIPe\Data\Fusion.txt"

Instant = load_list(Instant_fichier)
Fusion  = load_list(Fusion_fichier)

Instant_ia = np.array(Instant).reshape(-1, 1)
Fusion_ia  = np.array(Fusion)

model = SGDRegressor(max_iter=1000, eta0=0.001)
model.fit(Instant_ia, Fusion_ia)
a = model.coef_[0]
b = model.intercept_[0]
print(a, b)
print(model.score(Instant_ia, Fusion_ia))

def puissance(a, V):
    Pv = a * 2.256 * 10**(-13)
    P = Pv * (1/V) * 840
    return P

print(puissance(a, 2e-16))

plt.figure(figsize=(10, 5))
plt.plot(Instant, Fusion)
plt.plot(Instant, model.predict(Instant_ia), c='red')
plt.title("Nombre de fusion en fonction du temps")
plt.xlabel("Temps (s)")
plt.ylabel("Nombre de fusion")
plt.grid(True)
plt.show()