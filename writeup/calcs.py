import numpy as np

solarmass = 2e33
sigma_t = 6.652e-25
sigma = 5.670374e-5
c = 3e10
G = 6.67e-8
mp = 1.67262192e-24
me = 9.1093837e-28
ev_per_erg = 6.242e+11
k = 1.3807e-16

luminosity = (4 * np.pi * G * solarmass * mp * c) / sigma_t
print(f"Luminosity {luminosity} erg / s")

rs_sun = 2 * G * solarmass / c**2
flux_scale = luminosity / (4 * np.pi * rs_sun**2)
print(f"Flux {flux_scale} erg / s / cm^2")

temp_scale = k * (flux_scale / sigma)**0.25
print(f"Temp scale {temp_scale * ev_per_erg} eV")

scale_height = np.sqrt(4 * temp_scale / (mp * c**2))
print(f"Scale height scale {scale_height}")

corona_height = me * c**2 / 2 * ev_per_erg
print(f"Corona height scale {corona_height} / temp in eV")

travel_time = 2 * G * solarmass / c**2 * sigma_t / me
print(f"Travel scale {travel_time} cm^3 / g")

print(f"Optical depth {3.7e-3 * 250 * 2.4e12**(-1/8) * 2}")