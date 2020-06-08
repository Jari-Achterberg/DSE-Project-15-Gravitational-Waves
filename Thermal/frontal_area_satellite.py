import numpy as np

height = 30
width = 20
length = 20

alphas = np.arange(0, 0.5 * np.pi, 0.001 * np.pi)
betas = np.arange(0, 0.5 * np.pi, 0.001 * np.pi)


def rad_to_deg(rad):
    return rad * (180 / np.pi)


surface_areas = []
for alpha in alphas:
    for beta in betas:
        v_width = 20 * np.cos(alpha)
        v_height = 30 * np.cos(beta)
        v_area = v_width * v_height

        z_width = 20 * np.sin(alpha)
        z_height = 30 * np.cos(beta)
        z_area = z_width * z_height

        o_width = 20
        o_height = 20
        o_area = o_width * o_height * np.sin(beta)

        total_area = v_area + z_area + o_area

        surface_areas.append([total_area, rad_to_deg(alpha), rad_to_deg(beta)])

surface_areas = np.array(surface_areas)
max_area = max(surface_areas[:, 0])
max_values = np.where(surface_areas[:, 0] == max_area)
print(surface_areas[max_values])
