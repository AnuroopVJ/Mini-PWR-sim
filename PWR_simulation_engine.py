import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import scipy.constants

from data_parser import load_cross_section_data_from_json, get_energy_xs_arrays, linear_interpolate_xs

FUEL_POSITION = (0, 0)
FUEL_RADIUS = 300
WATER_DESNITY = 1.0 #g/cm^3
WATER_POSITION = (0, 0)
WATER_RADIUS = 600

avogadro = scipy.constants.Avogadro

def random_initial_vals():
    angle = np.random.uniform(0, 2 * np.pi)
    direction_vector = [np.cos(angle), np.sin(angle)]
    return direction_vector

def get_fission_yields_for_energy(product_distributions, MT_key, neutron_energy):
    results = {}
    if MT_key not in product_distributions:
        print(f"MT key {MT_key} not found.")
        return results
    subsections = product_distributions[MT_key].get("subsection", {})
    for subsec_key, subsec in subsections.items():
        yields = subsec.get("yields", {})
        energies = yields.get("Eint", [])
        yi = yields.get("yi", [])
        if not energies or not yi or len(energies) != len(yi):
            continue
        energies = np.array(energies)
        yi = np.array(yi)
        if neutron_energy <= energies[0]:
            interp_yi = yi[0]
        elif neutron_energy >= energies[-1]:
            interp_yi = yi[-1]
        else:
            interp_yi = np.interp(neutron_energy, energies, yi)
        results[subsec_key] = interp_yi
    return results

def check_neutron_intersection(h, k, r, direction):
    if direction[0] == 0:
        x = h
        y1 = k + math.sqrt(r**2 - (h - x)**2)
        y2 = k - math.sqrt(r**2 - (h - x)**2)
        return [(x, y1), (x, y2)]
    m = direction[1] / direction[0]
    A = 1 + m**2
    B = 2 * (m * (k - k) - h)
    C = h**2 + (k - k)**2 - r**2
    discriminant = B**2 - 4 * A * C
    if discriminant < 0:
        return []
    elif discriminant == 0:
        x = -B / (2 * A)
        y = m * x + k
        return [(x, y)]
    else:
        sqrt_discriminant = math.sqrt(discriminant)
        x1 = (-B + sqrt_discriminant) / (2 * A)
        x2 = (-B - sqrt_discriminant) / (2 * A)
        y1 = m * x1 + k
        y2 = m * x2 + k
        return [(x1, y1), (x2, y2)]

def sample_watt_spectrum(a=0.988e6, b=2.249e-6, E_max=15e6):
    E_peak = a / 2
    p_max = np.exp(-E_peak/a) * np.sinh(np.sqrt(b*E_peak))
    while True:
        E = np.random.uniform(0, E_max)
        p = np.exp(-E/a)*np.sinh(np.sqrt(b*E))
        if np.random.uniform(0, p_max) < p:
            return E

def calculate_scattering_angle():
    # Changed from -π/4 to π/4 to -π/2 to π/2 for more dramatic scattering
    return np.random.uniform(-np.pi/2, np.pi/2)

def update_direction(direction, angle_change):
    angle = np.arctan2(direction[1], direction[0]) + angle_change
    new_direction = [np.cos(angle), np.sin(angle)]
    return new_direction

def neutron_decision_maker(abs_cs, scat_cs, fission_cs):
    total_energy = abs_cs + scat_cs + fission_cs
    random_float = np.random.uniform(0, total_energy)
    if random_float <= abs_cs:
        return "absorb"
    elif random_float <= abs_cs + scat_cs:
        return "scatter"
    else:
        return "fission"

def get_cross_section_data(energy, data):
    fission_data = load_cross_section_data_from_json(data, 'MT_4')
    absorption_data = load_cross_section_data_from_json(data, 'MT_2')
    scattering_data = load_cross_section_data_from_json(data, 'MT_102', cs_type='product_distributions')
    energies_f, xs_f = get_energy_xs_arrays(fission_data)
    energies_a, xs_a = get_energy_xs_arrays(absorption_data)
    energies_s, xs_s = get_energy_xs_arrays(scattering_data)
    fission_xs = linear_interpolate_xs(energies_f, xs_f, energy)
    absorption_xs = linear_interpolate_xs(energies_a, xs_a, energy)
    scattering_xs = linear_interpolate_xs(energies_s, xs_s, energy)
    return fission_xs, absorption_xs, scattering_xs

def neutron_enters_fuel(prev_pos, new_pos, fuel_center, fuel_radius):
    # Check if either endpoint is in fuel
    if np.linalg.norm(np.array(prev_pos) - np.array(fuel_center)) <= fuel_radius:
        return True
    if np.linalg.norm(np.array(new_pos) - np.array(fuel_center)) <= fuel_radius:
        return True
    # Check for intersection with the fuel circle
    direction = np.array(new_pos) - np.array(prev_pos)
    if np.linalg.norm(direction) == 0:
        return False
    direction = direction / np.linalg.norm(direction)
    intersections = check_neutron_intersection(
        fuel_center[0], fuel_center[1], fuel_radius, direction
    )
    for pt in intersections:
        pt = np.array(pt)
        # Project intersection onto the segment
        t = np.dot(pt - np.array(prev_pos), direction)
        if 0 < t < np.linalg.norm(np.array(new_pos) - np.array(prev_pos)):
            return True
    return False

def distance_to_next_interaction(cross_section_fisssion, cross_section_absorption, cross_section_scattering):
    total_cross_section = cross_section_fisssion + cross_section_absorption + cross_section_scattering
    distance = -(1 / total_cross_section) * np.log(np.random.uniform(0, 1))
    distance = min(distance,10)  # Clamp to 10 cm max for visualization
    return distance

def calculate_distance_to_boundary(position, direction, radius):
    """Calculate distance to moderator/reactor boundary"""
    x, y = position
    dx, dy = direction
    a = dx**2 + dy**2
    b = 2 * (x*dx + y*dy)
    c = x**2 + y**2 - radius**2
    
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        return float('inf')
    
    d1 = (-b + np.sqrt(discriminant)) / (2*a)
    d2 = (-b - np.sqrt(discriminant)) / (2*a)
    
    if d1 >= 0 and d2 >= 0:
        return min(d1, d2)
    elif d1 >= 0:
        return d1
    elif d2 >= 0:
        return d2
    return float('inf')


def moderator_Water_logic(energy, position, direction):
    n = avogadro * WATER_DESNITY / 18.01528  # molecules/cm^3
    sigma_0 = 20e-24  # cm^2 (back to realistic value)
    E_0 = 0.025  # eV (thermal)
    
    # Modified cross section calculation
    sigma = sigma_0 * np.sqrt(E_0 / energy)  # More realistic energy dependence
    sigma = max(sigma, 1e-24)  # Minimum cross section
    lamda = 1 / (n * sigma)
    
    # Sample distance with minimum step size
    distance = np.random.exponential(lamda)
    distance = max(min(distance, 20), 0.1)  # Between 0.1 and 20 cm
    
    dist_to_boundary = calculate_distance_to_boundary(position, direction, WATER_RADIUS)
    
    if distance > dist_to_boundary:
        return dist_to_boundary, True
    return distance, False

def check_material_neutron(position):
    x, y = position
    r = np.sqrt(x**2 + y**2)
    if r <= FUEL_RADIUS:
        return 1
    elif r <= WATER_RADIUS:
        return 2
    else:
        return 0

# --- Simulation --- #

Energy_neutron = sample_watt_spectrum()
print(f"Sampled neutron energy: {Energy_neutron:.2e} eV")

# Position neutron just outside the fuel rod
position_neutron = [FUEL_RADIUS + 10, 280]  # much closer

# Initial direction: point toward the fuel center
direction_vector = np.array(FUEL_POSITION) - np.array(position_neutron)
direction_neutron = (direction_vector / np.linalg.norm(direction_vector)).tolist()

print(f"Initial neutron direction: {direction_neutron}")

energies_history = [Energy_neutron]
positions_history = [position_neutron.copy()]
steps = 50

while steps != 0:
    material = check_material_neutron(position_neutron)
    steps -= 1
    print(f"[steps]---------------------------{steps},{len(positions_history)}")
    fission_cs, abs_cs, scat_cs = get_cross_section_data(Energy_neutron, 'neutron_data.json')
    print("Cross sections (cm^2):", fission_cs, abs_cs, scat_cs)
    
    if material == 1:
        print("\n\nNeutron in FUEL ROD [FUEL ROD]\n\n")
        break
    if material == 2:
        print("Neutron in WATER [MODERATOR]")
        dist_to_interac, is_leakage = moderator_Water_logic(Energy_neutron, position_neutron, direction_neutron)
        print(f"Distance to next interaction: {dist_to_interac:.2f} cm")

        # Move neutron
        position_neutron[0] += dist_to_interac * direction_neutron[0]
        position_neutron[1] += dist_to_interac * direction_neutron[1]
        positions_history.append(position_neutron.copy())
        
        if is_leakage:
            print("NEUTRON LEAKED FROM MODERATOR")
            break
            
        # Only update direction and energy if not leaked
        angle_change = calculate_scattering_angle()
        direction_neutron = update_direction(direction_neutron, angle_change)
        
        # Update energy
        A = 1
        costheta = np.cos(angle_change)
        Eprime = Energy_neutron * ((A**2 + 1 + 2*A*costheta) / (A + 1)**2)
        Energy_neutron = Eprime
        energies_history.append(Energy_neutron)
        print(f"New neutron energy: {Energy_neutron:.2e} eV")


    else:
        print("NEUTRON EXITED THE REACTOR")
        position_neutron = [FUEL_RADIUS + 0.1, 0]  # much closer
        continue



# --- Plotting energy history ---
plt.figure(figsize=(8, 5))
plt.plot(range(len(energies_history)), energies_history, marker='o')
plt.yscale('log')
plt.xlabel('Moderation Step')
plt.ylabel('Neutron Energy (eV)')
plt.title('Neutron Energy Loss During Moderation')
plt.grid(True, which='both', ls='--')
plt.show()

# --- Plotting position history ---
positions_history = np.array(positions_history)

plt.figure(figsize=(6, 6))
plt.plot(positions_history[:, 0], positions_history[:, 1], marker='o', color='black', label='Neutron Path', linewidth=1.5, markersize=4)

# Draw fuel and moderator boundaries
fuel_circle = Circle(FUEL_POSITION, FUEL_RADIUS, color='orange', fill=False, linewidth=2, label='Fuel')
water_circle = Circle(WATER_POSITION, WATER_RADIUS, color='blue', fill=False, linewidth=2, label='Moderator')
plt.gca().add_patch(fuel_circle)
plt.gca().add_patch(water_circle)

plt.xlabel('x position (cm)')
plt.ylabel('y position (cm)')
plt.title('Neutron Path in Reactor')
plt.axis('equal')
plt.legend(loc='upper right')
plt.grid(True, which='both', ls='--', alpha=0.5)
plt.tight_layout()
plt.show()