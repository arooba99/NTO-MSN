import numpy as np

# Constants
epsilon_0 = 8.85e-12  # Vacuum permittivity in F/m

def calculate_interfacial_field(delta_q, epsilon_r, area):
    """
    Calculate the interfacial electric field.
    
    Parameters:
        delta_q: Charge transfer at the interface in Coulombs
        epsilon_r: Relative permittivity (dielectric constant)
        area: Interface area in square meters (m^2)
    
    Returns:
        Interfacial Electric field in V/m
    """
    E_int = delta_q / (epsilon_0 * epsilon_r * area)
    return E_int


def calculate_band_offset_fields(delta_CBO, delta_VBO, distance):
    """
    Calculate electric fields at band offsets.
    
    Parameters:
        delta_CBO: Conduction band offset
        delta_VBO: Valence band offset
        distance: Interfacial distance
    
    Returns:
        E_CBO, E_VBO in V/m
    """
    E_CBO = (delta_CBO) / distance
    E_VBO = (delta_VBO) / distance
    return E_CBO, E_VBO


# ======== Example Usage ========
if __name__ == "__main__":
    delta_q = 0.23 * 1.6e-19        # Charge transfer (Coulombs)
    epsilon_r = 1                   # Relative dielectric constant
    area = 2.56548328e-19              # Interface area in m²

    delta_CBO = 0.46          # GGA+U Conduction band offset
    delta_VBO = 0.53          # GGA+U Valence band offset
    # delta_CBO = 2.70          # HSE+U Conduction band offset
    # delta_VBO = 0.88          # HSE+U Valence band offset
    distance = 3.3e-10        # Interfacial distance in meters (3Å)

    E_int = calculate_interfacial_field(delta_q, epsilon_r, area)
    E_CBO, E_VBO = calculate_band_offset_fields(delta_CBO, delta_VBO, distance)

    print(f"Interfacial Electric Field (E_int): {E_int:.2e} V/m")
    print(f"CBO Electric Field (E_CBO): {E_CBO:.2e} V/m")
    print(f"VBO Electric Field (E_VBO): {E_VBO:.2e} V/m")
