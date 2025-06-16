import matplotlib.pyplot as plt

# === Heterostructure ===
VBM = -0.9018          # Valence Band Maximum (in eV)
CBM = 1.3176           # Conduction Band Minimum (in eV)
E_vacuum = 5.502       # Vacuum level (in eV)

# === NTO ===
# VBM = -4.6139          # Valence Band Maximum (in eV)
# CBM = -1.6730           # Conduction Band Minimum (in eV)
# E_vacuum = 1.793       # Vacuum level (in eV)

# === MSN ===
# VBM = 0.6732          # Valence Band Maximum (in eV)
# CBM = 3.2765           # Conduction Band Minimum (in eV)
# E_vacuum = 6.860       # Vacuum level (in eV)

# === Shifted to Vacuum Level ===
VBM_vac = VBM - E_vacuum
CBM_vac = CBM - E_vacuum

# Water redox levels (in eV vs vacuum)
water_reduction = -4.44  # H⁺/H₂
water_oxidation = -5.67  # O₂/H₂O

# === Plotting ===
plt.figure(figsize=(5, 6))
plt.bar(1, VBM_vac - CBM_vac, bottom=CBM_vac, width=0.4, color='skyblue', edgecolor='k', label='Band Gap')

# Redox and vacuum level lines
plt.axhline(y=water_reduction, color='red', linestyle='--', label='Water Reduction (-4.44 eV)')
plt.axhline(y=water_oxidation, color='green', linestyle='--', label='Water Oxidation (-5.67 eV)')
plt.axhline(y=0, color='gray', linestyle=':')

# Formatting
plt.xlim(0.5, 1.5)
plt.ylim(-2, -7)
plt.ylabel('Energy (eV vs Vacuum)', fontsize=12)
plt.xticks([])
plt.xlabel('')
plt.gca().invert_yaxis()  # Vacuum at top
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('band_alignment_HSE.png', dpi=300)