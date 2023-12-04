import sys
import uproot
import matplotlib.pyplot as plt
import numpy as np
import awkward as ak

# Check if the filename is provided as a command-line argument
if len(sys.argv) != 2:
    print("Usage: python script.py <filename.root>")
    sys.exit(1)

# Get the input filename from the command-line argument
filename = sys.argv[1]

# Open the ROOT file
fin = uproot.open(filename)

# Get the TTree
tin = fin["recodata"]

# Get the number of entries
nev = tin.num_entries

# Define the arrays to store the data
n = tin["n"].array()
x_ak = tin["x"].array()
y_ak = tin["y"].array()

# Flatten the nested arrays
x_ak_flat = ak.flatten(x_ak)
y_ak_flat = ak.flatten(y_ak)

# Convert flattened arrays to numpy arrays
x_np = ak.to_numpy(x_ak_flat)
y_np = ak.to_numpy(y_ak_flat)

# Apply random smearing using vectorized operations
x_smeared = np.random.uniform(x_np - 1.5, x_np + 1.5)
y_smeared = np.random.uniform(y_np - 1.5, y_np + 1.5)

# Create a 2D histogram
nbins = 396
xedges = np.linspace(-99, 99, nbins + 1)
yedges = np.linspace(-99, 99, nbins + 1)
hXY, _, _ = np.histogram2d(y_smeared, x_smeared, bins=[xedges, yedges])

# Plot the filled histogram
plt.imshow(hXY, extent=[-99, 99, -99, 99], origin='lower', cmap='viridis')
plt.colorbar(label='Counts')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.title('2D Histogram with Random Smearing')
plt.show()
