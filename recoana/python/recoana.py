import sys
import uproot

# Check if the filename is provided as a command-line argument
if len(sys.argv) != 2:
    print("Usage: python script.py <filename.root>")
    sys.exit(1)

# Get the input filename from the command-line argument
filename = sys.argv[1]

# Open the ROOT file
fin = uproot.open(filename)

# Get the TTree
tin = fin['recodata']

# Get the number of entries
nev = tin.num_entries

# Define the arrays to store the data
n = tin['n'].array()
x = tin['x'].array()
y = tin['y'].array()
t = tin['t'].array()

# Loop over the entries
for i in range(nev):
    # Access data for the i-th entry
    entry_n = n[i]
    entry_x = x[i]
    entry_y = y[i]
    entry_t = t[i]

    # Now you can use entry_n, entry_x, entry_y, entry_t in your analysis

    # Example: Print the values for each entry
    print(f"Entry {i + 1}: n={entry_n}, x={entry_x}, y={entry_y}, t={entry_t}")
