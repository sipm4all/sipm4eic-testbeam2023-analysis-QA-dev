# recoana: reconstructed data analysis

Reconstructed data are stored in `recodata.root` files with a plain ROOT TTree structure defined as follows

```
unsigned short n;
float x[65534];
float y[65534];
float t[65534];
auto recotree = new TTree("recodata", "recodata");
recotree->Branch("n", &n, "n/s");
recotree->Branch("x", &x, "x[n]/F");
recotree->Branch("y", &y, "y[n]/F");
recotree->Branch("t", &t, "t[n]/F");
```

The meaning of the branches is self-explanatory, but for the sake of being explicit

- `n` is the number of hits in the event
- `x` is the x-coordinate of the hit [mm]
- `y` is the y-coordinate of the hit [mm]
- `t` is the time of the hit [ns]

Examples are given to read the data both with ROOT and python scripts.
