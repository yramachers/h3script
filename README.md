# h3script
ROOT calculator script of Tritium beta decay spectrum

Simple Kurie plot script for Tritium; includes energy resolution smearing and neutrino mixing effects with mixings and mass splittings hard-coded, see betashape function in script.

- Usage: start ROOT
```
> .L h3script.C
> plotcurve(0.0)
for a zero lightest neutrino mass eigenstate or
> plotcurve(1.0e-5))
for a 10 meV lightest neutrino mass eigenstate
```

- Run plotcurve(neutrino mass in keV) function at the end
- edit hard-coded energy resolution value in plotcurve(mn)
- for plotting pseudo-data, set exposure by hand in the routine "plotCurve"
- Output: - integral spectrum fraction of range of interest
  and plot of derivative of Kurie curve to see shape changes clearly.
  
- Run tritiumKuriePlotter.C for overlaying various Kurie plots and/or pseudo-data on a single plot.
