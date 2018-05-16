## Functions provided 

### Models

```@docs
modelJ(theta) 
```

```@docs
randH(SiteEnergy, Edisorder, Jdisorder, modelJ, N)
```

```@docs
nondispersive_wavepacket(x0, λ)
```

```@docs
planewave(λ)
```

```@docs 
prepare_model()
```

### Propagators

```@docs
SiteEnergyFromDipoles(dipoles)
```

```@docs
DipolesFromDensity(dipoles,density,dampening)
```

```@docs
AdiabaticPropagation(dipoles,E,dampening)
```

```@docs
UnitaryPropagation(dipoles,E,psi,dt,dampening;slices::Int=1)
```

```@docs
TimeDependentPropagation(psi,H,dt;slices::Int=1,decompose::Bool=false,verbose::Bool=false) # propagate directly using full Hamiltonian=T+V
```

```@docs
TimeDependentPropagation(psi,H,dt,E) # propagate using eigenvalue
```

### Simulations

```@docs
Plot_S_psi_density_dipoles(S,psi,density,dipoles;title="",verbose::Bool=false)
```

```@docs
overlap(psia,psib)
```

```@docs
Decompose_H(H)
```

```@docs
Plot_H(H)
```

```@docs
outputpng()
```

```@docs
SCFthenUnitary(dampening, SCFcycles, Unitarycycles; PNG::Bool=false)
```

