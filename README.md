# 1D-TimeDependentPolaron

Polarons + Time-Dependent-Propagation in 1D. Some toy Julia codes.

Currently quick 1D code that treated the polarisation response as a load of dipoles (an N-long array of reals) centred on the sites. Time becomes the discretisation of the SCF cycle, and these dipoles have memory, and respond to the 'density' of the wavefunction. Currently the wavefunction is just the adiabatic solution from the Hamiltonian (i.e. eigvec(H)[:,1]^2).

A plot below shows the S-curve of the dipoles, generating the dimpled site energy surface (red) with the current electron density sitting in the middle of it.

![Screenshot](screenshot.png)

And you can watch some wavepacket propagation videos with lattice response here:

[![Lattice response to wavepacket](https://img.youtube.com/vi/3U_zD0kL460/0.jpg)](https://www.youtube.com/watch?v=3U_zD0kL460)

## Plan
- [x] 1D dipoles for polarisation of lattice
- [x] dipoles respond to adiabatic ground states (i.e. solve TISE for H)
- [x] evolve simulation in time; (adiabatic electronic structure; step dipoles
  to respond to lattice - 'AdiabaticPropagation')
- [x] Create wavepackets + plane waves
- [ ] figure out the realistic parameters of all the values set to 1 (starting with the dt=1 !)
- [x] add time-dependent Schr. equation for propogation of the wavefunction of interest ('UnitaryPropagation')
- [ ] ? surface hopping by overlap of this wavefn. with the adiabatic solution
