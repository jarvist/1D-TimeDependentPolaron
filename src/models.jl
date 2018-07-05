"""
Rose Teague & Jarvist Frost  04/07/2018
Script to set up the 1D models
Define functions:

- nondispersive_wavepacket
- planewave

"""


#-------------------------------------------------------------------------------
function nondispersive_wavepacket(x0, λ)
    #=
    Generates an N-unit 1D non-dispersive wavepacket
    Shamelessly copied from Wikipedia:
    https://en.wikipedia.org/wiki/Wave_packet
    ------------------
    Inputs
    ------------------
    x0 - center of wavepacket
    λ  - width of wavepacket
    ------------------
    Output
    ------------------
    ψ - wavepacket
    =#

    ψ=[ exp(-(x-x0)^2)*(cos(2π*(x-x0)/λ) - im * sin(2π*(x-x0)/λ) ) for x=1:N ]
    return ψ
end


#-------------------------------------------------------------------------------
function planewave(λ)
    #=
    Generates a N-unit 1D planewave with wavelength λ.
    ----------------
    Inputs
    ----------------
    λ - wavelength (r*a0)
    ----------------
    Outputs
    ----------------
    ψ - wavefunction
    =#

    k=2π/λ
    ψ=[ exp(-im* k*x) for x=1:N ]
    return ψ
end
