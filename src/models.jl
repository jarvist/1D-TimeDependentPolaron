# models.jl
# Specification of the 1D models

"""
    nondispersive_wavepacket(x0, λ)

Generates an N-unit 1D non-dispersive wavepacket
Shamelessly copied from Wikipedia:
https://en.wikipedia.org/wiki/Wave_packet

x0 - center of wavepacket
λ  - width of wavepacket
    
returns    
ψ - complex wavefunction
"""
function nondispersive_wavepacket(x0, λ)
    ψ=[ exp(-(x-x0)^2)*(cos(2π*(x-x0)/λ) - im * sin(2π*(x-x0)/λ) ) for x=1:N ]
    return ψ
end

"""
    planewave(λ)

Generates a N-unit 1D planewave with wavelength λ.

λ - wavelength (r*a0)

returns
ψ - wavefunction
"""
function planewave(λ)
    k=2π/λ
    ψ=[ exp(-im* k*x) for x=1:N ]
    return ψ
end

