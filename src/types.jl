# types.jl
# first stab at writing a tight-binding type to hold the model + data

struct TightBindingModel
    N::Int
    
    H::Array
    S::AbstractArray  #array view of diagonal
    J::AbstractArray  #array view of +1 diagonal

    psi::Array
    density::Array
    
    dipoles::Array
    dipole_timeconstant
    dipole_strength

    function TightBindingModel(N=20, dipole_timeconstant=0.2, dipole_strength=0.2) # Inner type constructor
        H=zeros(N,N)
        S=view(H,1:N+1:N*N) # Array view of diagonal 
        J=view(H,2:N+1:N*N) # Array view of (N-1) off diagonal # Not currently working

        psi=zeros(N)
        density=zeros(N)
 
        dipoles=zeros(N)

        new(N,H,S,J,psi,density,dipoles,dipole_timeconstant,dipole_strength)
    end
end


