# types.jl
# first stab at writing a tight-binding type to hold the model + data

struct TightBindingModel
    N::Int
    
    H::Array
    S::Array

    psi::Array
    
    dipoles::Array
end

function TightBindingModel(N=20) # Inner type constructor
    H=zeros(N,N)
    S=view(H,1:N:N*N) # array view; currently this gets obliterated into a copy at some point?

    psi=zeros(N)
    dipoles=zeros(N)

    return(TightBindingModel(N,H,S,psi,dipoles))
end


