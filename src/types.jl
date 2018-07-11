# types.jl
# first stab at writing a tight-binding type to hold the model + data

struct TightBindingModel
    N::Int
    
    H::Array
    S::AbstractArray  #array view of diagonal
    J::AbstractArray  #array view of +1 diagonal
    psi::Array
    
    dipoles::Array

    function TightBindingModel(N=20) # Inner type constructor
        H=zeros(N,N)
        S=view(H,1:N+1:N*N) # Array view of diagonal 
        J=view(H,2:N+1:N*N) # Array view of (N-1) off diagonal # Not currently working

        psi=zeros(N)
        dipoles=zeros(N)

        new(N,H,S,J,psi,dipoles)
    end
end


