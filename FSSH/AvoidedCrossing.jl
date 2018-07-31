A = 0.01
B = 1.6
C = 0.005
D = 1.0

#Diabatic interactions
function diabatic_potential(x)
    V = zeros(2,2)
    V[1] = (-sign(x))*A*(1-exp(-B*x*sign(x)))
    V[2] = C*exp(-D*x^2)
    V[4] = -V[1]
    V[3] = V[2]

    a = 1; b = -(V[1]+V[4]); c = det(V)
    lambdas = ( (-b + sqrt(disc(a,b,c)))/(2a), (-b - sqrt(disc(a,b,c)))/(2a) )
    return lambdas
end

function disc(a,b,c)
    return b^2 - 4*a*c
end
