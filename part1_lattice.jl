using LinearAlgebra, Random 
using CSV
using DataFrames
Random.seed!(1)

L = 6
N = L^2
t = 1.0

idx(x,y) = mod(x,L)*L + mod(y,L) + 1
pos = [(x,y) for x in 0:L-1 for y in 0:L-1]
bonds = [(idx(x,y), idx(x+1,y)) for x in 0:L-1 for y in 0:L-1]
append!(bonds, [(idx(x,y), idx(x,y+1)) for x in 0:L-1 for y in 0:L-1])

# H_sp[M] = kinetic - (U/2) * m_i . sigma_i   (Eq. A6, clean, half-filling)
function build_Hsp(M, U)
    H = zeros(ComplexF64, 2N, 2N)
    for i in 1:N
        mx,my,mz = M[i]
        Vi = -(U/2) .* [ mz          (mx-im*my) ;
                         (mx+im*my)  -mz        ]
        H[2i-1:2i, 2i-1:2i] .= Vi
    end
    for (i,j) in bonds
        H[2i-1:2i, 2j-1:2j] .= -t*I(2)
        H[2j-1:2j, 2i-1:2i] .= -t*I(2)
    end
    return Hermitian(H)
end

function set_L(Lval)
    global L, N, pos, bonds
    L = Lval; N = L^2
    pos = [(x,y) for x in 0:L-1 for y in 0:L-1]
    bonds = [(idx(x,y), idx(x+1,y)) for x in 0:L-1 for y in 0:L-1]
    append!(bonds, [(idx(x,y), idx(x,y+1)) for x in 0:L-1 for y in 0:L-1])
end