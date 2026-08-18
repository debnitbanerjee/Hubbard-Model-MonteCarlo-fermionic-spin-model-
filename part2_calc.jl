include("part1_lattice.jl")


log1pexp_stable(x) = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))

function F_eff(M, U, Tval)
    evals = eigvals(build_Hsp(M, U))
    beta = 1/Tval
    Fferm = -Tval * sum(log1pexp_stable(-beta*e) for e in evals)
    Fm = (U/4) * sum(mx^2+my^2+mz^2 for (mx,my,mz) in M)
    return Fferm + Fm
end

init_M(mmax=0.3) = [(mmax*(2rand()-1), mmax*(2rand()-1), mmax*(2rand()-1)) for _ in 1:N]

#metropolis 
function mc_sweep!(M, U, Tval, delta, Fold)
    beta = 1/Tval
    for _ in 1:N
        i = rand(1:N)
        mx,my,mz = M[i]
        Mtrial = copy(M)
        Mtrial[i] = (mx+delta*(2rand()-1), my+delta*(2rand()-1), mz+delta*(2rand()-1))
        Fnew = F_eff(Mtrial, U, Tval)
        if Fnew <= Fold || rand() < exp(-beta*(Fnew-Fold))
            M[i] = Mtrial[i]
            Fold = Fnew
        end
    end
    return Fold
end

# structure factor S(qx,qy) for one spin configuration
function structure_factor(M)
    qs = [2pi*p/L for p in 0:L-1]
    Sq = zeros(L,L)
    for (a,qx) in enumerate(qs), (b,qy) in enumerate(qs)
        A = ComplexF64[0,0,0]
        for i in 1:N
            x,y = pos[i]
            ph = exp(im*(qx*x+qy*y))
            A .+= collect(M[i]) .* ph
        end
        Sq[a,b] = sum(abs2, A)/N^2
    end
    return qs, Sq
end

# thermalize + measure <S(q)>, returns final config, q-grid, avg S(q)
function run_MC(U, Tval; nequil=1500, nmeas=4000, delta=0.3, M0=nothing)
    M = M0 === nothing ? init_M() : copy(M0)
    Fold = F_eff(M, U, Tval)
    for _ in 1:nequil
        Fold = mc_sweep!(M, U, Tval, delta, Fold)
    end
    qs = [2pi*p/L for p in 0:L-1]
    Sq_sum = zeros(L,L)
    for _ in 1:nmeas
        Fold = mc_sweep!(M, U, Tval, delta, Fold)
        _, Sq = structure_factor(M)
        Sq_sum .+= Sq
    end
    return M, qs, Sq_sum ./ nmeas
end