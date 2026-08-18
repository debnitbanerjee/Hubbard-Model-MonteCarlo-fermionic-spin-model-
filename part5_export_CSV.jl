include("part2_calc.jl")   
using CSV, DataFrames

outdir = raw"PASTE YOUR PATH"
isdir(outdir) || mkpath(outdir)


function Spi(qs, Sq)
    pi_idx = findfirst(q -> isapprox(q, pi; atol=1e-8), qs)
    return Sq[pi_idx, pi_idx]
end

function save_ST(U; Ts=0.02:0.02:1.0, nequil=1000, nmeas=3000)
    Svals = Float64[]
    M_prev = nothing
    for Tval in Ts
        M_prev, qs, Sq = run_MC(U, Tval; nequil=nequil, nmeas=nmeas, M0=M_prev)
        push!(Svals, Spi(qs, Sq))
    end
    df = DataFrame(T=collect(Ts), S_pipi=Svals)
    fname = joinpath(outdir, "S_vs_T_U$(U).csv")
    CSV.write(fname, df)
end

function Tc_of(Ts, Svals)
    d2 = diff(diff(Svals))
    for k in 1:length(d2)-1
        sign(d2[k]) != sign(d2[k+1]) && return Ts[k+2]
    end
    return Ts[argmin(diff(Svals)./diff(Ts))+1]
end

function save_TcU(Us; Ts=0.02:0.02:1.0, nequil=1000, nmeas=3000)
    Tcs = Float64[]
    for U in Us
        Svals = Float64[]
        M_prev = nothing
        for Tval in Ts
            M_prev, qs, Sq = run_MC(U, Tval; nequil=nequil, nmeas=nmeas, M0=M_prev)
            push!(Svals, Spi(qs, Sq))
        end
        push!(Tcs, Tc_of(collect(Ts), Svals))
    end
    df = DataFrame(U=Us, Tc=Tcs)
    fname = joinpath(outdir, "Tc_vs_U.csv")
    CSV.write(fname, df)
end


set_L(6)

for U in [1.0, 5.0, 8.0]
    save_ST(U)
end

save_TcU([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0])
