include("part2_calc.jl")
using Plots, Printf

# contour plot of <S(qx,qy)> at fixed (U,T)
function plot_Sq(U, Tval; nmeas=1500)
    _, qs, Sq = run_MC(U, Tval; nmeas=nmeas)
    Smax, imax = findmax(Sq); a,b = Tuple(imax)
    @printf("U=%.1f T=%.2f : peak S(q)=%.4f at q=(%.3f,%.3f)\n", U, Tval, Smax, qs[a], qs[b])
    contourf(qs, qs, Sq, xlabel="qx", ylabel="qy",
             title="⟨S(q)⟩, U=$U, T=$Tval", color=:viridis)
end

# S(pi,pi) vs T scan at fixed U
function T_scan(U; Ts=0.02:0.04:1.0, nmeas=2000)
    pi_idx = findfirst(q -> isapprox(q, pi; atol=1e-8), [2pi*p/L for p in 0:L-1])
    Svals = Float64[]
    M_prev = nothing
    for Tval in Ts
        M_prev, _, Sq = run_MC(U, Tval; nmeas=nmeas, M0=M_prev)
        push!(Svals, Sq[pi_idx, pi_idx])
        @printf("T=%.3f  S(π,π)=%.4f\n", Tval, Svals[end])
    end
    return collect(Ts), Svals
end

