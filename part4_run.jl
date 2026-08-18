include("part3_scans.jl")

function diagnose(U, Tval; nequil=1500, nmeas=4000)
    _, qs, Sq = run_MC(U, Tval; nequil=nequil, nmeas=nmeas)
    pi_idx = findfirst(q -> isapprox(q, pi; atol=1e-8), qs)
    Smax, imax = findmax(Sq)
    a, b = Tuple(imax)

    println("\n===== U=$U T=$Tval =====")
    println("S(π,π)   = ", round(Sq[pi_idx,pi_idx], digits=5))
    println("Global max = ", round(Smax, digits=5), " at q=($(round(qs[a],digits=3)), $(round(qs[b],digits=3)))")

    println("Full grid:")
    for a in 1:L, b in 1:L
        marker = (a==pi_idx && b==pi_idx) ? "  <-- (π,π)" : ""
        println("  q=($(round(qs[a],digits=2)),$(round(qs[b],digits=2)))  S=$(round(Sq[a,b],digits=4))$marker")
    end
    return qs, Sq
end

diagnose(1.0, 0.3)
diagnose(5.0, 0.3)
diagnose(8.0, 0.3)

# contours at U = 1, 5, 8

T0 = 0.3
p1 = plot_Sq(1.0, T0; nmeas=4000)
p2 = plot_Sq(5.0, T0; nmeas=4000)
p3 = plot_Sq(8.0, T0; nmeas=4000)
display(plot(p1, p2, p3, layout=(1,3), size=(1500,450)))