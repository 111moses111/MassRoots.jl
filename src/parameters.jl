# Certified parameters of the integral-mass test (Theorem 5 and §5.8 of the paper).

"""
    MassTestParameters

Parameters of a certified mass test, all exact:

* `S` – geometric scale (enlarged to `Q` for a modulated kernel), `ν = ⌈√S⌉`;
* `α` – kernel scale `S²`;
* `w` – radius of the neighbourhood of a zero used in the lower bound;
* `background` – upper bound `B₀` on the mass of a zero-free candidate block;
* `minmass` – lower bound `B₁` on the mass of a block that contains a zero;
* `threshold` – decision threshold `θ` between them;
* `p`, `atol = 2^-p` – absolute accuracy required for each value of `H`.

With every value of `H` known to within `atol`, the test `H(v) - H(u) > θ` is correct on
every candidate block.
"""
struct MassTestParameters
    S::BigInt
    ν::BigInt
    α::BigInt
    w::Rational{BigInt}
    background::Rational{BigInt}
    minmass::Rational{BigInt}
    threshold::Rational{BigInt}
    p::Int
    atol::Rational{BigInt}
end

"""
    mass_test_parameters(S, modulation = NoModulation())
    mass_test_parameters(prob::IntegerZeroProblem; modulation = NoModulation())

Parameters of Theorem 5 (`NoModulation`) or of §5.8 (`ShiftedPower`). A scale below 8 is
raised to 8, which keeps all three hypotheses valid.
"""
function mass_test_parameters(S::Integer, modulation = NoModulation())
    S = max(big(S), big(8))
    U = modulate_exact(modulation, big(1) // 1)
    V = modulate_exact(modulation, big(1) // 2)
    Q = max(S, big(cld(numerator(4U / V), denominator(4U / V))))
    ν = ceilsqrt(Q)
    α = Q^2
    w = 1 // (2Q * ν)                        # radius of the neighbourhood of a zero
    B1 = V / (Q * ν)                         # = V·2w, mass lower bound
    B0 = V / (2Q * ν)                        # background upper bound
    θ = 3V / (4Q * ν)
    p = max(0, ceillog2(16Q * ν / V))
    return MassTestParameters(Q, ν, α, w, B0, B1, θ, p, 1 // big(2)^p)
end

mass_test_parameters(prob::IntegerZeroProblem; modulation = NoModulation()) =
    mass_test_parameters(scale(prob), modulation)

function Base.show(io::IO, P::MassTestParameters)
    print(io, "MassTestParameters(S = ", P.S, ", α = S² (", bitlength(P.α), " bits), p = ",
          P.p, ", θ ≈ ", round(Float64(P.threshold); sigdigits = 4), ")")
end
