# Kernels g and the optional shifted-power modulation h.

"""
    GaussKernel()

`g(s) = exp(-s²)`. Allowed by Theorem 5; the simplest kernel for ball arithmetic.
"""
struct GaussKernel end

"""
    ErfcKernel()

`g(s) = erfc(s)`, the kernel used in §2 of the paper.
"""
struct ErfcKernel end

"""
    kernelvalue(g, s)

Value of the kernel `g` at `s`. Extensions specialise this for number types that need a
different implementation (e.g. `erfc` on Arb balls).
"""
kernelvalue(::GaussKernel, s) = exp(-(s * s))
kernelvalue(::ErfcKernel, s) = erfc(s)

"""
    NoModulation()

Identity modulation `h(q) = q`; equivalent to `ShiftedPower(1, 0)`.
"""
struct NoModulation end

"""
    ShiftedPower(β, ρ)

Background-corrected modulation `W(q) = (q + ρ)^β - ρ^β` (§2.3 and (5.2)), with an
integer `β ≥ 1` and a rational `ρ ≥ 0`.
"""
struct ShiftedPower
    β::Int
    ρ::Rational{BigInt}
    function ShiftedPower(β::Integer, ρ::Real = 0)
        β >= 1 || throw(ArgumentError("β must be ≥ 1"))
        r = Rational{BigInt}(ρ)
        r >= 0 || throw(ArgumentError("ρ must be ≥ 0"))
        return new(Int(β), r)
    end
end

modulate(::NoModulation, q) = q
function modulate(m::ShiftedPower, q)
    ρ = lift(q, m.ρ)
    return (q + ρ)^m.β - ρ^m.β
end

"`W(q)` evaluated exactly at a rational `q`."
modulate_exact(::NoModulation, q::Rational) = q
modulate_exact(m::ShiftedPower, q::Rational) = (q + m.ρ)^m.β - m.ρ^m.β
