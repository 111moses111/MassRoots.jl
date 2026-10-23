# Problems: a real analytic function on a search domain, given through the nonnegative φ
# used in the kernel (φ = f², or φ = f when f ≥ 0).

"""
    AbstractProblem

A search problem for the real zeros of a real analytic function `f`, given through the
nonnegative function `φ` used in the kernel: `φ = f²` in general, or `φ = f` when `f` is
already nonnegative. Every problem implements

* `phi(prob, x)` – the value `φ(x)` for any number type of `x` (real or complex: rigorous
  evaluators such as [`ArbEvaluator`](@ref) evaluate `φ` on complex balls);
* `domain(prob)` – the search interval `(a, b)` as exact rationals;
* `basepoint(prob)` – the base point `x₀` of `H(x) = ∫_{x₀}^{x} q`;
* `scalebits(prob)` – a bound on the bit size of intermediate arguments of `φ`, used to
  choose the working precision (default 0).
"""
abstract type AbstractProblem end

scalebits(::AbstractProblem) = 0
basepoint(prob::AbstractProblem) = domain(prob)[1]

"""
    IntegerZeroProblem <: AbstractProblem

A problem whose zeros on the search domain are integers from `candidates(prob)`, with a
scale `S = scale(prob)` satisfying the hypotheses of the mass-separation theorem
(Theorem 5 of the paper) on `I* = [first(candidates) - 1/2, last(candidates) + 1/2]`:

1. `|I*| ≤ S`;
2. `φ(x) ≥ 1/S` for every `x` in a candidate block that contains no zero;
3. `φ(d + s) ≤ S s^2` for every zero `d` and `|s| ≤ 1/2`.

Subtypes implement `phi`, `candidates`, `scale` and `isroot(prob, m)` (exact check of a
candidate). For your own problem use the constructor
`IntegerZeroProblem(φ, candidates, S; isroot, x0, scalebits)`.
"""
abstract type IntegerZeroProblem <: AbstractProblem end

function domain(prob::IntegerZeroProblem)
    c = candidates(prob)
    return (half(first(c) - 1), half(last(c)))
end

"""
    CustomIntegerProblem(φ, candidates, S; isroot, x0, scalebits = 0)

User-defined [`IntegerZeroProblem`](@ref). The caller is responsible for `S` satisfying
the three hypotheses; the search then inherits the guarantees of the paper.
"""
struct CustomIntegerProblem{F,R} <: IntegerZeroProblem
    φ::F
    cands::UnitRange{BigInt}
    S::BigInt
    isroot::R
    x0::Rational{BigInt}
    scalebits::Int
end

function IntegerZeroProblem(φ, cands::AbstractUnitRange, S::Integer; isroot,
                            x0 = half(first(cands) - 1), scalebits::Integer = 0)
    c = big(first(cands)):big(last(cands))
    return CustomIntegerProblem(φ, c, big(S), isroot, Rational{BigInt}(x0), Int(scalebits))
end

phi(p::CustomIntegerProblem, x) = p.φ(x)
candidates(p::CustomIntegerProblem) = p.cands
scale(p::CustomIntegerProblem) = p.S
isroot(p::CustomIntegerProblem, m::Integer) = p.isroot(m)
basepoint(p::CustomIntegerProblem) = p.x0
scalebits(p::CustomIntegerProblem) = p.scalebits

# --------------------------------------------------------------------------------------
"""
    FactorProblem(t)

Nontrivial divisor search for `t ≥ 4` (§3 and §5.6 of the paper):
`φ(x) = sin²(πt/x) + sin²(πx)`, candidates `2, …, t-1`, scale `S_t = 2t²`, base point 2.
"""
struct FactorProblem <: IntegerZeroProblem
    t::BigInt
    function FactorProblem(t::Integer)
        t >= 4 || throw(ArgumentError("FactorProblem needs t ≥ 4; t = 2, 3 are prime"))
        return new(big(t))
    end
end

function phi(p::FactorProblem, x)
    return sinpi2(lift(x, p.t) / x) + sinpi2(x)
end
candidates(p::FactorProblem) = big(2):(p.t - 1)
scale(p::FactorProblem) = 2 * p.t^2
isroot(p::FactorProblem, m::Integer) = 1 < m < p.t && iszero(rem(p.t, m))
basepoint(::FactorProblem) = big(2) // 1
scalebits(p::FactorProblem) = bitlength(p.t)

# --------------------------------------------------------------------------------------
"""
    QuadraticProblem(a, b, c)

Binary quadratic feasibility `a m² + b k = c` with `m, k ≥ 0` (§4 and §5.7):
`φ(x) = sin²(πx) + sin²(π (c - a x²)/b)`, candidates `0, …, N` with
`N = ⌊√⌊c/a⌋⌋`, scale `S = 10 D`, `D = b² + a²(2N+1)²`, base point 0.
"""
struct QuadraticProblem <: IntegerZeroProblem
    a::BigInt
    b::BigInt
    c::BigInt
    N::BigInt
    function QuadraticProblem(a::Integer, b::Integer, c::Integer)
        (a > 0 && b > 0 && c > 0) ||
            throw(ArgumentError("QuadraticProblem needs positive a, b, c"))
        A, B, C = big(a), big(b), big(c)
        return new(A, B, C, isqrt(div(C, A)))
    end
end

function phi(p::QuadraticProblem, x)
    y = (lift(x, p.c) - lift(x, p.a) * x * x) / lift(x, p.b)
    return sinpi2(x) + sinpi2(y)
end
candidates(p::QuadraticProblem) = big(0):p.N
scale(p::QuadraticProblem) = 10 * (p.b^2 + p.a^2 * (2p.N + 1)^2)
function isroot(p::QuadraticProblem, m::Integer)
    0 <= m <= p.N || return false
    r = p.c - p.a * big(m)^2
    return r >= 0 && iszero(rem(r, p.b))
end
basepoint(::QuadraticProblem) = big(0) // 1
scalebits(p::QuadraticProblem) = bitlength(p.c) + bitlength(p.a)

# --------------------------------------------------------------------------------------
"""
    GeneralProblem(f, a, b; square = true, x0 = a, scalebits = 0)

Any real-valued analytic function `f` of a real variable on `[a, b]`. By default the
problem uses `φ = f²`, which works whether or not `f` changes sign; pass `square = false`
to use `φ = f` when `f` is known to be nonnegative.
Zeros need not be integers. No mass gap is certified: [`localize`](@ref) takes the scale
`α` and the threshold from the caller.
"""
struct GeneralProblem{F} <: AbstractProblem
    f::F
    a::Rational{BigInt}
    b::Rational{BigInt}
    square::Bool
    x0::Rational{BigInt}
    scalebits::Int
end

function GeneralProblem(f, a::Real, b::Real; square::Bool = true, x0 = a,
                        scalebits::Integer = 0)
    A, B = rationalize_exact(a), rationalize_exact(b)
    A < B || throw(ArgumentError("need a < b"))
    return GeneralProblem(f, A, B, square, rationalize_exact(x0), Int(scalebits))
end

rationalize_exact(x::Integer) = big(x) // 1
rationalize_exact(x::Rational) = Rational{BigInt}(x)
rationalize_exact(x::AbstractFloat) = Rational{BigInt}(x)   # exact binary value

function phi(p::GeneralProblem, x)
    v = p.f(x)
    return p.square ? v * v : v
end
domain(p::GeneralProblem) = (p.a, p.b)
basepoint(p::GeneralProblem) = p.x0
scalebits(p::GeneralProblem) = p.scalebits
