# MassRoots.jl

Locating real zeros with integral mass tests, following
*M. Plistiev, "On Some Connections Between Computation of Integrals of Analytic Functions,
Root Finding, Integer Factorization and NP-Completeness"* (2026, submitted).

A nonnegative function `φ` becomes a kernel `q_α(x) = g(α φ(x))` with
`g(s) = exp(-s²)` or `g(s) = erfc(s)`. The kernel equals 1 at the zeros of `φ` and is
small elsewhere. The difference `H(v) - H(u)` of its primitive
`H(x) = ∫_{x₀}^{x} q_α` is the mass of `[u, v]`; a threshold on the mass tells whether
the interval contains a zero, and bisection locates it.

For problems whose zeros are integers and which have a certified scale `S`
(Theorem 5 of the paper), the choice `α = S²` and an absolute accuracy of
`p = ⌈log₂(32 S ⌈√S⌉)⌉` bits for each value of `H` make every test correct, and a
zero is found with at most `2 + ⌈log₂ M⌉` values of `H` for `M` candidates
(Theorem 6). Divisor search and the NP-complete quadratic problem `a m² + b k = c` are
built in.

**The package is built around one idea: the way `H` is computed is interchangeable.**
Problems, kernels and the search logic are fixed; an *evaluator* decides how `H(x)` is
computed. Every value of `H` goes through an `Oracle` that counts, times and logs it, so
different approaches can be compared on equal terms.

## Installation

```julia
pkg> add https://github.com/111moses111/MassRoots.jl
pkg> add Arblib        # optional: the rigorous evaluator
```

Julia 1.11 or newer.

## Quick start

```julia
using MassRoots, Arblib

r = find_divisor(1001, ArbEvaluator())     # LocateResult(:found, root = 7, queries = 12)
r.oracle.log                               # every value of H: point, value, error bound, time

factorize(588, ArbEvaluator())             # ([2, 2, 3, 7, 7], 44)

Q = QuadraticProblem(3, 7, 103)            # 3m² + 7k = 103
decide(Q, ArbEvaluator())                  # (true, …) – exactly two values of H
locate(Q, ArbEvaluator()).root             # 2

# H itself
H = HFunction(FactorProblem(91))           # certified α, Gaussian kernel
atol = mass_test_parameters(FactorProblem(91)).atol
evaluate(ArbEvaluator(), H, 89//2, atol)   # HValue: value, error bound, diagnostics
```

More in [`examples/`](examples): `quickstart.jl`, `custom_evaluator.jl`,
`compare_evaluators.jl`, `general_function.jl` (run with
`julia --project=examples examples/<file>.jl`).

## Building blocks

| Part | Types |
|---|---|
| Problems | `FactorProblem(t)`, `QuadraticProblem(a, b, c)`, `IntegerZeroProblem(φ, candidates, S; isroot)` for your own problem with integer zeros, `GeneralProblem(f, a, b; square)` for arbitrary zeros |
| Kernels | `GaussKernel()` (`exp(-s²)`), `ErfcKernel()`; optional modulation `ShiftedPower(β, ρ)` (§2.3, §5.8) |
| Parameters | `mass_test_parameters(problem)`: `S`, `α = S²`, `w`, threshold `θ = 3w/4`, `p`, `atol = 2⁻ᵖ`, all exact |
| Evaluators | see below; any subtype of `AbstractEvaluator` |
| Search | `locate`, `decide`, `find_divisor`, `factorize`, `localize` |

### Evaluators

| Evaluator | Number types | Rigorous | Notes |
|---|---|---|---|
| `QuadGKEvaluator(T)` | `Float64`, `BigFloat`, `Double64`, … | no | adaptive Gauss–Kronrod on a uniform pre-subdivision |
| `GaussLegendreEvaluator(T)` | same | no | fixed composite rule; cost known in advance |
| `ArbEvaluator()` | Arb balls (needs `using Arblib`) | yes | `acb_calc_integrate`; returns an enclosure of `H(x)` and tightens itself until the enclosure is within `atol` |
| `IdealEvaluator()` | — | — | **test oracle that knows the zeros**; checks the search logic and query counts, never a solver |

For `BigFloat` and Arb the working precision is chosen from the requested accuracy, the
size of `α`, the length of the interval and the size of the arguments of `φ`
(`MassRoots.required_precision`), so the input size is not limited. Override it with
the `precision` keyword.

The default pre-subdivision of the quadrature evaluators for certified problems has step
`2w = 1/(S⌈√S⌉)`. By the proof of Theorem 5 the kernel is at least 1/2 within `w` of
every zero, so no peak can fall between the nodes; the price is about `t⁴` segments for
divisor search. Without a pre-subdivision (`segment_length = nothing`) an adaptive rule
can step over a narrow peak and report no zero.

### Writing your own evaluator

```julia
struct MyEvaluator <: AbstractEvaluator end

function MassRoots.evaluate(ev::MyEvaluator, H::HFunction, x::Rational{BigInt},
                            atol::Rational{BigInt})
    # compute H(x) = ∫_{H.x0}^{x} integrand(H, ξ) dξ to absolute accuracy atol
    return HValue(value, error_bound; rigorous = false, info = (evals = n,))
end

find_divisor(1001, MyEvaluator())
```

`integrand(H, ξ)` works for any number type of `ξ`, real or complex. An evaluator may
also dispatch on the problem, e.g. `evaluate(::MyEvaluator, H::HFunction{FactorProblem}, x, atol)`,
to use the structure of `f_t`. See `examples/custom_evaluator.jl`.

### Your own problem

With integer zeros and a scale `S` satisfying the three hypotheses of Theorem 5 on
`I* = [first - 1/2, last + 1/2]` — `|I*| ≤ S`; `φ ≥ 1/S` on candidate blocks without
zeros; `φ(d + s) ≤ S s²` near every zero `d` for `|s| ≤ 1/2` — the guarantees carry over:

```julia
prob = IntegerZeroProblem(φ, lo:hi, S; isroot = m -> exact_check(m))
locate(prob, ArbEvaluator())
```

`φ` must accept the number types of the evaluators you use (complex balls for
`ArbEvaluator`). For zeros that are not integers, `localize` runs the same bisection
with `α`, the threshold and the resolution chosen by you; there is no certified gap
then (see `examples/general_function.jl`).

## Guarantees

A result is proved only for a certified problem evaluated with a rigorous evaluator: the
`Oracle` rejects any rigorous value whose error bound exceeds `atol`. With a
non-rigorous evaluator the answer can be wrong; the search reports `:inconsistent`
(a positive block with two negative halves) or `:unverified` (the final candidate fails
the exact check) when it notices.

## What it costs

The number of values of `H` grows like `log t`, but the cost of one value grows with
`t`: the integration interval has length about `t` and the peaks are narrower than
`t⁻²`. Divisor search on a 2-core machine (`examples/compare_evaluators.jl`; times
include all queries):

| t | queries | QuadGK Float64 | Gauss–Legendre Float64 | Arb, Gaussian kernel | Arb, erfc kernel |
|---|---|---|---|---|---|
| 21 = 3·7 | 7 | 13 s | 14 s | 0.06 s | 0.03 s |
| 91 = 7·13 | 8 | — | — | 0.3 s | 0.09 s |
| 323 = 17·19 | 10 | — | — | 1.5 s | 0.23 s |
| 1001 = 7·11·13 | 12 | — | — | 5.2 s | 1.7 s |
| 2021 = 43·47 | 13 | — | — | 11.9 s | 3.0 s |

The rigorous Arb integration is far cheaper here than the safe uniform subdivision,
and its time grows roughly linearly in `t`, i.e. exponentially in the bit length. How
the cost of computing `H` depends on the input is the open question of §6 of the paper,
and the evaluator interface exists to study it.

## Tests

```julia
pkg> test MassRoots
```

## License

MIT
