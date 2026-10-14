"""
    MassRoots

Localization of real zeros by integral mass tests.

A nonnegative function `φ` is turned into a kernel `q_α(x) = g(α φ(x))` that equals one
at the zeros of `φ` and is small elsewhere. Differences of the primitive
`H_α(x) = ∫_{x₀}^{x} q_α(ξ) dξ` measure the mass on an interval, and bisection with a
mass threshold locates a zero.

The package separates three things:

* **problems** – the function `φ`, its search domain and, for problems with integer
  zeros, the scale `S` that certifies the mass gap ([`FactorProblem`](@ref),
  [`QuadraticProblem`](@ref), [`IntegerZeroProblem`](@ref), [`GeneralProblem`](@ref));
* **evaluators** – interchangeable ways of computing `H_α(x)`
  ([`QuadGKEvaluator`](@ref), [`GaussLegendreEvaluator`](@ref), [`ArbEvaluator`](@ref),
  [`IdealEvaluator`](@ref), or your own subtype of [`AbstractEvaluator`](@ref));
* **search** – the bisection logic, which only ever sees values of `H_α` through a
  counting [`Oracle`](@ref) ([`locate`](@ref), [`decide`](@ref), [`find_divisor`](@ref),
  [`factorize`](@ref), [`localize`](@ref)).
"""
module MassRoots

using QuadGK
using SpecialFunctions: erfc

export AbstractProblem, IntegerZeroProblem, FactorProblem, QuadraticProblem, GeneralProblem
export phi, candidates, scale, isroot, basepoint, domain
export GaussKernel, ErfcKernel, NoModulation, ShiftedPower
export MassTestParameters, mass_test_parameters
export HFunction, integrand, HValue
export AbstractEvaluator, evaluate, QuadGKEvaluator, GaussLegendreEvaluator, ArbEvaluator,
       IdealEvaluator
export Oracle, query, nqueries
export LocateResult, locate, decide, find_divisor, factorize, localize

include("numbers.jl")
include("problems.jl")
include("kernels.jl")
include("parameters.jl")
include("hfunction.jl")
include("evaluators/interface.jl")
include("evaluators/quadgk.jl")
include("evaluators/gausslegendre.jl")
include("evaluators/ideal.jl")
include("evaluators/arb.jl")
include("oracle.jl")
include("search.jl")

end # module
