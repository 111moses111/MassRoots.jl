module MassRootsArblibExt

using MassRoots
using Arblib
import MassRoots: lift, liftto, tobigfloat, kernelvalue, arb_evaluate, required_precision,
                  integrand, ErfcKernel, ArbEvaluator, HFunction, HValue

# Exact constants enter at the precision of the ball they meet.
lift(x::Union{Arb,ArbRef}, v) = Arb(v; prec = precision(x))
lift(x::Union{Acb,AcbRef}, v) = Acb(v; prec = precision(x))
liftto(::Type{Arb}, v, prec::Integer) = Arb(v; prec)
liftto(::Type{Acb}, v, prec::Integer) = Acb(v; prec)

# BigFloat midpoint of a real ball, at the ball's precision.
tobigfloat(x::Union{Arb,ArbRef}) = BigFloat(Arblib.midref(x))
radius_big(x::Union{Arb,ArbRef}) = BigFloat(Arb(Arblib.radref(x)))

# erfc on balls through Arb's rigorous hypergeometric implementation.
function kernelvalue(::ErfcKernel, s::Union{Acb,AcbRef})
    r = Acb(; prec = precision(s))
    Arblib.hypgeom_erfc!(r, s)
    return r
end
function kernelvalue(::ErfcKernel, s::Union{Arb,ArbRef})
    r = Arb(; prec = precision(s))
    Arblib.hypgeom_erfc!(r, s)
    return r
end

function arb_evaluate(ev::ArbEvaluator, H::HFunction, x::Rational{BigInt},
                      atol::Rational{BigInt})
    x == H.x0 && return HValue(0, 0; rigorous = true, info = (precision = 0, attempts = 0))
    prec = something(ev.precision, required_precision(H, x, atol))
    # Arb treats the tolerance as a goal and can overshoot it by a small factor, so ask for
    # atol/64 and tighten further if needed.
    k = MassRoots.ceillog2(64 / atol)
    evl, dep = ev.eval_limit, ev.depth_limit
    local v
    for attempt in 1:max(1, ev.max_attempts)
        tol = Arblib.set_ui_2exp!(Mag(), one(UInt), -k)
        v = arb_integrate(H, x, prec, tol, something(ev.deg_limit, 0),
                          something(evl, 0), something(dep, 0), ev.verbose, attempt)
        v.radius <= atol && return v
        # Escalate: tolerance 2^-8 tighter, 16 more bits, limits ×4 (Arb's defaults are
        # eval_limit = 1000 prec + prec², depth_limit = 2 prec).
        evl = 4 * something(evl, 1000prec + prec^2)
        dep = 4 * something(dep, 2prec)
        k += 8
        prec += 16
    end
    return v          # still too wide: the Oracle reports it
end

function arb_integrate(H, x, prec, tol, deg, evl, dep, verbose, attempt)
    a = Acb(H.x0; prec)
    b = Acb(x; prec)
    opts = Arblib.calc_integrate_opt_struct(deg, evl, dep, 0, verbose)   # 0 = Arb default
    f = z -> integrand(H, Acb(z; prec))
    t = @elapsed res = Arblib.integrate(f, a, b; prec, rtol = 0.0, atol = tol, opts,
                                        warn_on_no_convergence = false)
    re = Arblib.realref(res)
    im = Arblib.imagref(res)
    # H(x) is real; the imaginary part is a ball around 0 and widens the bound.
    radius = radius_big(re) + abs(tobigfloat(im)) + radius_big(im)
    return HValue(tobigfloat(re), radius; rigorous = isfinite(radius),
                  info = (precision = prec, attempts = attempt, seconds = t))
end

end # module
