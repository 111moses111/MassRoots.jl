# The figure in the README: f, the kernel q_α and the primitive H for the double-zero
# example, with the intervals found by `localize`. Writes docs/assets/double_zero.svg
# (plain SVG, no plotting package needed).
#
#   julia --project=examples examples/readme_figure.jl

using MassRoots

f(x) = (x - 1)^2 * (x - 2)
α = 1000
prob = GeneralProblem(f, 0, 3)                         # φ = f²
H = HFunction(prob, GaussKernel(), α)
ev = QuadGKEvaluator(Float64; segment_length = 1 // 64)

xs = [k // 200 for k in 0:600]                         # grid on [0, 3]
fx = [f(Float64(x)) for x in xs]
qx = [Float64(integrand(H, Float64(x))) for x in xs]
Hx = [Float64(evaluate(ev, H, Rational{BigInt}(x), big(1) // 10^8).value) for x in xs]
ivals, _ = localize(prob, ev; α, threshold = 1e-3, resolution = 1 // 64)

# --- minimal SVG writer --------------------------------------------------------------------
W, PH, L, R, T, GAP = 760, 150, 70, 35, 40, 62          # width, panel height, margins
HT = T + 3 * PH + 2 * GAP + 52
X(x) = L + (W - L - R) * x / 3
panel_top(i) = T + (i - 1) * (PH + GAP)
Y(i, y, lo, hi) = panel_top(i) + PH * (hi - y) / (hi - lo)
fmt(v) = string(round(v; digits = 2))

io = IOBuffer()
println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$HT" viewBox="0 0 $W $HT" font-family="Helvetica, Arial, sans-serif" font-size="13">""")
println(io, """<rect width="$W" height="$HT" fill="white"/>""")

panels = [(fx, -2.2, 4.2, "#1f4e79", "f(x) = (x − 1)²(x − 2)"),
          (qx, 0.0, 1.1, "#c0504d", "kernel q(x) = exp(−(α f(x)²)²),  α = 1000"),
          (Hx, 0.0, 0.45, "#2e7d32", "primitive H(x) = ∫₀ˣ q;  shaded: intervals found by localize")]

for (i, (ys, lo, hi, col, title)) in enumerate(panels)
    top = panel_top(i)
    for (a, b) in ivals                                  # shaded intervals
        println(io, """<rect x="$(fmt(X(Float64(a))))" y="$top" width="$(fmt(X(Float64(b)) - X(Float64(a))))" height="$PH" fill="#f4b183" fill-opacity="0.35"/>""")
    end
    println(io, """<rect x="$L" y="$top" width="$(W - L - R)" height="$PH" fill="none" stroke="#444" stroke-width="1"/>""")
    if lo < 0 < hi                                       # zero line
        println(io, """<line x1="$L" x2="$(W - R)" y1="$(fmt(Y(i, 0, lo, hi)))" y2="$(fmt(Y(i, 0, lo, hi)))" stroke="#999" stroke-dasharray="4 3"/>""")
    end
    for xt in 0:0.5:3                                    # x ticks
        println(io, """<line x1="$(fmt(X(xt)))" x2="$(fmt(X(xt)))" y1="$(top + PH)" y2="$(top + PH + 5)" stroke="#444"/>""")
        println(io, """<text x="$(fmt(X(xt)))" y="$(top + PH + 18)" text-anchor="middle" fill="#333">$(xt)</text>""")
    end
    for yt in (lo, (lo + hi) / 2, hi)                    # y ticks
        println(io, """<text x="$(L - 8)" y="$(fmt(Y(i, yt, lo, hi) + 4))" text-anchor="end" fill="#333">$(fmt(yt))</text>""")
    end
    pts = join(("$(fmt(X(Float64(x)))),$(fmt(Y(i, clamp(y, lo, hi), lo, hi)))" for (x, y) in zip(xs, ys)), " ")
    println(io, """<polyline points="$pts" fill="none" stroke="$col" stroke-width="2"/>""")
    println(io, """<text x="$L" y="$(top - 8)" fill="#222">$title</text>""")
end
println(io, """<text x="$(W ÷ 2)" y="$(HT - 8)" text-anchor="middle" fill="#333">x   (double zero at x = 1 without a sign change, simple zero at x = 2)</text>""")
println(io, "</svg>")

out = joinpath(@__DIR__, "..", "docs", "assets", "double_zero.svg")
mkpath(dirname(out))
write(out, String(take!(io)))
println("wrote ", normpath(out), "; intervals: ", [(Float64(a), Float64(b)) for (a, b) in ivals])
