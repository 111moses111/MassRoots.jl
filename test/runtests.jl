using MassRoots
using MassRoots: half, ceillog2, ceilsqrt, bitlength, sinpi2
using Test
using Primes
using Arblib

primefactors(t) = sort!(reduce(vcat, [fill(big(p), k) for (p, k) in factor(t)]; init = BigInt[]))
querybound(M) = 2 + ceillog2(big(M))                 # Theorem 6: 2 + ⌈log₂ M⌉

@testset "MassRoots" begin

@testset "exact helpers" begin
    @test [ceillog2(m) for m in 1:9] == [0, 1, 2, 2, 3, 3, 3, 3, 4]
    @test ceillog2(1 // 3) == 0 && ceillog2(9 // 2) == 3
    @test [ceilsqrt(S) for S in (1, 2, 4, 5, 9, 10)] == [1, 2, 2, 3, 3, 4]
    @test half(3) == 7 // 2
end

@testset "parameters (Theorem 5, Propositions 1–2)" begin
    for t in (4, 15, 91, 1001, 2^20 - 1, big(2)^64 - 1, big(2)^200 + 1)
        P = mass_test_parameters(FactorProblem(t))
        n = bitlength(t)
        @test P.S == 2 * big(t)^2
        @test P.α == P.S^2
        @test P.ν == ceilsqrt(P.S)
        @test P.w == 1 // (2 * P.S * P.ν)
        @test P.background == P.w / 2 && P.minmass == P.w
        @test P.threshold == 3 * P.w / 4
        @test P.p <= 3n + 7                           # Proposition 1
        @test P.atol <= P.w / 16                      # (5.29)
        @test bitlength(P.α) <= 4n + 2                # Theorem 6
    end
    for (a, b, c) in ((1, 1, 1), (2, 3, 50), (3, 5, 1000), (2^20 - 1, 2^20 - 1, 2^20 - 1))
        P = mass_test_parameters(QuadraticProblem(a, b, c))
        @test P.p <= 3 * maximum(bitlength, (a, b, c)) + 15     # Proposition 2
    end
    P1 = mass_test_parameters(450)
    P2 = mass_test_parameters(450, ShiftedPower(1, 0))
    @test (P1.α, P1.threshold, P1.p) == (P2.α, P2.threshold, P2.p)   # modulation (1, 0) gives (5.18)
    P3 = mass_test_parameters(450, ShiftedPower(3, 1 // 2))
    @test P3.background < P3.threshold < P3.minmass
    @test P3.background + 2P3.atol < P3.threshold < P3.minmass - 2P3.atol
    @test mass_test_parameters(3).S == 8                          # raised to S ≥ 8
end

@testset "encodings (Theorems 1 and 3)" begin
    for t in (15, 91, 97, 360)
        P = FactorProblem(t)
        @test candidates(P) == 2:t-1
        for m in candidates(P)
            @test isroot(P, m) == (t % m == 0)
            m <= 40 || continue
            @test iszero(phi(P, big(m))) == (t % m == 0)
            @test iszero(phi(P, Float64(m))) == (t % m == 0)
        end
        @test phi(P, 2.5) > 0.5                     # half-integer: sin²(πx) = 1
    end
    for (a, b, c) in ((1, 1, 5), (2, 3, 50), (3, 7, 100), (5, 2, 7))
        Q = QuadraticProblem(a, b, c)
        sols = [m for m in 0:isqrt(c ÷ a) if (c - a * m^2) >= 0 && (c - a * m^2) % b == 0]
        @test [m for m in candidates(Q) if isroot(Q, m)] == sols
        for m in candidates(Q)
            @test iszero(phi(Q, big(m))) == (m in sols)
        end
    end
    @test_throws ArgumentError FactorProblem(3)
    @test_throws ArgumentError QuadraticProblem(0, 1, 1)
end

@testset "search logic with the ideal oracle (Theorem 6)" begin
    for t in 4:400, seed in 0:2
        r = find_divisor(t, IdealEvaluator(; seed))
        @test r.queries <= querybound(t - 2)
        if isprime(t)
            @test r.status === :none
        else
            @test r.status === :found && t % r.root == 0
        end
    end
    for t in (2, 3, 4, 12, 97, 360, 1001, 65536, 999_983, 2^5 * 3^3 * 7^2)
        f, q = factorize(t, IdealEvaluator())
        @test f == primefactors(t)
    end
    for a in 1:6, b in 1:6, c in 1:60
        Q = QuadraticProblem(a, b, c)
        expected = any(m -> isroot(Q, m), candidates(Q))
        ans, r = decide(Q, IdealEvaluator(; seed = c))
        @test ans == expected
        @test r.queries == 2                      # two nonadaptive queries
        l = locate(Q, IdealEvaluator(; seed = c))
        @test l.queries <= querybound(length(candidates(Q)))
        @test expected ? (l.status === :found) : (l.status === :none)
    end
    # a user-defined integer-zero problem behaves like the built-in one
    t = big(91)
    custom = IntegerZeroProblem(x -> sinpi2(MassRoots.lift(x, t) / x) + sinpi2(x),
                                2:90, 2t^2; isroot = m -> t % m == 0, x0 = 2)
    @test locate(custom, IdealEvaluator()).root == locate(FactorProblem(91), IdealEvaluator()).root
end

@testset "QuadGK evaluator" begin
    for t in (15, 21, 13)            # cost of the safe default subdivision grows like t⁴
        r = find_divisor(t, QuadGKEvaluator(Float64))
        @test isprime(t) ? r.status === :none : (r.status === :found && t % r.root == 0)
        @test r.queries <= querybound(t - 2)
    end
    r = find_divisor(6, QuadGKEvaluator(BigFloat))
    @test r.status === :found && 6 % r.root == 0
    r = find_divisor(15, QuadGKEvaluator(Float64); kernel = ErfcKernel())
    @test r.status === :found
    r = find_divisor(15, QuadGKEvaluator(Float64); modulation = ShiftedPower(2, 1))
    @test r.status === :found
    ans, _ = decide(QuadraticProblem(2, 3, 50), QuadGKEvaluator(Float64))
    @test ans                                     # 2·5² + 3·0 = 50
    ans, _ = decide(QuadraticProblem(3, 7, 2), QuadGKEvaluator(Float64))
    @test !ans
end

@testset "Gauss–Legendre evaluator" begin
    r = find_divisor(21, GaussLegendreEvaluator(Float64; degree = 6))
    @test r.status === :found && 21 % r.root == 0
end

@testset "Arb evaluator (rigorous)" begin
    for t in (15, 13, 21, 35)
        r = find_divisor(t, ArbEvaluator())
        @test isprime(t) ? r.status === :none : (r.status === :found && t % r.root == 0)
        @test all(e -> e.rigorous && e.radius <= r.params.atol, r.oracle.log)
    end
    r = find_divisor(15, ArbEvaluator(); kernel = ErfcKernel())
    @test r.status === :found
    ans, r = decide(QuadraticProblem(2, 3, 50), ArbEvaluator())
    @test ans && r.queries == 2
    # the rigorous ball contains the QuadGK value
    P = FactorProblem(15)
    H = HFunction(P)
    atol = mass_test_parameters(P).atol
    a = evaluate(ArbEvaluator(), H, 29 // 2, atol)
    q = evaluate(QuadGKEvaluator(Float64), H, 29 // 2, atol)
    @test abs(a.value - q.value) <= a.radius + 1e-12
end

@testset "general problem: example (2.6) of the paper" begin
    f(x) = sinpi(x / 2) + sinpi(4x)^2 + 0.1
    prob = GeneralProblem(f, 0, 4; square = true)
    ivals, o = localize(prob, QuadGKEvaluator(Float64; segment_length = 1 // 256);
                        α = 1000, threshold = 1e-4, resolution = 1 // 64)
    # reference zeros by sign changes on a fine grid
    xs = range(0, 4; length = 400_001)
    zs = [xs[i] for i in 1:length(xs)-1 if sign(f(xs[i])) != sign(f(xs[i+1]))]
    @test length(zs) == 14
    @test all(z -> any(iv -> iv[1] <= z <= iv[2], ivals), zs)
    @test length(ivals) == length(zs)              # one cluster per zero
end

end
