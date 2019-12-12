module TwoDimensionalTests

using TwoDimensional
using TwoDimensional: WeightedPoint
using Test, LinearAlgebra
import Base.MathConstants: φ

distance(a::Real, b::Real) = abs(a - b)

distance(a::NTuple{2,Real}, b::NTuple{2,Real}) =
    hypot(a[1] - b[1], a[2] - b[2])

distance(A::AffineTransform, B::AffineTransform) =
    max(abs(A.xx - B.xx), abs(A.xy - B.xy), abs(A.x - B.x),
        abs(A.yx - B.yx), abs(A.yy - B.yy), abs(A.y - B.y))

@testset "Points and bounding-boxes" begin
    @testset "Simple points" begin
        P1 = Point(1, 1.3)
        @test Point(P1) === P1
        @test eltype(P1) == Float64
        @test Point{eltype(P1)}(P1) === P1
        @test eltype(Point(1,2)) == Int
        @test Point(1,2) === Point(y=2, x=1)
        @test Point(CartesianIndex(7,8)) === Point(7,8)
        @test CartesianIndex(Point(2,3)) === CartesianIndex(2,3)
        @test convert(CartesianIndex, Point(2,3)) === CartesianIndex(2,3)
        @test convert(CartesianIndex{2}, Point(2,3)) === CartesianIndex(2,3)
        @test convert(Tuple, Point(2,3)) === (2,3)
        @test convert(Tuple{Float64,Int}, Point(2,3)) === (2.0,3)
        @test Point(Tuple(P1)...) === P1
        @test Point{Float32}(P1) === Float32.(P1)
        @test Tuple(Point(0.0,1)) === (0.0,1.0)
        @test nearest(Int, Point(1.2,-0.7)) === Point(1,-1)
        @test nearest(Point{Int}, Point(1.6,-0.7)) === Point(2,-1)
        @test round(Int, Point(1.2,-0.7)) === Point(1,-1)
        @test round(Point{Int}, Point(1.6,-0.7)) === Point(2,-1)
        @test hypot(P1) == hypot(P1.x, P1.y)
        @test atan(P1) == atan(P1.y, P1.x)
        @test 2*Point(3,4) === Point(6,8)
        @test Point(3,4)*3 === Point(9,12)
    end
    @testset "Weighted points" begin
        P2 = WeightedPoint(x=0.1, y=-6, w=0.1)
        @test WeightedPoint(P2.w, P2.x, P2.y) === P2
        @test WeightedPoint(P2) === P2
        @test eltype(P2) == Float64
        @test WeightedPoint{eltype(P2)}(P2) === P2
        @test WeightedPoint(Tuple(P2)...) === P2
        @test WeightedPoint{Float32}(P2) === Float32.(P2)
    end
end

@testset "AffineTransforms" begin
    tol = 1e-14
    I = AffineTransform()
    A = AffineTransform(1, 0, -3, 0.1, 1, +2)
    B = AffineTransform(-0.4,  0.1, -4.2, -0.3,  0.7,  1.1)
    C = AffineTransform( 2.3, -0.9, -6.1,  0.7, -3.1, -5.2)
    vectors = ((0.2,1.3), (-1,π), (-sqrt(2),3//4))
    scales = (2, 0.1, φ)
    angles = (-2π/11, π/7, 0.1)
    types = (BigFloat, Float64, Float32, Float16)
    @test_deprecated BigFloat(A)

    @testset "conversion" begin
        @test_throws ErrorException compose()
        for G in (I, A, B)
            @test eltype(G) == Float64
            for T in types
                @test typeof(AffineTransform{T}(G)) == AffineTransform{T}
                @test typeof(convert(AffineTransform{T}, G)) == AffineTransform{T}
                @test typeof(T(G)) == AffineTransform{T} # FIXME: deprecated
                @test eltype(AffineTransform{T}(G)) == T
            end
        end
    end

    @testset "identity" begin
        @test det(I) == 1
        @test distance(inv(I), I) ≤ 0
        for v in vectors
            @test distance(I(v), eltype(I).(v)) ≤ 0
        end
    end

    @testset "apply" begin
        for G in (I, A, B),
            v in vectors
            @test distance(G(v...), G(v)) ≤ 0
        end
    end

    @testset "composition" begin
        for G in (I, B, A)
            @test distance(G, compose(G)) ≤ 0
            for H in (A, B)
                @test distance(G*H, compose(G,H)) ≤ 0
                @test distance(G⋅H, compose(G,H)) ≤ 0
                @test distance(G∘H, compose(G,H)) ≤ 0
                for v in vectors
                    @test distance((G*H)(v), G(H(v))) ≤ tol
                end
            end
        end
        for T1 in types, T2 in types
            T = promote_type(T1, T2)
            @test eltype(T1(A)*T2(B)) == T
        end
        for v in vectors
            @test distance((A*B*C)(v), A(B(C(v)))) ≤ tol
        end
    end

    @testset "jacobian" begin
        for M in (I, B, A)
            @test jacobian(M) == abs(det(M))
        end
    end

    @testset "inverse" begin
        for M in (B, A)
            if det(M) == 0
                continue
            end
            @test distance(det(inv(M)), 1/det(M)) ≤ tol
            @test distance(M/M, M*inv(M)) ≤ tol
            @test distance(M\M, inv(M)*M) ≤ tol
            @test distance(M\M, I) ≤ tol
            @test distance(M/M, I) ≤ tol
            for v in vectors
                @test distance(M(inv(M)(v)), v) ≤ tol
                @test distance(inv(M)(M(v)), v) ≤ tol
                @test distance((M\M)(v), v) ≤ tol
                @test distance((M/M)(v), v) ≤ tol
            end
        end
        for T1 in types, T2 in types
            T = promote_type(T1, T2)
            @test eltype(T1(A)/T2(B)) == T
            @test eltype(T1(A)\T2(B)) == T
        end
    end

    @testset "scale" begin
        for M in (A, B, C),
            α in scales,
            v in vectors
            @test distance((α*M)(v), α.*M(v)) ≤ tol
            @test distance((M*α)(v), M(α.*v)) ≤ tol
        end
        for G in (A, B, C),
            α in scales,
            T in types
            @test eltype(T(α)*G) == eltype(G)
            @test eltype(G*T(α)) == eltype(G)
            H = T(G)
            @test eltype(α*H) == eltype(H)
            @test eltype(H*α) == eltype(H)
        end
    end

    @testset "translation" begin
        for M in (B, A),
            t in vectors,
            v in vectors
            @test distance(translate(t, M)(v), t .+ M(v)) ≤ tol
            @test distance(translate(t, M)(v), (t + M)(v)) ≤ tol
            @test distance(translate(M, t)(v), M(v .+ t)) ≤ tol
            @test distance(translate(M, t)(v), (M + t)(v)) ≤ tol
        end
        for G in (A, B, C),
            v in vectors,
            T in types
            @test eltype(T.(v) + G) == eltype(G)
            @test eltype(G + T.(v)) == eltype(G)
            H = T(G)
            @test eltype(v + H) == eltype(H)
            @test eltype(H + v) == eltype(H)
        end
    end

    @testset "rotation" begin
        for θ in angles,
            v in vectors
            R = rotate(+θ, I)
            Q = rotate(-θ, I)
            @test distance(R*Q, I) ≤ tol
            @test distance(Q*R, I) ≤ tol
            @test distance(rotate(θ, B)(v), (R*B)(v)) ≤ tol
            @test distance(rotate(B, θ)(v), (B*R)(v)) ≤ tol
        end
        for G in (A, B, C),
            θ in angles,
            T in types
            @test eltype(rotate(T(θ), G)) == eltype(G)
            @test eltype(rotate(G, T(θ))) == eltype(G)
            H = T(G)
            @test eltype(rotate(T(θ), H)) == eltype(H)
            @test eltype(rotate(H, T(θ))) == eltype(H)
        end
    end

    @testset "intercept" begin
        for M in (I, A, B)
            x, y = intercept(M)
            @test distance(M(x, y), (0,0)) ≤ tol
        end
    end

    @testset "show" begin
        for M in (I, A, B, C)
            @test occursin(r"AffineTransform{Float64}", string(M))
        end
    end
end

end # module
