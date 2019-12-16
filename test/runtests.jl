module TwoDimensionalTests

using TwoDimensional
using TwoDimensional: WeightedPoint, center, half, exterior, interior, area
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
        @test eltype(P1) == Float64
        @test Point(P1) === P1
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
        @test Point(Tuple(P1)) === P1
        @test Point{Float32}(P1) === Float32.(P1)
        @test convert(Point, P1) === P1
        @test convert(Point{Float32}, P1) === Float32.(P1)
        @test Tuple(Point(0.0,1)) === (0.0,1.0)
        @test nearest(Int, Point(1.2,-0.7)) === Point(1,-1)
        @test nearest(Point{Float64}, P1) === Point(map(round, Tuple(P1)))
        @test nearest(Point{Int}, Point(1.6,-0.7)) === Point(2,-1)
        @test round(Int, Point(1.2,-0.7)) === Point(1,-1)
        @test round(Point{Int}, Point(1.6,-0.7)) === Point(2,-1)
        @test hypot(P1) == hypot(P1.x, P1.y)
        @test atan(P1) == atan(P1.y, P1.x)
    end
    @testset "Weighted points" begin
        P2 = WeightedPoint(x=0.1, y=-6, w=0.1)
        @test eltype(P2) == Float64
        @test WeightedPoint(P2) === P2
        @test WeightedPoint{eltype(P2)}(P2) === P2
        @test WeightedPoint{Float64}(P2) === P2
        @test WeightedPoint(P2.w, P2.x, P2.y) === P2
        @test WeightedPoint(Tuple(P2)...) === P2
        @test WeightedPoint(Tuple(P2)) === P2
        @test WeightedPoint{Float32}(P2) === Float32.(P2)
        @test WeightedPoint(Point(3.1,4.2)) === WeightedPoint(w=1, x=3.1, y=4.2)
    end
    @testset "Bounding boxes" begin
        B = BoundingBox(2,3,4,5)
        @test eltype(B) === Int
        @test BoundingBox(B) === B
        @test BoundingBox{eltype(B)}(B) === B
        @test BoundingBox{Float32}(B) === BoundingBox{Float32}(B.xmin, B.xmax, B.ymin, B.ymax)
        @test BoundingBox(Point(2,3),Point(4,5)) === BoundingBox(2,4,3,5)
        @test BoundingBox(CartesianIndex(2,3),CartesianIndex(4,5)) === BoundingBox(2,4,3,5)
        @test BoundingBox(Tuple(B)...) === B
        @test BoundingBox(Tuple(B)) === B
        @test BoundingBox(rand(5,7)) === BoundingBox(1,5, 1,7)
        @test BoundingBox(-2:6,8:11) === BoundingBox(-2,6, 8,11)
        @test BoundingBox((2:4,-1:7)) === BoundingBox(2,4, -1,7)
        A = ones(7,8)
        A[1:2,:] .= 0
        A[4,:] .= 0
        A[end,:] .= 0
        @test BoundingBox(x -> x > 0, A) === BoundingBox(3,6, 1,8)
        A[:,1] .= 0
        A[:3:4] .= 0
        A[:,end-1:end] .= 0
        A[2,2] = 1
        @test BoundingBox(x -> x > 0, A) === BoundingBox(2,6, 2,6)
        @test BoundingBox{Float32}(nothing) === BoundingBox{Float32}(Inf,-Inf,Inf,-Inf)
        @test typemin(BoundingBox{Float64}) === BoundingBox(Inf,-Inf,Inf,-Inf)
        @test typemax(BoundingBox{Float64}) === BoundingBox(-Inf,Inf,-Inf,Inf)

        @test exterior(B) === B
        @test exterior(B + 0.1) === B + 1.0
        @test interior(B) === B
        @test interior(B + 0.1) === Float64.(B)
        @test area(BoundingBox(2,4,5,8)) == 6
        @test area(BoundingBox(2.0,4.0,5.0,8.0)) == 6.0

        C = Point(x = 0.5*(B.xmin + B.xmax), y = 0.5*(B.ymin + B.ymax))
        @test center(B) === C
        @test center(BoundingBox{Float64}(B)) === C
    end
    @testset "Arithmetic" begin
        @test 2*Point(3,4) === Point(6,8)
        @test Point(3,4)*3 === Point(9,12)
        @test 3\Point(3,4) === Point(3/3,4/3)
        @test Point(3,4)/3 === Point(3/3,4/3)
        @test -Point(2,-9) === Point(-2,9)
        @test Point(2,-9) + Point(3,7) === Point(5,-2)
        @test Point(2,-9) - Point(3,7) === Point(-1,-16)
        α = 3
        δ = 0.1
        t = Point(-0.2,0.7)
        B = BoundingBox(2,3,4,5)
        @test α*B === BoundingBox(xmin = α*B.xmin, xmax = α*B.xmax,
                                  ymin = α*B.ymin, ymax = α*B.ymax)
        @test B*α === BoundingBox(xmin = α*B.xmin, xmax = α*B.xmax,
                                  ymin = α*B.ymin, ymax = α*B.ymax)
        @test α\B === BoundingBox(xmin = (1/α)*B.xmin, xmax = (1/α)*B.xmax,
                                  ymin = (1/α)*B.ymin, ymax = (1/α)*B.ymax)
        @test B/α === BoundingBox(xmin = (1/α)*B.xmin, xmax = (1/α)*B.xmax,
                                  ymin = (1/α)*B.ymin, ymax = (1/α)*B.ymax)
        @test B + δ === BoundingBox(xmin = B.xmin - δ, xmax = B.xmax + δ,
                                    ymin = B.ymin - δ, ymax = B.ymax + δ)
        @test B - δ === BoundingBox(xmin = B.xmin + δ, xmax = B.xmax - δ,
                                    ymin = B.ymin + δ, ymax = B.ymax - δ)
        @test B + t === BoundingBox(xmin = B.xmin + t.x, xmax = B.xmax + t.x,
                                    ymin = B.ymin + t.y, ymax = B.ymax + t.y)
        @test B - t === BoundingBox(xmin = B.xmin - t.x, xmax = B.xmax - t.x,
                                    ymin = B.ymin - t.y, ymax = B.ymax - t.y)
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
