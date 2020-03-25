module TwoDimensionalTests

using TwoDimensional
using TwoDimensional: WeightedPoint, half, getaxisbounds
import TwoDimensional: distance
using Test, LinearAlgebra
import Base.MathConstants: φ

distance(a::Real, b::Real) = distance(promote(a, b)...)
distance(a::T, b::T) where {T<:Unsigned} = ifelse(a < b, b - a, a - b)
distance(a::T, b::T) where {T<:Real} = abs(a - b)

distance(A::NTuple{2,Real}, B::NTuple{2,Real}) =
    hypot(A[1] - B[1], A[2] - B[2])

distance(A::AffineTransform, B::AffineTransform) =
    max(abs(A.xx - B.xx), abs(A.xy - B.xy), abs(A.x - B.x),
        abs(A.yx - B.yx), abs(A.yy - B.yy), abs(A.y - B.y))

function randomize!(A::Matrix{Bool}, bits::Int)
    @inbounds @simd for i in eachindex(A)
        A[i] = (bits&1) == 1
        bits >>= 1
    end
    return A
end

# Naive implementation of the bounding box algorithm.
naiveboundingbox(A::AbstractMatrix{Bool}) = naiveboundingbox(identity, A)
function naiveboundingbox(f::Function, A::AbstractMatrix)
    I, J = axes(A)
    imin = jmin = typemax(Int)
    imax = jmax = typemin(Int)
    @inbounds for j in J, i in I
        if f(A[i,j])
            imin = min(imin, i)
            imax = max(imax, i)
            jmin = min(jmin, j)
            jmax = max(jmax, j)
        end
    end
    return BoundingBox(imin, imax, jmin, jmax)
end

@testset "Points and bounding-boxes" begin
    @testset "Simple points" begin
        types = (Int8, Int32, Int64, Float32, Float64)
        for T1 in types, T2 in types
            @test promote_type(Point{T1}, Point{T2}) ===
                Point{promote_type(T1,T2)}
        end
        P1 = Point(1, 1.3)
        @test eltype(P1) == Float64
        @test Point(P1) === P1
        @test Point{eltype(P1)}(P1) === P1
        @test eltype(Point(1,2)) == Int
        @test Point(1,2) === Point(y=2, x=1)
        @test Point(CartesianIndex(7,8)) === Point(7,8)
        @test CartesianIndex(Point(2,3)) === CartesianIndex(2,3)
        #@test convert(CartesianIndex, Point(2,3)) === CartesianIndex(2,3)
        #@test convert(CartesianIndex{2}, Point(2,3)) === CartesianIndex(2,3)
        #@test convert(Tuple, Point(2,3)) === (2,3)
        #@test convert(Tuple{Float64,Int}, Point(2,3)) === (2.0,3)
        @test Point(Tuple(P1)...) === P1
        @test Point(Tuple(P1)) === P1
        @test Point{Float32}(P1) === Float32.(P1)
        @test promote(P1) === (P1,)
        P2 = Point(3,11)
        P3 = Point{Float32}(-2.1,11.8)
        @test promote(P1, P2) === (P1,Float64.(P2))
        @test promote(P1, P2, P3) === (P1,Float64.(P2),Float64.(P3))
        @test convert(Point, P1) === P1
        @test convert(Point{Float32}, P1) === Float32.(P1)
        @test Tuple(Point(0.0,1)) === (0.0,1.0)
        for T in (Float32, Int16)
            @test zero(Point{T}) === Point(zero(T),zero(T))
            @test one(Point{T}) === Point(one(T),one(T))
            @test oneunit(Point{T}) === Point(oneunit(T),oneunit(T))
        end
        # round
        @test round(Point(1.2,-0.7)) === Point(1.0,-1.0)
        @test round(Point(1,-7)) === Point(1,-7)
        @test round(Float32, Point(1,-7)) === Point{Float32}(1,-7)
        @test round(Float32, Point(1.0,-7.6)) === Point{Float32}(1,-8)
        @test round(Int32, Point(1,-7)) === Point{Int32}(1,-7)
        @test round(Int, Point(1.2,-0.7)) === Point(1,-1)
        @test round(Int, Point(1,-4)) === Point(1,-4)
        @test round(Point{Float64}, P1) === Point(map(round, Tuple(P1)))
        @test round(Point{Int}, Point(1.6,-0.7)) === Point(2,-1)
        @test round(Int, Point(1.2,-0.7)) === Point(1,-1)
        @test round(Point{Int}, Point(1.6,-0.7)) === Point(2,-1)
        # floor
        @test floor(Point(-1,3)) === Point(-1,3)
        @test floor(Int, Point(-1,3)) === Point(-1,3)
        @test floor(Point{Int}, Point(-1,3)) === Point(-1,3)
        @test floor(Int16, Point(-1,3)) === Point{Int16}(-1,3)
        @test floor(Point{Int16}, Point(-1,3)) === Point{Int16}(-1,3)
        @test floor(Point(-1.7,3.2)) === Point(-2.0,3.0)
        @test floor(Int, Point(-1.7,3.2)) === Point(-2,3)
        @test floor(Float32, Point(-1,3)) === Point{Float32}(-1,3)
        @test floor(Point{Float32}, Point(-1.7,3.2)) === Point{Float32}(-2,3)
        # ceil
        @test ceil(Point(-1,3)) === Point(-1,3)
        @test ceil(Int, Point(-1,3)) === Point(-1,3)
        @test ceil(Point{Int}, Point(-1,3)) === Point(-1,3)
        @test ceil(Int16, Point(-1,3)) === Point{Int16}(-1,3)
        @test ceil(Point{Int16}, Point(-1,3)) === Point{Int16}(-1,3)
        @test ceil(Point(-1.7,3.2)) === Point(-1.0,4.0)
        @test ceil(Int, Point(-1.7,3.2)) === Point(-1,4)
        @test ceil(Float32, Point(-1,3)) === Point{Float32}(-1,3)
        @test ceil(Point{Float32}, Point(-1.7,3.2)) === Point{Float32}(-1,4)
        # other methods
        @test hypot(P1) == hypot(P1.x, P1.y)
        @test atan(P1) == atan(P1.y, P1.x)
        @test distance(P1, Point(0,0)) == hypot(P1)
        @test distance(Point(0x02,0x05), Point(0x00,0x00)) == hypot(0x02,0x05)
    end
    @testset "Weighted points" begin
        types = (Float32, Float64)
        for T1 in types, T2 in types
            @test promote_type(WeightedPoint{T1}, WeightedPoint{T2}) ===
                WeightedPoint{promote_type(T1,T2)}
        end
        P2 = WeightedPoint(x=0.1, y=-6, w=0.1)
        @test eltype(P2) == Float64
        @test WeightedPoint(P2) === P2
        @test WeightedPoint{eltype(P2)}(P2) === P2
        @test WeightedPoint{Float64}(P2) === P2
        @test WeightedPoint(P2.w, P2.x, P2.y) === P2
        @test WeightedPoint(Tuple(P2)...) === P2
        @test WeightedPoint(Tuple(P2)) === P2
        @test WeightedPoint{Float32}(P2) === Float32.(P2)
        @test WeightedPoint(Point(3.1,4.2)) ===
            WeightedPoint(w=1, x=3.1, y=4.2)
    end
    @testset "Bounding boxes" begin
        types = (Int8, Int32, Int64, Float32, Float64)
        for T1 in types, T2 in types
            @test promote_type(BoundingBox{T1}, BoundingBox{T2}) ===
                BoundingBox{promote_type(T1,T2)}
        end
        B = BoundingBox(2,3,4,5)
        δ = 0.1
        @test eltype(B) === Int
        @test BoundingBox(B) === B
        @test BoundingBox(2,3,4,5.0) === BoundingBox(2.0,3.0,4.0,5.0)
        @test BoundingBox{eltype(B)}(B) === B
        @test BoundingBox{Float32}(B) ===
            BoundingBox{Float32}(B.xmin, B.xmax, B.ymin, B.ymax)
        @test convert(BoundingBox, B) === B
        @test convert(BoundingBox{Int}, B) === B
        @test BoundingBox(Point(2,3),Point(4,5)) === BoundingBox(2,4,3,5)
        @test BoundingBox(CartesianIndex(2,3),CartesianIndex(4,5)) ===
            BoundingBox(2,4,3,5)
        @test first(B) === Point(B.xmin,B.ymin)
        @test last(B) === Point(B.xmax,B.ymax)
        @test isempty(typemax(BoundingBox{Int})) == false
        @test isempty(typemin(BoundingBox{Float64})) == true
        @test BoundingBox{Float32}(nothing) === typemin(BoundingBox{Float32})
        @test size(BoundingBox{Int32}(nothing)) === (0,0)
        Rx, Ry = 3:7, -6:-3
        B = BoundingBox(Rx,Ry)
        @test size(B) === (length(Rx),length(Ry))
        @test size(BoundingBox{Int32}(B)) === size(B)
        @test axes(B) === (Rx,Ry)
        @test axes(BoundingBox{Int32}(B)) === (Rx,Ry)
        for k in (1,2)
            @test size(B, k) === size(B)[k]
            @test axes(B, k) === axes(B)[k]
            @test size(BoundingBox{Int32}(B), k) === size(B)[k]
            @test axes(BoundingBox{Int32}(B), k) === axes(B)[k]
        end
        @test size(B,40) === 1
        @test_throws ErrorException size(B, 0)
        B1 = BoundingBox(2,3,4,5)
        B2 = BoundingBox(1.1,3.2,4.5,5.8)
        B3 = BoundingBox{Float32}(1.1,3.2,4.5,5.8)
        @test promote(B1, B2) === (Float64.(B1),B2)
        @test promote(B1, B2, B3) === (Float64.(B1),B2,Float64.(B3))
        @test BoundingBox(Tuple(B)...) === B
        @test BoundingBox(Tuple(B)) === B
        @test_deprecated BoundingBox(ones(5,7)) === BoundingBox(1,5, 1,7)
        @test BoundingBox(-2:6,8:11) === BoundingBox(-2,6, 8,11)
        @test BoundingBox((2:4,-1:7)) === BoundingBox(2,4, -1,7)
        @test CartesianIndices(BoundingBox(2:4,-1:7)) ===
            CartesianIndices((2:4,-1:7))
        @test getaxisbounds(9:13) === (9,13)
        @test getaxisbounds(13:9) === (13,12)
        @test getaxisbounds(13:-1:9) === (9,13)
        @test_throws ArgumentError getaxisbounds(13:-2:9)
        let A = ones(7,8)
            A[1:2,:] .= 0
            A[4,:] .= 0
            A[end,:] .= 0
            @test BoundingBox(x -> x > 0, A) === BoundingBox(3,6, 1,8)
            @test BoundingBox(A .> 0) === BoundingBox(3,6, 1,8)
            A[:,1] .= 0
            A[:3:4] .= 0
            A[:,end-1:end] .= 0
            A[2,2] = 1
            @test BoundingBox(x -> x > 0, A) === BoundingBox(2,6, 2,6)
            @test BoundingBox(A .> 0) === BoundingBox(2,6, 2,6)
        end
        @test BoundingBox{Float32}(nothing) ===
            BoundingBox{Float32}(Inf,-Inf,Inf,-Inf)
        @test typemin(BoundingBox{Float64}) === BoundingBox(Inf,-Inf,Inf,-Inf)
        @test typemax(BoundingBox{Float64}) === BoundingBox(-Inf,Inf,-Inf,Inf)
        # round
        @test round(BoundingBox(1.1,2.1,-3.6,7.7)) ===
            BoundingBox(1.0,2.0,-4.0,8.0)
        @test round(BoundingBox{Int32}, BoundingBox(1.1,2.1,-3.6,7.7)) ===
            BoundingBox{Int32}(1,2,-4,8)
        @test round(Int, BoundingBox(1,2,-3,7)) === BoundingBox(1,2,-3,7)
        @test round(Int32, BoundingBox(1.1,2.1,-3.6,7.7)) ===
            BoundingBox{Int32}(1,2,-4,8)
        @test round(Int32, BoundingBox(1,2,-4,8)) ===
            BoundingBox{Int32}(1,2,-4,8)
        @test round(Float32, BoundingBox{Int32}(1,2,-4,8)) ===
            BoundingBox{Float32}(1,2,-4,8)
        @test round(Float32, BoundingBox(1.1,2.7,-4.6,8.3)) ===
            BoundingBox{Float32}(1,3,-5,8)
        # exterior
        @test exterior(B) === B
        @test exterior(Int, B) === B
        @test exterior(Int16, B) === BoundingBox{Int16}(B)
        @test exterior(Float32, B) === BoundingBox{Float32}(B)
        @test exterior(BoundingBox{Int}, B) === B
        @test exterior(B + δ) === B + 1.0
        @test exterior(Int, B + δ) === B + 1
        @test exterior(BoundingBox{Int}, B) === BoundingBox{Int}(exterior(B))
        @test exterior(Float32, B + δ) ===
            BoundingBox{Float32}(exterior(B + δ))
        # interior
        @test interior(B) === B
        @test interior(Int, B) === B
        @test interior(Int16, B) === BoundingBox{Int16}(B)
        @test interior(Float32, B) === BoundingBox{Float32}(B)
        @test interior(BoundingBox{Int}, B) === B
        @test interior(B + δ) === Float64.(B)
        @test interior(Int, B - δ) === B - 1
        @test interior(BoundingBox{Int}, B) === BoundingBox{Int}(interior(B))
        @test interior(Float32, B + δ) === BoundingBox{Float32}(interior(B + δ))
        # other methods
        @test area(BoundingBox(2,4,5,8)) == 6
        @test area(BoundingBox(2.0,4.0,5.0,8.0)) == 6.0

        C = Point(x = 0.5*(B.xmin + B.xmax), y = 0.5*(B.ymin + B.ymax))
        @test center(B) === C
        @test center(BoundingBox{Float64}(B)) === C
        B1 = BoundingBox(-2:6, -7:-1)
        B2 = BoundingBox(1:8, -9:3)
        @test B1 ∪ B2 === BoundingBox(-2:8, -9:3)
        @test B1 ∩ B2 === BoundingBox(1:6, -7:-1)
        A = rand(7,8)
        X, Y = 2:4, 1:3
        B = BoundingBox(X, Y)
        @test A[X,Y] == A[B]
        @test view(A,X,Y) === view(A, B)
    end
    @testset "Bounding-box algorithm" begin
        if true
            # Exhaustively test all possibilities.
            A = Array{Bool,2}(undef, (5,4))
            n = (1 << length(A)) # number of possibilities
            for bits in 0:n-1
                randomize!(A, bits)
                @test BoundingBox(A) === naiveboundingbox(A)
            end
        else
            # Semi-exhaustive testing of the bounding box algorithm.
            A = zeros(Bool, (7,8))
            fill!(A,false)[1,3] = true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,true)[2:end-1,2:end-1] .= false
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[1,1] = true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[1,end] = true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[end,1] = true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[end,end] = true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[3,4] = true
            @test BoundingBox(A) === naiveboundingbox(A)
            A[2,7] = true
            @test BoundingBox(A) === naiveboundingbox(A)
            A[5,end] = true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[end,4] = true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[1,:] .= true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[:,1] .= true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[end,:] .= true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[:,end] .= true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[1,2:end-1] .= true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[2:end-1,1] .= true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[end,2:end-1] .= true
            @test BoundingBox(A) === naiveboundingbox(A)
            fill!(A,false)[2:end-1,end] .= true
            @test BoundingBox(A) === naiveboundingbox(A)
        end
    end
    @testset "Arithmetic" begin
        @test half(Float64) == 0.5
        @test half(Float32) == one(Float32)/2
        @test 2*Point(3,4) === Point(6,8)
        @test Point(3,4)*3 === Point(9,12)
        @test 3\Point(3,4) === Point(3/3,4/3)
        @test Point(3,4)/3 === Point(3/3,4/3)
        @test -Point(2,-9) === Point(-2,9)
        @test -BoundingBox(2,3,4,5) === BoundingBox(-3,-2,-5,-4)
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

    @testset "construction/conversion" begin
        @test_deprecated BigFloat(A) === BigFloat.(A)
        @test_deprecated Float32(A) === Float32.(A)
        for T1 in types, T2 in types
            @test promote_type(AffineTransform{T1}, AffineTransform{T2}) ===
                AffineTransform{promote_type(T1,T2)}
        end
        @test AffineTransform(A) === A
        @test AffineTransform{eltype(A)}(A) === A
        for G in (I, A, B)
            @test eltype(G) == Float64
            for T in types
                @test typeof(AffineTransform{T}(G)) == AffineTransform{T}
                @test typeof(convert(AffineTransform{T}, G)) ==
                    AffineTransform{T}
                @test eltype(AffineTransform{T}(G)) == T
                @test eltype(T.(G)) == T
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
            @test distance(G*v, G(v)) ≤ 0
            @test distance(G(Point(v)), Point(G(v))) ≤ 0
            @test distance(G*Point(v), Point(G(v))) ≤ 0
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
            @test eltype(T1.(A)*T2.(B)) == T
        end
        for v in vectors
            @test distance(compose(A,B,C)(v), A(B(C(v)))) ≤ tol
            @test distance((A*B*C)(v), A(B(C(v)))) ≤ tol
            @test distance(compose(A,C,B,C)(v), A(C(B(C(v))))) ≤ tol
            @test distance((A*C*B*C)(v), A(C(B(C(v))))) ≤ tol
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
            @test eltype(T1.(A)/T2.(B)) == T
            @test eltype(T1.(A)\T2.(B)) == T
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
            H = T.(G)
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
        for v in vectors
            @test distance(A - v, A + (-v[1], -v[2])) ≤ tol
            @test distance(A + Point(v), translate(A, v...)) ≤ tol
            @test distance(Point(v) + A, translate(v..., A)) ≤ tol
            @test distance(A - Point(v), A - v) ≤ tol
        end
        for G in (A, B, C),
            v in vectors,
            T in types
            @test eltype(T.(v) + G) == eltype(G)
            @test eltype(G + T.(v)) == eltype(G)
            H = T.(G)
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
            H = T.(G)
            @test eltype(rotate(T.(θ), H)) == eltype(H)
            @test eltype(rotate(H, T.(θ))) == eltype(H)
        end
    end

    @testset "intercept" begin
        for M in (I, A, B)
            x, y = intercept(M)
            @test distance(M(x, y), (0,0)) ≤ tol
            P = intercept(Point, M)
            @test distance(M*P, Point(0,0)) ≤ tol
        end
    end

    @testset "show" begin
        for M in (I, A, B, C)
            @test occursin(r"AffineTransform{Float64}", string(M))
        end
    end
end

end # module
