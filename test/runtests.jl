module TwoDimensionalTests

using TwoDimensional
using TwoDimensional: WeightedPoint, PointLike, compose, get_x, get_y, get_xy, factors_type, offsets_type
using TwoDimensional: half, get_axis_bounds
using Test, LinearAlgebra, Unitless
import Base.MathConstants: φ

# abs_diff(a, b) yields absolute difference.
# max_abs_diff(a, b) yields maximum absolute difference.
# sum_abs_diff(a, b) yields sum of absolute differences.
abs_diff(a::Number, b::Number) = abs_diff(promote(a, b)...)
abs_diff(a::T, b::T) where {T<:Unsigned} = a < b ? b - a : a - b
abs_diff(a::T, b::T) where {T<:Real} = abs(a - b)
for (pfx, op) in ((:max, :max), (:sum, :(+)))
    norm = Symbol(pfx,"_abs")
    @eval begin
        $norm(a::Number) = abs(a)
        $norm(a::Tuple) = mapreduce($norm, $op, a)
        function $norm(a::AbstractArray)
            init = $norm(zero(eltype(a)))
            return mapreduce($norm, $op, a, b; init = init)
        end
        $norm(a::AbstractPoint) = $norm(Tuple(a))
        $norm(a::BoundingBox) = $norm(Tuple(a))
        $norm(a::AffineTransform) = $norm(Tuple(a))
    end
    norm_diff = Symbol(norm,"_diff")
    @eval begin
        $norm_diff(a::Number, b::Number) = abs_diff(a, b)
        function $norm_diff(a::Tuple, b::Tuple)
            @assert length(a) == length(b)
            return mapreduce(abs_diff, $op, a, b)
        end
        function $norm_diff(a::AbstractArray, b::AbstractArray)
            @assert axes(a) == axes(b)
            init = abs_diff(zero(eltype(a)), zero(eltype(b)))
            return mapreduce(abs_diff, $op, a, b; init = init)
        end
        $norm_diff(a::AbstractPoint, b::AbstractPoint) = $norm_diff(Tuple(a), Tuple(b))
        $norm_diff(a::BoundingBox, b::BoundingBox) = $norm_diff(Tuple(a), Tuple(b))
        $norm_diff(a::AffineTransform, b::AffineTransform) = $norm_diff(Tuple(a), Tuple(b))
    end
end

# Default relative tolerance is eps.
relative_precision(x, y) = max(relative_precision(x), relative_precision(y))
relative_precision(x) = relative_precision(typeof(x))
relative_precision() = false
relative_precision(::Type{T}) where {T} = error("relative tolerance not defined for type $T")
relative_precision(::Type{T}) where {T<:AbstractFloat} = eps(T)
relative_precision(::Type{T}) where {T<:Real} = false
relative_precision(::Type{T}) where {T<:AbstractArray} = relative_precision(eltype(T))
relative_precision(::Type{T}) where {T<:AbstractPoint} = relative_precision(eltype(T))
relative_precision(::Type{T}) where {T<:BoundingBox} = relative_precision(eltype(T))
relative_precision(::Type{T}) where {T<:AffineTransform} = relative_precision(bare_type(T))
relative_precision(::Type{T}) where {T<:Number} = relative_precision(real_type(T))
relative_precision(::Type{<:NTuple{N,T}}) where {N,T} = relative_precision(T)
relative_precision(::Type{T}) where {T<:Tuple} = _relative_precision(false, (T.parameters)...)
_relative_precision(rtol) = rtol # end of recursion
_relative_precision(rtol, ::Type{T}, args...) where {T} =
    _relative_precision(max(rtol, relative_precision(T)), args...)

# This is very similar to `isapprox` but with our custom settings.
function ≈(a, b; atol = false, rtol = 4*relative_precision(a, b), norm=sum_abs)
    d = norm === sum_abs ? sum_abs_diff(a, b) :
        norm === max_abs ? max_abs_diff(a, b) : error("unknown norm")
    if iszero(rtol)
        return d ≤ atol
    else
        return d ≤ max(atol, rtol*max(norm(a), norm(b)))
    end
end

function unpack_bits!(A::Matrix{Bool}, bits::Unsigned)
    @assert length(A) ≤ 8*sizeof(bits)
    @inbounds @simd for i in eachindex(A)
        A[i] = (bits&1) == one(bits)
        bits >>= 1
    end
    return A
end

# Naive implementation of the bounding box algorithm.
naive_bounding_box(A::AbstractMatrix{Bool}) = naive_bounding_box(identity, A)
function naive_bounding_box(f::Function, A::AbstractMatrix)
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
    @testset "Miscellaneous" begin
        # Iterators.
        pnt = @inferred Point(2,3)
        @test pnt === @inferred Point(pnt...)
        @test pnt === @inferred Point(Tuple(pnt)...)
        @test pnt === @inferred Point(Tuple(pnt))
        x,y = pnt
        @test (x,y) === @inferred Tuple(pnt)
        @test (x,y) === (pnt.x, pnt.y)
        @test (x,y) === (pnt[1], pnt[2])
        @test_throws BoundsError pnt[0]
        @test_throws BoundsError pnt[3]

        wpnt = @inferred WeightedPoint(1.2,sqrt(2),3)
        @test wpnt === @inferred WeightedPoint(wpnt...)
        @test wpnt === @inferred WeightedPoint(Tuple(wpnt)...)
        @test wpnt === @inferred WeightedPoint(Tuple(wpnt))
        w,x,y = wpnt
        @test (w,x,y) === @inferred Tuple(wpnt)
        @test (w,x,y) === (wpnt.w, wpnt.x, wpnt.y)
        @test (w,x,y) === (wpnt[1], wpnt[2], wpnt[3])
        @test_throws BoundsError wpnt[0]
        @test_throws BoundsError wpnt[4]

        box = @inferred BoundingBox(1.2,sqrt(2),-3,11)
        @test box === @inferred BoundingBox(box...)
        @test box === @inferred BoundingBox(Tuple(box)...)
        @test box === @inferred BoundingBox(Tuple(box))
        xmin,xmax,ymin,ymax = box
        @test (xmin,xmax,ymin,ymax) === @inferred Tuple(box)
        @test (xmin,xmax,ymin,ymax) === (box.xmin, box.xmax, box.ymin, box.ymax)
        @test (xmin,xmax,ymin,ymax) === (box[1], box[2], box[3], box[4])
        @test_throws BoundsError box[0]
        @test_throws BoundsError box[5]
    end
    @testset "Simple points" begin
        types = (Int8, Int32, Int64, Float32, Float64)
        for T1 in types, T2 in types
            @test promote_type(Point{T1}, Point{T2}) === Point{promote_type(T1,T2)}
        end
        P1 = Point(1, 1.3)
        @test eltype(P1) == Float64
        @test Point(P1) === P1
        @test Point{eltype(P1)}(P1) === P1
        @test eltype(Point(1,2)) == Int
        @test Point(1,2) === Point(y=2, x=1)
        @test Point(CartesianIndex(7,8)) === Point(7,8)
        @test CartesianIndex(Point(2,3)) === CartesianIndex(2,3)
        @test Point(TwoDimensional.WeightedPoint(2,3,7)) === Point(3.0,7.0)
        #@test convert(CartesianIndex, Point(2,3)) === CartesianIndex(2,3)
        #@test convert(CartesianIndex{2}, Point(2,3)) === CartesianIndex(2,3)
        #@test convert(Tuple, Point(2,3)) === (2,3)
        #@test convert(Tuple{Float64,Int}, Point(2,3)) === (2.0,3)
        @test Point{Float32}(P1) === Float32.(P1)
        @test promote(P1) === (P1,)
        @test Point{Int}(x=Int16(8),y=Int16(9)) === Point(8,9)
        @test Point{Float64}(CartesianIndex(8,9)) === Point(8.0,9.0)
        @test Point{Float64}(2, 3) === Point(2.0, 3.0)
        @test Point{Float64}((8,9)) === Point(8.0,9.0)
        P2 = Point(3,11)
        P3 = Point{Float32}(-2.1,11.8)
        @test promote(P1, P2) === (P1,Float64.(P2))
        @test promote(P1, P2, P3) === (P1,Float64.(P2),Float64.(P3))
        @test convert(Point, P1) === P1
        @test convert(Point{Float32}, P1) === Float32.(P1)
        for T in (Float32, Int16)
            pnt = Point{T}(1,2)
            @test @inferred(zero(Point{T})) === Point(zero(T),zero(T))
            @test @inferred(one(Point{T})) === Point(one(T),one(T))
            @test @inferred(oneunit(Point{T})) === Point(oneunit(T),oneunit(T))
            @test @inferred(zero(pnt)) === Point(zero(T),zero(T))
            @test @inferred(one(pnt)) === Point(one(T),one(T))
            @test @inferred(oneunit(pnt)) === Point(oneunit(T),oneunit(T))
        end
        @test ntuple(i -> Point(3,5)[i], 2) == (3,5)
        # round
        @test round(Point(1.2,-0.7)) === Point(1.0,-1.0)
        @test round(Point(1,-7)) === Point(1,-7)
        if VERSION < v"1.11.0-beta1"
            @test_throws Exception round(Float32, Point(1,-7))
        else
            @test round(Float32, Point(1,-7)) === Point{Float32}(1,-7)
        end
        @test round(Int32, Point(1,-7)) === Point{Int32}(1,-7)
        @test round(Int, Point(1.2,-0.7)) === Point(1,-1)
        @test round(Int, Point(1,-4)) === Point(1,-4)
        @test round(Point, P1) === Point(map(round, Tuple(P1)))
        @test round(Point{Int}, Point(1.6,-0.7)) === Point(2,-1)
        @test round(Int, Point(1.2,-0.7)) === Point(1,-1)
        @test round(Point{Int}, Point(1.6,-0.7)) === Point(2,-1)
        @test round(Point(1.5, 2.5)) === Point(2.0, 2.0)
        @test round(Point(1.5, 2.5), RoundNearestTiesUp) === Point(2.0, 3.0)
        if VERSION < v"1.11.0-beta1"
            @test_throws Exception round(Float32, Point(1.5e0, 2.5e0), RoundNearestTiesUp)
        else
            @test round(Float32, Point(1.5e0, 2.5e0), RoundNearestTiesUp) === Point{Float32}(2, 3)
        end
        # floor
        @test floor(Point(-1,3)) === Point(-1,3)
        @test floor(Int, Point(-1,3)) === Point(-1,3)
        @test floor(Point{Int}, Point(-1,3)) === Point(-1,3)
        @test floor(Int16, Point(-1,3)) === Point{Int16}(-1,3)
        @test floor(Point{Int16}, Point(-1,3)) === Point{Int16}(-1,3)
        @test floor(Point(-1.7,3.2)) === Point(-2.0,3.0)
        @test floor(Int, Point(-1.7,3.2)) === Point(-2,3)
        # ceil
        @test ceil(Point(-1,3)) === Point(-1,3)
        @test ceil(Int, Point(-1,3)) === Point(-1,3)
        @test ceil(Point{Int}, Point(-1,3)) === Point(-1,3)
        @test ceil(Int16, Point(-1,3)) === Point{Int16}(-1,3)
        @test ceil(Point{Int16}, Point(-1,3)) === Point{Int16}(-1,3)
        @test ceil(Point(-1.7,3.2)) === Point(-1.0,4.0)
        @test ceil(Int, Point(-1.7,3.2)) === Point(-1,4)
        # other methods
        @test clamp(Point(-1.1, 6.3), BoundingBox(1,4,1,5)) === Point(1.0,5.0)
        @test hypot(P1) == hypot(P1.x, P1.y)
        @test atan(P1) == atan(P1.y, P1.x)
        @test distance(P1, Point(0,0)) == hypot(P1)
        @test distance(Point(0x02,0x05), Point(0x00,0x00)) == hypot(0x02,0x05)
    end
    @testset "Weighted points" begin
        types = (Float32, Float64)
        for T1 in types, T2 in types
            @test promote_type(WeightedPoint{T1}, WeightedPoint{T2}) === WeightedPoint{promote_type(T1,T2)}
        end
        P2 = WeightedPoint(x=0.1, y=-6, w=0.1)
        @test eltype(P2) == Float64
        @test WeightedPoint(P2) === P2
        @test WeightedPoint{eltype(P2)}(P2) === P2
        @test WeightedPoint{Float64}(P2) === P2
        @test WeightedPoint(P2.w, P2.x, P2.y) === P2
        @test WeightedPoint{Float32}(P2) === Float32.(P2)
        @test WeightedPoint(Point(3.1,4.2)) === WeightedPoint(w=1, x=3.1, y=4.2)
    end
    @testset "Bounding boxes" begin
        types = (Int8, Int32, Int64, Float32, Float64)
        for T1 in types, T2 in types
            @test promote_type(BoundingBox{T1}, BoundingBox{T2}) === BoundingBox{promote_type(T1,T2)}
        end
        B = BoundingBox(2,3,4,5)
        δ = 0.1
        @test eltype(B) === Int
        @test BoundingBox(B) === B
        @test BoundingBox{eltype(B)}(B) === B
        @test BoundingBox{Float32}(B) ===
            BoundingBox{Float32}(B.xmin, B.xmax, B.ymin, B.ymax)
        @test BoundingBox(2,3,4,5.0) === BoundingBox(2.0,3.0,4.0,5.0)
        @test convert(BoundingBox, B) === B
        @test convert(BoundingBox{eltype(B)}, B) === B
        @test convert(BoundingBox{Int16}, B) === BoundingBox{Int16}(B)
        # Construct a bounding-box from 4-tuple.
        @test BoundingBox((2,3,4,5)) === BoundingBox(2,3,4,5)
        @test BoundingBox{Float64}((2,3,4,5)) === BoundingBox(2.0,3.0,4.0,5.0)
        @test convert(BoundingBox, (2,3,4,5)) === BoundingBox(2,3,4,5)
        @test convert(BoundingBox{Float64}, (2,3,4,5)) ===
            BoundingBox(2.0,3.0,4.0,5.0)
        # Construct a bounding-box from 2 points.
        @test BoundingBox(Point(2,3),Point(4,5)) === BoundingBox(2,4,3,5)
        @test BoundingBox{Int16}(Point(2,3),Point(4,5)) ===
            BoundingBox{Int16}(2,4,3,5)
        @test BoundingBox((Point(2,3),Point(4,5))) === BoundingBox(2,4,3,5)
        @test BoundingBox{Int16}((Point(2,3),Point(4,5))) ===
            BoundingBox{Int16}(2,4,3,5)
        @test convert(BoundingBox,(Point(2,3),Point(4,5))) ===
            BoundingBox(2,4,3,5)
        @test convert(BoundingBox{Int16}, (Point(2,3),Point(4,5))) ===
            BoundingBox{Int16}(2,4,3,5)
        # Construct a bounding-box from 2 Cartesian indices.
        @test BoundingBox(CartesianIndex(2,3),CartesianIndex(4,5)) ===
            BoundingBox(2,4,3,5)
        @test BoundingBox{Int16}(CartesianIndex(2,3),CartesianIndex(4,5)) ===
            BoundingBox{Int16}(2,4,3,5)
        @test BoundingBox((CartesianIndex(2,3),CartesianIndex(4,5))) ===
            BoundingBox(2,4,3,5)
        @test BoundingBox{Int16}((CartesianIndex(2,3),CartesianIndex(4,5))) ===
            BoundingBox{Int16}(2,4,3,5)
        # Construct a bounding-box from 2 2-tuple.
        @test BoundingBox((2,3),(4,5)) === BoundingBox(2,4,3,5)
        @test BoundingBox{Int16}((2,3),(4,5)) === BoundingBox{Int16}(2,4,3,5)
        # Construct a bounding-box by keywords.
        @test BoundingBox(xmin=2,ymin=3,xmax=4,ymax=5) ===
            BoundingBox(2,4,3,5)
        @test BoundingBox{Float64}(xmin=2,ymin=3,xmax=4,ymax=5) ===
            BoundingBox(2.0,4.0,3.0,5.0)
        # Construct a bounding-box from CartesianIndices instance.
        R = CartesianIndices((2:4,3:5))
        @test BoundingBox(R) === BoundingBox(2,4,3,5)
        @test BoundingBox{Int16}(R) === BoundingBox{Int16}(2,4,3,5)
        @test convert(BoundingBox, R) === BoundingBox(2,4,3,5)
        @test convert(BoundingBox{Int16}, R) === BoundingBox{Int16}(2,4,3,5)
        # Construct a bounding-box from 2 ranges.
        R = (2:4, 3:5)
        @test BoundingBox(R...) === BoundingBox(2,4,3,5)
        @test BoundingBox{Int16}(R...) === BoundingBox{Int16}(2,4,3,5)
        @test BoundingBox(R) === BoundingBox(2,4,3,5)
        @test BoundingBox{Int16}(R) === BoundingBox{Int16}(2,4,3,5)
        @test convert(BoundingBox,R) === BoundingBox(2,4,3,5)
        @test convert(BoundingBox{Int16},R) === BoundingBox{Int16}(2,4,3,5)

        @test first(B) === Point(B.xmin,B.ymin)
        @test last(B) === Point(B.xmax,B.ymax)
        @test isempty(typemax(BoundingBox{Int})) == false
        @test isempty(typemin(BoundingBox{Float64})) == true
        @test BoundingBox{Float32}(nothing) === typemin(BoundingBox{Float32})
        @test size(BoundingBox{Int16}(nothing)) === (0,0)
        @test ntuple(i -> BoundingBox(3:7,5:8)[i], 4) == (3,7,5,8)
        Rx, Ry = 3:7, -6:-3
        B = BoundingBox(Rx,Ry)
        @test size(B) === (length(Rx),length(Ry))
        @test size(BoundingBox{Int16}(B)) === size(B)
        @test axes(B) === (Rx,Ry)
        @test axes(BoundingBox{Int16}(B)) === (Rx,Ry)
        for k in (1,2)
            @test size(B, k) === size(B)[k]
            @test axes(B, k) === axes(B)[k]
            @test size(BoundingBox{Int16}(B), k) === size(B)[k]
            @test axes(BoundingBox{Int16}(B), k) === axes(B)[k]
        end
        @test size(B,40) === 1
        @test_throws ErrorException size(B, 0)
        B1 = BoundingBox(2,3,4,5)
        B2 = BoundingBox(1.1,3.2,4.5,5.8)
        B3 = BoundingBox{Float32}(1.1,3.2,4.5,5.8)
        @test promote(B1, B2) === (Float64.(B1),B2)
        @test promote(B1, B2, B3) === (Float64.(B1),B2,Float64.(B3))
        @test_throws MethodError BoundingBox(ones(5,7))
        @test BoundingBox(-2:6,8:11) === BoundingBox(-2,6, 8,11)
        @test BoundingBox((2:4,-1:7)) === BoundingBox(2,4, -1,7)
        @test CartesianIndices(BoundingBox(2:4,-1:7)) ===
            CartesianIndices((2:4,-1:7))
        @test get_axis_bounds(9:13) === (9,13)
        @test get_axis_bounds(13:9) === (13,12)
        @test get_axis_bounds(13:-1:9) === (9,13)
        @test_throws ArgumentError get_axis_bounds(13:-2:9)
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
        @test round(BoundingBox{Int16}, BoundingBox(1.1,2.1,-3.6,7.7)) ===
            BoundingBox{Int16}(1,2,-4,8)
        @test round(Int, BoundingBox(1,2,-3,7)) === BoundingBox(1,2,-3,7)
        @test round(Int16, BoundingBox(1.1,2.1,-3.6,7.7)) ===
            BoundingBox{Int16}(1,2,-4,8)
        @test round(Int16, BoundingBox(1,2,-4,8)) ===
            BoundingBox{Int16}(1,2,-4,8)
        if VERSION < v"1.11.0-beta1"
            @test_throws Exception round(Float32, BoundingBox{Int16}(1,2,-4,8))
            @test_throws Exception round(Float32, BoundingBox(1.1,2.7,-4.6,8.3))
        else
            @test round(Float32, BoundingBox{Int16}(1,2,-4,8)) === BoundingBox{Float32}(1,2,-4,8)
            @test round(Float32, BoundingBox(1.1,2.7,-4.6,8.3)) === BoundingBox{Float32}(1,3,-5,8)
        end
        @test round(BoundingBox(1.5, 2.5, 3.5, 9.9)) === BoundingBox(2.0, 2.0, 4.0, 10.0)
        @test round(BoundingBox(1.5, 2.5, 3.5, 9.9), RoundNearestTiesUp) ===
            BoundingBox(2.0, 3.0, 4.0, 10.0)
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

        @test (Point(2,4) ∈ BoundingBox(1,2,3,4)) == true
        @test (Point(3,4) ∈ BoundingBox(1,2,3,4)) == false
        @test (Point(2,5) ∈ BoundingBox(1,2,3,4)) == false
        @test (Point(3,5) ∈ BoundingBox(1,2,3,4)) == false
        @test (WeightedPoint(1,2,4) ∈ BoundingBox(1,2,3,4)) == true
        @test (WeightedPoint(1,3,4) ∈ BoundingBox(1,2,3,4)) == false
        @test (WeightedPoint(1,2,5) ∈ BoundingBox(1,2,3,4)) == false
        @test (WeightedPoint(1,3,5) ∈ BoundingBox(1,2,3,4)) == false

        @test (BoundingBox(1,-2,3,4) ⊆ BoundingBox{Float32}(nothing)) == true
        @test (BoundingBox(1,-2,3,4) ⊆ BoundingBox(1,2,3,4)) == true
        @test (BoundingBox(1,2,3,4) ⊆ BoundingBox(1,2,3,4)) == true
        @test (BoundingBox(1,2,3,4) ⊆ BoundingBox(1,4,0,3)) == false
        @test (BoundingBox(0,2,3,4) ⊆ BoundingBox(1,2,3,5)) == false
        @test_deprecated (BoundingBox(0,2,3,4) ∈ BoundingBox(1,2,3,5)) == false

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
    @testset "Arithmetic" begin
        @test half(Float64) == 0.5
        @test half(Float32) == one(Float32)/2
        @test 2*Point(3,4) === Point(6,8)
        @test Point(3,4)*3 === Point(9,12)
        @test 3\Point(3,4) === Point(3/3,4/3)
        @test Point(3,4)/3 === Point(3/3,4/3)
        @test +Point(2,-9) === Point(2,-9)
        @test -Point(2,-9) === Point(-2,9)
        @test +BoundingBox(2,3,4,5) === BoundingBox(2,3,4,5)
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

@testset "Bounding-box algorithm" begin
    if true
        # Exhaustively test all possibilities.
        A = Array{Bool,2}(undef, (5,4))
        n = (one(UInt64) << length(A)) # number of possibilities
        for bits in 0:n-1
            unpack_bits!(A, bits)
            @test BoundingBox(A) === naive_bounding_box(A)
        end
    else
        # Semi-exhaustive testing of the bounding box algorithm.
        A = zeros(Bool, (7,8))
        fill!(A,false)[1,3] = true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,true)[2:end-1,2:end-1] .= false
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[1,1] = true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[1,end] = true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[end,1] = true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[end,end] = true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[3,4] = true
        @test BoundingBox(A) === naive_bounding_box(A)
        A[2,7] = true
        @test BoundingBox(A) === naive_bounding_box(A)
        A[5,end] = true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[end,4] = true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[1,:] .= true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[:,1] .= true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[end,:] .= true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[:,end] .= true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[1,2:end-1] .= true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[2:end-1,1] .= true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[end,2:end-1] .= true
        @test BoundingBox(A) === naive_bounding_box(A)
        fill!(A,false)[2:end-1,end] .= true
        @test BoundingBox(A) === naive_bounding_box(A)
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
        @testset "T = $T" for T in types
            @test_throws MethodError T(A)
            @test_throws MethodError T.(A)
        end
        for T1 in types, T2 in types
            @test promote_type(AffineTransform{T1}, AffineTransform{T2}) ===
                AffineTransform{promote_type(T1,T2)}
        end
        @test AffineTransform(A) === A
        @test AffineTransform{bare_type(A)}(A) === A
        for G in (I, A, B)
            for T in types
                H = @inferred AffineTransform{T}(G)
                @test typeof(H) <: AffineTransform{T}
                @test bare_type(H) === T
                @test real_type(H) === T
                @test floating_point_type(H) === T
                @test factors_type(H) === T
                @test offsets_type(H) === T
                @test (H === G) == (bare_type(H) === bare_type(G))
                @test Tuple(H) ≈ map(T, Tuple(G))
                H = @inferred convert(AffineTransform{T}, G)
                @test typeof(H) <: AffineTransform{T}
                @test bare_type(H) === T
                @test real_type(H) === T
                @test floating_point_type(H) === T
                @test factors_type(H) === T
                @test offsets_type(H) === T
                @test (H === G) == (bare_type(H) === bare_type(G))
                @test Tuple(H) ≈ map(T, Tuple(G))
                H = @inferred convert_bare_type(T, G)
                @test typeof(H) <: AffineTransform{T}
                @test bare_type(H) === T
                @test real_type(H) === T
                @test floating_point_type(H) === T
                @test factors_type(H) === T
                @test offsets_type(H) === T
                @test (H === G) == (bare_type(H) === bare_type(G))
                @test Tuple(H) ≈ map(T, Tuple(G))
                H = @inferred convert_real_type(T, G)
                @test typeof(H) <: AffineTransform{T}
                @test bare_type(H) === T
                @test real_type(H) === T
                @test floating_point_type(H) === T
                @test factors_type(H) === T
                @test offsets_type(H) === T
                @test (H === G) == (bare_type(H) === bare_type(G))
                @test Tuple(H) ≈ map(T, Tuple(G))
                H = @inferred convert_floating_point_type(T, G)
                @test typeof(H) <: AffineTransform{T}
                @test bare_type(H) === T
                @test real_type(H) === T
                @test floating_point_type(H) === T
                @test factors_type(H) === T
                @test offsets_type(H) === T
                @test (H === G) == (bare_type(H) === bare_type(G))
                @test Tuple(H) ≈ map(T, Tuple(G))
            end
        end
    end

    @testset "identity" begin
        @test isone(det(I))
        @test inv(I) ≈ I
        for v in vectors
            @test I(v) ≈ v
        end
    end

    @testset "apply" begin
        for G in (I, A, B),
            v in vectors
            @test G(v...) ≈ G(v)
            @test G*v ≈ G(v)
            @test G(Point(v)) ≈ Point(G(v))
            @test G*Point(v) ≈ Point(G(v))
        end
    end

    @testset "composition" begin
        for G in (I, B, A)
            @test G === @inferred compose(G)
            for H in (A, B)
                @test G*H ≈ compose(G,H)
                @test G⋅H ≈ compose(G,H)
                @test G∘H ≈ compose(G,H)
                for v in vectors
                    @test (G*H)(v) ≈ G(H(v))
                end
            end
        end
        for T1 in types, T2 in types
            T = promote_type(T1, T2)
            @test bare_type(AffineTransform{T1}(A)*AffineTransform{T2}(B)) == T
        end
        for v in vectors
            @test compose(A,B,C)(v) ≈ A(B(C(v)))
            @test (A*B*C)(v) ≈ A(B(C(v)))
            @test compose(A,C,B,C)(v) ≈ A(C(B(C(v))))
            @test (A*C*B*C)(v) ≈ A(C(B(C(v))))
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
            @test det(inv(M)) ≈ 1/det(M)
            @test M/M ≈ M*inv(M)
            @test M\M ≈ inv(M)*M
            @test M\M ≈ I
            @test M/M ≈ I
            for v in vectors
                @test M(inv(M)(v)) ≈ v
                @test inv(M)(M(v)) ≈ v
                @test (M\M)(v) ≈ v
                @test (M/M)(v) ≈ v
            end
        end
        for T1 in types, T2 in types
            T = promote_type(T1, T2)
            @test bare_type(AffineTransform{T1}(A)/AffineTransform{T2}(B)) == T
            @test bare_type(AffineTransform{T1}(A)\AffineTransform{T2}(B)) == T
        end
    end

    @testset "scale" begin
        for M in (A, B, C), α in scales, v in vectors
            @test (α*M)(v) ≈ α.*M(v)
            @test (M*α)(v) ≈ M(α.*v)
            for T in types
                @test bare_type(T(α)*M) === promote_type(T, bare_type(M))
                @test bare_type(M*T(α)) === promote_type(T, bare_type(M))
                H = @inferred AffineTransform{T}(M)
                @test bare_type(T(α)*H) === bare_type(H)
                @test bare_type(H*T(α)) === bare_type(H)
            end
        end
    end

    @testset "translation" begin
        for M in (B, A), t in vectors, v in vectors
            @test translate(t, M)(v) ≈ t .+ M(v)
            @test translate(t, M)(v) ≈ (t + M)(v)
            @test translate(M, t)(v) ≈ M(v .+ t)
            @test translate(M, t)(v) ≈ (M + t)(v)
        end
        for v in vectors
            @test A - v ≈ A + (-v[1], -v[2])
            @test A + Point(v) ≈ translate(A, v...)
            @test Point(v) + A ≈ translate(v..., A)
            @test A - Point(v) ≈ A - v
        end
        for G in (A, B, C), v in vectors, T in types
            @test bare_type(T.(v) + G) === promote_type(T, bare_type(G))
            @test bare_type(G + T.(v)) === promote_type(T, bare_type(G))
            H = @inferred AffineTransform{T}(G)
            @test bare_type(T.(v) + H) === T
            @test bare_type(H + T.(v)) === T
        end
    end

    @testset "rotation" begin
        for θ in angles
            R = @inferred rotate(+θ, I)
            Q = @inferred rotate(-θ, I)
            @test R*Q ≈ I
            @test Q*R ≈ I
            for G in (A, B, C), v in vectors
                @test rotate(θ, G)(v) ≈ (R*G)(v)
                @test rotate(G, θ)(v) ≈ (G*R)(v)
            end
        end
        for G in (A, B, C), θ in angles, T in types
            @test bare_type(rotate(T(θ), G)) === promote_type(T, bare_type(G))
            @test bare_type(rotate(G, T(θ))) === promote_type(T, bare_type(G))
            H = @inferred AffineTransform{T}(G)
            @test bare_type(rotate(T(θ), H)) === T
            @test bare_type(rotate(H, T(θ))) === T
        end
    end

    @testset "intercept" begin
        for M in (I, A, B)
            x, y = intercept(M)
            @test M(x, y) ≈ (0,0) atol=16*eps(Float64)
            P = intercept(Point, M)
            @test M*P ≈ Point(0,0) atol=16*eps(Float64)
        end
    end

    @testset "show" begin
        for M in (I, A, B, C)
            @test occursin(r"\bAffineTransform\b", string(M))
        end
    end
end

end # module
