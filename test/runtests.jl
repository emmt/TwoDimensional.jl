module TwoDimensionalTests

using TwoDimensional
using TwoDimensional: PointLike, compose, get_x, get_y, get_xy, factors_type, offsets_type
using TwoDimensional: get_axis_bounds
using Test, LinearAlgebra, TypeUtils
import Base.MathConstants: φ

struct NamedPoint{T} <: AbstractPoint{T}
    name::String
    x::T
    y::T
end

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

if VERSION < v"1.2"
    # `mapreduce` with multiple iterators requires Julia 1.2 or later.
    Base.mapreduce(f, op, a, b; kwds...) = reduce(op, map(f, a, b); kwds...)
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
    return BoundingBox(imin:imax, jmin:jmax)
end

@testset "TwoDimensional" begin
    # Constructors.
    @testset "Points ($T)" for T in (Int16, Int, Float32)
        @assert !(T === Float64) # this is assumed by the tests
        @test TwoDimensional.Point2D{Float32} === Point{Float32}
        @test TwoDimensional.AbstractPoint2D{Float32} === AbstractPoint{Float32}
        xy = 1, 2
        pnt = @inferred Point(map(T, xy)...)
        pnt_f64 = @inferred Point{Float64}(xy)
        pnt_misc = @inferred NamedPoint{T}("yet another point", xy...)
        @test pnt isa Point{T}
        @test typeof(pnt) === Point{T}
        @test coord_type(pnt) === coord_type(typeof(pnt)) === T
        @test eltype(pnt) === eltype(typeof(pnt)) === T
        @test eltype(pnt_misc) === eltype(typeof(pnt_misc)) === T
        @test Tuple(pnt) === (pnt...,)
        @test pnt === Point(pnt...,)
        @test length(pnt) === length(getfield(pnt, 1)) == 2
        @test Set(Base.propertynames(pnt)) == Set((:x, :y))
        @test firstindex(pnt) === 1
        @test lastindex(pnt) === 2
        @test Base.IteratorSize(pnt) === Base.IteratorSize(typeof(pnt)) === Base.HasLength()
        @test Base.IteratorEltype(pnt) === Base.IteratorEltype(typeof(pnt)) === Base.HasEltype()
        x, y = pnt
        r, θ = hypot(x, y), atan(y, x)
        @test (x, y) === map(T, xy)
        @test (x, y) === @inferred Tuple(pnt)
        @test (x, y) === (pnt.x, pnt.y)
        @test (x, y) === (pnt[1], pnt[2])
        @test (x, y) === (first(pnt), last(pnt))
        @test_throws BoundsError pnt[0]
        @test_throws BoundsError pnt[3]
        @test_throws KeyError pnt.vec
        @test occursin(r"^Point\b", string(pnt))
        @test pnt === @inferred Point(map(T, xy))
        @test pnt === @inferred Point(pnt...)
        @test pnt === @inferred Point(Tuple(pnt)...)
        @test pnt === @inferred Point(Tuple(pnt))
        @test pnt === @inferred Point{T}(xy...)
        @test pnt === @inferred Point{T}(xy)
        @test pnt === @inferred Point{T}(Int8(x), Int16(y))
        @test pnt === @inferred Point{T}((Int8(x), Int16(y)))
        @test pnt === @inferred Point{T}(pnt_f64)
        @test pnt === @inferred Point{T}(pnt_misc)
        @test pnt === @inferred Point(pnt_misc)
        @test pnt === @inferred Point{T}(; x = x, y = y)
        @test pnt === @inferred Point(; x = T(x), y = T(y))
        @test pnt === @inferred Point{T}(CartesianIndex(xy))
        if T <: Int
            @test pnt === @inferred Point(CartesianIndex(xy))
            @test CartesianIndex(xy) === CartesianIndex(pnt)
            @test CartesianIndex(xy) === CartesianIndex(pnt_misc)
        end
        # Conversions.
        @test pnt === @inferred convert(Point, pnt)
        @test pnt === @inferred convert(Point{T}, pnt)
        @test pnt === @inferred convert(Point{T}, pnt_f64)
        @test pnt === @inferred convert(Point{T}, xy)
        @test pnt === @inferred convert(Point{T}, (Int8(x), Int16(y)))
        @test pnt === @inferred convert(Point{T}, CartesianIndex(xy))
        @test convert_coord_type(coord_type(pnt_f64), pnt) === pnt_f64
        @test convert_coord_type(coord_type(pnt), pnt_f64) === pnt
        if T <: AbstractFloat
            @test @inferred(float(pnt)) === pnt
            @test @inferred(float(typeof(pnt))) === typeof(pnt)
        elseif T <: Integer
            @test @inferred(float(pnt)) === @inferred(Point(float(pnt.x), float(pnt.y)))
            @test @inferred(float(typeof(pnt))) === Point{float(T)}
        end
        # Functions specific to points.
        @test hypot(pnt) ≈ r
        @test norm(pnt) ≈ r
        @test abs(pnt) ≈ r
        @test abs2(pnt) ≈ r^2
        @test abs2(pnt) ≈ dot(pnt, pnt)
        @test atan(pnt) ≈ θ
        @test Point(; r = r, θ = θ) ≈ pnt
        let a = Point(-1, 2), b = Point(3.0, 7.0)
            @test a ⋅ b ≈ a.x*b.x + a.y*b.y
            @test a ⋅ b === dot(a, b)
            @test a ⋅ b === b ⋅ a
            @test a * b ≈ a.x*b.y - a.y*b.x
            @test a * b === cross(a, b)
            @test a * b === -(b * a)
        end
    end

    @testset "Rectangles ($T)" for T in (Int16, Int, Float32)
        @assert !(T === Float64) # this is assumed by the tests
        @test TwoDimensional.Rectangle2D{Float32} === Rectangle{Float32}
        start, stop = (-1, 2), (3, 4)
        @assert start[1] < stop[1] && start[2] < stop[2] # this is assumed by the tests
        rec = @inferred Rectangle(T.(start), T.(stop))
        rec_f64 = @inferred Rectangle{Float64}(start, stop) # same with other type
        @test rec isa Rectangle{T}
        @test typeof(rec) === Rectangle{T}
        @test coord_type(rec) === coord_type(typeof(rec)) === T
        @test eltype(rec) === eltype(typeof(rec)) === Point{T}
        @test Tuple(rec) === (rec...,)
        @test rec === Rectangle(rec...,)
        @test length(rec) === length(getfield(rec, 1)) == 2
        @test Set(Base.propertynames(rec)) == Set((:x0, :x1, :y0, :y1, :start, :stop))
        @test firstindex(rec) === 1
        @test lastindex(rec) === 2
        @test Base.IteratorSize(rec) === Base.IteratorSize(typeof(rec)) === Base.HasLength()
        @test Base.IteratorEltype(rec) === Base.IteratorEltype(typeof(rec)) === Base.HasEltype()
        (x0, y0), (x1, y1) = (rec.x0, rec.y0), (rec.x1, rec.y1)
        @test x0 ≤ x1 && y0 ≤ y1
        @test rec === Rectangle((x0, y0), (x1, y1))
        @test rec === Rectangle((x1, y0), (x0, y1))
        @test rec === Rectangle((x0, y1), (x1, y0))
        @test rec === Rectangle((x1, y1), (x0, y0))
        @test (Point(x0, y0), Point(x1, y1)) === @inferred Tuple(rec)
        @test (Point(x0, y0), Point(x1, y1)) === (rec.start, rec.stop)
        @test (Point(x0, y0), Point(x1, y1)) === (rec[1], rec[2])
        @test (Point(x0, y0), Point(x1, y1)) === (first(rec), last(rec))
        @test_throws BoundsError rec[0]
        @test_throws BoundsError rec[3]
        @test_throws KeyError rec.vec
        @test occursin(r"^Rectangle\b", string(rec))
        @test rec === @inferred Rectangle(rec...)
        @test rec === @inferred Rectangle(Tuple(rec)...)
        @test rec === @inferred Rectangle(Tuple(rec))
        @testset "start = $start, stop = $stop" for (start, stop) in (
            ((x0,y0), (x1,y1)),
            (Int8.((x0,y0)), Float64.((x1,y1))),
            ((Int8(x0),Int16(y0)), (Float64(x1),Float32(y1))),
            (Point(x0,y0), Point(x1,y1)),
            (Point{Int8}(x0,y0), Point{Float64}(x1,y1)),
            (CartesianIndex(Int(x0),Int(y0)), CartesianIndex(Int(x1),Int(y1))),
            )
            @test rec === @inferred Rectangle{T}(start, stop)
            @test rec === @inferred Rectangle{T}((start, stop))
            @test rec === @inferred Rectangle{T}(; x0 = start[1], y0 = start[2], x1 = stop[1], y1 = stop[2])
            @test rec === @inferred Rectangle{T}(; start = start, stop = stop)
            if start isa Point && stop isa Point
                @test rec === @inferred Rectangle(Point{T}(start), Point{T}(stop))
                @test rec === @inferred Rectangle(; start = Point{T}(start), stop = Point{T}(stop))
            end
            if start isa CartesianIndex && stop isa CartesianIndex
                @test rec === @inferred Rectangle(T.(Tuple(start)), T.(Tuple(stop)))
                @test rec === @inferred Rectangle(; start = T.(Tuple(start)), stop = T.(Tuple(stop)))
                @test rec === @inferred Rectangle(; x0 = T(start[1]), y0 = T(start[2]), x1 = T(stop[1]), y1 = T(stop[2]))
                if T <: Int
                    @test rec === @inferred Rectangle(start, stop)
                end
            end
        end
        # Other constructors.
        @test rec === @inferred Rectangle(BoundingBox(rec.start, rec.stop))
        @test rec === @inferred Rectangle{T}(BoundingBox{Float64}(rec.start, rec.stop))
        pnt = Point(3.1,2.7)
        small_rec = @inferred(Rectangle(pnt))
        @test Tuple(small_rec) === (pnt, pnt)
        @test coord_type(small_rec) === coord_type(pnt)
        pnt = NamedPoint{Float32}("some point", 3.1, -2.7)
        small_rec = @inferred(Rectangle(pnt))
        @test Tuple(small_rec) === (Point(pnt), Point(pnt))
        @test coord_type(small_rec) === coord_type(pnt)
        # Conversions.
        @test rec === @inferred convert(Rectangle, rec)
        @test rec === @inferred convert(Rectangle{T}, rec)
        @test rec === @inferred convert(Rectangle{T}, rec_f64)
        @test rec === @inferred convert(Rectangle{T}, Tuple(rec_f64))
        @test convert_coord_type(coord_type(rec_f64), rec) === rec_f64
        @test convert_coord_type(coord_type(rec), rec_f64) === rec
        if T <: AbstractFloat
            @test @inferred(float(rec)) === rec
            @test @inferred(float(typeof(rec))) === typeof(rec)
        elseif T <: Integer
            @test @inferred(float(rec)) === @inferred(Rectangle(float(rec.start), float(rec.stop)))
            @test @inferred(float(typeof(rec))) === Rectangle{float(T)}
        end
        # Functions specific to rectangles.
        @test area(rec) == (rec.x1 - rec.x0)*(rec.y1 - rec.y0)
    end

    @testset "Circles ($T)" for T in (Int16, Int, Float32)
        @assert !(T === Float64) # this is assumed by the tests
        @test TwoDimensional.Circle2D{Float32} === Circle{Float32}
        c, r = (-1, 2), 3
        circ = @inferred Circle(T.(c), T(r))
        circ_f64 = @inferred Circle{Float64}(c, r) # same with other type
        @test circ isa Circle{T}
        @test typeof(circ) === Circle{T}
        @test coord_type(circ) === coord_type(typeof(circ)) === T
        @test eltype(circ) === Union{Point{T}, T}
        @test Tuple(circ) === (circ...,)
        @test circ === Circle(circ...,)
        @test length(circ) == 2
        @test Set(Base.propertynames(circ)) == Set((:center, :radius, :diameter,))
        @test Base.IteratorSize(circ) === Base.IteratorSize(typeof(circ)) === Base.HasLength()
        @test Base.IteratorEltype(circ) === Base.IteratorEltype(typeof(circ)) === Base.HasEltype()
        c, r = circ
        (x, y), r′ = circ
        @test r′ === r && c.x === x && c.y === y
        @test circ.diameter == r + r
        @test circ === @inferred Circle(c, r)
        @test circ === @inferred Circle((x, y), r)
        @test circ === @inferred Circle((c, r))
        @test circ === @inferred Circle(((x, y), r))
        @test (Point(x, y), r) === @inferred Tuple(circ)
        @test (Point(x, y), r) === (circ...,)
        @test (Point(x, y), r) === (center(circ), radius(circ))
        @test (Point(x, y), r) === (circ.center, circ.radius)
        @test circ.diameter === diameter(circ) ≈ 2*circ.radius
        @test_throws MethodError circ[1]
        @test_throws KeyError circ.vec
        @test occursin(r"^Circle\b", string(circ))
        @test circ === @inferred Circle(circ...)
        @test circ === @inferred Circle(Tuple(circ)...)
        @test circ === @inferred Circle(Tuple(circ))
        @test circ === @inferred Circle{T}(c, r)
        @test circ === @inferred Circle(Point(c), r)
        @test circ === @inferred Circle{T}(Point(c), Float64(r))
        @test circ === @inferred Circle{T}(NamedPoint("center", x, y), r)
        @test circ === @inferred Circle(NamedPoint("center", x, y), r)
        if T === Int
            let c = CartesianIndex(x, y)
                @test circ === @inferred Circle(c, r)
                @test circ === @inferred Circle((c, r))
                @test circ === @inferred Circle(center=c, radius=r)
            end
        end
        # Conversions.
        @test circ === @inferred convert(Circle, circ)
        @test circ === @inferred convert(Circle{T}, circ)
        @test circ === @inferred convert(Circle{T}, circ_f64)
        @test circ === @inferred convert(Circle{T}, Tuple(circ_f64))
        @test convert_coord_type(coord_type(circ_f64), circ) === circ_f64
        @test convert_coord_type(coord_type(circ), circ_f64) === circ
        if T <: AbstractFloat
            @test @inferred(float(circ)) === circ
            @test @inferred(float(typeof(circ))) === typeof(circ)
        elseif T <: Integer
            @test @inferred(float(circ)) === @inferred(Circle(float(circ.center), float(circ.radius)))
            @test @inferred(float(typeof(circ))) === Circle{float(T)}
        end
        # Functions specific to circles.
        @test area(circ) ≈ π*r*r
    end

    @testset "Polygons ($T)" for T in (Int16, Int, Float32)
        @assert !(T === Float64) # this is assumed by the tests
        @test TwoDimensional.Polygon2D{Float32,Vector{Point{Float32}}} === Polygon{Float32,Vector{Point{Float32}}}
        pnts = ((-1, 2), (3, 4), (5, -6))
        poly = @inferred Polygon{T}(pnts)
        poly_f64 = @inferred Polygon{Float64}(pnts) # same with other type
        @test poly isa Polygon{T}
        @test typeof(poly) <: Polygon{T}
        @test coord_type(poly) === coord_type(typeof(poly)) === T
        @test eltype(poly) === eltype(typeof(poly)) === Point{T}
        @test vec(poly) == collect(map(Point{T}, pnts))
        @test vec(poly) === TwoDimensional.elements(poly)
        @test vec(poly) === TwoDimensional.vertices(poly)
        @test poly === Polygon(vec(poly))
        @test length(poly) === length(pnts)
        @test firstindex(poly) === 1
        @test lastindex(poly) === length(pnts)
        @test eachindex(poly) == firstindex(poly):lastindex(poly)
        @test size(poly) === (length(poly),)
        @test axes(poly) == (firstindex(poly):lastindex(poly),)
        @test keys(poly) === eachindex(poly)
        @test values(poly) === vec(poly)
        @test Base.IteratorSize(poly) === Base.IteratorSize(typeof(poly)) === Base.HasLength()
        @test Base.IteratorEltype(poly) === Base.IteratorEltype(typeof(poly)) === Base.HasEltype()
        @test_throws BoundsError poly[firstindex(poly) - 1]
        @test_throws BoundsError poly[lastindex(poly) + 1]
        @test Set(Base.propertynames(poly)) == Set((:vertices,))
        @test poly.vertices === vec(poly)
        @test_throws KeyError poly.non_existing_property
        @test occursin(r"^Polygon\b", string(poly))
        let arr = @inferred(collect(poly)),
            tup = (vec(poly)...,),
            arr_xy = map(p -> (p.x, p.y), vec(poly)),
            arr_xy_mix = map(p -> (p.x, Float64(p.y)), vec(poly))
            # Same values but not same object.
            @test arr !== vec(poly)
            @test arr == vec(poly)
            @test poly === @inferred Polygon(vec(poly))
            @test poly === @inferred Polygon{T}(vec(poly))
            @test poly ==  @inferred Polygon(arr)
            @test poly ==  @inferred Polygon(arr...)
            @test poly ==  @inferred Polygon((arr...,))
            @test poly ==  @inferred Polygon{T}(arr)
            @test poly ==  @inferred Polygon{T}(arr...)
            @test poly ==  @inferred Polygon{T}((arr...,))
            @test poly ==  @inferred Polygon(arr_xy)
            @test poly ==  @inferred Polygon(arr_xy...)
            @test poly ==  @inferred Polygon((arr_xy...,))
            @test poly ==  @inferred Polygon{T}(arr_xy)
            @test poly ==  @inferred Polygon{T}(arr_xy...)
            @test poly ==  @inferred Polygon{T}((arr_xy...,))
            @test poly ==  @inferred Polygon{T}(arr_xy_mix)
            @test poly ==  @inferred Polygon{T}(arr_xy_mix...)
            @test poly ==  @inferred Polygon{T}((arr_xy_mix...,))
        end
        if T <: Int
            let arr = map(CartesianIndex, vec(poly)), tup = (arr...,)
                @test poly == @inferred Polygon(arr)
                @test poly == @inferred Polygon{T}(arr)
                @test poly == @inferred Polygon(tup...)
                @test poly == @inferred Polygon{T}(tup...)
                @test poly == @inferred Polygon(tup)
                @test poly == @inferred Polygon{T}(tup)
            end
        end
        # Bounding-box of a polygon and conversely.
        box  = @inferred BoundingBox(poly)
        @test (box.xmin, box.xmax) === extrema(map(p -> p.x, vec(poly)))
        @test (box.ymin, box.ymax) === extrema(map(p -> p.y, vec(poly)))
        p = @inferred Polygon(box)
        @test p isa Polygon{T}
        @test length(p) == 4
        pnts = (Point(box.xmin, box.ymin), Point(box.xmax, box.ymin),
                Point(box.xmax, box.ymax), Point(box.xmin, box.ymax),)
        @test p[1] ∈ pnts
        @test p[2] ∈ pnts
        @test p[3] ∈ pnts
        @test p[4] ∈ pnts
    end

    @testset "BoundingBoxes ($T)" for T in (Int16, Int, Float32)
        @assert !(T === Float64) # this is assumed by the tests
        @test TwoDimensional.BoundingBox2D{Float32} === BoundingBox{Float32}
        start, stop = (-1, 2), (3, 4)
        @assert start[1] < stop[1] && start[2] < stop[2] # this is assumed by the tests
        box = @inferred BoundingBox(T.(start), T.(stop))
        box_f64 = @inferred BoundingBox{Float64}(start, stop)
        @test box isa BoundingBox{T}
        @test typeof(box) === BoundingBox{T}
        @test coord_type(box) === coord_type(typeof(box)) === T
        @test eltype(box) === eltype(typeof(box)) === Point{T}
        @test Tuple(box) === (box...,)
        @test box === BoundingBox(box...,)
        @test length(box) === length(getfield(box, 1)) == 2
        @test Set(Base.propertynames(box)) == Set((:xmin, :xmax, :ymin, :ymax, :start, :stop))
        @test firstindex(box) === 1
        @test lastindex(box) === 2
        @test Base.IteratorSize(box) === Base.IteratorSize(typeof(box)) === Base.HasLength()
        @test Base.IteratorEltype(box) === Base.IteratorEltype(typeof(box)) === Base.HasEltype()
        (xmin, ymin), (xmax, ymax) = (box.xmin, box.ymin), (box.xmax, box.ymax)
        @test false === isempty(BoundingBox((xmin, ymin), (xmax, ymax)))
        @test true  === isempty(BoundingBox((xmax, ymin), (xmin, ymax)))
        @test true  === isempty(BoundingBox((xmin, ymax), (xmax, ymin)))
        @test true  === isempty(BoundingBox((xmax, ymax), (xmin, ymin)))
        @test (Point(xmin, ymin), Point(xmax, ymax)) === @inferred Tuple(box)
        @test (Point(xmin, ymin), Point(xmax, ymax)) === (box.start, box.stop)
        @test (Point(xmin, ymin), Point(xmax, ymax)) === (box[1], box[2])
        @test (Point(xmin, ymin), Point(xmax, ymax)) === (first(box), last(box))
        @test_throws BoundsError box[0]
        @test_throws BoundsError box[3]
        @test_throws KeyError box.vec
        @test occursin(r"^BoundingBox\b", string(box))
        @test box === @inferred BoundingBox(box...)
        @test box === @inferred BoundingBox(Tuple(box)...)
        @test box === @inferred BoundingBox(Tuple(box))
        @testset "start = $start, stop = $stop" for (start, stop) in (
            ((xmin,ymin), (xmax,ymax)),
            (Int8.((xmin,ymin)), Float64.((xmax,ymax))),
            ((Int8(xmin),Int16(ymin)), (Float64(xmax),Float32(ymax))),
            (Point(xmin,ymin), Point(xmax,ymax)),
            (Point{Int8}(xmin,ymin), Point{Float64}(xmax,ymax)),
            (CartesianIndex(Int(xmin),Int(ymin)), CartesianIndex(Int(xmax),Int(ymax))),
            )
            @test box === @inferred BoundingBox{T}(start, stop)
            @test box === @inferred BoundingBox{T}((start, stop))
            @test box === @inferred BoundingBox{T}(; xmin = start[1], ymin = start[2], xmax = stop[1], ymax = stop[2])
            @test box === @inferred BoundingBox{T}(; start = start, stop = stop)
            if start isa Point && stop isa Point
                @test box === @inferred BoundingBox(Point{T}(start), Point{T}(stop))
                @test box === @inferred BoundingBox(; start = Point{T}(start), stop = Point{T}(stop))
            elseif start isa CartesianIndex && stop isa CartesianIndex
                @test box === @inferred BoundingBox(T.(Tuple(start)), T.(Tuple(stop)))
                @test box === @inferred BoundingBox(; start = T.(Tuple(start)), stop = T.(Tuple(stop)))
                @test box === @inferred BoundingBox(; xmin = T(start[1]), ymin = T(start[2]), xmax = T(stop[1]), ymax = T(stop[2]))
                if T <: Int
                    @test box === @inferred BoundingBox(start, stop)
                end
            end
        end
        # Other constructors.
        xrng, yrng = Int(xmin):Int(xmax), Int(ymin):Int(ymax)
        @test box === @inferred BoundingBox{T}(xrng, yrng)
        @test box === @inferred BoundingBox{T}((xrng, yrng))
        @test box === @inferred BoundingBox{T}(; xmin = first(xrng), xmax = last(xrng), ymin = first(yrng), ymax = last(yrng))
        @test box === @inferred BoundingBox(Rectangle(box.start, box.stop))
        @test @inferred(BoundingBox(box.start)) === @inferred(BoundingBox(box.start, box.start))
        pnt = Point(3.1,2.7)
        smallest_box = @inferred(BoundingBox(pnt))
        @test Tuple(smallest_box) === (pnt, pnt)
        @test !isempty(smallest_box)
        # Conversions.
        @test box === @inferred convert(BoundingBox, box)
        @test box === @inferred convert(BoundingBox{T}, box)
        @test box === @inferred convert(BoundingBox{T}, box_f64)
        @test box === @inferred convert(BoundingBox{T}, Tuple(box_f64))
        @test convert_coord_type(coord_type(box_f64), box) === box_f64
        @test convert_coord_type(coord_type(box), box_f64) === box
        if T <: AbstractFloat
            @test @inferred(float(box)) === box
            @test @inferred(float(typeof(box))) === typeof(box)
        elseif T <: Integer
            @test @inferred(float(box)) === @inferred(BoundingBox(float(box.start), float(box.stop)))
            @test @inferred(float(typeof(box))) === BoundingBox{float(T)}
        end
        # Functions specific to bounding-boxes.
    end

    @testset "Coordinate type" begin
        let pnt_i16 = Point{Int16}(-1,2),
            rec_i32 = Rectangle{Int32}(pnt_i16, 2*pnt_i16),
            box_f32 = BoundingBox{Float32}(rec_i32),
            pnt_f64 = Point(1.0, 3.0)
            @test Int16   === @inferred coord_type(pnt_i16)
            @test Int16   === @inferred coord_type(typeof(pnt_i16))
            @test Int32   === @inferred coord_type(rec_i32)
            @test Int32   === @inferred coord_type(typeof(rec_i32))
            @test Float32 === @inferred coord_type(box_f32)
            @test Float32 === @inferred coord_type(typeof(box_f32))
            @test Float64 === @inferred coord_type(pnt_f64)
            @test Float64 === @inferred coord_type(typeof(pnt_f64))
            @test Int32   === @inferred coord_type(typeof(pnt_i16), typeof(rec_i32))
            @test Float32 === @inferred coord_type(typeof(pnt_i16), typeof(rec_i32), typeof(box_f32))
            @test Float64 === @inferred coord_type(typeof(pnt_i16), typeof(rec_i32), typeof(box_f32), typeof(pnt_f64))
            @test @inferred(promote_coord_type(pnt_i16)) === pnt_i16
            @test @inferred(promote_coord_type(rec_i32)) === rec_i32
            @test @inferred(promote_coord_type(box_f32)) === box_f32
            @test @inferred(promote_coord_type(pnt_f64)) === pnt_f64
            if VERSION < v"1.8" # inference does not wokr here for Julia versions older than 1,8
                @test promote_coord_type(pnt_i16, rec_i32) === (Point{Int32}(pnt_i16), rec_i32)
                @test promote_coord_type(pnt_i16, rec_i32, box_f32) === (Point{Float32}(pnt_i16), Rectangle{Float32}(rec_i32), box_f32)
                @test promote_coord_type(pnt_i16, rec_i32, box_f32, pnt_f64) === (Point{Float64}(pnt_i16), Rectangle{Float64}(rec_i32), BoundingBox{Float64}(box_f32), pnt_f64)
            else
                @test @inferred(promote_coord_type(pnt_i16, rec_i32)) === (Point{Int32}(pnt_i16), rec_i32)
                @test @inferred(promote_coord_type(pnt_i16, rec_i32, box_f32)) === (Point{Float32}(pnt_i16), Rectangle{Float32}(rec_i32), box_f32)
                @test @inferred(promote_coord_type(pnt_i16, rec_i32, box_f32, pnt_f64)) === (Point{Float64}(pnt_i16), Rectangle{Float64}(rec_i32), BoundingBox{Float64}(box_f32), pnt_f64)
            end

            @test @inferred(float(typeof(pnt_i16))) === Point{float(coord_type(pnt_i16))}
            @test @inferred(float(pnt_i16))         === Point{float(coord_type(pnt_i16))}(pnt_i16)
            @test @inferred(float(typeof(rec_i32))) === Rectangle{float(coord_type(rec_i32))}
            @test @inferred(float(rec_i32))         === Rectangle{float(coord_type(rec_i32))}(rec_i32)
            @test @inferred(float(typeof(box_f32))) === typeof(box_f32)
            @test @inferred(float(box_f32))         === box_f32
            @test @inferred(float(typeof(pnt_f64))) === typeof(pnt_f64)
            @test @inferred(float(pnt_f64))          === pnt_f64

            @test @inferred(convert_bare_type(Int16, typeof(pnt_i16)))   === typeof(pnt_i16)
            @test @inferred(convert_bare_type(Int16, pnt_i16))           === pnt_i16
            @test @inferred(convert_bare_type(Int, typeof(pnt_i16)))     === Point{Int}
            @test @inferred(convert_bare_type(Int, pnt_i16))             === Point{Int}(pnt_i16)
            @test @inferred(convert_bare_type(Float64, typeof(pnt_i16))) === Point{Float64}
            @test @inferred(convert_bare_type(Float64, pnt_i16))         === Point{Float64}(pnt_i16)
            @test @inferred(convert_bare_type(Int, typeof(rec_i32)))     === Rectangle{Int}
            @test @inferred(convert_bare_type(Int, rec_i32))             === Rectangle{Int}(rec_i32)
            @test @inferred(convert_bare_type(Float32, typeof(rec_i32))) === Rectangle{Float32}
            @test @inferred(convert_bare_type(Float32, rec_i32))         === Rectangle{Float32}(rec_i32)
            @test @inferred(convert_bare_type(Float64, typeof(pnt_f64))) === typeof(pnt_f64)
            @test @inferred(convert_bare_type(Float64, pnt_f64))         === pnt_f64
            @test @inferred(convert_bare_type(Float32, typeof(pnt_f64))) === Point{Float32}
            @test @inferred(convert_bare_type(Float32, pnt_f64))         === Point{Float32}(pnt_f64)
            @test @inferred(convert_bare_type(Float64, typeof(box_f32))) === BoundingBox{Float64}
            @test @inferred(convert_bare_type(Float64, box_f32))         === BoundingBox{Float64}(box_f32)
            @test @inferred(convert_bare_type(Float32, typeof(box_f32))) === typeof(box_f32)
            @test @inferred(convert_bare_type(Float32, box_f32))         === box_f32

            @test @inferred(convert_real_type(Int16, typeof(pnt_i16)))   === typeof(pnt_i16)
            @test @inferred(convert_real_type(Int16, pnt_i16))           === pnt_i16
            @test @inferred(convert_real_type(Int, typeof(pnt_i16)))     === Point{Int}
            @test @inferred(convert_real_type(Int, pnt_i16))             === Point{Int}(pnt_i16)
            @test @inferred(convert_real_type(Float64, typeof(pnt_i16))) === Point{Float64}
            @test @inferred(convert_real_type(Float64, pnt_i16))         === Point{Float64}(pnt_i16)
            @test @inferred(convert_real_type(Int, typeof(rec_i32)))     === Rectangle{Int}
            @test @inferred(convert_real_type(Int, rec_i32))             === Rectangle{Int}(rec_i32)
            @test @inferred(convert_real_type(Float32, typeof(rec_i32))) === Rectangle{Float32}
            @test @inferred(convert_real_type(Float32, rec_i32))         === Rectangle{Float32}(rec_i32)
            @test @inferred(convert_real_type(Float64, typeof(pnt_f64))) === typeof(pnt_f64)
            @test @inferred(convert_real_type(Float64, pnt_f64))         === pnt_f64
            @test @inferred(convert_real_type(Float32, typeof(pnt_f64))) === Point{Float32}
            @test @inferred(convert_real_type(Float32, pnt_f64))         === Point{Float32}(pnt_f64)
            @test @inferred(convert_real_type(Float64, typeof(box_f32))) === BoundingBox{Float64}
            @test @inferred(convert_real_type(Float64, box_f32))         === BoundingBox{Float64}(box_f32)
            @test @inferred(convert_real_type(Float32, typeof(box_f32))) === typeof(box_f32)
            @test @inferred(convert_real_type(Float32, box_f32))         === box_f32

            @test @inferred(convert_floating_point_type(Float64, typeof(pnt_i16))) === Point{Float64}
            @test @inferred(convert_floating_point_type(Float64, pnt_i16))         === Point{Float64}(pnt_i16)
            @test @inferred(convert_floating_point_type(Float32, typeof(rec_i32))) === Rectangle{Float32}
            @test @inferred(convert_floating_point_type(Float32, rec_i32))         === Rectangle{Float32}(rec_i32)
            @test @inferred(convert_floating_point_type(Float64, typeof(pnt_f64))) === typeof(pnt_f64)
            @test @inferred(convert_floating_point_type(Float64, pnt_f64))         === pnt_f64
            @test @inferred(convert_floating_point_type(Float32, typeof(pnt_f64))) === Point{Float32}
            @test @inferred(convert_floating_point_type(Float32, pnt_f64))         === Point{Float32}(pnt_f64)
            @test @inferred(convert_floating_point_type(Float64, typeof(box_f32))) === BoundingBox{Float64}
            @test @inferred(convert_floating_point_type(Float64, box_f32))         === BoundingBox{Float64}(box_f32)
            @test @inferred(convert_floating_point_type(Float32, typeof(box_f32))) === typeof(box_f32)
            @test @inferred(convert_floating_point_type(Float32, box_f32))         === box_f32
        end
    end

    # Common mathematical operations on geometric objects.
    @testset "Arithmetic ($(typeof(obj)))" for obj in (
        @inferred(Point(-1, 3)),
        @inferred(Rectangle((-1.0f0, 2.0f0), (3.0f0, 4.0f0))),
        @inferred(Circle((-1.0f0, 2.0f0), 3.0f0)),
        @inferred(BoundingBox((-1.0, 2.0), (3.0, 4.0))),
        @inferred(Polygon((-1.0, 2.0), (3.0, 4.0), (5.0, -6.0))),
        )

        # The following tests require that the bounding-box be not empty.
        if obj isa BoundingBox
            @test !isempty(obj)
        end

        # Unary plus does nothing.
        @test +obj === obj

        # Unary minus negates coordinates and should as multiplying by -1.
        if obj isa Polygon
            # Vertices are the same but are stored in another object.
            @test -obj == -one(coord_type(obj))*obj
        else
            @test -obj === -one(coord_type(obj))*obj
        end
        if obj isa Point
            @test -obj === @inferred Point(-obj.x, -obj.y)
        elseif obj isa Rectangle
            @test -obj === @inferred Rectangle(-first(obj), -last(obj))
            @test -obj === @inferred Rectangle((-obj.x0, -obj.y0), (-obj.x1, -obj.y1))
        elseif obj isa Polygon
            @test -obj == @inferred Polygon(map(-, vec(obj)))
        elseif obj isa BoundingBox
            @test !isempty(obj)
            @test -obj === @inferred BoundingBox(-last(obj), -first(obj))
            @test -obj === @inferred BoundingBox((-obj.xmax, -obj.ymax), (-obj.xmin, -obj.ymin))
        end

        # Multiplicative identity.
        @test isone(one(obj))
        @test one(obj) === one(typeof(obj)) === one(coord_type(obj))
        if obj isa Polygon
            # Vertices are the same but are stored in another object.
            @test one(obj)*obj == obj
            @test obj*one(obj) == obj
            @test one(obj)\obj == @inferred float(obj)
            @test obj/one(obj) == @inferred float(obj)
            @test (-one(obj))*obj == -obj
            @test obj*(-one(obj)) == -obj
            @test (-one(obj))\obj == -float(obj)
            @test obj/(-one(obj)) == -float(obj)
        else
            @test one(obj)*obj === obj
            @test obj*one(obj) === obj
            @test one(obj)\obj === @inferred float(obj)
            @test obj/one(obj) === @inferred float(obj)
            @test (-one(obj))*obj === -obj
            @test obj*(-one(obj)) === -obj
            @test (-one(obj))\obj === -float(obj)
            @test obj/(-one(obj)) === -float(obj)
        end

        # Multiplying/dividing a geometric object by a number amounts to
        # scaling around origin.
        α = 7/2
        if obj isa Point
            @test α*obj === Point(α*obj.x, α*obj.y)
            @test α\obj === Point(α\obj.x, α\obj.y)
        elseif obj isa Rectangle
            @test α*obj === @inferred Rectangle(α*first(obj), α*last(obj))
            @test α*obj === @inferred Rectangle((α*obj.x0, α*obj.y0), (α*obj.x1, α*obj.y1))
            @test α\obj === @inferred Rectangle(α\first(obj), α\last(obj))
            @test α\obj === @inferred Rectangle((α\obj.x0, α\obj.y0), (α\obj.x1, α\obj.y1))
        elseif obj isa Polygon
            @test α*obj == @inferred Polygon(map(vertex -> α*vertex, vec(obj)))
            @test α\obj == @inferred Polygon(map(vertex -> α\vertex, vec(obj)))
        elseif obj isa BoundingBox
            @test α*obj === @inferred BoundingBox(α*first(obj), α*last(obj))
            @test α*obj === @inferred BoundingBox((α*obj.xmin, α*obj.ymin), (α*obj.xmax, α*obj.ymax), )
            @test α\obj === @inferred BoundingBox(α\first(obj), α\last(obj))
            @test α\obj === @inferred BoundingBox((α\obj.xmin, α\obj.ymin), (α\obj.xmax, α\obj.ymax), )
            @test (-α)*obj === @inferred BoundingBox(-α*last(obj), -α*first(obj))
            @test (-α)*obj === @inferred BoundingBox((-α*obj.xmax, -α*obj.ymax), (-α*obj.xmin, -α*obj.ymin))
            @test (-α)\obj === @inferred BoundingBox(-α\last(obj), -α\first(obj))
            @test (-α)\obj === @inferred BoundingBox((-α\obj.xmax, -α\obj.ymax), (-α\obj.xmin, -α\obj.ymin))
        end
        if obj isa Polygon
            # Vertices are the same but are stored in another object.
            @test obj*α == α*obj
            @test obj/α == α\obj
        else
            @test obj*α === α*obj
            @test obj/α === α\obj
        end

        # Additive identity.
        @test iszero(zero(obj))
        @test zero(obj) === zero(typeof(obj))
        if obj isa Polygon
            # Vertices are the same but are stored in another object.
            @test obj + zero(obj) == obj
            @test zero(obj) + obj == obj
            @test obj - zero(obj) == obj
            @test zero(obj) - obj == -obj
        else
            @test obj + zero(obj) === obj
            @test zero(obj) + obj === obj
            @test obj - zero(obj) === obj
            @test zero(obj) - obj === -obj
        end

        # Adding a point to a geometric object amounts to shifting the object.
        pnt = Point{Int8}(7,-3)
        if obj isa Polygon
            # Vertices are the same but are stored in another object.
            @test obj + pnt == pnt + obj
            @test obj - pnt == obj + (-pnt)
            @test obj - pnt == -pnt + obj
        else
            @test obj + pnt === pnt + obj
            @test obj - pnt === obj + (-pnt)
            @test obj - pnt === -pnt + obj
        end
        if obj isa Point
            @test obj + pnt === Point(obj.x + pnt.x, obj.y + pnt.y)
            @test obj - pnt === Point(obj.x - pnt.x, obj.y - pnt.y)
            @test pnt - obj === Point(pnt.x - obj.x, pnt.y - obj.y)
        elseif obj isa Rectangle
            @test obj + pnt === @inferred Rectangle(first(obj) + pnt, last(obj) + pnt)
            @test obj + pnt === @inferred Rectangle((obj.x0 + pnt.x, obj.y0 + pnt.y), (obj.x1 + pnt.x, obj.y1 + pnt.y))
            @test obj - pnt === @inferred Rectangle(first(obj) - pnt, last(obj) - pnt)
            @test obj - pnt === @inferred Rectangle((obj.x0 - pnt.x, obj.y0 - pnt.y), (obj.x1 - pnt.x, obj.y1 - pnt.y))
            @test pnt - obj === @inferred Rectangle(pnt - first(obj), pnt - last(obj))
            @test pnt - obj === @inferred Rectangle((pnt.x - obj.x0, pnt.y - obj.y0), (pnt.x - obj.x1, pnt.y - obj.y1))
        elseif obj isa Polygon
            @test obj + pnt == @inferred Polygon(map(vertex -> vertex + pnt, vec(obj)))
            @test pnt + obj == @inferred Polygon(map(vertex -> pnt + vertex, vec(obj)))
            @test obj - pnt == @inferred Polygon(map(vertex -> vertex - pnt, vec(obj)))
            @test pnt - obj == @inferred Polygon(map(vertex -> pnt - vertex, vec(obj)))
        elseif obj isa BoundingBox
            @test obj + pnt === @inferred BoundingBox(first(obj) + pnt, last(obj) + pnt)
            @test obj + pnt === @inferred BoundingBox((obj.xmin + pnt.x, obj.ymin + pnt.y), (obj.xmax + pnt.x, obj.ymax + pnt.y))
            @test obj - pnt === @inferred BoundingBox(first(obj) - pnt, last(obj) - pnt)
            @test obj - pnt === @inferred BoundingBox((obj.xmin - pnt.x, obj.ymin - pnt.y), (obj.xmax - pnt.x, obj.ymax - pnt.y))
            @test pnt - obj === @inferred BoundingBox(pnt - last(obj), pnt - first(obj))
            @test pnt - obj === @inferred BoundingBox((pnt.x - obj.xmax, pnt.y - obj.ymax), (pnt.x - obj.xmin, pnt.y - obj.ymin))
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
            @test @inferred(float(A)) === A
            @test @inferred(float(typeof(A))) === typeof(A)
            for T1 in types, T2 in types
                T = promote_type(T1,T2)
                @test promote_type(AffineTransform{T1,T1,T1}, AffineTransform{T2,T2,T2}) ===
                    AffineTransform{T,T,T}
            end
            @testset "T = $T" for T in types
                @test_throws MethodError T(A)
                if VERSION ≥ v"1.3"
                    # FIXME: This makes old Julia versions panic on illegal instruction...
                    @test_throws MethodError T.(A)
                end
            end
            @test AffineTransform(A) === A
            @test AffineTransform{bare_type(A)}(A) === A
            for G in (I, A, B)
                for T in types
                    H = @inferred AffineTransform{T}(G)
                    @test typeof(H) <: AffineTransform{T}
                    @test @inferred(float(H)) === H
                    @test @inferred(float(typeof(H))) === typeof(H)
                    @test @inferred(eltype(H)) === Union{offsets_type(H),factors_type(H)}
                    @test @inferred(eltype(typeof(H))) === Union{offsets_type(typeof(H)),factors_type(typeof(H))}
                    @test bare_type(H) === T
                    @test real_type(H) === T
                    @test floating_point_type(H) === T
                    @test factors_type(H) === T
                    @test offsets_type(H) === T
                    @test (H === G) == (bare_type(H) === bare_type(G))
                    @test Tuple(H) ≈ map(T, Tuple(G))
                    H = @inferred convert(AffineTransform{T}, G)
                    @test typeof(H) <: AffineTransform{T}
                    @test @inferred(float(H)) === H
                    @test @inferred(float(typeof(H))) === typeof(H)
                    @test @inferred(eltype(H)) === Union{offsets_type(H),factors_type(H)}
                    @test @inferred(eltype(typeof(H))) === Union{offsets_type(typeof(H)),factors_type(typeof(H))}
                    @test bare_type(H) === T
                    @test real_type(H) === T
                    @test floating_point_type(H) === T
                    @test factors_type(H) === T
                    @test offsets_type(H) === T
                    @test (H === G) == (bare_type(H) === bare_type(G))
                    @test Tuple(H) ≈ map(T, Tuple(G))
                    H = @inferred convert_bare_type(T, G)
                    @test typeof(H) <: AffineTransform{T}
                    @test @inferred(float(H)) === H
                    @test @inferred(float(typeof(H))) === typeof(H)
                    @test @inferred(eltype(H)) === Union{offsets_type(H),factors_type(H)}
                    @test @inferred(eltype(typeof(H))) === Union{offsets_type(typeof(H)),factors_type(typeof(H))}
                    @test bare_type(H) === T
                    @test real_type(H) === T
                    @test floating_point_type(H) === T
                    @test factors_type(H) === T
                    @test offsets_type(H) === T
                    @test (H === G) == (bare_type(H) === bare_type(G))
                    @test Tuple(H) ≈ map(T, Tuple(G))
                    H = @inferred convert_real_type(T, G)
                    @test typeof(H) <: AffineTransform{T}
                    @test @inferred(float(H)) === H
                    @test @inferred(float(typeof(H))) === typeof(H)
                    @test @inferred(eltype(H)) === Union{offsets_type(H),factors_type(H)}
                    @test @inferred(eltype(typeof(H))) === Union{offsets_type(typeof(H)),factors_type(typeof(H))}
                    @test bare_type(H) === T
                    @test real_type(H) === T
                    @test floating_point_type(H) === T
                    @test factors_type(H) === T
                    @test offsets_type(H) === T
                    @test (H === G) == (bare_type(H) === bare_type(G))
                    @test Tuple(H) ≈ map(T, Tuple(G))
                    H = @inferred convert_floating_point_type(T, G)
                    @test typeof(H) <: AffineTransform{T}
                    @test @inferred(float(H)) === H
                    @test @inferred(float(typeof(H))) === typeof(H)
                    @test @inferred(eltype(H)) === Union{offsets_type(H),factors_type(H)}
                    @test @inferred(eltype(typeof(H))) === Union{offsets_type(typeof(H)),factors_type(typeof(H))}
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
            @test A\B ≈ inv(A) ∘ B
            @test B\A ≈ inv(B) ∘ A
            @test A/B ≈ A ∘ inv(B)
            @test B/A ≈ B ∘ inv(A)

            @test A\C ≈ inv(A) ∘ C
            @test C\A ≈ inv(C) ∘ A
            @test A/C ≈ A ∘ inv(C)
            @test C/A ≈ C ∘ inv(A)

            @test B\C ≈ inv(B) ∘ C
            @test C\B ≈ inv(C) ∘ B
            @test B/C ≈ B ∘ inv(C)
            @test C/B ≈ C ∘ inv(B)
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

        @testset "ldiv" begin
            for M in (I, A, B)
                x, y = @inferred TwoDimensional.ldiv(M)
                @test M(x, y) ≈ (0,0) atol=16*eps(Float64)
                b = Point(1.0, -2.0)
                c = @inferred TwoDimensional.ldiv(M, b)
                @test c isa Point
                @test M(c) ≈ b atol=16*eps(Float64)
                @test M\b === c
                b = Point(3.0, -1.0)
                c = @inferred TwoDimensional.ldiv(M, b)
                @test c isa Point
                @test M(c) ≈ b atol=16*eps(Float64)
                @test M\b === c
            end
        end

        @testset "show" begin
            for M in (I, A, B, C)
                @test occursin(r"\bAffineTransform\b", string(M))
            end
        end
    end
end

@testset "Bounding-box algorithm" begin
    if true # false for debugging to avoid too many errors, true for production
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

end # module
