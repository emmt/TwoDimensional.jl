"""
    circ = Circle{T}(center::Point, radius)
    circ = Circle{T}((x, y), radius)
    circ = Circle{T}(; center=..., radius=...)

construct a circle of given `center` and `radius` or `diameter`. The center may
be specified by its coordinates `(x,y)` along the Cartesian axes. Parameter `T`
is the type used to store coordinates, it may be omitted.

Circles have the following properties reflecting the keywords accepted by their
constructor (`diameter` is provided for convenience):

    circ.center   -> center::Point{T}
    circ.radius   -> radius::T
    circ.diameter -> 2*circ.radius

Circles can be iterated to retrieve their center and their radius (in that order):

    center, radius = circ
    (x0, y0), radius = circ

See also [`Point`](@ref), [`BoundingBox`](@ref).

"""
Circle{T}(center::PointLike, radius) where {T} = Circle{T}(Point(center), radius)
Circle(center::PointLike, radius) = Circle(Point(center), radius)
function Circle(center::Point, radius)
    T = promote_type(coord_type(center), typeof(radius))
    return Circle{T}(center, radius)
end

# Unpack arguments provided by a tuple.
Circle((center, radius)::Tuple{PointLike,Number}) = Circle(center, radius)
Circle{T}((center, radius)::Tuple{PointLike,Number}) where {T} = Circle{T}(center, radius)

# Keyword-only constructors. NOTE Having `diameter` as a possible keyword is not a good idea
# because dividing it by 2 enforces that coordinate type be at least `Float64`.
@inline Circle(; center, radius) = Circle(center, radius)
@inline Circle{T}(; center, radius) where {T} = Circle{T}(center, radius)

# Convert/copy constructors.
Circle(circ::Circle) = circ
Circle{T}(circ::Circle{T}) where {T} = circ
Circle{T}(circ::Circle) where {T} = Circle{T}(elements(circ)...)
Base.convert(::Type{T}, circ::T) where {T<:Circle} = circ
Base.convert(::Type{T}, obj::CircleLike) where {T<:Circle} = T(obj)

# Make circles indexable iterators.
Base.length(        circ::Circle) = 2
Base.eltype(        circ::Circle) = eltype(typeof(circ))
Base.IteratorSize(  circ::Circle) = IteratorSize(typeof(circ))
Base.IteratorEltype(circ::Circle) = IteratorEltype(typeof(circ))
Base.eltype(        ::Type{<:Circle{T}}) where {T} = Union{Point{T},T}
Base.IteratorSize(  ::Type{<:Circle}) = HasLength()
Base.IteratorEltype(::Type{<:Circle}) = HasEltype()
@inline Base.iterate(circ::Circle, i::Int = 1) =
    i == 1 ? (center(circ), 2) :
    i == 2 ? (radius(circ), 3) : nothing

# Properties of circles.
Base.propertynames(::Circle) = (:center, :radius, :diameter)
Base.getproperty(circ::Circle, key::Symbol) =
    key === :center   ? center(circ)   :
    key === :radius   ? radius(circ)   :
    key === :diameter ? diameter(circ) : throw(KeyError(key))

Base.show(io::IO, circ::Circle{T}) where {T} =
    print(io, "Circle{", T,
          "}(center = (", circ.center.x, ", ", circ.center.y,
          "), radius = ", circ.radius, ")")
