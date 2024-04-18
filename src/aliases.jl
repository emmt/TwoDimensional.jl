"""
    AbstractPoint2D{T}

is an alias for [`TwoDimensional.AbstractPoint{T}`](@ref AbstractPoint).

"""
const AbstractPoint2D{T} = AbstractPoint{T}

"""
    Point2D{T}

is an alias for [`TwoDimensional.Point{T}`](@ref Point).

"""
const Point2D{T} = Point{T}

"""
    BoundingBox2D{T}

is an alias for [`TwoDimensional.BoundingBox{T}`](@ref BoundingBox).

"""
const BoundingBox2D{T} = BoundingBox{T}

"""
    AffineTransform2D{T,R,S}

is an alias for [`TwoDimensional.AffineTransform{T}`](@ref AffineTransform).

"""
const AffineTransform2D{T,R,S} = AffineTransform{T,R,S}
