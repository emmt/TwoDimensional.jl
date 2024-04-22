"""
    TwoDimensional.AbstractPoint2D{T}

is an alias for [`TwoDimensional.AbstractPoint{T}`](@ref TwoDimensional.AbstractPoint).

"""
const AbstractPoint2D{T} = AbstractPoint{T}

"""
    TwoDimensional.Point2D{T}

is an alias for [`TwoDimensional.Point{T}`](@ref TwoDimensional.Point).

"""
const Point2D{T} = Point{T}

"""
    TwoDimensional.Rectangle2D{T}

is an alias for [`TwoDimensional.Rectangle{T}`](@ref TwoDimensional.Rectangle).

"""
const Rectangle2D{T} = Point{T}

"""
    TwoDimensional.BoundingBox2D{T}

is an alias for [`TwoDimensional.BoundingBox{T}`](@ref TwoDimensional.BoundingBox).

"""
const BoundingBox2D{T} = BoundingBox{T}

"""
    TwoDimensional.AffineTransform2D{T,R,S}

is an alias for [`TwoDimensional.AffineTransform{T,R,S}`](@ref
TwoDimensional.AffineTransform).

"""
const AffineTransform2D{T,R,S} = AffineTransform{T,R,S}
