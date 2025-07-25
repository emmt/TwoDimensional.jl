"""
    TwoDimensional.AbstractPoint2D{T}

is an alias for [`TwoDimensional.AbstractPoint{T}`](@ref TwoDimensional.AbstractPoint).

""" AbstractPoint2D
@public AbstractPoint2D
const AbstractPoint2D{T} = AbstractPoint{T}

"""
    TwoDimensional.Point2D{T}

is an alias for [`TwoDimensional.Point{T}`](@ref TwoDimensional.Point).

""" Point2D
@public Point2D
const Point2D{T} = Point{T}

"""
    TwoDimensional.Rectangle2D{T}

is an alias for [`TwoDimensional.Rectangle{T}`](@ref TwoDimensional.Rectangle).

""" Rectangle2D
@public Rectangle2D
const Rectangle2D{T} = Rectangle{T}

"""
    TwoDimensional.Circle2D{T}

is an alias for [`TwoDimensional.Circle{T}`](@ref TwoDimensional.Circle).

""" Circle2D
@public Circle2D
const Circle2D{T} = Circle{T}

"""
    TwoDimensional.Polygon2D{T}

is an alias for [`TwoDimensional.Polygon{T}`](@ref TwoDimensional.Polygon).

""" Polygon2D
@public Polygon2D
const Polygon2D{T} = Polygon{T}

"""
    TwoDimensional.BoundingBox2D{T}

is an alias for [`TwoDimensional.BoundingBox{T}`](@ref TwoDimensional.BoundingBox).

""" BoundingBox2D
@public BoundingBox2D
const BoundingBox2D{T} = BoundingBox{T}

"""
    TwoDimensional.MaskElement2D{T,S}

is an alias for [`TwoDimensional.MaskElement{T,S}`](@ref
TwoDimensional.MaskElement). `T` is the coordinate type, `S` is the type of the
elementary geometrical object defining the shape of the mask.

""" MaskElement2D
@public MaskElement2D
const MaskElement2D{T,S} = MaskElement{T,S}

"""
    TwoDimensional.AffineTransform2D{T,R,S}

is an alias for [`TwoDimensional.AffineTransform{T,R,S}`](@ref
TwoDimensional.AffineTransform).

""" AffineTransform2D
@public AffineTransform2D
const AffineTransform2D{T,R,S} = AffineTransform{T,R,S}
