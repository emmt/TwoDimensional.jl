var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "TwoDimensional is a Julia package which provides useful types and methods to define and manipulate 2-dimensional points, bounding boxes and affine coordinate transforms.Other related packages:CoordinateTransformationsThe source code of TwoDimensional is available on GitHub."
},

{
    "location": "#Exported-Types-1",
    "page": "Introduction",
    "title": "Exported Types",
    "category": "section",
    "text": "using TwoDimensionalgives you types AffineTransform{T}, Point{T} and BoundingBox{T} parameterized by the type T of their components (T must be floating point for AffineTransform{T}).To avoid conflicts with other packages, you may use/import TwoDimensional.Suffixed which gives you types AffineTransform2D{T}, Point2D{T} and BoundingBox2D{T} instead, that is with suffix 2D.You can also fine tune what you want.  For instance:using TwoDimensional: AffineTransform, Point2D"
},

{
    "location": "#Table-of-contents-1",
    "page": "Introduction",
    "title": "Table of contents",
    "category": "section",
    "text": "Pages = [\"install.md\", \"AffineTransforms.md\",\n         \"Points.md\", \"BoundingBoxes.md\", \"reference.md\"]"
},

{
    "location": "#Index-1",
    "page": "Introduction",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "install/#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "install/#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "TwoDimensional is not yet an offical Julia package but it is easy to install it from Julia as explained here."
},

{
    "location": "install/#Using-the-package-manager-1",
    "page": "Installation",
    "title": "Using the package manager",
    "category": "section",
    "text": "At the REPL of Julia, hit the ] key to switch to the package manager REPL (you should get a ... pkg> prompt) and type:pkg> add https://github.com/emmt/TwoDimensional.jlwhere pkg> represents the package manager prompt and https protocol has been assumed; if ssh is more suitable for you, then type:pkg> add git@github.com:emmt/TwoDimensional.jlinstead.  To check whether the TwoDimensional package works correctly, type:pkg> test TwoDimensionalLater, to update to the last version (and run tests), you can type:pkg> update TwoDimensional\npkg> build TwoDimensional\npkg> test TwoDimensionalIf something goes wrong, it may be because you already have an old version of TwoDimensional.  Uninstall TwoDimensional as follows:pkg> rm TwoDimensional\npkg> gc\npkg> add https://github.com/emmt/TwoDimensional.jlbefore re-installing.To revert to Julia\'s REPL, hit the Backspace key at the ... pkg> prompt."
},

{
    "location": "install/#Installation-in-scripts-1",
    "page": "Installation",
    "title": "Installation in scripts",
    "category": "section",
    "text": "To install TwoDimensional in a Julia script, write:if VERSION >= v\"0.7.0-\"; using Pkg; end\nPkg.add(PackageSpec(url=\"https://github.com/emmt/TwoDimensional.jl\", rev=\"master\"));or with url=\"git@github.com:emmt/TwoDimensional.jl\" if you want to use ssh.This also works from the Julia REPL."
},

{
    "location": "AffineTransforms/#",
    "page": "Affine Coordinate Transforms",
    "title": "Affine Coordinate Transforms",
    "category": "page",
    "text": ""
},

{
    "location": "AffineTransforms/#Affine-Coordinate-Transforms-1",
    "page": "Affine Coordinate Transforms",
    "title": "Affine Coordinate Transforms",
    "category": "section",
    "text": "TwoDimensional provides types and methods to deal with 2-dimensional affine coordinate transforms.An affine 2D transform C is defined by 6 real coefficients, Cxx, Cxy, Cx, Cyx, Cyy and Cy.  Such a transform maps (x,y) as (xp,yp) given by:xp = Cxx*x + Cxy*y + Cx\nyp = Cyx*x + Cyy*y + CyThe immutable type AffineTransform{T} is used to store an affine 2D transform with coefficients of floating-point type T, it can be created by:I = AffineTransform{T}() # yields the identity with type T\nC = AffineTransform{T}(Cxx, Cxy, Cx, Cyx, Cyy, Cy)If the parameter T is omitted, it is guessed from the types of the coefficients."
},

{
    "location": "AffineTransforms/#Operations-with-affine-2D-transforms-1",
    "page": "Affine Coordinate Transforms",
    "title": "Operations with affine 2D transforms",
    "category": "section",
    "text": "Many operations are available to manage or apply affine transforms."
},

{
    "location": "AffineTransforms/#Perform-a-coordinate-transformation-1",
    "page": "Affine Coordinate Transforms",
    "title": "Perform a coordinate transformation",
    "category": "section",
    "text": "There are several possibilities to apply an affine transform A to coordinates (x,y):(xp, yp) = A(x,y)       # apply affine transform A to coordinates (x,y)\n(xp, yp) = A*(x,y)      # idem\n(xp, yp) = A(v)         # idem, with v = (x,y)\n(xp, yp) = A*v          # idemwhere x, y, xp and yp are reals while v = (x,y) is a 2-tuple of real coordinates.  Coordinates may also be specified as Point instances:A(Point(x,y)) -> Point(xp, yp)\nA*Point(x,y)  -> Point(xp, yp)"
},

{
    "location": "AffineTransforms/#Reciprocal-transform-1",
    "page": "Affine Coordinate Transforms",
    "title": "Reciprocal transform",
    "category": "section",
    "text": "The reciprocal transform of A is given by:B = inv(A)"
},

{
    "location": "AffineTransforms/#Compose-coordinate-transformations-1",
    "page": "Affine Coordinate Transforms",
    "title": "Compose coordinate transformations",
    "category": "section",
    "text": "To compose 2 (or more) transforms A and B, do one of:C = compose(A, B, ...)\nC = A∘B\nC = A*B\nC = A⋅Ball these statements yields an object C which applies B then A.  Note that ∘ and ⋅ can be typed by \\\\circ<tab> and \\\\cdot<tab>.Left and right divisions respectively write:R = A/B   # right division, same as: R = compose(A, inv(B))\nL = A\\B   # left division, same as: L = compose(inv(A), B)"
},

{
    "location": "AffineTransforms/#Translation-of-coordinates-1",
    "page": "Affine Coordinate Transforms",
    "title": "Translation of coordinates",
    "category": "section",
    "text": "An operator B which which applies a translations by (x,y) after a coordinate transform A (possibly the identity) can be created by one of the following statements:B = translate(x, y, A)\nB = translate(v, A)\nB = v + Awhere v is a 2-tuple of coordinates, i.e. v = (x,y) or a Point, i.e. v = Point(x,y).To perform the translation before the coordinate transform A, do:B = translate(A, x, y)\nB = translate(A, v)\nB = A + v"
},

{
    "location": "AffineTransforms/#Rotation-of-coordinates-1",
    "page": "Affine Coordinate Transforms",
    "title": "Rotation of coordinates",
    "category": "section",
    "text": "B = rotate(θ, A)   # B = apply A then rotate by angle θ\nC = rotate(A, θ)   # C = rotate by angle θ then apply A"
},

{
    "location": "AffineTransforms/#Scaling-of-coordinates-1",
    "page": "Affine Coordinate Transforms",
    "title": "Scaling of coordinates",
    "category": "section",
    "text": "B = scale(ρ, A)    # B = apply A then scale by ρ\nB = ρ*A            # idem\nC = scale(A, ρ)    # C = scale by ρ then apply A\nC = A*ρ            # idem"
},

{
    "location": "AffineTransforms/#Miscellaneous-1",
    "page": "Affine Coordinate Transforms",
    "title": "Miscellaneous",
    "category": "section",
    "text": "det(A) returns the determinant of the linear part of the affine transform A.jacobian(A) returns the Jacobian of the affine transform A, that is the absolute value of the determinant of its linear part."
},

{
    "location": "AffineTransforms/#Type-conversion-1",
    "page": "Affine Coordinate Transforms",
    "title": "Type conversion",
    "category": "section",
    "text": "As a general rule, the floating-point type T of an AffineTransform{T} is imposed for all operations and for the result.  The floating-point type of the composition of several coordinate transforms is the promoted type of the transforms which have been composed.The type of the coefficients of the affine transform A  is given by:eltype(A)To convert the floating-point type of the coefficients of A to be T, do one of:B = T.(A)\nB = AffineTransform{T}(A)\nB = convert(AffineTransform{T}, A)"
},

{
    "location": "Points/#",
    "page": "Points",
    "title": "Points",
    "category": "page",
    "text": ""
},

{
    "location": "Points/#Points-1",
    "page": "Points",
    "title": "Points",
    "category": "section",
    "text": "An object whose type is derived from AbstractPoint{T} (or AbstractPoint2D{T} if TwoDimensional.Suffixed is used instead of TwoDimensional) has 2-D coordinates: its abscissa and ordinate respectively named x and y.  The parameter T is the type of the coordinates and can be retrieved by the eltype method."
},

{
    "location": "Points/#Aliases-1",
    "page": "Points",
    "title": "Aliases",
    "category": "section",
    "text": "Call:using TwoDimensional.Suffixedinstead of:using TwoDimensionalto have Point2D, WeightedPoint2D and AbstractPoint2D provided as respective aliases to TwoDimensional.Point TwoDimensional.WeightedPoint and TwoDimensional.AbstractPoint."
},

{
    "location": "Points/#Construction-1",
    "page": "Points",
    "title": "Construction",
    "category": "section",
    "text": "The most simple concrete type is Point{T} constructed by:Point(x,y)where (x,y) are the coordinates of the point.  Weighted points of type WeightedPoint{T} associate a weight and coordinates:WeightedPoint(w,x,y)where (x,y) are the coordinates of the point and w its weight (nonnegative by convention).Coordinates and weights can also be specified by keywords:Point(x=xval, y=yval)orWeightedPoint(w=wgt, x=xval, y=yval)There are no default values for keywords w, x and y so they must all be specified."
},

{
    "location": "Points/#Conversion-1",
    "page": "Points",
    "title": "Conversion",
    "category": "section",
    "text": "Simple points can be constructed from a 2-tuple of coordinates or from an instance of 2-dimensional CartesianIndex:v = (x,y)\nI = Cartesianindex(x,y)\nPoint(v)    # yields Point(x,y)\nPoint(I)    # yields Point(x,y)and reciprocally:P = Point(x, y)\nTuple(P)          # yields (P.x, P.y)\nCartesianIndex(P) # yields Cartesianindex(P.x, P.y)Coordinate type conversion, say to type T, is done by:P = Point(x, y)\nPoint{T}(P)\nconvert(Point{T}, P)\nT.(P)The latter form involves broadcasting rules and may be a bit slower."
},

{
    "location": "Points/#Operations-on-Points-1",
    "page": "Points",
    "title": "Operations on Points",
    "category": "section",
    "text": "The addition (resp. subtraction) of two points adds (resp. subtracts) their coordinates:Point(x1,y1) + Point(x2,y2)   # yields Point(x1+x2,y1+y2)\nPoint(x1,y1) - Point(x2,y2)   # yields Point(x1-x2,y1-y2)Unary minus of a point negates its coordinates:-Point(x,y)   # yields Point(-x,-y)A point may be multiplied or divided by a scalar, say α, to scale its coordinates:α*Point(x,y)  # yields Point(α*x,α*y)\nPoint(x,y)*α  # yields Point(α*x,α*y)\nα\\Point(x,y)  # yields Point((1/α)*x,(1/α)*y)\nPoint(x,y)/α  # yields Point((1/α)*x,(1/α)*y)Taking the hypothenuse or the arctangent of a point P = Point(x,y) yield its distance to the origin O = Point(0,0) and the angle between OP and the abscissa axis:hypot(Point(x,y))  # yields hypot(x, y)\natan(Point(x,y))   # yields atan(y, x)The distance between two points is given by the distance method:distance(Point(x1,y1),Point(x2,y2)  # yields hypot(x1-x2,y1-y2)The nearest point to an instance obj of Point is given by:round([T,] obj)which rounds the coordinates of obj to the nearest integer.  Optional argument T is to specify the type of the result or the type of the coordinates of the result."
},

{
    "location": "BoundingBoxes/#",
    "page": "Bounding-Boxes",
    "title": "Bounding-Boxes",
    "category": "page",
    "text": ""
},

{
    "location": "BoundingBoxes/#Bounding-Boxes-1",
    "page": "Bounding-Boxes",
    "title": "Bounding-Boxes",
    "category": "section",
    "text": "2-D bounding boxes are build by:BoundingBox(xmin,xmax,ymin,ymax)to represent a 2D rectangular box whose sides are aligned with the coordinate axes and containing points of coordinates (x,y) such that xmin ≤ x ≤ xmax and ymin ≤ y ≤ ymax.The type of the coordinates, say T, can be explicitly specified:BoundingBox{T}(xmin,xmax,ymin,ymax)The type of the coordinates and can be retrieved by the eltype method."
},

{
    "location": "BoundingBoxes/#Aliases-1",
    "page": "Bounding-Boxes",
    "title": "Aliases",
    "category": "section",
    "text": "Call:using TwoDimensional.Suffixedinstead of:using TwoDimensionalto have BoundingBox2D provided as an alias to TwoDimensional.BoundingBox."
},

{
    "location": "BoundingBoxes/#Construction-1",
    "page": "Bounding-Boxes",
    "title": "Construction",
    "category": "section",
    "text": "The coordinates of a bounding box can be specified by keywords:BoundingBox(xmin=x0, ymin=y0, xmax=x1, ymax=y1)There are no default values for keywords xmin, xmax, ymin and ymax so all must be specified.A bounding box can be constructed from a 4-tuple of coordinates and conversely:BoundingBox((x0,x1,y0,y1))      # yields BoundingBox(x0,x1,y0,y1)\nTuple(BoundingBox(x0,x1,y0,y1)) # yields (x0,x1,y0,y1)A bounding box can be constructed from its first and last points (i.e. at the lower-left and upper right opposite corners) specified as instances of Point or of CartesianIndex:BoundingBox(Point(x0,y0), Point(x1,y1))\nBoundingBox(CartesianIndex(x0,y0), CartesianIndex(x1,y1))both yield the same result as:BoundingBox(x0,x1,y0,y1)Conversely, methods first(B) and last(B) respectively yield the lower left and upper right corners of the bounding-box B (as a Point instance):first(BoundingBox(x0,x1,y0,y1)) # yields Point(x0,y0)\nlast(BoundingBox(x0,x1,y0,y1))  # yields Point(x1,y1)Integer-valued unit-ranges can be specified to define a bounding box.  For example:BoundingBox(x0:x1,y0:y1)    # 2 unit-range\nBoundingBox((x0:x1,y0:y1))  # a 2-tuple of unit rangeThis makes possible writing:BoundingBox(axes(A))to get the bounding box corresponding to all indices of array A which is also given by:BoundingBox(A)Conversely:axes(BoundingBox(x0,x1,y0,y1))yields the axes of a bounding-box with integer coordinates, that is (x0:x1,y0:y1).  To get the k-th axis of a bounding-box B, call axes(B,k).To loop over the Cartesian indices edfined by a bounding-box B with integer coordinates, you can just write:for I in CartesianIndices(B)\n   ...\nendA bounding box may also be constructed by applying a predicate function to the elements of a 2-dimensional array:BoundingBox(f, A)yields the bounding box of all integer coordinates (x,y) such that f(A[x,y]) yields true."
},

{
    "location": "BoundingBoxes/#Conversion-1",
    "page": "Bounding-Boxes",
    "title": "Conversion",
    "category": "section",
    "text": "Coordinate type conversion, say to type T, is done by:B = BoundingBox(x0,x1,y0,y1)\nBoundingBox{T}(B)\nconvert(BoundingBox{T}, B)\nT.(B)The latter form involves broadcasting rules and may be a bit slower."
},

{
    "location": "BoundingBoxes/#Union-and-Intersection-of-Bounding-Boxes-1",
    "page": "Bounding-Boxes",
    "title": "Union and Intersection of Bounding-Boxes",
    "category": "section",
    "text": "The union of bounding boxes b1, b2, ... is given by one of:B1 ∪ B2 ∪ ...\nunion(B1, B2, ...)wich both yield the smallest bounding box containing the bounding boxes B1, B2, ...The intersection of bounding boxes B1, B2, ... is given by one of:B1 ∩ B2 ∩ ...\nintersect(B1, B2, ...)wich both yield the largest bounding box contained into the bounding boxes B1, B2, ...The maximal or minimal bounding-box with coordinates of type T that can be constructed are respectively given by typemax(BoundingBox{T}) and typemin(BoundingBox{T}).  These can be useful to initiate a shrinking ar a growing bounding-box.  The call:BoundingBox{T}(nothing)yields the same result as typemin(BoundingBox{T}).The method isempty(B) yields whether a bounding-box B is empty or not."
},

{
    "location": "BoundingBoxes/#Interrior,-Exterior,-Nearest,-etc.-1",
    "page": "Bounding-Boxes",
    "title": "Interrior, Exterior, Nearest, etc.",
    "category": "section",
    "text": "Given the bounding-box B, interior(B) and exterior(B) respectively yield the largest interior and smallest exterior bounding boxes with integer bounds.round(B) or round(T,B) yield a bounding box whose limits are those of the bounding-box B rounded to the nearest integer values.center(B) yields the Point whose coordinates are the geometrical center of the bounding-box B.area(B) yields the area of a bounding-box B."
},

{
    "location": "BoundingBoxes/#Arithmetic-and-Basic-Methods-1",
    "page": "Bounding-Boxes",
    "title": "Arithmetic and Basic Methods",
    "category": "section",
    "text": "Adding or subtracting a scalar δ to a bounding box B can be used to add or remove a margin δ to the bounding box B:BoundingBox(x0,x1,y0,y1) + δ # yields BoundingBox(x0-δ,x1+δ,y0-δ,y1+δ)Adding or subtracting a point P to a bounding box B can be used to shift the bounding box B:BoundingBox(x0,x1,y0,y1) + Point(x,y) # yields BoundingBox(x0+x,x1+x,y0+y,y1+y)eltype(B) yields the type of the coordinates of a bounding-box B.Basic methods size(B[,k]) and axes(B[,k]) can be applied to an integer-valued bounding-box B."
},

{
    "location": "reference/#",
    "page": "Reference",
    "title": "Reference",
    "category": "page",
    "text": ""
},

{
    "location": "reference/#Reference-1",
    "page": "Reference",
    "title": "Reference",
    "category": "section",
    "text": "The following provides detailled documentation about types and methods provided by the TwoDimensional package.  This information is also available from the REPL by typing ? followed by the name of a method or a type."
},

{
    "location": "reference/#TwoDimensional.AffineTransforms.AffineTransform",
    "page": "Reference",
    "title": "TwoDimensional.AffineTransforms.AffineTransform",
    "category": "type",
    "text": "Affine 2D Transforms\n\nAn affine 2D transform A is defined by 6 real coefficients, Axx, Axy, Ax, Ayx, Ayy and Ay.  Such a transform maps (x,y) as (xp,yp) given by:\n\nxp = Axx*x + Axy*y + Ax\nyp = Ayx*x + Ayy*y + Ay\n\nThe immutable type AffineTransform is used to store an affine 2D transform A, it can be created by:\n\nI = AffineTransform{T}() # yields the identity with type T\nA = AffineTransform{T}(Axx, Axy, Ax, Ayx, Ayy, Ay)\n\nThe parameter T above is used to specify the floating-point type for the coefficients; if omitted it is guessed from the types of the coefficients.\n\nOperations with affine 2D transforms\n\nMany operations are available to manage or apply affine transforms:\n\n(xp, yp) = A(x,y)       # apply affine transform A to coordinates (x,y)\n(xp, yp) = A*(x,y)      # idem\n(xp, yp) = A(v)         # idem, with v = (x,y)\n(xp, yp) = A*v          # idem\n\nA(Point(x,y)) -> Point(xp, yp)\nA*Point(x,y)  -> Point(xp, yp)\n\nC = compose(A, B, ...)  # compose 2 (or more) transforms, C = apply B then A\nC = A∘B                 # idem\nC = A*B                 # idem\nC = A⋅B                 # idem\n\nB = translate(x, y, A)  # B = apply A then translate by (x,y)\nB = translate(v, A)     # idem with v = (x,y)\nB = v + A               # idem\n\nB = translate(A, x, y)  # B = translate by (x,y) then apply A\nB = translate(A, v)     # idem with v = (x,y)\nB = A + v               # idem\n\nB = rotate(θ, A)   # B = apply A then rotate by angle θ\nC = rotate(A, θ)   # C = rotate by angle θ then apply A\n\nB = scale(ρ, A)    # B = apply A then scale by ρ\nB = ρ*A            # idem\nC = scale(A, ρ)    # C = scale by ρ then apply A\nC = A*ρ            # idem\n\nB = inv(A)         # reciprocal coordinate transform\nC = A/B            # right division, same as: C = compose(A, inv(B))\nC = A\\B            # left division, same as: C = compose(inv(A), B)\n\n\"∘\" and \"⋅\" can be typed by \\circ<tab> and \\cdot<tab>.\n\nType conversion\n\nAs a general rule, the floating-point type T of an AffineTransform{T} is imposed for all operations and for the result.  The floating-point type of the composition of several coordinate transforms is the promoted type of the transforms which have been composed.\n\nCalling eltype(A) yields floating-point type of the coefficients of the 2D affine transform A.  To convert the floating-point type of the coefficients of A to be T, do one of:\n\nB = T.(A)\nB = AffineTransform{T}(A)\nB = convert(AffineTransform{T}, A)\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.AffineTransforms.scale",
    "page": "Reference",
    "title": "TwoDimensional.AffineTransforms.scale",
    "category": "function",
    "text": "Scaling an affine transform\n\nThere are two ways to combine a scaling by a factor ρ with an affine transform A.  Left-scaling as in:\n\nB = scale(ρ, A)\n\nresults in scaling the output of the transform; while right-scaling as in:\n\nC = scale(A, ρ)\n\nresults in scaling the input of the transform.  The above examples yield transforms which behave as:\n\nB(v) = ρ.*A(v)\nC(v) = A(ρ.*v)\n\nwhere v is any 2-element tuple.\n\nThe same results can be obtained with the * operator:\n\nB = ρ*A    # same as: B = scale(ρ, A)\nC = A*ρ    # same as: B = scale(A, ρ)\n\nSee also: AffineTransform, rotate, translate.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.AffineTransforms.rotate",
    "page": "Reference",
    "title": "TwoDimensional.AffineTransforms.rotate",
    "category": "function",
    "text": "Rotating an affine transform\n\nThere are two ways to combine a rotation by angle θ (in radians counterclockwise) with an affine transform A.  Left-rotating as in:\n\nB = rotate(θ, A)\n\nresults in rotating the output of the transform; while right-rotating as in:\n\nC = rotate(A, θ)\n\nresults in rotating the input of the transform.  The above examples are similar to:\n\nB = R∘A\nC = A∘R\n\nwhere R implements rotation by angle θ around (0,0).\n\nSee also: AffineTransform, scale, translate.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.AffineTransforms.translate",
    "page": "Reference",
    "title": "TwoDimensional.AffineTransforms.translate",
    "category": "function",
    "text": "Translating an affine transform\n\nAffine transforms can be letf- or right-translated.\n\ntranslate(x, y, A)\n\nor\n\ntranslate((x,y), A)\n\nyield an affine transform which translate the output of affine transform A by offsets x and y.\n\ntranslate(A, x, y)\n\nor\n\ntranslate(A, (x,y))\n\nyield an affine transform which translate the input of affine transform A by offsets x and y.\n\nThe same results can be obtained with the + operator:\n\nB = (x,y) + A    # same as: B = translate((x,y), A)\nB = A + (x,y)    # same as: B = translate(A, (x,y))\n\nSee also: AffineTransform, rotate, scale.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.AffineTransforms.jacobian",
    "page": "Reference",
    "title": "TwoDimensional.AffineTransforms.jacobian",
    "category": "function",
    "text": "jacobian(A) returns the Jacobian of the affine transform A, that is the absolute value of the determinant of its linear part.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.AffineTransforms.intercept",
    "page": "Reference",
    "title": "TwoDimensional.AffineTransforms.intercept",
    "category": "function",
    "text": "intercept(A) returns the tuple (x,y) such that A(x,y) = (0,0).\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.AffineTransforms.compose",
    "page": "Reference",
    "title": "TwoDimensional.AffineTransforms.compose",
    "category": "function",
    "text": "compose(A,B) yields the affine transform which combines the two affine transforms A and B, that is the affine transform which applies B and then A.  Composition is accessible via: A∘B, A*B or A⋅B (\"∘\" and \"⋅\" can be typed by \\circ<tab> and \\cdot<tab>).\n\nIt is possible to compose more than two affine transforms.  For instance, compose(A,B,C) yields the affine transform which applies C then B, then A.\n\n\n\n\n\n"
},

{
    "location": "reference/#Affine-2D-Coordinate-Transforms-1",
    "page": "Reference",
    "title": "Affine 2D Coordinate Transforms",
    "category": "section",
    "text": "AffineTransform\nscale\nrotate\ntranslate\njacobian\nintercept\ncompose"
},

{
    "location": "reference/#TwoDimensional.AbstractPoint",
    "page": "Reference",
    "title": "TwoDimensional.AbstractPoint",
    "category": "type",
    "text": "Any object whose type is derived from AbstractPoint{T} has at least 2 fields: x its abscissa and y its ordinate, both of type T.\n\nSee also: Point, WeightedPoint.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.Point",
    "page": "Reference",
    "title": "TwoDimensional.Point",
    "category": "type",
    "text": "Point(x,y)\n\nyields an instance of a 2D point of coordinates (x,y).\n\nA point may be multiplied or divided by a scalar to scale its coordinates.  The addition (resp. subtraction) of two points adds (resp. subtracts) their coordinates.\n\nCoordinates can be specified by keywords:\n\nPoint(x=xval, y=yval)\n\nThere are no default values for keywords x and y so both must be specified.\n\nSee also: WeightedPoint, AbstractPoint.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.WeightedPoint",
    "page": "Reference",
    "title": "TwoDimensional.WeightedPoint",
    "category": "type",
    "text": "A WeightedPoint{T} has just 3 fields: w its weight, x its abscissa and y its ordinate, all of type T.  By convention w ≥ 0 but this is not checked for efficiency reasons.\n\nSee also: Point, AbstractPoint.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.distance",
    "page": "Reference",
    "title": "TwoDimensional.distance",
    "category": "function",
    "text": "distance(A, B)\n\nyields the Euclidean distance between the 2 points A and B.\n\n\n\n\n\n"
},

{
    "location": "reference/#Points-1",
    "page": "Reference",
    "title": "Points",
    "category": "section",
    "text": "AbstractPoint\nPoint\nWeightedPoint\ndistance"
},

{
    "location": "reference/#TwoDimensional.BoundingBox",
    "page": "Reference",
    "title": "TwoDimensional.BoundingBox",
    "category": "type",
    "text": "BoundingBox(xmin,xmax,ymin,ymax) yields an instance of a 2D rectangular bounding-box whose sides are aligned with the coordinate axes and containing points of coordinates (x,y) such that xmin ≤ x ≤ xmax and ymin ≤ y ≤ ymax.  The box is empty if xmin > xmax or ymin > ymax.\n\nA bounding-box can be constructed from the first and last points (i.e. at the lower-left and upper right opposite corners) of the box:\n\nBoundingBox(P0::Point, P1::Point)\nBoundingBox(I0::CartesianIndex{2}, I1::CartesianIndex{2})\n\nCoordinates can be specified by keywords:\n\nBoundingBox(xmin=x0, ymin=y0, xmax=x1, ymax=y1)\n\nThere are no default values for keywords xmin, xmax, ymin and ymax so all must be specified.\n\nSee also Point, interior, exterior.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.area",
    "page": "Reference",
    "title": "TwoDimensional.area",
    "category": "function",
    "text": "area(B) yields the area of the bounding-box B.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.center",
    "page": "Reference",
    "title": "TwoDimensional.center",
    "category": "function",
    "text": "center(B::BoundingBox) -> c::Point\n\nyields the central point of the bounding-box B.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.interior",
    "page": "Reference",
    "title": "TwoDimensional.interior",
    "category": "function",
    "text": "interior([T,] B)\n\nyields the largest bounding -box with integer valued bounds which is contained by the bounding-box B.  Optional argument T is to specify the type of the result or of the coordinates of the result which is the same as B by default.\n\nSee also: exterior, round.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.exterior",
    "page": "Reference",
    "title": "TwoDimensional.exterior",
    "category": "function",
    "text": "exterior([T,] B)\n\nyields the smallest bounding-box with integer valued bounds which contains the bounding-box B.  Optional argument T is to specify the type of the result or of the coordinates of the result which is the same as B by default.\n\nSee also: interior, round.\n\n\n\n\n\n"
},

{
    "location": "reference/#Bounding-Boxes-1",
    "page": "Reference",
    "title": "Bounding-Boxes",
    "category": "section",
    "text": "BoundingBox\narea\ncenter\ninterior\nexterior"
},

{
    "location": "reference/#Base.round-Tuple{Point}",
    "page": "Reference",
    "title": "Base.round",
    "category": "method",
    "text": "round([T,] obj)\n\nyields the object that is the nearest to obj by rounding its coordinates to the nearest integer.  Argument T can be the type of the result (a point or a bounding-box) or the type of the coordinates of the result.\n\nSee also: interior, exterior.\n\n\n\n\n\n"
},

{
    "location": "reference/#Methods-1",
    "page": "Reference",
    "title": "Methods",
    "category": "section",
    "text": "round(::Point)"
},

{
    "location": "reference/#TwoDimensional.AbstractPoint2D",
    "page": "Reference",
    "title": "TwoDimensional.AbstractPoint2D",
    "category": "type",
    "text": "AffineTransform{T} is an alias for TwoDimensional.AbstractPoint{T}.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.Point2D",
    "page": "Reference",
    "title": "TwoDimensional.Point2D",
    "category": "type",
    "text": "Point2D{T} is an alias for TwoDimensional.Point{T}.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.WeightedPoint2D",
    "page": "Reference",
    "title": "TwoDimensional.WeightedPoint2D",
    "category": "type",
    "text": "WeightedPoint2D{T} is an alias for TwoDimensional.WeightedPoint{T}.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.BoundingBox2D",
    "page": "Reference",
    "title": "TwoDimensional.BoundingBox2D",
    "category": "type",
    "text": "BoundingBox2D{T} is an alias for TwoDimensional.BoundingBox{T}.\n\n\n\n\n\n"
},

{
    "location": "reference/#TwoDimensional.AffineTransform2D",
    "page": "Reference",
    "title": "TwoDimensional.AffineTransform2D",
    "category": "type",
    "text": "AffineTransform2D{T} is an alias for TwoDimensional.AffineTransform{T}.\n\n\n\n\n\n"
},

{
    "location": "reference/#Aliases-1",
    "page": "Reference",
    "title": "Aliases",
    "category": "section",
    "text": "AbstractPoint2D\nPoint2D\nWeightedPoint2D\nBoundingBox2D\nAffineTransform2D"
},

]}
