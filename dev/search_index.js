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
    "location": "#Table-of-contents-1",
    "page": "Introduction",
    "title": "Table of contents",
    "category": "section",
    "text": "Pages = [\"install.md\", \"AffineTransforms.md\"]"
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

]}
