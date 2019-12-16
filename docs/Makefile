SOURCES = \
   src/AffineTransforms.md \
   src/BoundingBoxes.md \
   src/index.md \
   src/install.md \
   src/Points.md \
   src/reference.md

default: build/index.html

clean:
	rm -f *~ */*~ */*/*~

build/index.html: $(SOURCES)
	julia --color=yes make.jl

.PHOONY: default clean
