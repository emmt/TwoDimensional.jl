SOURCES = \
    Project.toml \
    make.jl \
    src/xforms.md \
    src/boxes.md \
    src/index.md \
    src/install.md \
    src/points.md \
    src/reference.md

default: build/index.html

clean:
	rm -f *~ */*~ */*/*~

build/index.html: $(SOURCES)
	julia --color=yes make.jl

.PHOONY: default clean
