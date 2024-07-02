#! /bin/bash
exe=$(realpath "$0")
dir=$(dirname "$exe")
cd "$dir/.."
#pkg=$PWD
#cd "$dir"
#echo "{$exe}{$dir}{$pkg}"
exec julia --color=yes --project="docs/" -e \
    'using Pkg;
     Pkg.develop(PackageSpec(path=pwd()));
     Pkg.instantiate();
     cd("docs");
     include("make.jl");'
