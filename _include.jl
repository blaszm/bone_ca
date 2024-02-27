using Ferrite, FerriteGmsh
using CoherentStructures, Distances # kommt bald raus
using BlockArrays, SparseArrays, LinearAlgebra, Tensors
using TimerOutputs, Printf, HDF5
using KrylovMethods
using Distributed

# careful with order!!

#include("mesh/rvemesh_constr.jl")
include("mesh/mesh.jl"); # mesh tools micro system
include("mesh/geometry_check.jl")
include("macroscale/material_macro.jl") # material parameters macro system
include("material/material.jl"); # material parameters micro system
include("material/microproblem.jl"); # construction of micro system
include("macroscale/tangent_macro.jl")
include("macroscale/assemble_macro.jl")
include("macroscale/elmt_macro.jl")
include("macroscale/elmtair_macro.jl") # element routine air
include("macroscale/mesh_macro.jl") # mesh tools macro system
include("macroscale/model_macro.jl")

include("material/applypbc.jl") # periodic boundary conditions helper # kommt bald raus

include("material/assemble.jl") # assembly micro system
include("material/elmtcb.jl") # element routine cortical bone
include("material/elmtbm.jl") # element routine bone marrow
include("material/solve_RVE.jl") # micro model solver
include("material/states.jl") # output state and flux quantities from states
include("material/volaver.jl") # calculate volume averages of quantities
include("output/writeinfo.jl") # output info
include("output/writeoutput.jl") # output Paraview
include("output/hdffiles.jl") # read, write and copy hdf-files