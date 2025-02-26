push!(LOAD_PATH, "modules/")
#push!(DEPOT_PATH, JULIA_DEPOT_PATH)
using Pkg

using Gurobi

import Data
import Parameters
import Formulations
import Heuristics
import Solution
import UpperBounds



# Read the parameters from command line
params = Parameters.readInputParameters(ARGS)

open(params.outputfile,"a") do f
    write(f,"$(params.instName)")
end


# Read instance data
if occursin(".col", params.instName)
    inst = Data.readDataCol(params.instName,params)
elseif occursin(".clq", params.instName)
    inst = Data.readDataClq(params.instName,params)
end


if params.method == "exact"
    cor,colors,order = Heuristics.coloringHeuristics(inst,params, params.problem == "cgcol")
    zeta = UpperBounds.newZeta(inst)
    deltaTwo = UpperBounds.deltaTwo(inst) + 1
    degSequence = UpperBounds.degSequence(inst, params)
    if params.form == "std"
        Formulations.stdFormulation(inst, params, colors, order, min(zeta,deltaTwo), degSequence)
    elseif params.form == "rep"
        Formulations.repFormulation(inst, params, colors, order, min(zeta,deltaTwo), degSequence)
    end
end
