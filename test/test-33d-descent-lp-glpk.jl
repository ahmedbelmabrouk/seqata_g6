# Création d'une struct encapsulant la spécification d'un test LP
specs = [
    LpSolverSpec(:lp, :glpk, LpTimingSolver, GLPK.Optimizer),
    #    LpSolverSpec(:lp1, :glpk, Lp1TimingSolver, GLPK.Optimizer),
    #    LpSolverSpec(:lp2, :glpk, Lp2TimingSolver, GLPK.Optimizer),
    #    LpSolverSpec(:lp3, :glpk, Lp3TimingSolver, GLPK.Optimizer),
    #    LpSolverSpec(:lp4, :glpk, Lp4TimingSolver, GLPK.Optimizer),
]

# ===========
# Préparartion des arguments
Args.set(:itermax, 300)

# ===========
inst = instance_build_mini10()
@test inst.nb_planes == 10

for spec in specs
    local test_name = "Test descente LP $(spec.algo):$(spec.external_lp_solver)"
    lg1("$test_name ... ")
    @testset "$test_name" begin
        Log.pushlevel!(0) # on passe temporairement à 0
        test_one_lp_descent(spec)
        Log.poplevel!()
    end
    ln1(" => fait.")
end
