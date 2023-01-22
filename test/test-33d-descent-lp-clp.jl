
# Création d'une struct encapsulant la spécification d'un test LP
specs = [
    LpSolverSpec(:lp, :clp, LpTimingSolver, Clp.Optimizer),
    #    LpSolverSpec(:lp1, :clp, Lp1TimingSolver, Clp.Optimizer),
    #    LpSolverSpec(:lp2, :clp, Lp2TimingSolver, Clp.Optimizer),
    #    LpSolverSpec(:lp3, :clp, Lp3TimingSolver, Clp.Optimizer),
    #    LpSolverSpec(:lp4, :clp, Lp4TimingSolver, Clp.Optimizer),
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
