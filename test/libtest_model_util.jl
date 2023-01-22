
# Ces usings sont nécessaires pour utilisation de tests isolés !
@ms using JuMP
try
    @ms using CPLEX
catch
end
try
    @ms using GLPK
catch
end
try
    @ms using Cbc
catch
end
try
    @ms using Clp
catch
end

# Création d'une struct encapsulant la spécification d'un test LP
struct LpSolverSpec
    algo::Union{Symbol,Nothing}
    external_lp_solver::Union{Symbol,Nothing}
    lp_timing_solver::Union{DataType,Nothing}
    backend_model::Union{DataType,Nothing}
end

function test_one_lp_descent(spec)

    # ===========
    ln2("Création sol avec algo ")
    lg2(spec.algo, ":", spec.external_lp_solver)
    lg2(" (", spec.lp_timing_solver, "... ")
    Args.set(:external_lp_solver, spec.external_lp_solver)
    sol = Solution(inst, algo = spec.algo)
    model = sol.solver.model
    initial_sort!(sol, presort = :shuffle)

    @test isa(sol.solver, spec.lp_timing_solver)
    # On vérifie par exemple que le nom du solver externe "Clp" effectif
    # correspondant bien au symbole :clp du solveur demandé
    @test Symbol(lowercase(solver_name(model))) == spec.external_lp_solver
    lg2("ok (avec sol.cost=$(sol.cost))\n")


    # ===========
    lg2("Création DescentSolver nb_cons_reject_max=$(Args.get(:itermax))")
    lg2(" (cost=700<=850.0?)... ")
    sv = DescentSolver(inst)
    # println("\nconstruction faite")
    sv.do_save_bestsol = false
    sv.nb_cons_reject_max = Args.get(:itermax) # Devra être écrasé
    sv.durationmax = 1.0                       # Devra être écrasé

    Log.pushlevel!(0) # on passe temporairement à 0
    solve!(sv, startsol = sol, nb_cons_reject_max = 200, durationmax = 3.0)
    Log.poplevel!()

    @test sv.durationmax == 3.0
    @test sv.nb_cons_reject_max == 200
    @test isa(sv.bestsol.solver, spec.lp_timing_solver)
    # @test sv.bestsol.cost <= 850.0
    # @test sv.bestsol.cost <= 900.0
    # @test sv.bestsol.cost >= 700.0
    @test 700 <= sv.bestsol.cost <= 1_000.0
    ln2("ok ($(sv.bestsol.cost))")

end
