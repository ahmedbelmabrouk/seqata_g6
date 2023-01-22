# EXEMPLE DE LIGNE DE COMMANDE TESTÉE :
# ./bin/run.jl tim -i data/ampl/01.ampl -p 3,4,5,6,7,8,9,1,10,2 -x clp
#

# ===========
# Construction de l'instance mini
inst = instance_build_mini10()
@test inst.nb_planes == 10

# ==========
# Ordre des avions à tester
planes_str = "p3,p4,p5,p6,p7,p8,p9,p1,p10,p2"
# planes = matchall(r"[\w]+", planes_str)
# planes = collect( (m.match for m = eachmatch(r"[\w]+", planes_str) ))
planes = [m.match for m in eachmatch(r"[\w]+", planes_str)]

test_name = "Création et résolution du timing d'une solution pour algo :earliest"
lg1(test_name, "... ")
@testset "$test_name" begin

    sol = Solution(inst, algo = :earliest, update = false)
    @test sol.solver == nothing

    Seqata.init_solver(sol, sol.timing_algo_solver)
    @test isa(sol.solver, EarliestTimingSolver)

    # Random.shuffle!(sol.planes); solve!(sol)
    shuffle!(sol, do_update = true)
    @test sol.cost > 1000.0  # cost >> 700 si mélange presque surement

    # On impose l'ordre des avions comme souhaité
    set_from_names!(sol, planes)

    # On doit retrouver cet ordre dans la solution
    @test join(get_names(sol), ",") == planes_str

    # On résout le sous-problème de timing de cette solution
    solve!(sol)
    # @test sol.cost == 700.0 # avec :lp
    @test sol.cost == 2830.0 # avec :earliest

end
ln1(" fait.")


# ====================================================================
test_name = "Création et résolution du timing d'une solution pour algo :lp avec :cplex"
lg1(test_name, "... ")
@testset "$test_name" begin

    Args.set("timing_algo_solver", :lp)
    Args.set("external_lp_solver", :cplex)
    # Args.set("planes", planes_str) # idem que l'option -p

    # sol = Solution(inst, update=false, algo=:lp)
    sol = Solution(inst, update = false)   # alpo=:lp par défaut
    @test sol.solver == nothing

    Seqata.init_solver(sol, sol.timing_algo_solver)
    @test isa(sol.solver, LpTimingSolver)

    @test Symbol(lowercase(solver_name(sol.solver.model))) == :cplex

    # On mélange les avions puis on met à jour la solution
    # print("Mélange et évaluation de la solution (cos>1000 ?)...")
    shuffle!(sol, do_update = true)
    @test sol.cost > 1000.0  # cost >> 700 car mélange presque surement

    # On impose l'ordre des avions comme souhaité
    set_from_names!(sol, planes)

    # On doit retrouver cet ordre dans la solution
    @test join(get_names(sol), ",") == planes_str

    # On résout le sous-problème de timing de cette solution
    solve!(sol)
    @test sol.cost == 700.0
end
ln1(" fait.")


# ====================================================================
test_name = "Création et résolution du timing d'une solution pour algo :lp avec :clp"
lg1(test_name, "... ")
@testset "$test_name" begin
    Args.set("external_lp_solver", :clp)
    sol = Solution(inst, update = false, algo = :lp)
    @test sol.solver == nothing

    Seqata.init_solver(sol, sol.timing_algo_solver)
    @test isa(sol.solver, LpTimingSolver)

    shuffle!(sol, do_update = true)
    @test sol.cost > 1000.0  # cost >> 700 si mélange presque surement

    # On impose l'ordre des avions comme souhaité
    set_from_names!(sol, planes)

    # On doit retrouver cet ordre dans la solution
    @test join(get_names(sol), ",") == planes_str

    # On résout le sous-problème de timing de cette solution
    solve!(sol)
    @test sol.cost == 700.0

end
ln1(" fait.")
