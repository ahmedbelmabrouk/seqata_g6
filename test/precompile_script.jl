#!/bin/sh
#= La ligne shell suivante est en commentaire multiligne pour julia
exec julia --project --color=yes --startup-file=no --depwarn=no -- "$0" "$@"
=#
# Ce fichier est destininer à exécuter le maximum de fonction possible de
# cette application tout en prennant le moins de temos possible.
# Son but est d'être utilise pour une précompilation manulle
# (par exemple par appel au script bin/buils_seqta_sysimg.jl)
# Il est constitué à partir du version simplifié de runtests.jl.

# Le realpath est nécessaire si l'un des répertoires parents est un lien symbolique :
this_appdir = dirname(dirname(realpath(@__FILE__())))
using Pkg
Pkg.activate(dirname(dirname(realpath(@__FILE__()))))
Pkg.instantiate()

ENV["JULIA_USING_ALL"] = 1  # for loading all package

pushfirst!(ARGS, "none") # On précise l'action par défaut
include("$this_appdir/src/Seqata.jl")
using .Seqata
using .Log  # pour pouvoir taper ln1() au lieu de Log.ln1()

println("BEGIN all tests at ", ms(), "sec")
using Test
using Random
using Printf
using Dates


function test_loggins()
    lg1("test_loggins ($(ms()))... ")
    # Activation des méthodes de logging
    lg0();lg1();lg2();lg3();lg4();lg5();
    ln0();ln1();ln2();ln3();ln4();ln5();
    Log.pushlevel!(Log.level() + 1)
    Log.poplevel!()

    Args.get(:loglevel)
    Args.set(:loglevel, 0)
    ln1(" fait ($(ms()))")
end


function test_solution_stp_solver()
    ln1("test_solution_stp_solver ($(ms()))... ")
    inst = instance_build_mini10()

    # construction et résolution (update=true par défaut)
    sol = Solution(inst, update=true, algo=:earliest)
    shuffle!(sol)
    ok = is_feasable(sol)
    nbviols, violtxt = get_viol_description(sol)
    # set_from_names!(sol, shuffle!(sol.planes)) # JE VEUX LES plane.name
    set_from_names!(sol, getfield.(shuffle!(sol.planes), :name) )

    # Exécution des algos de TimingSolcer autonomes
    algos = [:dp, :faye]
    for algo in algos
        lg1("test de $(algo) ($(ms()))... ")
        shuffle!(sol)
        try
            sol = Solution(inst, update=true, algo=algo)
            ln1("fait.")
        catch
            ln1()
            println(join(stacktrace(), "\n\n"))
            println("problème avec $(algo) => précompilation ignorée")
        end
    end

    # Exécution des algos de TimingSolcer dépendants d'un solver externe
    xalgos = [:lp] # un seul algo lp pour Seqata mais lp2 lp3... pour Alap
    xsolvers = [:cplex, :gurobi, :clp, :cbc, :tulip, :glpk]
    for algo in xalgos, xsolver in xsolvers
        lg1("test de $(algo)/$(xsolver) ($(ms()))... ")
        try
            Args.set("external_lp_solver", xsolver)
            sol = Solution(inst, algo=algo)
            ln1("fait.")
        catch
            ln1()
            println(join(stacktrace(), "\n\n"))
            println("problème avec $(algo)/$(xsolver) => précompilation ignorée")
        end
    end


    # Relecture d'un fichier solution
    solfile = "$this_appdir/test/data/alp_01.sol"
    sol = Solution(inst, solfile)

    # ln1(" fait ($(ms()))")
    ln1("test_solution_stp_solver ($(ms()))... ait ($(ms()))")
end

function test_explore_mini()
    lg1("test_explore_mini ($(ms()))... ")

    inst = instance_build_mini10()
    Log.pushlevel!(0) # on passe temporairement à 0
    solver = ExploreSolver(inst)
    initial_sort!(solver.cursol, presort = :shuffle)
    copy!(solver.bestsol, solver.cursol)

    itermax = 10
    solve!(solver, itermax)
    Log.poplevel!() # remise au loglevel initial
    ln1(" fait ($(ms()))")
end


function test_descent_mini()
    lg1("test_descent_mini ($(ms()))... ")
    inst = instance_build_mini10()
    Log.pushlevel!(0) # on passe temporairement à 0
    solver = DescentSolver(inst)
    initial_sort!(solver.cursol, presort = :shuffle)
    copy!(solver.bestsol, solver.cursol)

    solve!(solver, nb_cons_reject_max=2)
    Log.poplevel!() # remise au loglevel initial
    ln1(" fait ($(ms()))")
end


function test_mutation()
    lg1("test_mutation ($(ms()))... ")
    inst = instance_build_mini10()
    Log.pushlevel!(0) # on passe temporairement à 0

    muts,label = generate_nbh(n, "D2+2T2+2T2g1+P3")

    Log.poplevel!() # remise au loglevel initial
    ln1(" fait ($(ms()))")
end



function main()
    methods = [
        test_loggins;
        test_solution_stp_solver;
        test_explore_mini;
        test_descent_mini;
    ]
    # ln0("Liste des méthodes de test à précompiler")
    # for m in methods
    #     ln0("m=$(m)")
    # end

    for m in methods
        ln0("method=$(m) => ($ms()s)")
        m()
    end
end

main()
println("END all tests at ", ms(), "sec")
