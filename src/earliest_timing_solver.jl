export EarliestTimingSolver, symbol, to_s, solve!

"""
    EarliestTimingSolver

Résoud du Sous-Problème de Timing en plaçant les avions le plus tôt possible.

Ce solveur est largement sous-optimum, mais il est rapide à calculer et fourni
à coup sûr une solution valide si l'ordre des avions donnés par la solutiom 
courante est réalisable.

Il est assoié à l'option `----timing-algo-solver earliest` (alias `-t earliest`).
"""
mutable struct EarliestTimingSolver
    # inst contient l'instance initiale et ses données
    inst::Instance
    nb_calls::Int
    nb_infeasable::Int

    # Le constructeur
    function EarliestTimingSolver(inst::Instance)
        this = new()
        this.inst = inst
        this.nb_calls = 0
        this.nb_infeasable = 0
        return this
    end
end

# Permettre de retrouver le nom de notre XxxxTimingSolver à partir de l'objet
function symbol(sv::EarliestTimingSolver)
    return :earliest
end

function Base.show(io::IO, sv::EarliestTimingSolver)
    print(io, to_s(sv))
end

# Fonction d'affichage, to string
function to_s(sv::EarliestTimingSolver)
    buf = IOBuffer()
    print(buf, '\n')
    print(buf, "  nb_calls : ", sv.nb_calls)
    print(buf, "  nb_infeasable : ", sv.nb_infeasable)
    String(take!(buf))
end


#= Résoud le sous-problème de timing (dates d'atterrissage) des avions à
 permutation fixée par l'ordre des avions de l'objet sol, puis met à jour
 la solution sol.
=#
function solve!(sv::EarliestTimingSolver, sol::Solution)
    sv.nb_calls += 1

    # On sous-traite le positionnement des dates x[] à l'une des méthodes de la
    # classe Solution.
    solve_to_earliest!(sol::Solution, do_update_cost = true)

    return nothing
end

# function test_earliest_timing_solver()
#     println("BEGIN EarliestTimingSolver")
#     inst = Instance("$(Args.appli)/data/ampl/01.ampl")
#     sv = EarliestTimingSolver(inst)
#     sol = Solution(inst)
#     solve!(sv, sol)
#     println(to_s_long(sol))
#     println("END EarliestTimingSolver")
# end
