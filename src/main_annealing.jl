@ms include("annealing_solver.jl")

export main_annealing

function main_annealing()
    ln1("="^70)
    ln1("Début de l'action annealing")
    inst = Instance(Args.get(:infile))

    # Construction de la solution initiale :
    sol = Solution(inst)
    ln1("Solution correspondant à l'ordre de l'instance")
    ln1(to_s(sol))

    # ON POURRAIT AUSSI REPARTIR DE LA SOLUTION DU GLOUTON INTELLIGENT 
    initial_sort!(sol)
    ln1("Solution initiale envoyée au solver")
    ln1(to_s(sol))

    if Args.get(:nb_cons_no_improv_max) != 0   # 0 pour automatique
        my_improv_max = Args.get(:nb_cons_no_improv_max)
    else # 0 pour automatique
        if Args.get(:itermax) != 0
            my_improv_max = Args.get(:itermax)
        else # 0 pour automatique
            my_improv_max = 1_000_000_000 # infini
        end
    end


    # Choix des options pour le solver
    user_opts = Dict(
        # :startsol           => nothing,  # nothing pour auto à partir de l'instance
        :startsol => sol,  # nothing pour auto à partir de l'instance
        :step_size => inst.nb_planes,   # à renommer en step_size
        # :temp_init          => -1.0, # -1.0 pour automatique
        :temp_init => nothing, # nothing pour automatique
        :temp_init_rate => 0.30,  # valeur standard : 0.8
        :temp_mini => 0.000_001,
        # :temp_coef          => 0.999_95,
        :temp_coef => 0.95,
        :nb_cons_reject_max => 1_000_000_000, # infini
        # :nb_cons_no_improv_max => 500*inst.size*inst.size,
        :nb_cons_no_improv_max => my_improv_max,
    )

    sv = AnnealingSolver(inst, user_opts)
    ln1(get_stats(sv))

    ms_start = ms() # nb secondes depuis démarrage avec précision à la ms
    solve!(sv)
    ms_stop = ms()

    bestsol = sv.bestsol
    print_sol(bestsol)

    nb_calls = bestsol.solver.nb_calls
    nb_infeasable = bestsol.solver.nb_infeasable
    nb_sec = round(ms_stop - ms_start, digits = 3)
    nb_call_per_sec = round(nb_calls / nb_sec, digits = 3)
    println("Performance: ")
    println("  nb_calls=$nb_calls")
    println("  nb_infeasable=$nb_infeasable")
    println("  nb_sec=$nb_sec")
    println("  => nb_call_per_sec = $nb_call_per_sec call/sec")

    ln1("Fin de l'action annealing")
end
