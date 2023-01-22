@ms include("stupid_solver.jl")

export main_carlo

"""
    main_carlo()
Méthode principale d'exécution de l'action `carlo`.

Cette méthode exploite les arguments, construit des variables `cursol` et
`bestsol` de type Solution, crée une variable solveur de type `StupidSolver` et
appelle sa méthodes solve!(0 pour un nombre `itermax` d'itération.
"""
function main_carlo()
    println("="^70)
    println("Début de l'action carlo")

    # On peut afficher l'ensemble des options disponibles du programme (cf fichier args.jl)
    # Args.show_args()    # Déjà fait dans Args.jl

    # On construit l'instance dont le chemin est passé en paramètre
    inst = Instance(Args.get(:infile))

    # Si l'option itermax vaut 0 (i.e automatique) on choisit une valeur "pertinente"
    itermax_default = 500 * inst.nb_planes
    itermax = Args.get(:itermax) == 0 ? itermax_default : Args.get(:itermax)

    # Le timingSolver sous-jacent sera choisi automatiquement en fonction de l'option
    # --timing-algo-solver (alias -t )
    # (Cette option serait accessible par : Args.get(:timing_algo_solver) )
    cursol = Solution(inst, update = false)
    bestsol = Solution(inst, update = false)

    @show cursol
    @show bestsol
    @show itermax
    @show Log.level()

    ms_start = ms() # en seconde depuis le démarrage avec précision à la ms

    # Cas particulier : si itermax==1 : on tente le meilleur coup possible en triant
    # les avions sur leur target (ou sur l'ordre passé en paramètre)
    # IMPORTANT : voyez la méthode initial_sort! dans la classe/fichier solution.jl
    if itermax == 1
        if Args.get(:presort) == :none  # :none pour automatique
            initial_sort!(cursol, presort = :target)
        else
            initial_sort!(cursol, presort = Args.get(:presort))
        end
        # solve!(cursol)   # résolution du STP inutil car déjà fait par initial_sort
        copy!(bestsol, cursol)
    else
        sv = StupidSolver(inst)
        solve!(sv, itermax=itermax)
        copy!(bestsol, sv.bestsol)
    end

    ms_stop = ms()

    # le print_sol final n'est exécuté que si lg1() retourne true
    lg1() && print_sol(bestsol)

    nb_calls = bestsol.solver.nb_calls
    nb_infeasable = bestsol.solver.nb_infeasable
    nb_sec = round(ms_stop - ms_start, digits = 3)
    nb_call_per_sec = round(nb_calls / nb_sec, digits = 3)
    println("Performance: ")
    println("  nb_calls=$nb_calls")
    println("  nb_infeasable=$nb_infeasable")
    println("  nb_sec=$nb_sec")
    println("  => nb_call_per_sec = $nb_call_per_sec call/sec")

    println("Fin de l'action carlo")
end
