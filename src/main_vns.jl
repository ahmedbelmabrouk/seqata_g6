@ms include("vns_solver.jl")

export main_vns

"""
main_vns()
Méthode principale d'exécution de l'action `vns`.

Cette méthode lit l'instance et créer un objet sv::VnsSolver.
exploite les paramètres `--duration` et `--itermax`,
sous-traite la résolution à la méthode solve! de l'objet sv.
Elle récupère alors le meilleur résultat obtenu, l'enregistre et
affiche des statistiques de la résolution.
"""
function main_vns()
    # Résolution de l'action
    println("="^70)
    println("Début de la méthode main_vns")
    inst = Instance(Args.get(:infile))

    sv = VnsSolver(inst)

    # duration = 3600.0 # secondes
    duration = 1.0*Args.get(:duration)
    itermax = Args.get(:itermax)

    # Voir aussi option startsol de la méthode solve!
    ms_start = ms() # nb secondes depuis démarrage avec précision à la ms
    solve!(sv, durationmax = duration)
    ms_stop = ms()

    bestsol = sv.bestsol
    write(bestsol)
    print_sol(bestsol)

    nb_calls = bestsol.solver.nb_calls
    nb_infeasable = bestsol.solver.nb_infeasable
    nb_sec = round(ms_stop - ms_start, digits = 3)
    nb_call_per_sec = round(nb_calls / nb_sec, digits = 3)
    println("Performances: ")
    println("  nb_calls=$nb_calls")
    println("  nb_infeasable=$nb_infeasable")
    println("  nb_sec=$nb_sec")
    println("  => nb_call_per_sec = $nb_call_per_sec call/sec")

    println("Fin de la méthode main_vns")
end
