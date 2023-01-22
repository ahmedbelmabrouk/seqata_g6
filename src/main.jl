# Pour des raison d'efficacité de chargement, on ne charge que le fichiers 
# nécessaire
action = Args.get(:action)
# @show action

action == :annealing  && @ms include("main_annealing.jl")
action == :carlo      && @ms include("main_carlo.jl")
action == :descent    && @ms include("main_descent.jl")
action == :dmip       && @ms include("main_dmip.jl")
action == :explore    && @ms include("main_explore.jl")
action == :mip        && @ms include("main_mip.jl")
action == :stats      && @ms include("main_stats.jl")
action == :steepest   && @ms include("main_steepest.jl")
action == :test       && @ms include("main_test.jl")
action == :timing     && @ms include("main_timing.jl")
action == :validate   && @ms include("main_validate.jl")
action == :vns        && @ms include("main_vns.jl")

"""
    main()
Méthode principale de lancement des actions.

L'exécution de chaque action est sous-traitée à une méthode de la forme
`main_action`.
Par exemple l'action `carlo` de la ligne de commande (CLI) entrainera l'appelle à
à la méthode [`main_carlo()`](@ref).

"""
function main()
    # @ms Args.parse_commandline(ARGS) # Déjà fait dans le module
    ln1("main() BEGIN") # ln1 n'est utilisable que si level est connu
    lg1() && Args.show_args()

    # date1= now() # en secondes entières
    time1 = time() # en secondes, précision microsecondes

    action = Args.get(:action)
    println("DANS main() : action=$action")

    DONE = false
    if action == :validate
        main_validate()
    elseif action == :timing
        main_timing()
    elseif action == :carlo
        main_carlo()
    elseif action == :explore
        main_explore()
    elseif action == :descent
        main_descent()
    elseif action == :steepest
        main_steepest()
    elseif action == :vns
        main_vns()
    elseif action == :annealing
        main_annealing()
    elseif action == :mip
        main_mip()
    elseif action == :dmip
        main_dmip()
    elseif action == :stats
        main_stats()
    elseif action == :test
        try
            main_test()
        catch e
            println(join(stacktrace(), "\n\n"))
            println("\nERREUR: relancer main_test")
        end
    elseif action == :none
        println("Aucune action indiquée")
        # println("Actions possibles : ", join(actions(), ","))
        println(Args.get_syntaxe())
        # println(Args.get_actions())
        # exit(1)
    else
        println("Erreur : action $(action) non implémentée (dans main.jl)")
        println(Args.get_syntaxe())
        exit(1)
    end
    if DONE
        # heure de fin du traitement
        time2 = time()
        sec = round((time() - time1), digits = 3) # on veut limiter la précision à la ms
        ln1("Durée totale du main 1000*(time2-time1) : $(sec)s")
        ln1("main() END")
    end
end
