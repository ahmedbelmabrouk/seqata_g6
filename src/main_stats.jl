export main_stats

"""
    main_stats()
Méthode principale d'exécution de l'action `stats`.

Cette méthode lit l'instance.
Si ne niveau de verbosité est suffisant, réaffiche l'instance dans
différent format (dont AMPL).
Dans tous les cas des statistiques sur l'instance sont affichées.
"""
function main_stats()
    println("="^70)
    println("Début de l'action stats\n")

    # Hack provisoire et pas très propre pour pouvoir traiter plusieurs 
    # fichiers d'entrée
    # - infile : contient en principe le chemin de l'instance
    # - outfile : contient en principe le chemin de la solution en sortie
    # - files : contient les chemins éventuels supplémentaire
    # Ici on accumule tous ces chemins comme les instance à analyser

    files = String[Args.get(:infile)]
    if Args.get(:solfile) != nothing
        push!(files, Args.get(:solfile))
    end
    if length(Args.get(:files)) != 0
        append!(files, Args.get(:files))
    end

    # @show(files)
    for file in files
        inst = Instance(file)
        if lg4()
            println("="^70)
            println("Regéneration de l'instance au format ampl")
            println(to_s_long(inst, format = "ampl"))
            println("="^70)
            println("Regéneration de l'instance au format alp")
            println(to_s_long(inst, format = "alp"))
            println("="^70)
            println("Regéneration de l'instance au format alpx")
            println(to_s_long(inst, format = "alpx"))
        end
        println("="^70)
        println(to_s_stats(inst, verbose=lg3()))
    end

    println("Fin de l'action stats")
end
