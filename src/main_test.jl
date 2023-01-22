# @ms include("xxxx.jl")

export main_test

# Cette action est destinée à la mise au point de code (nouveaux solveurs,
# nouveaux voisinages, ...) 
# Une fois la fonctionnalité opérationnelle, il faut en faire en test unitaire
# (i.e. en créant un nouveau fichier dans tests/text_xxx.jl)
"""
    main_test()
Méthode principale d'exécution de l'action `test`.

Cette action est destinée à la mise au point de code (nouveaux solveurs,
nouveaux voisinages, ...) 
Une fois la fonctionnalité opérationnelle, il faut en faire en nouveau 
test unitaire (i.e. en créant un nouveau fichier dans tests/text_xxx.jl)
de façon à assurer la non régression de la fonctionnalité.
"""
function main_test()

    if Args.get(:infile) == "NO_INFILE"
        Args.set(:infile, "data/alpx/01.alpx")
    end
    println("="^70)
    println("Début de l'action test")
    inst = Instance(Args.get(:infile))

    sol = Solution(inst)
    println(to_s(sol)) # OK
    solve!(sol)
    println(to_s(sol)) # OK
    Random.shuffle!(sol.planes)
    println(to_s(sol)) # OK
    solve!(sol)
    println(to_s(sol)) # OK

    println("Fin de l'action test")
end
