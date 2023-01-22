export main_timing

"""
    main_timing()
Méthode principale d'exécution de l'action `timing`.

Cette méthode lit une instance, lit une liste de noms d'avion
puis résoud le sous-problème de timing (STP).
Le choix de l'algorithme de résolution du STP est défini par l'option
`--timing-algo-solver` alias `-t` qui est exploité directement par le
constructeur de l'objet `Solution` (TODO FAIT LIEN IEN API DE SOLUTION).

La solution résolue est alors affichée sur la sortie standard et 
enregistrée.

# Exemple

    ./bin/run.jl tim -t lp -i data/01.alp  -p 3,4,5,6,7,8,9,1,10,2
    => 700.0

"""
function main_timing()
    println("Début de l'action timing")

    # @error "main_timing: désolé la méthode main_timing n'est pas implantée"
    # Résolution de l'action
    println("="^70)
    println("Résolution du timing pour $(Args.get(:infile))\n")
    inst = Instance(Args.get(:infile))
    # args.names=[1,6,8,4,12,9,11,3,10,2,19,24,20,7,5,50,15,23,18,14,13,25,17,26,
    #          43,16,44,27,32,28,22,33,29,47,34,48,49,46,38,35,45,39,31,40,36,30,42,41,37,21]
    # args.names = %w(3 4 5 6 8 9 7 1 10 2)
    # puts inst.to_s
    if inst.nb_planes != length(Args.get(:planes))
        println("\nERREUR taille (=$(inst.nb_planes) de l'instance différente ")
        println(" de la taille (=$(length(Args.get(:planes))) de la solution")
        println()
        exit(1)
    end
    sol = Solution(inst, update=false)
    set_from_names!(sol, Args.get(:planes))

    solve!(sol, do_update_cost=true)
    print_sol(sol)

    lg1("Création du fichier \"$(guess_solname(sol))\"... ")
    write(sol)


    println("Fin de l'action timing")
end
