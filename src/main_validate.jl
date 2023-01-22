export main_validate

"""
    main_validate()
Méthode principale d'exécution de l'action `validate`.

Cette méthode lit une instance et une solution associée, 
Elle vérifie la conformité de cette solution vus-à-vs des contraintes
et affiche les différents viols éventuels.
"""
function main_validate()
    ln2("="^70)
    ln2("main_validate: BEGIN")

    # Résolution de l'action
    ln1("main_validate: fichier d'instance : $(Args.get(:infile))")
    inst = Instance(Args.get(:infile))

    ln1("main_validate: fichier de la solution : $(Args.get(:solfile))")
    sol = Solution(inst, Args.get(:solfile))

    ln1("main_validate: examen de la solution (par get_viol_description)")
    nbviols, violtxt = get_viol_description(sol)

    if nbviols == 0
        msg = "Solution correcte de coût : $(sol.cost)"
        println(to_sc(msg, :GREEN))
        ln3(to_s(sol))
    else
        msg = "Solution incorrecte : il y a $nbviols erreurs !"
        println(to_sc(msg, :RED))
        ln1(violtxt)
        ln3(to_s(sol))
    end

    ln1("main_validate: END")
    ln2("="^70)
end
