# Chargement des packages de l'application.
# Deux cas sont possibles :
# - Soit on essaie de charger les seuls packages indispensables à l'action 
#   demandée par la ligne de commande (avec USING_ALL == false)
# - soit on ne connait pas l'action à exécuter (e.g (i.e en mode test ou 
#   en mode interactif) et on charge alors tous les packages, ce qui est plus
#   long (avec USING_ALL == true).
#
# En effet les packages tels que JuMP, CPLEX, Cbc, ... sont assez lents à charger
# et ne sont pas toujours indispensables selon les options.

# Fichier de chargements CONDITIONNEL des packages stictement nécessaires
# en fontion de l'action demandée sur la ligne de commande.

MS_LOG && ln1("===== USING USING_ALL==$(USING_ALL) ($(ms()))")

if USING_ALL

    # JuMP et les solveurs mathématiques externes
    @ms using JuMP
    const MOI = JuMP.MathOptInterface

    # Not all of these solvers are to be available
    try @ms using CPLEX catch end
    try @ms using GLPK catch end
    try @ms using Cbc catch end
    try @ms using Clp catch end
    try @ms using Gurobi catch end
    try @ms using Tulip catch end

    # Une fonctions utile pour la création d'un modèle multi-solveur
    @ms include("lp_model_util.jl")

    # Les timingSolver
    # @ms include("dynprog_timing_solver.jl")
    @ms include("earliest_timing_solver.jl")
    # @ms include("faye_timing_solver.jl")
    @ms include("lp_timing_solver.jl")

    # Les fichiers de lancement de chaque actions (avec leurs propre usings)

    @ms include("main_annealing.jl")
    @ms include("main_carlo.jl")
    @ms include("main_descent.jl")
    @ms include("main_dmip.jl")
    @ms include("main_explore.jl")
    @ms include("main_mip.jl")
    @ms include("main_stats.jl")
    @ms include("main_steepest.jl")
    @ms include("main_test.jl")
    @ms include("main_timing.jl")
    @ms include("main_validate.jl")
    @ms include("main_vns.jl")
    
    @ms include("interactive.jl") # ajout le 24/12/2021 pour précompilation

else
    # Doit-on charger JuMP ?
    # oui si un des cas suivant est vrai
    # - on utilise une approche frontale du problème complet (mip ou dmpi)
    # - un des solveurs du sous pb de timing utilise la PL
    #   (e.g. LpTimingSolver, Lp2TimingSolver, ...)
    #
    _tas = Args.get("timing_algo_solver")  # e.g.  :lp, :earliest, :dp, ...
    _xsv = Args.get("external_lp_solver") # e.g.  :cplex, :glpk, ...
    _need_jump = false

    if Args.get("action") in [:mip, :dmip]
        _need_jump = true
    end

    if occursin(r"^lp", String(_tas)) # cat on peut avoit lp1, lp2 ...
        _need_jump = true
    end

    if MS_LOG && ln1()
        println("   action=$(Args.get("action"))")
        println("   _need_jump=$_need_jump")
        println("   external_lp_solver: _xsv=$(_xsv))")
        println("   timing_algo_solver:  _tas=$(_tas))")
        println("   USING_ALL=$(USING_ALL)")
    end

    if _need_jump
        MS_LOG && ln1("===== USING in _need_jump")
        @ms using JuMP
        const MOI = JuMP.MathOptInterface

        @ms include("lp_model_util.jl")

        # Le solveur externe à charger dépend de l'option :external_lp_solver
        if _xsv == :glpk
            @ms using GLPK
        elseif _xsv == :cplex
            @ms using CPLEX
        elseif _xsv in [:clp, :cbc]
            @ms using Cbc
            @ms using Clp
        elseif _xsv == :gurobi
            @ms using Gurobi
        elseif _xsv == :tulip
            @ms using Tulip
        else
            error("external_lp_solver inconnu : $_xsv")
        end
    end

    # Chargement des TimingSolvers (pour la brique STP)
    # _tas == :dp       &&  @ms include("dynprog_timing_solver.jl")
    _tas == :earliest &&  @ms include("earliest_timing_solver.jl")
    # _tas == :faye     &&  @ms include("faye_timing_solver.jl")
    _tas == :lp       &&  @ms include("lp_timing_solver.jl")

end

if isinteractive()
    # Quelques définition ou macros utiles seulement en mode interactifs
    @ms include("interactive.jl")
end

# CHARGEMENT FINAL DE LA MÉTHODE METHODE PRINCIPALE MAIL !
# 
@ms include("main.jl")

MS_LOG && ln1("===== USINGS END ($(ms()))")
