export new_lp_model

# Quelques fonctions utilitaires pour les manipulations de solvers PL ou MIP

# Création d'un modèle de solver externe (:cplex, :glpk, ...)
# Options possibles :
#  solver : symbole :cplex, :glpk, :clp, :cbc
#     :cbc : permet de faire du MIP
#     :clp : resteint aux problème LP mais plus efficace de cbc pour les
#            relaxations.
#  mode : (:lp ou :mip) utile pour choisir entre :clp ou :cbc
#  log_level (expérimental) : entier 0 pour silence, sinon 1, 2, ...
#
# Autre solver à supporter plus tard Gurobi (concurrent de cplex)
#
# DEPENDANCE EXTERNE : module Args (pour Args.get("external_lp_solver"))
#
function new_lp_model(;solver = :auto,
                     mode = :lp,
                     log_level = 0, # EXPLOITATION PARTIELLE
                     time_limit = 0,
                     )
    if !(mode in [:mip, :lp])
        error("ERROR: unknown mode $(mode). Should be :cplex, :glpk, ...")
    end
    if solver == :auto
        solver = Args.get("external_lp_solver")
    end
    if solver in [:clp, :cbc]
        # solver = (mode == :mip ? :cbc : :clp )
        if mode==:mip
            solver = :cbc
        else
            solver = :clp
        end
    end
    if ! (solver in [:cplex, :clp, :cbc, :glpk, :gurobi, :tulip])
        error("ERROR: unknown solver $(solver). Should be :cplex, :glpk, ...")
    end

    if solver == :glpk
        # voir https://github.com/JuliaOpt/GLPK.jl pour les options
        model = JuMP.Model(GLPK.Optimizer)
        # CPX_PARAM_TILIM
    elseif solver == :cplex
        model = JuMP.Model(CPLEX.Optimizer)
        if time_limit != 0
            set_optimizer_attribute(model, "CPX_PARAM_TILIM", time_limit)
        end
    elseif solver == :clp
        # voir https://github.com/JuliaOpt/Clp.jl pour les options
        model = JuMP.Model(Clp.Optimizer)
        # set_optimizer_attribute(model, "LogLevel", 0) # (MajusCule pour Clp)
    elseif solver == :cbc
        # Voir https://github.com/JuliaOpt/Cbc.jl pour les options
        model = JuMP.Model(Cbc.Optimizer)
        # set_optimizer_attribute(model, "logLevel", 0) # (minusCule pour Cbc)
    elseif solver == :tulip
        # voir https://ds4dm.github.io/Tulip.jl/stable/reference/options/ pour les options
        # voir https://github.com/ds4dm/Tulip.jl pour les sources
        model = JuMP.Model(Tulip.Optimizer)
    elseif solver == :gurobi
        # Pour les options, voir :
        # - https://github.com/JuliaOpt/Gurobi.jl
        # - https://www.gurobi.com/documentation/current/refman/parameters.html
        global GRB_ENV
        if !@isdefined(GRB_ENV)
            GRB_ENV = Gurobi.Env() # affiche warning pénible (academic)
        end
        model = JuMP.Model(() -> Gurobi.Optimizer(GRB_ENV))
        if time_limit != 0
            set_optimizer_attribute(model, "TimeLimit", time_limit)
        end
    else
        error("LpTimingSolver: unknown external_lp_solver inconnu : $solver")
    end
    if log_level <= 3
        JuMP.set_silent(model) # JuMP.unset_silent(model) : réautorise l'affichage
        # MOI.set(model, MOI.Silent(), true) # Variante de bas niveau
    end
    return model
end
