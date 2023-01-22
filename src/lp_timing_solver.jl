export LpTimingSolver, symbol, solve!

"""
    LpTimingSolver

Résoud du Sous-Problème de Timing par programmation linéaire.

Ce solveur résoud le sous-problème de timing consistant à trouver les dates
optimales d'atterrissage des avions à ordre fixé.
Par rapport aux autres solvers (e.g DescentSolver, AnnealingSolver, ...), il
ne contient pas d'attribut bestsol

VERSION -t lp4 renommé en Lp : modèle diam version sans réoptimisation (coldrestart)
  - modèle diam simplifié par rapport à lp3 : sans réoptimisation (coldrestart)
  - pas de réoptimisation : on recrée le modèle à chaque nouvelle permu
    d'avion dans le solver
  - seules les contraintes de séparation nécessaires à la permu sont créées
  - gestion de l'option --assume_trineq true|false (true par défaut, cf lp1)
  - contraintes de coût simple (diam) : une seule variable de coût par avion
    plus un contrainte par segment :
       cost[i] >= tout_segment[i]
"""
mutable struct LpTimingSolver
    inst::Instance
    # Les attributs spécifiques au modèle
    model::Model  # Le modèle MIP
    x         # vecteur des variables d'atterrissage
    cost      # variable du coût de la solution
    costs     # variables du coût de chaque avion

    nb_calls::Int    # POUR FAIRE VOS MESURES DE PERFORMANCE !
    nb_infeasable::Int

    # Le constructeur
    function LpTimingSolver(inst::Instance)
        this = new()

        this.inst = inst

        # Création et configuration du modèle selon le solveur externe sélectionné
        this.model = new_lp_model() # SERA REGÉNÉRÉ DANS CHAQUE solve!()

        this.nb_calls = 0
        this.nb_infeasable = 0

        return this
    end
end

# Permettre de retrouver le nom de notre XxxxTimingSolver à partir de l'objet
function symbol(sv::LpTimingSolver)
    return :lp
end

function solve!(sv::LpTimingSolver, sol::Solution)

    #error("\n\nMéthode solve!(sv::LpTimingSolver, ...) non implantée: AU BOULOT :-)\n\n")

    #par analogie au fichier mip_discret_solver.jl
    sv.model = new_lp_model()
    sv.nb_calls += 1
    planes = sol.planes
    n = sv.inst.nb_planes
    sep_mat = sv.inst.sep_mat
    model = sv.model
    lbmin = lb_min(sv.inst)
    ubmax = ub_max(sv.inst)

    #
    # 1. Création du modèle spécifiquement pour cet ordre d'avion de cette solution
    #
    @variable(model,x[p in 1:n] >=0 )
    sv.x=x
    @variable(model,y[p in 1:n] >=0)
    @variable(model,z[p in 1:n] >=0)

    @expression(model,costs[p in 1:n], planes[p].ep *y[p] + planes[p].tp * z[p])

    @constraint(model,[p in 1:n],y[p]>= planes[p].target - x[planes[p].id])  # 1 
    @constraint(model,[p in 1:n],y[p]>= 0) # 2  => 1 et 2 coulent de l'expression de Max
    @constraint(model,[p in 1:n],z[p]>= x[planes[p].id] - planes[p].target) # 3
    @constraint(model,[p in 1:n],z[p]>= 0) #4 => 3 et 4 coulent de l'expression de Max
 
    @expression(model, total_cost, sum(costs[p] for p in 1:n))

    sv.costs = costs
    sv.cost = total_cost

    @objective(model,Min,total_cost)

    @constraint(model, [p in 1:n], planes[p].lb <= x[planes[p].id] )
    @constraint(model, [p in 1:n], x[planes[p].id] <= planes[p].ub)
    for p1 in 1:n-1
        for p2 in p1+1 :n
        @constraint(model, x[planes[p2].id] >= x[planes[p1].id] + sep_mat[planes[p1].kind,planes[p2 ].kind])
    end
end

    # 2. résolution du problème à permu d'avion fixée
    #
    JuMP.optimize!(sv.model)

    # 3. Test de la validité du résultat et mise à jour de la solution
    if JuMP.termination_status(sv.model) == MOI.OPTIMAL
        # tout va bien, on peut exploiter le résultat

        # 4. Extraction des valeurs des variables d'atterrissage
        #
        # ATTENTION : les tableaux x et costs sont dans l'ordre de
        # l'instance et non pas de la solution !
        for (i, p) in enumerate(sol.planes)
            sol.x[i] = round(Int, value(sv.x[p.id]))
        end
        # Mise à jour des coûts (par avion et global) de la solution à partir
        # des dates d'atterrissage. On sous-traite cette mise à jour à la méthode
        # update_costs de la "classe" Solution
        # update_costs!(sol) # diam : serait couteux car recalcule les pénalités
        update_costs!(sol, add_viol_penality=false) # diam => gain de 6% sur alp13

    else
        # La solution du solver est invalide : on utilise le placement au plus
        # tôt de façon à disposer malgré tout d'un coût pénalisé afin de pouvoir
        # continuer la recherche heuristique de solutions.
        sv.nb_infeasable += 1
        solve_to_earliest!(sol)
    end
end
