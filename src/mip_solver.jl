export MipSolver, solve
# unexport : buildMipModel

"""
    MipSolver

Résoud le problème frontale de manière exacte en PLNE.

Exploite l'hypothèse de sujet Seqata, à savoir
que la fonction de coût est linéaire en deux morceaux.

Description du modèle utilisé

### Les données

```
- P (planes) : ensemble des avions p (d'attributs p.lb, p.target, p.ub, ...)
  i est l'indice de p dans les variables de décision
- cost[p,t] : coût de pénalité de l'avion p s'il atterrit à la date t
```

### Les variables de décision

```
x[i in 1:n] : 
date d'atterrissage effective de l'avion i

costs[i in 1:n] : 
coût de chaque avion i compte tenu de sa date d'atterrissage effective

b[i,j] : boolean vrai si l'avion i atterrit avant l'avion j
```
### L'objectif

```
Minimiser total_cost = sum(costs[i]  for i in 1:n
```

### Les contraintes

    AU BOULOT !
    ...

"""
mutable struct MipSolver
    inst::Instance
    bestsol::Solution     # meilleure Solution rencontrée
    # Les attribut spécifiques au modèle
    model::Model  # Le modèle MIP
    x         # vecteur des variables d'atterrissage
    b         # vecteur des variables de précédence (before)
    cost      # variable du cout de la solution
    costs     # variables du cout de chaque avion

    # Le constructeur
    function MipSolver(inst::Instance)
        ln2("MipSolver : constructeur avec $(Args.get("external_lp_solver"))")
        this = new()
        this.inst = inst
        this.bestsol = Solution(inst)
        solver = Args.get("external_lp_solver")
        # Création et configuration du modèle selon le solveur interne sélectionné
        this.model = new_lp_model(mode=:mip, log_level=1)
        return this
    end
end

# Création du modèle frontal pour le problème Seqata (linéaire en deux morceaux)
#
function buildMipModel(sv::MipSolver)

    #error("\n\nMéthode buildMipModel(MipSolver) non implantée : AU BOULOT :-)\n\n")
    # ...

    #
    # 1. Création du modèle spécifiquement pour cet ordre d'avion de cette solution
    planes = sv.inst.planes
    n = sv.inst.nb_planes
    model = sv.model

    # Quelques raccourcis car utilisés un peu partout
    lbmin = lb_min(sv.inst)
    ubmax = ub_max(sv.inst)



    @variable(sv.model, b[1:n,1:n],Bin)
    @variable(sv.model, x[1:n] >= 0 ) #Heure d'atterrissage
    @variable(sv.model, c[1:n] >=0 ) #Coût,float
    
    sv.x=x
    sv.costs = c
    @expression(
        sv.model,
        total_cost,
        sum(c[i] for i in 1:n)
    )
    @objective(sv.model, Min, total_cost)

    sv.cost = total_cost

    @constraint(sv.model, [p in planes], c[p.id]>=p.ep*(p.target-x[p.id])) #Pour le coût
    @constraint(sv.model, [p in planes], c[p.id]>=p.tp*(x[p.id]-p.target))

    @constraint(sv.model, [p in planes], x[p.id]>=p.lb) #Pour être dans le bon intervalle
    @constraint(sv.model, [p in planes], x[p.id]<=p.ub)
    M= sum(get_sep(sv.inst,planes[i],planes[j]) for i in 1:n , j in 1:n)
    println("====> ",M)
    @constraint(sv.model,[ i in 1:n , j in 1:n , i!=j], x[planes[j].id]-x[planes[i].id]  + M* (1-b[i,j]) >= get_sep(sv.inst,planes[i],planes[j])) #Entre deux avions
    for i in 1:n
        for j in 1:n 
            if j!=i
                @constraint(sv.model, b[i,j]+ b[j,i] == 1)
                for k in 1:n
                    if k!= j && k!=i
                        @constraint(sv.model, b[i,j]+ b[j,k] - b[i,k] <= 1) #contrainte d'inegalité triangulaire sur l'ordre 
                    end
                end
            end
        end
    end
    #@constraint(sv.model,[i in 1:n, j in 1:n, k in 1:n; i!=j!=k], b[i,j]+ b[j,k] - b[i,k] <= 1) #contrainte d'inegalité triangulaire sur l'ordre 
end

# Résoud le problème complet : calcul l'ordre, le timing (dates d'atterrissage)
# et le cout total de la solution optimale du problème (coûts linéaires)
# Cette fonction crée le modèle PLNE, le résoud puis met l'objet solution à
# jour.
#
function solve!(sv::MipSolver)
    ln2("BEGIN solve!(MipSolver)")
    ln2( "="^60 )

    lg2("Construction du modèle ($(ms())) ... ")
    buildMipModel(sv)
    ln2("fait ($(ms())).")

    lg2("Lancement de la résolution ($(ms())) ... ")
    optimize!(sv.model)
    ln2("fait ($(ms())).")

    lg2("Test de validité de la solution ($(ms())) ... ")
    if JuMP.termination_status(sv.model) != MOI.OPTIMAL
        print("ERREUR : pas de solution pour :\n    ")
        @show JuMP.termination_status(sv.model)
        println(to_s(sv.bestsol))
        exit(1)
    end
    ln2("fait ($(ms())).")

    lg2("Exploitation des résultats (mise à jour de la solution)($(ms())) ... ")
    # Extraction des valeurs entières des variables d'atterrissage

    # Il reste maintenant à mettre à jour notre objet solution à partir du
    # résultat de la résolution MIP
    for (i, p) in enumerate(sv.inst.planes)
        sv.bestsol.planes[i] = p
        sv.bestsol.x[i] = round(Int, value(sv.x[p.id]))  # car value() retourne Float64 !
        sv.bestsol.costs[i] = value(sv.costs[p.id])
    end
    sv.bestsol.cost = round(value(sv.cost), digits=Args.get(:cost_precision))

    # On trie juste la solution par date d'atterrissage croissante des avions
    # pour améliorer la présentation de la solution
    sort!(sv.bestsol)
    ln2("fait. ($(ms()))")

    ln2("END solve!(MipSolver)")
end
