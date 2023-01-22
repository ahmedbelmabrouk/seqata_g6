export Solution, current_algo_symbol
export get_names, set_from_names!, sort!
export to_s, to_s_long
export get_viol_penality, get_viol_description, is_feasable
export solve!, disturb!, update_costs!
export solve_to_earliest!
export swap!, shift!, permu!, initial_sort!
export guess_solname, print_sol
# unexport Random.shuffle!
# unexport Base.write
# unexport Base.show
# unexport Base.copy!
# unexport init_solver, hash

"""
    Solution

Représente une solution au problème Seqata avec ses attributs et méthodes.

Les principaux attributs sont les suivants :
- `inst::Instance` : l'instance associée
- `planes::Vector{Plane}` : lliste de avions dans l'ordre de la solution
- x et costs : vecteurs des variables de décisions associées à chaque avion
- cost : le coût global
- timing_algo_solver : le solveur utilisé par la solution pour résoudre le
  Sous-Problème de Timing (STP) à ordre des avions fixé).
  Ce solveur paut être choisi via l'option `--timing_algo_solver` de la ligne
  de commande et est utilisé par la méthode `solve!`

et les méthodes pour la manipuler dont :
- solve!, to_s, to_s_long, get_viol_penality, ...
"""
mutable struct Solution
    inst::Instance
    planes::Vector{Plane}
    x::Vector{Int}
    costs::Vector{Float64}
    cost::Float64
    timing_algo_solver     # :earliest, :lp, :dp, ...
    # lpTimingSolver       # solver lp si c'est :lp qui est choisi
    # earliestTimingSolver # solver si c'est :earliest qui est choisi

    # l'attribut contiendra l'objet solver dont le type dépend du symbole
    # contenu dans l'attribut timing_algo_solver précédent
    solver # Exemple : EarliestTimingSolver (type fonction du symbol timing_algo_solver)

    # Le constructeur peut prendre :
    # - soit une instance qui sert à initialiser la solution
    # - soit un autre objet solution pour en faire une copie
    #
    function Solution(inst::Instance; update::Bool = true, algo = :ARGS)
        this = new()
        this.inst = inst
        this.planes = copy(inst.planes)
        # @show inst.nb_planes
        this.x = zeros(Int, inst.nb_planes)
        this.costs = zeros(inst.nb_planes)
        this.cost = typemax(Float64) # valeur Inf pour un Float64)

        # @show algo
        # @show Args.t_algos()
        if algo == :ARGS
            this.timing_algo_solver = Args.get("timing_algo_solver")
        else
            if String(algo) in Args.t_algos()
                this.timing_algo_solver = algo
            else
                error("timing_algo_solver \"$(algo)\" inconnu !")
            end
        end
        this.solver = nothing    # sera initialisé lors du premier solve!()

        ln3("Solution: Construction d'une solution avec $(this.timing_algo_solver)")
        if update
            init_solver(this, this.timing_algo_solver)
            solve!(this)
        end
        return this
    end
    # ATTENTION la référence au timingSolver est copiées superficiellement !
    function Solution(sol::Solution)
        this = new()
        # On recopie une solution existante
        this.inst = sol.inst
        this.planes = copy(sol.planes)
        this.x = copy(sol.x)
        this.costs = copy(sol.costs)
        this.cost = sol.cost
        this.timing_algo_solver = sol.timing_algo_solver
        this.solver = sol.solver
        # solve!(this)  # ceci est en principe inutile
        return this
    end
end

# Ce constructeur lit une solution préenregistrée dans un fichier
function Solution(inst::Instance, filename::String; update::Bool = true, algo = :ARGS)
    this = Solution(inst, update = update, algo = :ARGS)
    readsol(this, filename)
end

# Recopie dans l'objet sol le contenu de l'objet other
function Base.copy!(sol::Solution, other::Solution)
    sol.inst = other.inst
    copyto!(sol.planes, other.planes)
    copyto!(sol.x, other.x)
    copyto!(sol.costs, other.costs)
    sol.cost = other.cost  # ceci maintient la cohérence du coût avec la solution
    sol.timing_algo_solver = other.timing_algo_solver
    sol.solver = other.solver # ATTENTION : COPIE SUPERFICIELLE
    return sol   # pratique pour enchainer les appels de méthode (en POO !)
end

# Crée et affecte l'attribut solver pour la résolution des dates d'atterrissage.
# Le type du solver créée est mémorisé dans l'attribut symbole
# sol.timing_algo_solver (par exemple :ealiest, :lp, :dp, ...)
#
function init_solver(sol::Solution, algo::Symbol)
    if algo == :earliest
        # sol.solver = nothing # avant création d'un TimingSolver indépendant
        sol.solver = EarliestTimingSolver(sol.inst)
    elseif algo == :lp
        sol.solver = LpTimingSolver(sol.inst)
    elseif algo == :dp
        sol.solver = DynprogTimingSolver(sol.inst)
        # error("init_solver:  DynprogTimingSolver non implanté pour : $(algo)")
    elseif algo == :faye
        sol.solver = FayeTimingSolver(sol.inst)
        # error("init_solver:  FayeTimingSolver non implanté pour : $(algo)")
    else
        error("init_solver:  algo invalide : $(algo)")
    end
    sol.timing_algo_solver = algo
end

# Return le symbol de l'algo à partir du TimingSolver mémorisé dans la
# variable sol.solver
# Si l'objet n'est pas un solver, le symbol :earliest est retourné
#
function current_algo_symbol(sol::Solution)
    if sol.solver != nothing
        return symbol(sol.solver)
    else
        return :earliest
    end
end

# Calcul du hachage d'une solution.
# Deux solutions identiques doivent donner le même hachage.
# Deux solutions sont considérées comme identique ssi elles ont les avions
# dans le même ordre, et avec les mêmes dates d'atterrissage.
# HYP: les deux solutions sont à jour (par solve!)
#
function hash(sol::Solution)
    # on construit une chaine de caractères représentative de la solution
    # (ce qui revient plus ou moins à sérialiser l'objet sol)
    io = IOBuffer()
    for p in sol.planes
        print(io, p.name, ":", round(Int, sol.x[p.id]), "-")
    end
    return Base.hash(String(take!(io)))
end

# Retourne le tableaux des noms d'avions dans l'ordre de la solution.
function get_names(sol::Solution)
    # Retourne le tableaux des noms d'avions dans l'ordre de la solution
    # RUBY: @planes.map{| plane |  plane.name}
    # Autre méthode par les compréhensions : [p.name for p in sol.planes]
    map(p::Plane -> p.name, sol.planes)
end

# Affecte la solution en fonction de l'ordre des noms d'avion définit
# en paramètre
# Exemple
#    set_from_names!(mysol, [3, 4, 5, 6, 8, 9, 7, 1, 10, 2])
#
# ATTENTION : couteux O(n^2), mais pratique pour affichage et tests.
#
function set_from_names!(sol::Solution, names)
    if length(names) != length(sol.planes)
        error(
            "\nERREUR : taille de names $(length(names)) incorrecte " *
            "attendue : $(length(sol.planes))",
        )
    end
    # TODO : tester l'opération dot-parenthesis (16/07/2018)
    # sol.planes .= get_plane_from_name.(sol.inst, names)
    for i = 1:length(names)
        sol.planes[i] = get_plane_from_name(sol.inst, names[i])
    end
    return sol
end

# Trie les avions de l'objet sol selon leur date d'atterrissage croissante.
# La solution n'est pas recalculée, mais les vecteurs x et planes sont
# triés de manière cohérente.
# .
function sort!(sol::Solution)
    # Pour cela, on récupère l'ordre correspondant au tri des dates d'atterrissage
    # (le résultat est la permutation correspondant au tri à effectuer)
    permu = sortperm(sol.x)

    # On applique cette permutation aux deux tableaux de notre solution pour les
    # trier dans le même ordre
    sol.x = sol.x[permu]
    sol.planes = sol.planes[permu]
    sol.costs = sol.costs[permu]
end

function Base.show(io::IO, sol::Solution)
    Base.write(io, to_s(sol))
end

# Convertit une Solution en chaine
# e.g  "cost=123.456  :[p2,p5,p1,...,p10]"
#
function to_s(sol::Solution)
    # prc = Args.get(:cost_precision) + 2 # cmt le 30/10/2018 car 7 par défaut
    prc = Args.get(:cost_precision)
    cost = round(sol.cost, digits = prc)

    txt =
        "cost=" *
        rpad(cost, 8) *
        " :" *
        # "[" * sol.planes.map {|plane| plane.name }.join(",") << "]"
        "[" *
        join(map(p -> p.name, sol.planes), ",") *
        "]"
    return txt
end

function to_s_long(sol::Solution)
    io = IOBuffer()
    prc = Args.get(:cost_precision)
    cost = round(sol.cost, digits = prc)
    viol_cost = get_viol_penality(sol)

    println(io, "# ALP solution version 1.0")
    println(io, "name $(sol.inst.name)  # nb_planes=$(length(sol.planes))")
    if viol_cost != 0.0
        println(io, "# Solution **invalide** ")
        println(io, "viols_cost=$(viol_cost)")
    else
        println(io, "# Solution valide")
    end
    println(io, "timestamp $(Dates.now())")
    println(io, "cost $(cost)")
    println(io, "order [", join([p.name for p in sol.planes], ","), "]")
    println(io)
    println(io, "#       name     t   dt  cost        # comments")
    for i = 1:length(sol.planes)
        p = sol.planes[i]
        print(io, @sprintf("landing %4s ", p.name))
        print(io, @sprintf(" %4s ", sol.x[i]))
        print(io, @sprintf("%4s  ", sol.x[i] - p.target))

        # Ceci est horrible, mais la chaine de formatage de @sprintf ne peut
        # pas contenir de variable !
        #   print(io, @sprintf("% 7.$(prc)f ", sol.costs[i]) ) # CECI PLANTE
        if prc <= 5
            print(io, @sprintf("%9.5f ", sol.costs[i]))
        else
            print(io, @sprintf("%11.7f ", sol.costs[i]))
        end

        print(io, "  # ")
        print(io, @sprintf("E=%-4d ", p.lb))
        print(io, @sprintf("T=%-4d ", p.target))
        print(io, @sprintf("L=%-4d ", p.ub))
        print(io, "  ")
        print(io, @sprintf(" ep=%3.2f  ", p.ep))
        print(io, @sprintf(" tp=%3.2f  ", p.tp))

        # On affiche les temps d'écart et temps attendu entre les trois derniers
        # prédessesseurs pour vérifier le respect des temps d'écart.
        print(io, " sep ")
        isbad = false
        # sep_string = ""
        for j = i-1:-1:max(i - 3, 1)
            p_prev = sol.planes[j]
            sep_real = sol.x[i] - sol.x[j]
            sep_th = get_sep(sol.inst, p_prev, p)
            # sep_string *= " $sep_real($sep_th)"
            if sep_real < sep_th
                isbad = true
            end
            print(io, " $sep_real($sep_th)")
        end
        if isbad
            # print(io, to_sc(" VIOL "), sep_string)
            println(io, to_sc(" VIOL "))
        else
            # print(io, " -ok- ", sep_string)
            println(io, " -ok- ")
        end
    end
    println(io)
    String(take!(io))
end

# Vérifie la faisabilité de la solution par rapport aux dates d'atterrissage
# réelle de la solution.
# Retourne un coût de pénalité (nul si solution valide) si des contraintes sont
# violées. Ce coût prend en compte :
# - l'amplitude des viols de la contrainte de séparation,
# - l'amplitude des viols des bornes lb et ub de chaque avion.
#
# Cette méthode suppose que l'information sol.x de la date d'atterrissage
# réelle de chaque avion est à jour.
#
# Les viols sont détectés même si l'inégalité triangulaire n'est pas respectée
# dans l'instance (instance incorrecte par hypothèse)
#
# Si level est suffisant : affiche message si solution non faisable
#
function get_viol_penality(sol::Solution)
    total_penality = 0.0
    unit_penality = 1000.0
    SEP_MAX = sep_max(sol.inst)

    ### ASSERT : p1=sol.planes[1] est bien positionné sinon il

    # On commence par prendre en compte la pénalité de placement du premier avion
    # On calcule l'écart (viol sur les bornes) de la date d'atterissage x[1]
    # cet écart est positif sauf en cas de viol.
    p1 = sol.planes[1]
    dt = min(sol.x[1] - p1.lb, p1.ub - sol.x[1])
    if dt < 0
        total_penality += -dt * unit_penality
    end

    for i2 = 2:length(sol.planes)
        p2 = sol.planes[i2]

        # Étape 1 : calcul du viol des temps de séparation pour p2

        # parcours les prédecesseurs en partant de de p2
        for i1 = i2-1:-1:1
            if sol.x[i2] > sol.x[i1] + SEP_MAX
                # tous les temps d'écarts des prédécesseurs sont respectés
                break
            end
            p1 = sol.planes[i1]
            # s_val = sol.inst.sep_mat[p1.kind, p2.kind]
            s_val = get_sep(sol.inst, p1, p2)
            if sol.x[i2] < sol.x[i1] + s_val
                total_penality += unit_penality * (sol.x[i1] + s_val - sol.x[i2])
                # il faudrait pousser j
                feasible = false
                if lg4() >= 4
                    println("ERREUR écart entre type des avions $(p1.name)->$(p2.name)")
                    println("   x1=$(sol.x[i1]) pour $(p1)")
                    println("   x2=$(sol.x[i2]) pour $(p2)")
                    println("   écart réel=$(sol.x[i2]-sol.x[i1]) au lieu de $(s_val)")
                end
            end
        end

        # Étape 2 : calcul du viol de la Upper Bound
        if sol.x[i2] > p2.ub
            total_penality += unit_penality * (sol.x[i2] - p2.ub)
            lg4() && println("UB infaisable pour l'avion $(p2.name)")
        end

        # Étape 3 : calcul du viol de la Lower Bound
        if sol.x[i2] < p2.lb
            total_penality += unit_penality * (p2.lb - sol.x[i2])
            lg4() && println("LB infaisable pour l'avion $(p2.name)")
        end
    end
    if lg4() >= 4 && total_penality == 0.0
        println("solution faisable.")
    end
    feasible = (total_penality == 0.0)
    if !feasible && lg4()
        println("get_viol_penality: Infaisable : total_penality=", total_penality)
    end
    return total_penality
end

# Vérifie la faisabilité dates d'atterrissage x de la solution par rapport
# aux contraintes de l'instance.
# Si pas de viol, retourne une chaine vide
# Si viol, retourne un texte avec une ligne par contrainte violée.
# Les contraintes testés sont les suivante :
# - la date d'atterissage doit entre dans l'intervale permi pour l'avion
# - tous les temps de séparation doivent être respertés
#
# Cette méthode suppose que l'information sol.x de la date d'atterrissage
# réelle de chaque avion est à jour.
#
# Les viols sont détectés même si l'inégalité triangulaire n'est pas respectée
# dans l'instance (instance incorrecte par hypothèse)
#
function get_viol_description(sol::Solution)
    io = IOBuffer()

    nbviols = 0
    for i = 1:length(sol.planes)
        p = sol.planes[i]
        x = sol.x[i]
        lb = p.lb
        ub = p.ub
        if sol.x[i] - p.lb < 0
            println(io, p.name, " atterrit trop tôt  x=$x < lb=$lb")
            nbviols += 1
        end
        if p.ub - sol.x[i] < 0
            println(io, p.name, " atterrit trop tard ub=$ub < x=$x")
            nbviols += 1
        end
    end

    SEP_MAX = sep_max(sol.inst)
    for i2 = 2:length(sol.planes)
        p2 = sol.planes[i2]
        x2 = sol.x[i2]
        # parcours les prédecesseurs en partant de de p2
        for i1 = i2-1:-1:1
            p1 = sol.planes[i1]
            x1 = sol.x[i1]
            #
            # On pourrait ne pas explorer tout si l'on se contentait de détecter
            # une invalidité. Mais on veut afficher **tous** les viols
            #   if x2 > x1 + SEP_MAX
            #       # tous les temps d'écarts des prédécesseurs sont respectés
            #       # (s'il n'y a aucun e viol avant x1)
            #       break
            #   end
            #
            sep = get_sep(sol.inst, p1, p2)
            if x2 < x1 + sep
                print(io, "$(p1.name)->$(p2.name) : écart insuffisant ")
                print(io, "x[$(p1.name)]=$x1->x[$(p2.name)]=$x2 => ")
                println(
                    io,
                    "sep($(p1.name),$(p2.name))=$(sol.x[i2]-sol.x[i1]) ",
                    "au lieu de ",
                    get_sep(sol.inst, p1, p2),
                )
                nbviols += 1
            end
        end
    end
    violtxt = String(take!(io))
    return (nbviols, violtxt)
end

# Calcul la faisabilité intrinsèque de la solution en fonction de l'ordre
# imposé pour les avions.
#
# Ne modifie pas la solution
#
# N'utilise pas l'attribut x[i] de la solution associée au placement
# de chaque avion.
#
# Principe :
# Si on peut placer les avions au plus tôt en respectant les temps de
# séparation sans violer la borne sup du placement des avions suivants,
# alors la solution est faisable.
#
function is_feasable(sol::Solution; param_bidon = nothing)
    # ASSUME_TRINEQ = Args.get(:assume_trineq)  # Projet Alap
    ASSUME_TRINEQ = false # L'inégalité triangulaire n'est pas présupposée

    n = length(sol.planes)
    # x : vecteur des placements au plus tôt des avions
    x = Vector{Int}(undef, n)

    i1 = 1
    p1 = sol.planes[i1]
    x[i1] = p1.lb
    # @assert (x[i1] <= p1.ub), "Instance infaisable pour avion $(p1) !"
    @assert x[i1] <= p1.ub

    for i2 = 2:length(sol.planes)
        p2 = sol.planes[i2]
        x[i2] = max(p2.lb, x[i1] + get_sep(sol.inst, p1, p2))
        if !ASSUME_TRINEQ
            # inégalité triangulaire non vérifiée : on doit vérifier
            # les prédécesseurs antérieurs à p1
            for i_prev = i1-1:-1:1
                p_prev - sol.planes[i_prev]
                this_sep = x[i2] - x[i_prev]
                if this_sep >= MAX_SEP
                    # c'est bon : p2 est forcément suffisamment loin devant
                    # tous ces prédécesseurs
                    break
                end
                # on retarde p2 pour satisfaire cette séparation
                x[i2] = max(x[i2], x[i_prev] + get_sep(sol.inst, p_prev, p2))
            end
        end
        # maintenant que p2 est placé au plus tôt en respectant les contraintes
        # de séparation, on vérifie qu'il n'est pas hors borne supérieure
        if x[i2] > p2.ub
            # if lg1() println("viol x[$(i2)]=$(x[i2]) mais $(p2).ub=$(p2.ub)") end
            return false
        end
    end
    # Tous les avions ont pu être placés au plus tôt.
    # La solution est donc timing-faisable
    return true
end

"""
    solve!(sol::Solution; do_update_cost::Bool = false)

Résoud le sous-problème de timind de la solution à partir de l'ordre des
avions défini par l'attribut `planes`.

Le choix de l'algorithme utilisé est mémorisé dans l'attribut `timing_algo_solver`
et peut être modifié via la ligne d ecommande.

Le "!" dans le nom `solve!`` est une convention indiquant que l'objet
est modifié.

`do_update_cost` (false par défaut) : si vrai, met à jours les coûts.

REMARQUE :

- Cette méthode peut être appelée par un solveur global qui modifie
  la solution courante. Cependant cet appel est souvent effectué par une méthode
  interne à la classe Solution (i.e. dans ce fichier, e.g xxx)
- l'option `do_update_cost`` n'est utile que si les coûts individuels des avions
  ne sont pas déjà faits par le TimingSolver utilisé.
  Mais il est souhaitable que le TimingSolver se charge lui-même de mettre
  à jour les coûts individuels des avions.
- Si aucun TimingSolver n'a été spécifié lors de la création de la solution
  alors on utilise la résolution au plus tôt (sous-optimale).

TODO: prévoir paramètre from_index=1 pour permettre un recalcul partiel.
Par exemple cela permettra de recalculer plus efficacement des dates d'atterrissage
d'une solution du problème avec 500 avions si l'on a seulement permuté les deux
derniers avions !

"""
function solve!(sol::Solution; do_update_cost::Bool = false)
    if sol.solver == nothing
        init_solver(sol, sol.timing_algo_solver)
    end
    solve!(sol.solver, sol)
    if do_update_cost
        update_costs!(sol)
    end
    return sol
end

# Mélange la solution et mets à jour (par défaut) les timings
#
# do_update (true par défaut) : si vrai , résoud le pb de timing
#     REMARQUE : l'option do_update ne semble plus utilisée
function Random.shuffle!(sol::Solution; do_update::Bool = true)
    Random.shuffle!(sol.planes)
    if do_update
        solve!(sol) # NEW 03/04/2019
    end
    return sol
end

# Remue légèrement la solution et mets à jour (par défaut) les timings
function disturb!(
    sol::Solution;
    idx_first = 1,
    idx_last = -1,
    shift_max = 1,
    nb_shift = 10,
    do_update::Bool = true,
)

    # VÉRIFICATION DU DOMAINE DES PARAMÈTRES
    #
    if idx_last == -1
        idx_last = length(sol.planes)
    end
    @assert 1 < idx_last <= length(sol.planes)
    @assert 1 <= idx_first < idx_last

    # shift_max == -1 => on prend le maximum possible
    if shift_max == -1
        shift_max = idx_last - idx_first
    end
    @assert shift_max >= 1 "shift_max doit être stricement positif (ou -1 pour max)"

    # APPLICATION DES MOUVEMENTS (de type shift)
    #
    for i = 1:nb_shift
        dist = rand(1:shift_max)
        idx1 = rand(idx_first:idx_last)
        idx2_min = max(idx1 - dist, idx_first)
        idx2_max = min(idx1 + dist, idx_last)
        idx2 = rand(idx2_min:idx2_max)
        shift!(sol, idx1, idx2, do_update = false)
    end

    # RECALCUL ÉVENTUELLE DES COÛTS DE LA SOLUTION
    #
    if do_update
        # do_update_cost est utile ici par la solution obtenue peut-être infaisable
        # et il faut pouvoir évaluer le coût de cette solution
        solve!(sol, do_update_cost = true) # mise à jour après le dernier shift
    end
end

"""
    update_costs!(sol::Solution; add_viol_penality = true)

Met à jour les coûts de chaque avion (sol.costs) et le coût global
(sol.cost) à partir des dates d'atterrissage supposées connues (sol.x).

Si `add_viol_penality` est true, les contraintes sont vérifiées et un coût
supplémentaire est ajouté.

PRECONDITION : le Vector s.costs est dans le même ordre que s.planes
"""
function update_costs!(sol::Solution; add_viol_penality = true)


    # Mise à jour des coûts individuels de chaque avion
    for i = 1:length(sol.planes)
        p = sol.planes[i]
        sol.costs[i] = get_cost(p, sol.x[i])
    end

    # Mise à jour du coût global de la solution
    #
    sol.cost = sum(sol.costs)
    # @show Args.get(:cost_precision)
    # @show sol.cost
    sol.cost = round(sol.cost, digits = Args.get(:cost_precision)) # 30/10/2018
    # @show sol.cost

    if lg4() >= 4
        println("get_viol_penality=$(get_viol_penality(sol))")
    end
    if add_viol_penality
        sol.cost += get_viol_penality(sol)
    end
    return sol
end

# Met à jour la date d'atterrissage de chaque avion au plus tôt, de façon à
# respecter les contraintes de précédence au mieux.
# Principe :
#   Les contraintes de séparation sont **toujours respectées**.
#   Au besoin on viole la borne d'atterrissage des avions successeurs (p.ub)
#
function solve_to_earliest!(sol::Solution; do_update_cost = true)
    # println("="^70)
    # println("=== BEGIN solve_to_earliest! pour $(to_s(s))")
    # println("")
    SEP_MAX = sep_max(sol.inst)
    # Méthode : on parcourt tous les successeurs (plane2) et on s'assure qu'ils
    # sont suffisamment éloignés de **tous** leurs prédécesseurs (les plane1).
    for i2 = 1:length(sol.planes)
        plane2 = sol.planes[i2]
        # @show plane2
        sol.x[i2] = plane2.lb # date au plus tot
        # parcours les prédecesseurs de plane2 du plus proche au plus éloigné
        for i1 = i2-1:-1:1
            # @show sol.planes[i1]
            if (sol.x[i2] - sol.x[i1]) >= SEP_MAX
                # i2 est suffisamment en avance quels que soit le type de
                # l'avions i1 qui le précède
                break
            end
            plane1 = sol.planes[i1]
            if sol.x[i2] < sol.x[i1] + sol.inst.sep_mat[plane1.kind, plane2.kind]
                # il faut pousser plane2
                sol.x[i2] = sol.x[i1] + sol.inst.sep_mat[plane1.kind, plane2.kind]
            end
        end
        if sol.x[i2] > plane2.ub
            if lg4() >= 4
                println("pb infaisable pour l'avion $(plane2.name)")
            end
            # TODO: On pourrait ajouter un coût de pénalité pour les
            # contraintes de bornes non respectées
        end
    end
    if do_update_cost
        update_costs!(sol)
    end
    return sol
end

# swap! : Permute deux avions dans la solution
# idx1, idx2 : les indices des avions à permuter
#   si un indice n'est pas défini : il est choisi aléatoiremenet.
#   retourne les indices des avions permutées
# do_update : recalcule le coût (true par défaut)
#
function swap!(sol::Solution, idx1 = -1, idx2 = -1; do_update = true)
    # RUBY: idx1 = rand(@planes.size)  if not idx1
    if idx1 == -1
        idx1 = rand(1:length(sol.planes))
    end
    if idx2 == -1
        idx2 = rand(1:length(sol.planes))
    end
    # Il se peut que idx1 == idx2 mais tant pis
    sol.planes[idx1], sol.planes[idx2] = sol.planes[idx2], sol.planes[idx1]
    if idx1 != idx2 && do_update
        solve!(sol)
    end
    return (idx1, idx2)
end

# shift! : Déplace avions dans la solution
# idx1, idx2 : les indices origine et cible de l'avion à déplacer
#   si un indice n'est pas défini : il est choisi aléatoiremenet.
#   retourne les indices du déplacement effectué.
# do_update : recalcule le coût (true par défaut).
#
function shift!(sol::Solution, idx1 = -1, idx2 = -1; do_update = true)
    if idx1 == -1
        idx1 = rand(1:length(sol.planes))
    end
    if idx2 == -1
        idx2 = rand(1:length(sol.planes))
    end
    # Il se peut que idx1 == idx2 mais tant pis
    if idx1 == idx2
        do_update && solve!(sol)
        return (idx1, idx1)
    elseif abs(idx1 - idx2) == 1
        return swap!(sol, idx1, idx2, do_update = do_update)
    end
    # println("AVANT SHIFT $idx1->$idx2 : ", to_s(sol))
    # 03/05/2019 : shift! est définie dans le fichier array_util.jl
    shift!(sol.planes, idx1, idx2)

    if do_update
        solve!(sol)
    end
    # println("APRES SHIFT $idx1->$idx2 : ", to_s(sol))
    return (idx1, idx2)
end

# permu! : Permute des avions dans la solution.
# indices1 : tableau des indices d'avions à permuter
# indices2 : tableau des nouveaux indices des avions permutés
# HYPOTHESE : indices2 est une permutation valide de indices1
#
# Les avions de la solution ont pour indice 0:nb_planes-1
# Par exemple, si indices1 vaut [0, 2, 7, 12]
# et si indices2 vaut [12, 7, 0, 2],
# L'opération consiste à permuter les avions positionnés en :
#    0, 2, 7, 12
# de façon à ce que leur nouvelle position soit en :
#   12, 7, 0, 2
#
# Exemple
# On veut transformer la liste d'avion :
#    @plane = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
#               0   1*  2   3   4   5   6*  7*  8    9
# en :
#    @plane = [10, 70, 30, 40, 50, 60, 80, 20, 90, 100]
#               0   1*  2   3   4   5   6*  7*  8    9
#
# Pour cela l'appel doit être d ela forme :
#  mysol.permu! [1, 6, 7], [7, 1, 6]
#
# do_update : racalcule le coût (true par défaut)
#
# HYPOTHESE :
# TODO: effectuer une mise à jour différencielle du coût
#
function permu!(sol::Solution, indices1, indices2; do_update = true)
    sol.planes[indices1] = sol.planes[indices2]
    # sol.planes[[indices1...]] = sol.planes[[indices2...]] # avec des tuples iiii
    if do_update
        solve!(sol)
    end
    return indices2
end

# Effectue un tri de la solution courante selon le critère passer en paramètre.
# Si aucun critère presort n'est imposé en argumant, celui spécifié par
# l'option --presort est utilisé
# Le timing de la solution estmis à jour.
#
function initial_sort!(sol::Solution; presort = :ARGS)

    # On part de l'instance en triant les avions par ordre de lb croissant
    if presort == :ARGS
        presort = Args.get("presort")
    end
    # @show presort
    if presort == :inst
        # le champs id correcpond au numéro d'ordre d el'avon dans l'instance
        Base.sort!(sol.planes, by = p -> p.id)
    elseif presort == :rinst
        Base.sort!(sol.planes, by = p -> p.id, rev = true)
    elseif presort == :target
        Base.sort!(sol.planes, by = p -> p.target)
    elseif presort == :rtarget
        Base.sort!(sol.planes, by = p -> p.target, rev = true)
    elseif presort == :lb
        Base.sort!(sol.planes, by = p -> p.lb)
    elseif presort == :rlb
        Base.sort!(sol.planes, by = p -> p.lb, rev = true)
    elseif presort == :ub
        Base.sort!(sol.planes, by = p -> p.ub)
    elseif presort == :rub
        Base.sort!(sol.planes, by = p -> p.ub, rev = true)
    elseif presort == :shuffle
        # Base.shuffle!(sol.planes)
        Random.shuffle!(sol.planes)
    elseif presort == :none
        # rien à faire
    else
        error("Valeur incorrecte pour presort : ", presort)
    end

    # Autre solution pour trier de type d'avions (A GARDER COMME EXEMPLES) :
    # On passe la fonction qui compare deux éléments :
    # Base.sort!(sol.planes, lt=(p1,p2)->p1.lb < p2.lb) # Marche aussi (comparator)
    # Base.sort!(sol.planes, rev=true,lt=(p1,p2)->p1.lb < p2.lb) # ok aussi THE WORST CASE

    solve!(sol, do_update_cost = true)
    return nothing
end

# Enregistre la solution dans un fichier.
# Par défaut, le nom est construit automatiquement à partir des caractéristiques
# de l'instance et du coût de la solution (arrondi à deux décimales)
#
function Base.write(sol::Solution, filename = "")
    dir = Args.get(:outdir) # e.g "_tmp"
    if !isdir(dir)
        dir = "."
    end
    cost = round(sol.cost, digits = 2) # cost intervient dans le nom du fichier
    if filename == ""
        filename = guess_solname(sol)
    end
    # solpath = dir*"/"*filename
    solpath = "$dir/$filename"
    open(solpath, "w") do fh
        # Base.write(fh, to_s(sol))
        # Base.write(fh, "\n")
        Base.write(fh, to_s_long(sol))
    end
end

# Retourne un nom de la forme :
#   alp_01_p10=700.0.sol              # VERSION SIMPLIFIÉE POUR SEQATA
#   alp_13_p500_k201=6440.8959533.sol # VERSION ALAP pour 100 segments
#
function guess_solname(sol::Solution)
    prc = Args.get(:cost_precision)
    cost = round(sol.cost, digits = prc)
    return "$(sol.inst.name)=$(cost).sol"
end

# Fonction d'impression de haut niveau
function print_sol(sol::Solution, msg::AbstractString = "")
    if msg != ""
        println(msg)
    end
    println("\nmeilleure solution trouvée :")
    println("="^70)
    println(to_s(sol))
    println("="^70)
    println(to_s_long(sol))
    println("="^70)
end

