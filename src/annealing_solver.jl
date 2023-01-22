export AnnealingSolver
export finished, get_stats, guess_temp_init
export solve

mutable struct AnnealingSolver
    inst::Instance
    opts::Dict

    temp_init::Float64 # température courante
    temp_mini::Float64 # température mini avant arrêt
    temp_coef::Float64 # coefficiant de refroidissement
    temp::Float64      # température courante

    nb_test::Int     # Nombre total de voisins testés
    nb_move::Int     # Nombre de voisins acceptés (améliorant ou non)
    nb_reject::Int   # Nombre de voisins refusés
    nb_steps::Int    # Nombre de paliers à température constante effectués
    step_size::Int   # Nombre d'itérations à température constante
    iter_in_step::Int # Nombre d'itér effectués dans le palier courant

    # Nombre maxi de refus consécutifs
    nb_cons_reject        # Nombre de refus consécutifs
    nb_cons_reject_max    # Nombre maxi de refus consécutifs

    # Nombre tests non améliorants consécutifs
    nb_cons_no_improv     # Nombre de tests non améliorants
    nb_cons_no_improv_max # Nombre de maxi tests non améliorants

    cursol::Solution        # Solution courante
    bestsol::Solution       # meilleure Solution rencontrée
    testsol::Solution       # nouvelle solution courante potentielle

    function AnnealingSolver(inst::Instance, user_opts = Dict())
        # ...
        error("\n\nConstructeur de AnnealingSolver non implanté : AU BOULOT :-)\n\n")
        # ...
        this = new()
        this.inst = inst

        this.opts = Dict(
            :startsol => nothing,    # nothing pour auto à partir de l'instance
            :step_size => 1,
            :temp_init => -1.0, # -1.0 pour automatique
            :temp_init_rate => 0.75,  # valeur standard : 0.8
            :temp_mini => 0.000_001,
            :temp_coef => 0.999_95,
            :nb_cons_reject_max => 1_000_000_000, # infini
            :nb_cons_no_improv_max => 5000 * inst.nb_planes,
        )
        nb_options = length(this.opts)
        merge!(this.opts, user_opts)

        if length(this.opts) != nb_options
            error("Au moins une option du recuit est inconnue dans :\n$(keys(user_opts))")
        end

        if user_opts[:startsol] == nothing
            this.cursol = Solution(inst)
        else
            this.cursol = user_opts[:startsol]
        end

        # On calcule éventuellement la température initiale automatiquement
        if this.opts[:temp_init] == nothing
            this.opts[:temp_init] =
                guess_temp_init(this.cursol, this.opts[:temp_init_rate], 100)
        end

        # À POURSUIVRE (AU BOULOT)
        # ...

        this.bestsol = Solution(this.cursol)
        this.testsol = Solution(this.cursol)

        return this
    end
end
# stop : retourne true ssi l'état justifie l'arrêt de l'algorithme
# (dommage qu'on ne puisse pas l'appeler stop? comme en ruby !)
# On pourra utiliser d'autres critères sans toucher au programme principal
#
function finished(sv::AnnealingSolver)
    # AU BOULOT !
    return false
end

function get_stats(sv::AnnealingSolver)
    txt = "
    Paramètres de l'objet AnnealingSolver :
    step_size=         $(sv.step_size)
    temp_init=         $(sv.temp_init)
    temp_init_rate=    $(sv.opts[:temp_init_rate])
    temp_mini=         $(sv.temp_mini)
    temp_coef=         $(sv.temp_coef)
    nb_cons_reject_max=$(sv.nb_cons_reject_max)
    Etat de l'objet AnnealingSolver :
    nb_steps=$(sv.nb_steps) step_size=$(sv.step_size)
    nb_cons_reject=$(sv.nb_cons_reject) nb_cons_reject_max=$(sv.nb_cons_reject_max)
    nb_cons_no_improv=$(sv.nb_cons_no_improv) nb_cons_no_improv_max=$(sv.nb_cons_no_improv_max)
    nb_test=$(sv.nb_test)
    nb_move=$(sv.nb_move)
    nb_reject=$(sv.nb_reject)
    temp=$(sv.temp) temp_init=$(sv.opts[:temp_init])
    testsol.cost=$(sv.testsol.cost)
    cursol.cost=$(sv.cursol.cost)
    bestsol.cost=$(sv.bestsol.cost)
    sv.testsol.solver.nb_infeasable=$(sv.testsol.solver.nb_infeasable)
    "
    return replace(txt, r"^ {4}" => "")
end

# Calcul d'une température initiale de manière à avoir un
# taux d'acceptation TAUX en démarrage
#
# arguments :
#   - taux_cible : pourcentage représentant le taux d'acceptation cible(e.g. 0.8)
#   - nb_degrad_max : nbre de degradation à accepter pour le calcul de la moyenne
#
# Principe :
#   On lance une suite de mutations (succession de mouvement systématiquement
#   acceptés). On relève le nombre et la moyenne des mouvements conduisant à une
#   dégradation du coût de la solution.
#
#   degrad : dégradation moyenne du coût pour deux mutations consécutives de coût
#       croissant
#
#   La probabilité standard d'acceptation d'une mauvaise solution est :
#       p = e^{ -degrad/T } = 0.8    =>    T = t_init = -degrad / ln(p)
#
#   avec :
#       p = taux_cible = proba(t_init)
#       degrad = moyenne des dégradations de l'énergie
#       T = t_init = la température initiale à calculer
#
# Exemple :
#   On va lancer des mutations jusqu'à avoir 1000 dégradations.
#   Si par exemple le coût des voisins forme une suite de la forme :
#
#       990, 1010, 990, 1010, 990,...
#
#   On devra faire 2000 mutations pour obtenir 1000 dégradations de valeur 20,
#   d'où t_init = -degrad / ln(proba)
#       proba = 0.8   =>  t_init = degrad * 4.5
#       proba = 0.37  =>  t_init = degrad
#
# ATTENTION :
#  Cette fonction n'est **pas** une méthode de AnnealingSolver.
#  Elle a juste besoin d'une solution et du type de mouvement à effectuer.
#  Ici, on suppose que le seul mouvement possible est swap!(sol::Solution)
#  Mais il faudra pouvoir paramétrer cette méthode pour des voisinages différents.
#
function guess_temp_init(sol::Solution, taux_cible = 0.8, nb_degrad_max = 1000)
    # A COMPLÉTER EVENTUELLEMENT
    t_init = 0    # stupide : pour faire une descente pure !
    # Initialisations diverses et calculs savants !
    # ...
    return t_init
end

function solve!(sv::AnnealingSolver)
    println("BEGIN solve!(AnnealingSolver)")

    error("\n\nMéthode solve!(sv::AnnealingSolver) non implantée: AU BOULOT :-)\n\n")
    # ...

    lg2() && println(get_stats(sv))
    println("END solve!(AnnealingSolver)")
end
