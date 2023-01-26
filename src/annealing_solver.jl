export AnnealingSolver
export finished, solve! ,get_stats, guess_temp_init
#export solve!
"""
anealing solver

"""
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
    prevsol::Solution       # Prend la valeur de solution courante avant changement
    AnnealingSolver() = new() # Constructeur par défaut
    
end

function AnnealingSolver(inst::Instance, user_opts = Dict())
    # ...
    #error("\n\nConstructeur de AnnealingSolver non implanté : AU BOULOT :-)\n\n")
    # ...
    this =AnnealingSolver()
    this.inst = inst
    this.nb_test = 0
    this.nb_move = 0

    this.opts = Dict(
        :startsol => nothing,    # nothing pour auto à partir de l'instance
        :step_size => 1,
        :temp_init => -1.0, # -1.0 pour automatique
        :temp_init_rate => 0.75,  # valeur standard : 0.8
        :temp_mini => 0.000_001,
        :temp_coef => 0.999_95,
        :nb_cons_reject_max => 1_000_000_000, # infini
        :nb_cons_no_improv_max => 5000 * inst.nb_planes)
        nb_options = length(this.opts) 
    #print(typeof(this.opts[:startsol]))
    merge!(this.opts, user_opts)

    if length(this.opts) != nb_options
        error("Au moins une option du recuit est inconnue dans :\n$(keys(user_opts))")
    end

    if user_opts[:startsol] == nothing
        # Pas de solution initiale => on en crée une
        this.cursol = Solution(inst)
        if Args.get(:presort) == :none
            # initial_sort!(this.cursol, presort=:shuffle) # proto
             initial_sort!(this.cursol, presort = :target) # diam
        else
            initial_sort!(this.cursol, presort = Args.get(:presort))
        end
    else
        this.cursol = startsol
        if lg2()
            println("Dans SteepestSolver : this.cursol = this.opts[:startsol] ")
            println("this.cursol", to_s(this.cursol))
        end
    end
    
    this.bestsol = Solution(this.cursol)
    this.testsol = Solution(this.cursol)
    this.prevsol = Solution(this.cursol)
    return this
end

# stop : retourne true ssi l'état justifie l'arrêt de l'algorithme
# (dommage qu'on ne puisse pas l'appeler stop? comme en ruby !)
# On pourra utiliser d'autres critères sans toucher au programme principal

function finished(sv::AnnealingSolver)
    #sv.duration = time_ns() / 1_000_000_000 - sv.starttime
    #too_long = sv.duration >= sv.durationmax
    #print("nb test ===>",sv.nb_test)
    #too_long = sv.nb_test >= 10 000 
    too_long =false
    other = false 
    stop = too_long || other
    if stop
        if lg1()
            println("\nSTOP car :")
            println("     sv.nb_test=$(sv.nb_test)")
            #println("     sv.nb_move=$(nb_move)")
            #println("     (à compléter par vos autres critères d'arrêt)")
            println(get_stats(sv))
        end
        return true
    else
        return false
    end
end

function get_stats(sv::AnnealingSolver)
    txt = "
    
    
    nb_test=$(sv.nb_test)
    nb_move=$(sv.nb_move)
    
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
    t_init = 50    # stupide : pour faire une descente pure !
    # Initialisations diverses et calculs savants !
    # ...
    return t_init
end
function new_temp(t,c = 0.01)
    return (1-c)t
end

function solve!(sv::AnnealingSolver,startsol::Union{Nothing,Solution} = nothing)
    println("BEGIN solve!(AnnealingSolver)")

    #error("\n\nMéthode solve!(sv::AnnealingSolver) non implantée: AU BOULOT :-)\n\n")
    # ...
    """
    if durationmax != 0.0
        sv.durationmax = Float64(durationmax)
    end
    """
    #Phase de Construction / initialisation 
    if startsol != nothing
        sv.cursol = startsol
        copy!(sv.bestsol, sv.cursol) # on réinitialise bestsol à cursol
        copy!(sv.testsol, sv.cursol)
        if lg2()
            println("Dans steepestSolver : sv.cursol = sv.opts[:startsol] ")
            println("sv.cursol : ", to_s(sv.cursol))
        end
    else
        # on garde la dernière solution sv.cursol
    end
    t=guess_temp_init(sv.cursol)
    #sv.starttime = time_ns() / 1_000_000_000

    if lg3()
        println("Début de solve : get_stats(sv)=\n", get_stats(sv))
    end

    ln1("\niter <nb_test> =<nb_move>+<nb_reject> <movedesc> => bestcost=...")
    #Phase d'amelioration 
    while !finished(sv) && t > 0.001
        #DEBUT====> VND=variable neighborhood descent <====
                 #initialisation
        #neighborhood_structures=[["P3","p4","p5","p6","p7","p8"],["T3","t4","t5","t6","t7","t8"],["s3","T4","S4","d4","D12.7","2T4","2T4g12"]]
        neighborhood_structures=["P3","T4","S4","d4","D12.7","2T4","2T4g12"] #une série de voisinages aleatoire
        #neighborhood_structures=["P3","p4","p5","p6","p7","p8"] #une série de voisinages croissants et non redondants
        #neighborhood_structures=["T3","t4","t5","t6","t7","t8"] #une série de voisinages croissants et non redondants
                 #amelioration
        improv=true
        copy!(sv.prevsol,sv.cursol)
        while  improv
            improv=false
            k=1
            #voisins=generate_nbh(sv.inst.nb_planes,neighborhood_structures[k])[1] 
            while k < length(neighborhood_structures)
                prevcost=sv.cursol.cost
                voisins=generate_nbh(sv.inst.nb_planes,neighborhood_structures[k])[1]
                #Exploration complete de Voisinage :on cherche le meilleur Voisinage
                
                for v in voisins
                    tempsol = Solution(sv.cursol)
                    permu!(tempsol,v.indices1,v.indices2)
                    sv.nb_test += 1
                    degrad = tempsol.cost - prevcost
                    ln4("degrad=$(degrad)")
                    if degrad < 0 # Ce voisin est meilleur : on le considere
                        tempsol.cost < sv.testsol.cost
                        copy!(sv.testsol,tempsol)
                    end
                    lg3("+")
                end
                #Move or NOT 
                if !isequal(sv.cursol.planes,sv.testsol.planes)
                    copy!(sv.cursol,sv.testsol) #On prend le meilleur candidat que l'on ait trouve
                    sv.nb_move+=1
                    k=1
                    improv=true
                    println("=============NEW SOL IMPROVED===============")
                else
                    k+=1 #on incremente k
                    println("=======We didn't move===: ",k,"=====")
                end
                
            end
        end
        #FIN====> VND=variable neighborhood descent <====
        
        #Critere d'acceptation 
        delta = sv.cursol.cost - sv.prevsol.cost
        println(">> DELTA= ",delta)
        if (delta < 0) 
            # Mise à jour de sv.bestsol
            copy!(sv.bestsol, sv.cursol)
            if lg1()
                msg = string("\niter ", sv.nb_test, ", move ", sv.nb_move)
                if lg2()
                    # affiche coût + ordre des avions
                    msg *= string(" => ", to_s(sv.bestsol))
                else
                    # affiche seulement le coût
                    msg *= string(" => ", sv.bestsol.cost)
                end
                print(msg)
            end
        elseif ((x=rand()) < (y=exp(-delta/t)))
            # on garde la dernière solution sv.cursol
            println("la valeur de x: ",x, "et y: ",y)
            if lg1()
                msg = string("\niter ", sv.nb_test, ", move ", sv.nb_move)
                if lg2()
                    # affiche coût + ordre des avions
                    msg *= string(" => ", to_s(sv.cursol))
                else
                    # affiche seulement le coût
                    msg *= string(" => ", sv.cursol.cost)
                end
                print(msg)
            end
        else
            #on n'accepte pas la nouvelle solution sv.cursol
            copy!(sv.cursol, sv.prevsol)
        end
        #Mise a jour de T
        t=new_temp(t)
        println("\n La nouvelle temperature : ",t)

    end # fin while !finished
    #lg2() && println(get_stats(sv))
    ln2("END solve!(AnnealingSolver)")
end
