export SteepestSolver
export finished, solve!, sample_two_shifts
export record_bestsol, get_stats

"""
SteepestSolver

Résoud le problème global par la statégie de descente profonde.

"""
mutable struct SteepestSolver
    inst::Instance
    nb_test::Int          # Nombre total de voisins testé
    nb_move::Int          # nombre de voisins acceptés 

    duration::Float64     # durée réelle (mesurée) de l'exécution
    durationmax::Float64  # durée max de l'exécution (--duration)
    starttime::Float64    # heure de début d'une résolution

    cursol::Solution      # Solution courante
    bestsol::Solution     # meilleure Solution rencontrée
    testsol::Solution     # nouvelle solution potentielle

    bestiter::Int
    do_save_bestsol::Bool
    SteepestSolver() = new() # Constructeur par défaut
end

function SteepestSolver(inst::Instance; startsol::Union{Nothing,Solution} = nothing)
    ln3("Début constructeur de SteepestSolver")

    this = SteepestSolver()
    this.inst = inst
    this.nb_test = 0
    this.nb_move = 0

    this.durationmax = 1.0 * 366 * 24 * 3600   # soit 1 année par défaut !
    this.duration = 0.0 # juste pour initialisation
    this.starttime = 0.0 # juste pour initialisation

    if startsol == nothing
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
    this.do_save_bestsol = true
    return this
end

# Retourne true ssi l'état justifie l'arrêt de l'algorithme
#
function finished(sv::SteepestSolver)
    sv.duration = time_ns() / 1_000_000_000 - sv.starttime
    too_long = sv.duration >= sv.durationmax
    other = false
    stop = too_long || other
    if stop
        if lg1()
            println("\nSTOP car :")
            println("     sv.duration=$(sv.duration)")
            println("     sv.durationmax=$(sv.durationmax)")
            println("     (à compléter par vos autres critères d'arrêt)")
            println(get_stats(sv))
        end
        return true
    else
        return false
    end
end

function solve!(
    sv::SteepestSolver;
    startsol::Union{Nothing,Solution} = nothing,
    durationmax::Float64 = 0.0,
    fb::Bool=false,
    nbh::String="AUTO"
)
    ln2("BEGIN solve!(SteepestSolver)")

    if durationmax != 0.0
        sv.durationmax = Float64(durationmax)
    end

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

    sv.starttime = time_ns() / 1_000_000_000

    if lg3()
        println("Début de solve : get_stats(sv)=\n", get_stats(sv))
    end

    ln1("\niter <nb_test> =<nb_move>+<nb_reject> <movedesc> => bestcost=...")

    nom_voisinage=""
    if isequal(nbh,"AUTO")
        nom_voisinage="p2"
    else
        nom_voisinage=nbh
    end
    #print("====> le nom du voisinage est: ",nom_voisinage)
    while !finished(sv)
        movefb=false
        #parametres
        prevcost=sv.cursol.cost
        #nous faisons appel à la fonction generate_nbh
        voisins=generate_nbh(sv.inst.nb_planes,nom_voisinage)[1]
        shuffle!(voisins) #on parcours aleatoirement les voisins
        for v in voisins
            tempsol = Solution(sv.cursol)
            permu!(tempsol,v.indices1,v.indices2)
            sv.nb_test += 1
            degrad = tempsol.cost - prevcost
            ln4("degrad=$(degrad)")
            if degrad < 0
                # Ce voisin est meilleur : on le considere
                
                if fb 
                    sv.nb_move+=1
                    copy!(sv.cursol,tempsol)
                    movefb=true
                    break
                elseif tempsol.cost < sv.testsol.cost
                    copy!(sv.testsol,tempsol)
                end
                lg3("+")
            end
        end
        if !fb && !isequal(sv.cursol.planes,sv.testsol.planes)
            copy!(sv.cursol,sv.testsol) #On prend le meilleur candidat que l'on ait trouve
            sv.nb_move+=1
        elseif !movefb
            break #Pas de meilleur candidat
        end
        # mise a jour éventuelle de la meilleure solution
        if sv.cursol.cost < sv.bestsol.cost
            # La sauvegarde dans bestsol n'est utile que si on ne fait une descente pure
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
        end
    end # fin while !finished

    ln2("END solve!(SteepestSolver)")
end

function record_bestsol(sv::SteepestSolver; movemsg = "")
    copy!(sv.bestsol, sv.cursol)
    println("À COMPLÉTER POUR SEQATA !")
end
function get_stats(sv::SteepestSolver)
    txt = """
    ==Etat de l'objet SteepestSolver==
    (à compléter !)
    """
end

# END TYPE SteepestSolver
