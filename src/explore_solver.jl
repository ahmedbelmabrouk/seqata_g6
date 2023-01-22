export ExploreSolver, solve

"""
    ExploreSolver

Résoud le problème global par statégie d/ exploration aveugle.

Un voisin est tiré aléatoirement (voisinage large) et est systématiquement 
accepté.
Par défaut, la solution initiale est choisie aléatoirement
"""
mutable struct ExploreSolver
    inst::Instance
    nb_test::Int          # Nombre total de voisins testés
    nb_move::Int          # nombre de voisins acceptés (améliorants ou non)
    nb_degrad::Int        # nombre de mouvements dégradants
    nb_improve::Int       # Nombre de mouvements améliorants

    cursol::Solution      # Solution courante
    bestsol::Solution     # meilleure Solution rencontrée

    function ExploreSolver(inst::Instance)
        this = new()
        this.inst = inst
        this.nb_test = 0
        this.nb_move = 0
        this.nb_degrad = 0
        this.nb_improve = 0

        this.cursol = Solution(inst)
        ln1("Solution correspondant à l'ordre de l'instance")
        ln1(to_s(this.cursol))

        if Args.get(:presort) == :none
            initial_sort!(this.cursol, presort = :shuffle)
        else
            initial_sort!(this.cursol, presort = Args.get(:presort))
        end
        ln1("Solution initiale envoyée au solver")
        ln1(to_s(this.cursol))

        this.bestsol = Solution(this.cursol)

        return this
    end
end

function solve!(sv::ExploreSolver, itermax_max::Int)
    ln2("BEGIN solve!(ExploreSolver, itermax_max=$itermax_max)")
    itermax = 1 # car on veut faire une seule itération si on passe itermax_max=1

    lg1("iter <nb_move>=<nb_improve>+<nb_degrade> => <bestcost>")

    while itermax <= itermax_max
        prevcost = sv.cursol.cost

        # On pourrait paramétrer la distance maxi du voisin
        shift_max = -1   # -1 pour pas de limite !
        # shift_max = 2
        # shift_max = 1
        if shift_max == -1
            i1 = i2 = -1
        else
            i1 = rand(1:length(sv.inst.planes))
            i2 = i1
            while i2 == i1
                i2_min = max(i1 - shift_max, 1)
                i2_max = min(i1 + shift_max, length(sv.inst.planes))
                i2 = rand(i2_min:i2_max)
            end
        end
        # @show i1, i2

        swap!(sv.cursol, i1, i2)
        # println("APRES SWAP: ", to_s(sv.cursol))
        sv.nb_move += 1
        degrad = sv.cursol.cost - prevcost
        ln4("degrad=$(degrad)")
        if degrad < 0
            # Ce voisin est meilleur : on l'accepte
            lg3("+")
            sv.nb_improve += 1
            # mise a jour éventuelle de la meilleure solution
            if sv.cursol.cost < sv.bestsol.cost
                # La sauvegarde dans bestsol n'est utile que si on ne fait une descente pure
                copy!(sv.bestsol, sv.cursol)
                if lg1()
                    msg =
                        string("\niter ", sv.nb_move, "=", sv.nb_improve, "+", sv.nb_degrad)
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
        else # degrad < 0
            # Ce voisin est plus mauvais : on l'accepte aussi (car exploration) !!
            sv.nb_degrad += 1
            if Log.level() in 3:3      # idem à : if lg3() && !lg4()
                print("-")
            end
            if lg4()
                msg = string(
                    "\n     ",
                    sv.nb_move,
                    ":",
                    sv.nb_improve,
                    "+/",
                    sv.nb_degrad,
                    "- cursol=",
                    to_s(sv.cursol),
                )
                print(msg)
            end
        end
        itermax += 1
    end  # while itermax
    ln2("\nEND solve!(ExploreSolver)")
end
