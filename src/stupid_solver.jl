export StupidSolver
export solve_one!

mutable struct StupidSolver
    inst::Instance
    cursol::Solution      # Solution courante
    bestsol::Solution     # meilleure Solution rencontrée
    StupidSolver() = new() # Constructeur par défaut
end

function StupidSolver(inst::Instance)
    ln2("Début constructeur de StupidSolver")
    this = StupidSolver() # Appel le constructeur par défaut
    this.inst = inst
    this.cursol = Solution(this.inst)
    this.bestsol = Solution(this.cursol)
    return this
end


function solve_one!(sv::StupidSolver)
    ln4("BEGIN solve_one!(StupidSolver)")

    #
    # Phase 1 : on va définir un ordre "intelligent" pour les avions
    #

    # La liste des avions restant à positionner (i.e tous ceux de l'instance !)
    planes_todo = copy(sv.inst.planes)

    # La liste des avions ordonnés pour la solution à construire (initialisé à vide)
    planes_done = Vector{Plane}()
    # Puisqu'on en connait sa taille finale on réserve la place nécessaire
    sizehint!(planes_done, sv.inst.nb_planes)

    # On va choisir intelligemment le prochain avion à insérer
    # 
    # Quelques fonctions julia très évolués pour les échantillonnages :
    #   https://juliastats.org/StatsBase.jl/stable/sampling/#StatsBase.sample!
    #
    while !isempty(planes_todo)
        # Quel est le meilleur avion candidat ?
        p_idx = rand(1:length(planes_todo))

        # On l'a trouvé : on le récupère !
        plane = planes_todo[p_idx]

        # On met à jour nos deux listes d'avions
        push!(planes_done, plane)
        deleteat!(planes_todo, p_idx)

        # On pourrait faire plus simplement
        #   plane = popat!(planes_todo, p_idx)
        #   push!(planes_done, plane)
    end

    #
    # Phase 2 : on utilise la liste d'avions pour imposer l'ordre de la solution
    #

    # On recopie notre nouvelle liste d'avion dans celle de la solution
    # Attention de ne pas changer la référence de l'attribut bestsol.planes par 
    # une simple affectation (même si ça pourrait marcher dans certain cas).
    copy!(sv.cursol.planes, planes_done)

    # On demane à l'objet solution de demander à son XxxxTimingSolver de resoudre 
    # le Sous-Problème de timing
    solve!(sv.cursol)

    ln4("END solve_one!(StupidSolver)")
end

function solve!(sv::StupidSolver; itermax::Int = typemax(Int))
    # ln2("BEGIN solve!(StupidSolver)")

    for i = 1:itermax
        solve_one!(sv)

        # lg(i, " ")
        if sv.cursol.cost < sv.bestsol.cost
            # On mémorise ce nouveau record
            copy!(sv.bestsol, sv.cursol)

            # On enregistre cette nouvelle solution dans un fichier
            write(sv.bestsol)

            if lg1()
                ln3("\n") # saut de ligne pour séparer l'affichage des "." précédents
                # Ce code n'est exécuté que si --level vaut au moins 1
                print(i, ":", sv.bestsol.cost)
                if lg2()
                    # Ce code n'est exécuter que si --level vaut au moins 2
                    print(to_s(sv.bestsol))
                end
                println()
            end
        else
            lg3(".")
        end
    end
    ln2("END solve!(StupidSolver)")
end
