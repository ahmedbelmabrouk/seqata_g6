export Instance
export sep_max, ub_max, lb_min, etp_max, get_sep
export to_s_long, to_s_alp, to_s_alpx, to_s_stats
export get_plane_from_name, get_inequality_viols
# export write  # ce fichier surcharge la méthode Base.write

"""
    Instance

 Encapsule tout ce qui concerne des données d'une instance
 (aucune intelligence ni méthode de résolution).

 Les principaux attributs sont :
 - name: nom de l'instance extrait du fichier
 - planes: un vercteur de Plane (dans l'ordre du fichier d'instance)
 - nb_planes:  nombre d'avion (redondant avec length(planes) )
 - nb_kinds:  nombre de type d'avion
 - sep_mat::Matrix{Int} : matrice des temps de séparation indicée 
   par les numéros d'avion. Ne pas utiliser mais préférer les accesseurs
   get_sep(Instance, Plane, Plane) dont les paramètres sont directement 
   deux objets de type Plane
"""
mutable struct Instance
    planes::Array{Plane}
    name::String
    nb_kinds::Int
    nb_planes::Int
    freeze_time::Int     # UNUSED FOR THIS PROBLEM
    sep_mat::Matrix{Int} # tableau d'éléments de type Int et de dimension 2

    function Instance()
        this = new()
        return this
    end

    function Instance(nb_planes::Int)
        println("ERREUR: génération d'instance aléatoire non implémentée !")
        # voir generate(inst::Instance, nb_planes)
        exit(1)
    end

    # Lit une instance en fonction du paramètre format
    # - Si format vaut alp, ampl ou orlib : ce format est utilisé
    # - Si format vaut AUTO alors l'argument de la ligne de commande
    #   est utilisé :
    #   - si l'argument impose un format, celui-ci est utilisé
    #   - sinon (AUTO) le format est deviné à partir de l'instance
    function Instance(infile::AbstractString; format = "AUTO")
        this = new()
        this.name = "NO_NAME"

        if !isfile(infile)
            println("\nERREUR : fichier \"$infile\" inexistant ou illisible !\n")
            exit(1)
        end
        read_alp(this, infile)
        return this
    end
end
# Quelques méthodes auxiliaires simplifiant le codage
sep_max(inst::Instance) = maximum(inst.sep_mat)
ub_max(inst::Instance) = maximum(p -> p.ub, inst.planes)
lb_min(inst::Instance) = minimum(p -> p.lb, inst.planes)

# on veut "01" à partir de "alp_01_10" mais aussi "01" à partir de "xx_yy-zzz-01-10avions"
shortname(inst::Instance) = match(r"^\D*?(\d+)", inst.name).captures[1] 

# Plus grande valeur de pénalité (d'avance ou de retard)
etp_max(inst::Instance) = maximum(p -> max(p.ep, p.tp), inst.planes)

# Version de haut niveau pour accéder au temps de séparation entre avion
# (les méthodes de bas niveau peuvent accéder directement au tableau inst.sep_mat)
function get_sep(inst::Instance, p1::Plane, p2::Plane)
    inst.sep_mat[p1.kind, p2.kind]
end

# Enregistre l'instance dans un fichier au format demandé
function Base.write(inst::Instance, filename::AbstractString; format = "alpx")
    fh = open(filename, "w")
    println(fh, to_s_long(inst))
    close(fh)
end

# Quelques fonctions raccourcis
to_s_alp(inst::Instance) = to_s_long(inst)
# to_s_alp(inst::Instance) = to_s_long(inst, format="alp")
# to_s_alpx(inst::Instance) = to_s_long(inst, format="alpx")

# génère une chaine au format demandé (ampl, alp, alpx)
# Pour le projet SEQATA, seul le format alp est supporté

function to_s_long(inst::Instance; format::String = "alp")
    if format != "alp"
        error("Format \"$format\" inconnu : seule le forme :alp es supporté!")
    end
    io = IOBuffer()
    println(io, "# ALP instance version 1.0")
    println(io)
    println(io, "name ", inst.name)
    println(io, "nb_planes ", inst.nb_planes)
    println(io, "nb_kinds ", inst.nb_kinds)
    println(io, "freeze_time ", inst.freeze_time)
    println(io)
    # println(io, "#    name  kind   at     E     T     L    ep    tp ")
    println(io, to_s_alp_plane_header())
    for p in inst.planes
        println(io, to_s_long(p))
    end
    println(io)

    println(io, "# Separation time between aircraft kinds")
    for k1 = 1:inst.nb_kinds, k2 = 1:inst.nb_kinds
        println(io, "sep ", k1, " ", k2, " ", inst.sep_mat[k1, k2])
    end
    String(take!(io))
end

# génère une chaine de statistiques sur l'instance
function to_s_stats(inst::Instance; verbose=false)
    io = IOBuffer()
    print(io, "Statistiques sur l'instance $(inst.name)\n")
    println(io)
    println(io, "  name: ", inst.name)
    println(
        io,
        "  nb_planes: ",
        inst.nb_planes,
        " (",
        inst.planes[1].name,
        "..",
        inst.planes[end].name,
        ")",
    )
    println(io, "  nb_kinds: ", inst.nb_kinds)
    println(io, "  freeze_time: ", inst.freeze_time)

    # println(io, "="^70)
    println(io)
    println(io, "  Caractéristiques des avions :\n")
    println(
        io,
        "  at:     ",
        minimum(p -> p.at, inst.planes),
        "..",
        maximum(p -> p.at, inst.planes),
    )
    println(
        io,
        "  lb:     ",
        minimum(p -> p.lb, inst.planes),
        "..",
        maximum(p -> p.lb, inst.planes),
    )
    println(
        io,
        "  target: ",
        minimum(p -> p.target, inst.planes),
        "..",
        maximum(p -> p.target, inst.planes),
    )
    println(
        io,
        "  ub:     ",
        minimum(p -> p.ub, inst.planes),
        "..",
        maximum(p -> p.ub, inst.planes),
    )
    println(
        io,
        "  ub-lb:  ",
        minimum(p -> (p.ub - p.lb), inst.planes),
        "..",
        maximum(p -> (p.ub - p.lb), inst.planes),
    )
    println(
        io,
        "  ep:     ",
        minimum(p -> p.ep, inst.planes),
        "..",
        maximum(p -> p.ep, inst.planes),
    )
    println(
        io,
        "  tp:     ",
        minimum(p -> p.tp, inst.planes),
        "..",
        maximum(p -> p.tp, inst.planes),
    )

    viols = get_inequality_viols(inst)
    println(io)
    println(io, "  Test de get_inequality_viols : nb_viols=$(length(viols))")
    if verbose
        for viol in viols
            (k1, k2, k3, sep12, sep23, sep13) = viol
            print(io, "  $k1->$k3=$sep13 >= ")
            println(io, "$k1->$k2=$(sep12) + $k2->$k3=$(sep23))")
        end
    end
    String(take!(io))
end

# Retourne l'avion à partir de son nom
#
function get_plane_from_name(inst::Instance, name::AbstractString)
    # TODO : utiliser plutot une fonction de recherche de Julia
    for plane in inst.planes
        if plane.name == name
            return plane
        end
    end
    # Cas particulier où les noms d'avion sont de la forme "p1" alors que l'on
    # recherche le nom "1"
    for plane in inst.planes
        if plane.name == "p$name"
            return plane
        end
    end
    error("\nAucun avion de nom : $(name)\n")
end

# Retourne la liste des viols de l'hypothese de l'inégalité triangulaire
#
#    sep(k1,k3) > sep(k1,k2) + sep(k2,k3)  => viol existe
#
# En effet, l'inégalité triangulaire impose que pour toute paire de types
# d'avions (k1,k3), on ne puisse pas insérer un autre avion de type k2 qui
# permette de raccourcir l'écart des temps d'atterrissage entre un avion de type
# k1 et un avion de type k3.
#
# Un viol est représenté par un tuple de six entiers
# (k1, k2, k3, sep12, sep23, sep13)
#
function get_inequality_viols(inst::Instance)
    viols = Vector{Tuple{Int,Int,Int,Int,Int,Int}}()
    if inst.nb_kinds < 2
        # pas de problème d'inégalité triangulaire car moins de 3 types d'avions !
        return viols
    end
    for k1 = 1:inst.nb_kinds, k3 = 1:inst.nb_kinds, k2 = 1:inst.nb_kinds
        @inbounds sep13 = inst.sep_mat[k1, k3]
        @inbounds sep12 = inst.sep_mat[k1, k2]
        @inbounds sep23 = inst.sep_mat[k2, k3]
        if sep13 > sep12 + sep23
            viol = (k1, k2, k3, sep12, sep23, sep13)
            push!(viols, viol)
        end
    end
    viols
end

# END TYPE Instance
