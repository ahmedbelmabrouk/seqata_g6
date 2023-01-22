export Plane, get_cost_basic, get_cost
export get_cost_basic # seqata
export to_s, to_s_alp, to_s_long

# # Quelques variables globale pour émuler les variables  de class en POO
# current_plane_id = 0

"""
    Plane
 Encapsule les données d'un avion.

 - id: numéro interne commençant en 1, Utilisé pour l'indexation des vecteurs
 - name: nom arbitraire, par exemple le numéro commençant en "p1", "p2", ...
   Mais pour l'instant, c'est un entier pour être conforme à ampl
 - kind: le type de l'avion (car kind est un mot clé réservé en julia !)
 - at: heure d'apparition de l'avion dans l'oeil du radar (inutilisé)
 - lb, target, ub: les heures mini, souhaitées et maxi d'atterrissage
 - ep, tp : heures d'atterrissage au lus tot et ou plus tard

 REMARQUE : le type Plane est IMMUTABLE, il est lié à un avions de
 de l'instance et est défini une fois pour toute.

"""
struct Plane
    id::Int
    name::AbstractString
    kind::Int
    at::Int # appearing time : UNUSED IN THIS PROBLEM
    lb::Int # lowest bound = lowest time
    target::Int
    ub::Int # upper bound = upper time
    ep::Float64 # earliness penalty
    tp::Float64 # tardiness penalty

    # Précalul des coûts (de 1 à p.ub)
    # costs::Vector{Float64}
    costs::Tuple{Vararg{Float64}}

    # Constucteur complet avec paramètres passés par association.
    # 
    function Plane(;
            id, name, kind,
            at,
            lb, target, ub, 
            ep, tp,
        )

        @assert lb <= target <= ub

        # On précalcule les coûts réels en fonction de la date d'atterrissage.
        # Lors de la lecture du coût (par get_cost), les valeurs hors borne lb:ub 
        # seront remplacées par une pénalité très élevée (dépendant du viol).
        # 
        # TODO : optimiser la taille du tableaux costs pour éviter de mémoriser 
        # les valeurs en dessous de ub :
        # tmp_costs = fill(BIG_COST, ub-lb+1) # pour optimiser la mémoire utilisée
        # 
        tmp_costs = Vector{Float64}(undef, ub)
        @inbounds for t in lb:ub
            tmp_costs[t] = t < target ? ep * (target - t) : tp * (t - target)
        end

        return new(
            id, name, kind,
            at,
            lb, target, ub, 
            ep, tp,
            Tuple{Vararg{Float64}}(tmp_costs),
        )
    end
end

# Méthode Julia pour convertir tout objet en string (Merci Matthias)
function Base.show(io::IO, p::Plane)
    Base.write(io, to_s(p))
end


"""
get_cost_basic(p::Plane, t::Int)

Calcul du coût d'un avion pour une fonction en V sans pénalité si hors borne.
"""
function get_cost_basic(p::Plane, t::Int)
    return t < p.target ? p.ep * (p.target - t) : p.tp * (t - p.target)
end

"""
    get_cost(p::Plane, t::Int; violcost::Float64 = 1000.0)

Retourne le coût de l'avion, éventuellement pénalisé si hors bornes.
"""
@inline function get_cost(p::Plane, t::Int; violcost::Float64 = 1000.0)
    if (t in p.lb:p.ub)
        @inbounds return p.costs[t]
    else
        return violcost * max(p.lb - t, t - p.ub)
    end
end

# return simplement le name (e.g. "p1")
function to_s(p::Plane)
    p.name
end
# return e.g. : "[p1,p2,p3,..,p10]"
function to_s(planes::Vector{Plane})
    # string("[", join( [p.name for p in planes], "," ), "]")
    # string("[", join( (p->p.name).(planes), "," ), "]")
    string("[", join(getfield.(planes, :name), ","), "]")
end

"""
    to_s_alp(p::Plane)

Retourne la représentation String de l'avion conforme au format d'instance alp
"""
to_s_alp(p::Plane) = to_s_long(p)
# to_s_alp(p::Plane) = to_s_long(p, format = "alp")
# to_s_alpx(p::Plane) = to_s_long(p, format = "alpx")

function to_s_long(p::Plane)
    io = IOBuffer()
    print(io, "plane ")
    print(io, lpad(p.name, 3), " ")
    print(io, lpad(p.kind, 4), " ")
    print(io, lpad(p.at, 5), " ")
    print(io, lpad(p.lb, 5), " ")
    print(io, lpad(p.target, 5), " ")
    print(io, lpad(p.ub, 5), "    ")
    print(io, lpad(p.ep, 4), " ")
    print(io, lpad(p.tp, 4), " ")
    String(take!(io))
end
# Retourne un commentaire décrivant une ligne au forme alp ou alpx
# Attention : pour le projet Seqata, seul le format alp existe)
# soit: #    name  kind   at     E     T     L    ep    tp
# soit: #    name  kind   at     E     T     L    dt1 cost1   dt2 cost2 ...
# Il n'y a pas de return final
#
function to_s_alp_plane_header()
    io = IOBuffer()
    print(io, "#    name  kind   at     E     T     L")
    print(io, "    ep    tp ")
    String(take!(io))
end

# Affiche les éléments définis ( != -1 ) du tableau précalculé costs 
# - d'une part les costs précalculés en chaque date
# - d'autre part les coûts calculés sur demande pour des dates arbitraires
#   (et mémoïsés)
#
function to_s_costs(p::Plane)
    io = IOBuffer()
    print(io, p.name, "=>costs[]= ")
    for t = 1:length(p.costs)
        p.costs[t] <= -1.0 && continue
        print(io, " ", t, ":", p.costs[t])
    end
    String(take!(io))
end

#####################################################################
# COMPLÉMENT POUR DÉBUT D"EXTENSION POUR COÛT DES AVIONS MULTI-PENTES.
# N'EST PAS UTILISÉ POUR LES ÉLÈVES POUR L'INSTANT
#
# Un breackpoint est un point de cassure de la fonction de coût d'un avion.
# Un timecost est un couple lié à un breakpoint de l'avion et associant :
# - la date absolue du breackpoint (entre lb et ub inclues)
# - le coût réel de l'avion s'il atterrit à cette date
#
# Le type BreakPoint encapsule :
# - la date d'un breakpoint d'avion,
# - la pente du coût à partir de cette date
#
# TODO: simplifier pour Seqata en créant un nouveau type BreakPoint immutable
# avec les attributs suivants :
#   time::Int, cost::Float64, slope::Float64
# 
# BreakPoint doit être mutable car les champs time et slope sont mis à jour 
# dans FayeTimingSolver
mutable struct BreakPoint
    time::Int
    slope::Float64
end

# return les timecosts en fonction des caratéristiques de l'avion.
# En principe, il y a trois timecosts associés aux dates lb, target et ub
# 
# Mais on peut avoir seulement deux timecosts dans le cas dégénéré pour lequel
# le target correspond à une des extrémités (lb ou ub)
# Le coûts sont supposés linéaire "en V assymétrique"
# 
function get_timecosts_from_eptp(p::Plane)

    prc = Args.get(:cost_precision)
    timecosts = Vector{Tuple{Int,Float64}}()

    cost_lb = round(p.ep * (p.target - p.lb), digits = prc)
    tc_lb = (p.lb, cost_lb)
    push!(timecosts, tc_lb)

    if p.lb != p.target && p.target != p.ub
        # Si target est une des extémités, on n'aura que deux timecosts
        tc_target = (p.target, 0.0)
        push!(timecosts, tc_target)
    end

    cost_ub = round(p.tp * (p.ub - p.target), digits = prc)
    tc_ub = (p.ub, cost_ub)
    push!(timecosts, tc_ub)

    return timecosts
end

function get_bpoints(p::Plane)
    # n = length(p.timecosts) # ORO de Alap
    n = length(get_timecosts_from_eptp(p)) # NEW pour Seqata mais INEFFICACE À REVOIR
    bpoints = Vector{BreakPoint}(undef, n)
    timecosts = get_timecosts_from_eptp(p)
    for i = 1:n-1
        time1, cost1 = timecosts[i]
        time2, cost2 = timecosts[i+1]
        bpoints[i] = BreakPoint(time1, (cost1 - cost2) / (time1 - time2))
    end
    bpoints[n] = BreakPoint(timecosts[n][1], 1000000.0)
    return bpoints
end


