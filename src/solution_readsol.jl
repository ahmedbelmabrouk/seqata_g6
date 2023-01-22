export readsol
export parse_solfile # n'est utilisé que dans les tests

# readsol
#
# Liste des paramètres :
# - sol: l'objet Solution préconstruit et à mettre jour par cette méthode,
# - filename: chemin du fichier solution à charger
#
# Lit un fichier solution au format alp, et effectue les opérations suivantes :
# - extrait tous les paramètres possibles du fichier (sauf les commentaires)
# - met à jour les attributs de l'objet Solution en fonction des informations
#   présentes dans le fichier
# Cette méthode **ne vérifie pas** la validité de la solution vis-à-vis des
# contraintes : seule la cohérence de la solution lue par rapport à son instance
# est assurée.
#
# Les informations suivantes sont lues mais non exploitées :
# - name (e.g. "alp_01_p10_k3")
# - timestamp (e.g "2019-07-08T16:57:22.343")
#
# L'ordre des avions de la solution est défini de la manière suivante selon que
# le fichier de solution à lire est plus ou moins complet :
# 1. Seule la ligne "order" est définie (pas de ligne "landing" dans le fichier)
#    Dans ce cas seul le champ planes de la solution est mis à jour dans l'ordre
#    des noms de la ligne "order".
# 2. Seules les lignes "landing" sont présentes (la ligne "order" est absente)
#    Dans ce cas l'ordre des avions est donné par l'ordre des ligne "landing".
# 3. La ligne "order" **et** les lignes landing sont définies.
#    Dans ce cas l'ordre des avions est défini par la ligne "order" et les lignes
#    "landing" peuvent être dans un ordre quelconque
#
# Les vérifications de cohérence suivantes sont faites :
# - si les informations "order" et/ou "landing" sont présentes, elles doivent
#   être complètes (i.e couvrir tous les avions de l'instance, et seulement
#   ceux-ci).
# - si les informations de date (t) sont fournies pour tous les avions, alors
#   l'attribut x de la solution sera mis à jour.
#   Sinon seules les informations sur l'ordre des avions est utile, et n'a de sens
#   que si la ligne "order" est abscente.
# - si les informations dt et cost sont présentes sur les lignes "landing"
#   alors elle doivent être correctes.
# - si l'information globale "costs" est présente, elle doit être correcte.
#
# - **ne vérifie pas** la validité de la solution vis-à-vis des contraintes :
#   seule la cohérence de la solution lue par rapport à son instance est
#   vérifiée.
#
function readsol(sol::Solution, filename::String)
    DEBUG = false
    ln4("readsol BEGIN $filename")

    data = parse_solfile(sol.inst, filename)
    # Affichage des données lues (pour vérification)
    #
    if DEBUG
        # @show data   # ok, c'est simple mais pas beau !
        println.("  data.", keys(data), " => ", values(data))
        for lan in data.landings
            @show lan
            println("===  $(lan.plane.name) : x=$(lan.time) = cost=$(lan.cost)")
        end
    end

    # On remplit les attributs de la solution (x,costs et cost) à partir du
    # du fichier lue
    if data.has_times
        sol.cost = data.cost
    end
    # data.name != "UNDEF"   (sol.name = data.name)
    # data.has_times         (sol.cost = data.cost)

    for (i, landing) in enumerate(data.landings)
        plane = landing.plane
        sol.planes[i] = landing.plane

        # Les heures et coûts ne sont affectés que si définis dans le fichier.
        if data.has_times
            sol.x[i] = landing.time
        end
        if data.has_times
            sol.costs[i] = landing.cost
        end
    end

    if DEBUG
        println("Rappel de l'instance traitée :")
        println(to_s_long(sol.inst))

        println("Affichage de la solution lue :")
        println(to_s_long(sol))
    end

    ln4("readsol END")
    sol
end

# parse_alp_solfile : analyse un fichier représentant une solution au format .alp
# Retourne une structure représentant les informations brutes du ficher de la
# solution :
# inst: objet de l'instance
# filename: nom du fichier solution
#
function parse_solfile(inst::Instance, filename::String)

    ln4("parse_solfile BEGIN $filename")
    DEBUG = false
    ATOL = 10^-5 #  tolérance absolue pour la précision des coûts

    ok = true # restera true si aucune erreur dans le fichier
    errio = IOBuffer()

    # Initialisation des paramètres à lire
    name::String = "UNDEF" # nom de l'instance
    timestamp::String = "UNDEF"
    cost::Float64 = 0.0
    order_names::Vector{String} = String[] # tableau vide

    # Les paramètres d'atterrissage de chaque avion
    landing_names = String[]   # name: nom de chaque avion (idem Vector{String}())
    landing_times = Int[]      # time: date atterissage (=x dans objet Solution)
    landing_dts = Int[]        # dt: écarts par rapport à la date target
    landing_costs = Float64[]  # cost: coût de pénalité de l'avions

    # PRINCIPE DE L'ANALYSE
    # - on lit chaque ligne séparément
    # - on supprime les commentaires (commençant par "#")
    # - on ignore les lignes vides ou avec espaces
    # - pour chaque ligne
    #   - on extrait la clé (premier mot)
    #   - on exploite les valeurs (les mots suivants)
    #
    lines = readlines(filename)
    lid = 0  # current line number
    while lid < length(lines)
        line = lines[lid+=1]  # autoincrément du numéro de la ligne lue

        # Suppression de tous les commentaires.
        # On peut utiliser un syntaxe fonctionnelle : strip(replace(xxx))
        # ou une syntaxe de style pipe : replace(xxx) |> strip)
        # line = strip(replace(line => r"#.*$" => ""))
        line = replace(line, r"#.*$" => "") |> strip

        # puis on ignore les lignes vides
        if line == ""
            continue
        end
        key, val = extract_key_val(line)
        # @show line
        # @show key,val
        if key == "name"
            # e.g name alp_01_p10_k3
            name = val
            if name != inst.name
                # @warn("ATTENTION name sol:$(name) != inst:$(inst.name)")
                ln2("ATTENTION name sol:$(name) != inst:$(inst.name)")
            end
            continue
        end
        if key == "timestamp"
            # e.g timestamp 2019-07-08T16:57:22.343
            timestamp = val
            continue
        end
        if key == "cost"
            # e.g cost 700.0
            cost = parse(Float64, val)
            continue
        end
        if key == "order"
            # e.g order [p7,p1,p3,...] OU BIEN   order p7 p1 p3...
            # À partir de julia-1.3 on peut faire
            # order_names = map(range->val[range], findall(r"[\w]+", val))
            order_names = collect((m.match for m in eachmatch(r"[\w]+", val)))
            # order_names = [ "p7", "p1", "p3", ...]
            continue
        end
        if key == "landing"
            # landing   p7   135   -4  120.0000000
            words = split(val, r"\s+")
            # @show words
            if length(words) < 1
                # error("A langing line should have at least two fields name and t")
                # error("$lid: A langing line should have at least the field \"name\"")
                msg = "$lid: A langing line should have at least the field \"name\""
                println(errio, msg)
                ok = false
            end
            # pid += 1
            push!(landing_names, words[1])
            if length(words) >= 2
                push!(landing_times, parse(Int, words[2]))
            else
                push!(landing_times, missing)
            end
            if length(words) >= 3
                push!(landing_dts, parse(Int, words[3]))
            else
                push!(landing_dts, missing)
            end
            if length(words) >= 4
                push!(landing_costs, parse(Float64, words[4]))
            else
                push!(landing_costs, missing)
            end
            continue
        end
    end

    # Création de quelques variables auxiliaires pour lisibilité ultérieure
    has_order = length(order_names) != 0
    has_landings = length(landing_names) != 0
    # @show any(x->x==missing, landing_ts)

    # has_times vaut true si aucune baleur ne vaut missing
    has_times = has_landings && all(!ismissing, landing_times)

    # le nom des avions dans l'ordre de la solution
    if has_order
        names = order_names
    elseif has_landings
        names = landing_names
    else
        # println("ERREUR : manque la liste des avions dans le fichier $filename !")
        # exit(1)
        msg = "ERREUR : manque la liste des avions dans le fichier $filename !"
        println(errio, msg)
        ok = false
    end

    # On contruit le tableau des atterrissages dans l'ordre d'etterrissage réel
    landings = []
    for pid = 1:length(names)
        landing = (
            plane = get_plane_from_name(inst, landing_names[pid]),
            time = landing_times[pid],
            dt = landing_dts[pid],
            cost = landing_costs[pid],
        )
        push!(landings, landing)
    end

    ###########################################################################
    # Vérification des informations de la solution lue par rapport à son instance
    ###########################################################################

    # inst_names = [p.name for p in inst.planes]
    inst_names = getfield.(inst.planes, :name)

    # Si la ligne "order" est définie, elle doit être complète
    # Pour ceci, on compare le tableau trié des noms d'avions de l'instance
    # et ceux de la ligne order lue (variable landing_names)
    if has_order && (sort(order_names) != sort(inst_names))
        msg =
            "ERREUR order : les noms de la ligne order ne correspondent " *
            "pas à l'instance."
        # println.(errio, "  ", sort(order_names), " =?= ", sort(inst_names))
        println(errio, msg)
        ok = false
    end

    # Si les lignes "landing" sont définies, elle doivent être complètes
    # Pour ceci, on compare le tableau trié des noms d'avions de l'instance
    # et ceux des lignes "landing" (variable landing_names)
    if has_landings && (sort(landing_names) != sort(inst_names))
        msg =
            "ERREUR landings : Les noms des lignes landing ne correspondent " *
            "pas à l'instance."
        # println.(errio, "  ", sort(landing_names), " =?= ", sort(inst_names))
        println(errio, msg)
        ok = false
    end

    # si les "landing" et les heures d'atterrissages sont définies, on vérifie
    # l'exactitude des infos dt et cost (si présentes).
    if has_times
        for landing in landings
            # ATTENTION il faut déclarer ce "cost" de cet atterrissage local pour
            # ne pas écraser le "cost" total précédent
            local cost
            plane, time, dt, cost = landing  # on dispatch le tuple landing

            if !ismissing(dt) && (dt != time - plane.target)
                msg =
                    "ERREUR $(plane.name) : Retard prétendu incorrect : " *
                    "prétendu=$dt  attendu=$(time-plane.target)"
                println(errio, msg)
                ok = false
            end
            if !ismissing(cost) && !isapprox(cost, get_cost(plane, time), atol = ATOL)
                msg =
                    "ERREUR $(plane.name) : Coût prétendu incorrect : " *
                    "prétendu=$cost  attendu=$(get_cost(plane, time))"
                println(errio, msg)
                ok = false
            end
        end
    end

    sumcosts = mapreduce(land -> land.cost, +, landings)
    # @show cost
    # @show has_times
    # @show sumcosts

    # si les coût individuels des avions sont indiqués, alors le coût globale
    # doit être indiqué et être correct (totélence 5 décimales)
    if has_times && (cost != 0.0)
        sumcosts = mapreduce(land -> land.cost, +, landings)
        if !isapprox(cost, sumcosts, atol = ATOL)
            msg =
                "ERREUR $(name) : coût global incorrect (atol=$ATOL): " *
                "prétendu=$cost attendu=$sumcosts"
            println(errio, msg)
            ok = false
        end
    end

    # S'il y a des d'erreurs d'analyse : on les affiche et on arête
    if !ok
        errmsg = String(take!(errio))
        println(stderr, "Des erreurs ont été rencontrées dans le fichier\n")
        println(stderr, "   $filename\n")
        println(stderr, errmsg)
        error("On arrête !")
        # exit(1)
    end

    result = (
        name = name,
        timestamp = timestamp,
        cost = cost,
        order_names = order_names,
        names = names,
        landings = landings,
        landing_names = landing_names,
        landing_times = landing_times,
        # landing_dts     = landing_dts, # non exporté 
        landing_costs = landing_costs,
        has_order = has_order,
        has_landings = has_landings,
        has_times = has_times,
    )
    DEBUG && @show result
    return result
end
