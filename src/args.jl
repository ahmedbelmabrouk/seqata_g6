"""
    Args
Module de gestion des arguments passés à la ligne de commande.
Ce module est spécifique au module Seqata, mais suffisamment générique
pour être adapté à une autre application Julia.
"""
module Args

using ArgParse
using DataStructures: OrderedDict
using Random
using Statistics
using Crayons

using ..Log # pour accès direct aux méthodes ln1(...)

export get_syntaxe, get_usage, parse_commandline
export to_s_dict, get_args, show_args, get, set
export actions_desc, action
export actions, x_solvers, t_algos, presorts
export reset_seed
if !@isdefined(args)
    # global const APPDIR = dirname(dirname(@__FILE__()))
    global args = nothing
    global settings = nothing
    include("abbrev_dict.jl")
end
if @isdefined(args)
    return
end

# Accesseur à un paramètre de la ligne de commande.
# Sans argument, cette méthode retourne le dict qui peut être modifier par
# l'appelant
#
function get(k = "")
    global args
    if args == nothing
        error("Dans args.jl/get : args n'est pas encore défini !")
    end
    # # On s'assure que k est une String, même si un symbole était fourni.
    # k = String(k)
    # On s'assure que k est un Symbol, même si une String était fournie.
    k = Symbol(k)
    if k == ""
        return args
    else
        return args[k]
    end
end

# Permet de lire ou de modifier dynamiquement un paramètre du programme global.
#
# e.g
#   # on diminue provisoirement ce paramètre pour une descente d'un GRASP
#   Args.set(:nb_cons_no_improv_max, 100)
#
function set(k, v)
    global args
    if args == nothing
        error("Dans args.jl/set : args n'est pas encore défini !")
    end
    k = Symbol(k)
    args[k] = v
end

function show_args()
    global args
    if args == nothing
        error("Dans args.jl/show_args() : args n'est pas encore défini !")
    end
    println("Args : état des paramètres")
    show_dict(args)
end

function actions_desc()
    OrderedDict{String,String}(
        "explore"   => "résolution par exploration pure",
        "descent"   => "résolution par une descente avec choix aléatoire du voisin",
        "steepest"  => "résolution par steepest descent",
        "vns"       => "résolution par méthode à voisinage variable",
        "annealing" => "résolution par recuit simulé (sim. annealing)",
        "recuit"    => "résolution par recuit simulé (sim. annealing)",
        "carlo"     => "résolution par méthode de Monte-Carlo",

        "mip"       => "résolution exacte par solver MIP (CPLEX, GLPK, CLP, ...)",
        "dmip"      => "résolution exacte par solver DiscretMip (CPLEX, GLPK, CLP, ...)
                   (générique mais lent, avec une variable binaire par date d'atterrissage)",
        "timing" => "résoud le timing d'une solution à partir d'une instance et d'une liste
                    de noms.
                    --infile : nom du fichier d'instance traitée
                    --planes : liste de noms d'avions imposant l'ordre d'une solution.
                    exemple :  \"[3,4,5,6,8,9,7,1,10,2]\"
                    ou bien :  \" 3 4 5 6 8 9 7 1 10 2 \"",
        "validate" =>
            "vérifie une solution à partir d'une instance et d'une solution complète.
            --infile : nom du fichier d'instance
            --solfile : le fichier solution à valider (au format alp).
            ",
        "stats" =>
            "statistiques sur l'instance (dont test de l'inégalité triangulaire, ...)",
        "test" =>
            "exécute la procédure de test (pour mise au point d'une nouvelle action) ",
        "help" => "Syntaxe de ce programme et description des actions disponibles",
    )
end

# function actions()
#     actions = collect(keys(actions_desc()))
# end
actions() = collect(keys(actions_desc()))
x_solvers() = ["cplex", "cbc", "clp", "glpk", "gurobi", "tulip"]
t_algos() = ["dp", "faye", "lp", "earliest"]
presorts() =
    ["none", "inst", "rinst", "target", "rtarget", "lb", "rlb", "ub", "rub", "shuffle"]
# informats() = ["AUTO", "alp", "alpx", "orlib", "ampl"]
# outformats() = ["AUTO", "alp", "alpx", "ampl"]
# reductors() = ["mean", "median"]

function get_syntaxe(error_msg = "")
    exe = basename(Base.PROGRAM_FILE) # myapp.jl
    # exe = Base.PROGRAM_FILE         # /home/.../.../.../myapp.jl

    if error_msg != ""
        error_msg = "\n    $(error_msg)\n"
    end
    txt = "$(error_msg)
    Syntaxe : $(exe) [action] [--key1 val1] [--flag2] [--key3 val3] ...]\n
    Actions autorisées : $(join(actions(), ","))
    Aide en ligne : $(exe) --help    (ou -h)
    "
end

function get_usage(error_msg = "")
    global settings
    io = IOBuffer()
    println(io, "="^70)
    println(io, "### Liste des actions disponibles")
    println(io, get_actions_dest())
    println(io, "="^70)
    println(io, "### Résumé de la syntaxe")
    println(io, usage_string(settings))
    println(
        io,
        """
Utiliser -h ou --help pour la liste détaillée des options
""",
    )
    String(take!(io))
end

function get_actions_dest()
    io = IOBuffer()
    println(
        io,
        """
Une action est toute abréviation non ambiguë de la liste suivante :
""",
    )
    for (key, val) in actions_desc()
        # println(io, string(key, "\n  $val\n"))
        println(io, key, "\n  ", val)
    end
    String(take!(io))
end

function to_s_dict(args)
    io = IOBuffer()
    for (k, v) in args
        # @show typeof(k)
        if typeof(v) <: Dict
            println(io, " $(rpad(k, 10))  => Dict:")
            # for (local,local) in variable_parent => FONCTIONNE !
            for (k2, v2) in v
                # println(io, "     $k2 => $v2")
                println(io, "     $(rpad(k2, 10))  =>  $v2")
            end
        else
            # println("$k => $v")
            # @show k, v
            if v == nothing
                println(io, " $(rpad(k, 10))  =>  nothing")
            else
                println(io, " $(rpad(k, 10))  =>  $v")
            end
        end
    end
    String(take!(io))
end
function show_dict(args)
    println(to_s_dict(args))
end

# réinitialise la seed du générateur aléatoire avec la valeur
# passée en paramètre.
# Cette méthode peut (doit) être appelée par exemple avant chaque lancement
# d'un test de performance multiple pour rendre le résultat déterministe.
#
# ATTENTION : si ne paramètre vaut 0, le germe n'est jamais le même.
# Cette méthode est donc surtout utile si l'on vaut un comportement déterministe
# à différent endroit du programme.
#
# - si seed == -1 => on utilise la valeur de Args.get[:seed]
# - si seed == 0 => aléatoire (idem si Args.get[:seed] == 0 )
# - si seed > 0 : impose ce germe au générateur aléatoire
#
# Exemple :
# reset_seed()
# => réinitialise le germe avec la valeur passé en paramètre de l'application
#
# reset_seed(seed=Args.get(:seed)+i)
# => initialise le germe avec une valeur dépendante de l'itération i
#
# REMARQUE : cette méthode n'utilise pas le module Log car est doit
# être utilisable avant d'avoir configuré le Log.level() !
function reset_seed(; seed = -1, level = 2)
    if seed == -1
        seed = Args.get(:seed)
    end
    if seed == 0
        used_seed = time_ns()
    else
        used_seed = seed
    end
    if level >= 2
        println("Reset random generator with seed=$used_seed")
    end
    Random.seed!(used_seed)
end

function parse_commandline(inargs = ["none"])
    global args, settings

    # Si l'argument est une String : on le casse en tableau de String
    # Pratique surtout pour les tests mais inutile en fonctionnement normal)
    # e.g parse_commandline("val -p p3,p1,p6,...,p10")
    if isa(inargs, String)
        inargs = split(inargs)
    end

    settings = ArgParseSettings()
    # settings.prog = "seqata.jl"
    settings.description = "Validateur/Solver pour le projet Sequata"
    settings.version = "Seqata-0.3"
    settings.add_version = true
    settings.allow_ambiguous_opts = true
    settings.autofix_names = true
    @add_arg_table! settings begin
        "action"
        help = "action à exécuter.
                Actions autorisées : [$(join(actions(), ", "))]"
        # help = "action à exécuter INCOMPLET....."
        # required = true
        default = "none" # Ne pas imposer d'action facilite les tests unitaires

        "files"
        # ATTENTION : SUBACTIONS OU FICHIERS DOIVENT ETRE CONSECUTIFS !
        # SINON ERROR : TOO MANY ARGUMENTS !
        # (mais pas de subactions pour projet Seqata !)
        nargs = '*'
        help = "arguments supplémentaires (e.g liste de fichiers...).
                L'utilisation est réservée à certaines actions"
        required = false
        # default = Any["popo"] # Vector vide de Any
        default = Any[] # Vector vide de Any

        "--infile", "-i"
        help = "Nom du fichier d'instance."
        # required = true  # 20/06/2018 inhibé pour permettre action help
        default = "NO_INFILE"

        "--outdir", "-d"
        help = "Répertoire d'enregistrement des fichiers de sortie
                (pour les solutions ou la génération d'instance).
                Le répertoire doit exister, sinon \".\" est utilisé."
        default = "_tmp"
        # range_tester = (v->isdir(v))

        "--solfile", "-s" #, "--sol", "--solution", "-S"
        help = "Non du fichier d'une solution initiale ou à valider."
        default = "NO_SOLFILE"

        "--planes", "-p"
        help =
            "Liste de noms d'avions imposant l'ordre d'une solution. " *
            "e.g \"[3,1,4,2...]\" ou \"3,1,4,2...\" ou \"3 1 4 2...\""

        "--force", "-f"
        help = "Force l'écrasement d'un fichier existant"
        action = :store_true

        "--itermax", "-n"
        help = "Nombre d'itérations (sa fonction dépend de l'action : 0 pour AUTO). "
        arg_type = Int
        default = 0    # 0 pour AUTOMATIC

        "--itermax2", "-m"
        help = "Nombre d'itérations pour l'algo secondaire : 0 pour AUTO). "
        arg_type = Int
        default = 0    # 0 pour AUTOMATIC

        "--nb-cons-no-improv-max", "--nim"
        help = "Nombre de non-améliorations consecutives maxi (pour arrêt du recuit/tabou)  (0 pour AUTO)."
        arg_type = Int
        default = 0    # 0 pour AUTOMATIC

        "--nbh"
        help = "choix du voisinage à utiliser (selon les solveurs utilisés).
                **PAS DE VÉRIF POUR L'INSTANT !**
                "
        arg_type = String
        default = "AUTO"    # 19/11/2021

        "--rcl-size", "-k"
        help = "Longueur de la liste de candidat pour le glouton randommisé : 0 pour AUTO)"
        arg_type = Int
        default = 0    # 0 pour AUTOMATIC

        "--presort"
        help =
            "tri initial des avions selon la date d'atterrissage cible. " *
            "Tri possible :  [$(join(presorts(), ", "))] " *
            "(none pour aucun tri ou automatique selon les solveurs)."
        # default = "target"
        default = "none"
        range_tester = (v -> v in presorts())

        "--external-lp-solver", "-x"
        help = "Solver externe : [$(join(x_solvers(), ", "))]"
        default = "cbc"
        range_tester = (v -> v in x_solvers())

        "--timing-algo-solver", "-t"
        help =
            "Algo de résolution pour le sous-problème de timing." *
            "(possible : [$(join(t_algos(), ", "))])"
        default = "lp"
        range_tester = (v -> v in t_algos())

        "--cost-precision", "--prc"
        help =
            "nombre de décimales pour l'affichage des coûts." *
            "Utilisé aussi en interne pour nomaliser le coûts de chaque " *
            "avions en chaque date possible."
        arg_type = Int
        range_tester = (v -> v >= 0) # peut être nul (mais non testé !)
        # default = 5
        default = 7 # modif le 30/10/2018

        "--move-to-first-best", "--fb"
        help = "(steepest, vns ou tabou) se déplacer dès la première amélioration absolue."
        action = :store_true

        "--loglevel", "--level", "-L"
        help = "Niveau de verbosité (pour debug)"
        arg_type = Int
        default = 2

        "--param", "-T"
        help = "Paramètre libre (pour mise au point d'une variante du code)"
        arg_type = String
        default = ""

        "--seed"
        help = "Graine pour le générateur aléatoire"
        arg_type = Int
        default = 0
        range_tester =
            v -> begin
                if v < 0
                    println("\n\nLe germe aléatoire doit être positif ou nul\n")
                    return false
                end
                return true
            end

        "--duration", "--dur"
        help = "Durée de chaque test en secondes pour l'action perf.
                0 pour automatique (1s pour action perf, infini pour taboo).\n
                utilisé avec les actions perf, report et descent, taboo
                ATTENTION : PAS ENCORE EXPLOITÉE POUR RECUIT, NI EXPLORE !"
        arg_type = Int
        default = 0


    end

    global args
    # args = ArgParse.parse_args(settings)
    args = ArgParse.parse_args(inargs, settings, as_symbols = true)
    args[:external_lp_solver] = Symbol(args[:external_lp_solver])
    args[:timing_algo_solver] = Symbol(args[:timing_algo_solver])
    args[:presort] = Symbol(args[:presort])
    args[:nbh] = Symbol(args[:nbh])
    if false
        show_args()
    end

    #
    # ON FAIT QUELQUES VERIFICATIONS ET POST-TRAITEMENTS
    #

    # l'option loglevel est transférée à la globale du module Log
    Log.level(args[:loglevel])

    # on extrait la vraie action à partir de son abréviation
    #
    if args[:action] != "none"
        abbrev = args[:action]
        abbrevs = abbrevsdict(actions()) # voir fichier abbrev.jl
        # @show abbrev
        if !haskey(abbrevs, abbrev)
            msg = "ERREUR : action \"$(abbrev)\" inconnue (dans args.jl)"
            # println(Crayon(foreground = :red, background = :light_yellow, bold = true))
            # println(Crayon(foreground = :red))
            # println("\n$(msg)")
            # println(get_syntaxe())
            # println(Crayon(reset = true))
            # exit(0) # ferait planter les tests unitaires
            throw(ArgumentError(msg))
        end
        # On récupère le nom complet de l'action
        # (si on a passé l'action "val" on remplace par "validate")
        args[:action] = abbrevs[abbrev]
    end
    # e.g. On transforme la chaine "validate" en symbole :validate
    args[:action] = Symbol(args[:action])

    # si action Help : on affiche la liste des actions dispos, puis on quitte
    #
    if args[:action] == :help
        println("="^70)
        println("### Valeur des paramètres après analyse :")
        show_args()
        println(get_usage())
        exit(0)
    end

    # L'affection des paramètres infile et solfile n'est à faire que si
    # on n'est pas dans l'action none.
    if args[:action] != :none

        # Si on a pas passé --infile, on récupère le premier argument sans clé
        # et on le supprime de la liste
        #
        if args[:infile] == "NO_INFILE" && length(args[:files]) >= 1
            args[:infile] = popfirst!(args[:files])
        end

        # Vérification de l'existence du fichier d'instance à traiter
        #
        if args[:infile] != "NO_INFILE" && !isfile(args[:infile])
            # println("ERREUR : fichier inexistant ou illisible : $(args[:infile])")
            # exit(1)
            msg = "ERREUR : fichier inexistant ou illisible : $(args[:infile])"
            throw(ArgumentError(msg))
        end

        # Si on a pas passé --solfile, on récupère le premier argument sans clé
        # et on le supprime de la liste
        #
        if args[:solfile] == "NO_SOLFILE" && length(args[:files]) >= 1
            args[:solfile] = popfirst!(args[:files])
        end

        # si action est timing : on doit extraire la liste des avions à valider
        #
        if args[:action] == :timing
            # planes = matchall(r"[\w]+", args[:planes])
            if typeof(args[:planes]) != String
                msg = "\nERREUR : manque paramètre --planes p1,p2,..."
                throw(ArgumentError(msg))
            end
            planes = collect((m.match for m in eachmatch(r"[\w]+", args[:planes])))
            if length(planes) < 3
                msg = "\nERREUR : pas assez d'avions : $(length(planes)) (au moins trois !)"
                throw(ArgumentError(msg))
            end
            args[:planes] = planes
            println("planes=[", join(args[:planes],","), "]")
        end

        # si action est validate : on extrait le nom du fichier solution à valider
        #
        if args[:action] == :validate
            # Le paramètre solfile doit exister (via option --solfile ou sans clé)
            if args[:solfile] != "NO_SOLFILE" && !isfile(args[:solfile])
                msg = "ERREUR : fichier inexistant ou illisible : $(args[:solfile])"
                throw(ArgumentError(msg))
            end
        end
    end

    # Vérification du répertoire de sortie
    #
    # @show Args.args
    # show_args()
    if !isdir(args[:outdir])
        # println("ERREUR : répertoire inexistant ou illisible : $(args[:outdir])")
        # exit(1)
        if Log.level() >= 3
            println("L3 WARNING : répertoire inexistant ou illisible : $(args[:outdir])")
            println("L3 => on utilisera le répertoire courant : (\".\")")
        end
        args[:outdir] = "."
    end

    # L'action :recuit n'est qu'un alias pour l'action :annealing
    if args[:action] == :recuit
        args[:action] = :annealing
    end

    Args.reset_seed(level = 1)

    return args
end
end # module

