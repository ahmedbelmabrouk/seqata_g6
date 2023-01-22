#!/bin/sh
#= Les lignes **shell** suivantes sont des commentaires multilignes ignorées
#  par julia mais elle permet de lancer julia avec des options arbitraires.
#
#### DÉBUT DU CODE EN SHELL UNIX (e.g bash, ....)
#
if test -z $SEQATA_SYSIMG; then
    export SEQATA_SYSIMG="/tmp/julia_sysimg_seqata.so"
fi

echo "SEQATA_SYSIMG=$SEQATA_SYSIMG"
echo "JULIA_USING_ALL=$JULIA_USING_ALL"
echo "JULIA_MS_LOG=$JULIA_MS_LOG"
if test -e "$SEQATA_SYSIMG"; then
    _sysimg=-J"$SEQATA_SYSIMG"
    echo "ATTENTION : on utilise l'image SEQATA_SYSIMG=$SEQATA_SYSIMG"
    echo "VÉRIFIER QU'ELLE N'EST PAS PÉRIMÉE."
    echo "VOUS POUVEZ REGÉNÉRER L'IMAGE PAR :"
    echo "   rm \"$SEQATA_SYSIMG\""
    echo "   ./bin/build_sysimg.jl"
else
    _sysimg=""
fi
exec julia --project --color=yes --startup-file=no $_sysimg -- "$0" "$@"
# Autre options possibles pour julia
#  --depwarn=no  : pour supprimer d'éventuels warnings
#  --startup-file=no : évite de charger le fichier utilisateur
=#

#### DÉBUT DU CODE JULIA

# Le but de ce script est de tester les performances de la brique STP pour
# tout TimingSolver disponible.
# Ce script fait partie du projet Seqata corrigée mais il doit pouvoir tester la
# performance des brique STP des projet d'élève.
# On suppose que le projet des éleves et conforme au prototype fourni (même API
# pour les solveurs).
#
# Utilisation :# ./bin/showd_stp_perfs.jl [projdir]
# - si le paramètre projdir est passé : il doit correspondre au projet d'élève à mesurer.
# - sinon on var mesure le projet courant.
#
# TODO (04/02/2022)
# - prévoir des option pour
#   --duration alias --dur pour la durée d'un test (e.g 30s)
#   --use_descent true|false pour tester via la descent ou non (défaut=false)
#

# Le realpath est nécessaire si l'un des répertoires parents est un lien symbolique :
appdir = dirname(dirname(realpath(@__FILE__())))


##############################################################
# LECTURE PRECOCE DES ARGUMENTS (VANT LECTURE PAR SEQATA)
##############################################################
function show_usage()
    println("show_stp_perf.jl [-h] [projetdir] [durmax [true|false]")
    println("  projetdir: répertoire du projet seqata à traiter")
    println("             Si ce répertoire est passé, il doit l'être en premier")
    println("  durmax: entier de 1 à 120 : durée d'un test pour un algo d etiming")
    println("  use_descent: (true|false) : si true le test est effectué via une descente")
    println("Exemple:")
    println("  show_stp_perf.jl . 2")
    println("  => mesure du répertoire courant (d'élève) avec durée 2sec et mesure STP directe")
    println("  show_stp_perf.jl . 2 true")
    println("  => mesure du répertoire courant (d'élève) avec durée 2sec et mesure STP via descente")
    println("  show_stp_perf.jl")
    println("  => mesure projet diam pendant 5s directe i.e. (sans descente)")
    println("  show_stp_perf.jl 10")
    println("  => mesure projet diam pendant 10s directe i.e. (sans descente)")
end
mutable struct LocalArgs
    durmax::Int
    use_descent::Bool
    hack_pour_eleves::Bool
    projdir::String
    LocalArgs() = new()
end
function show_fields(args::LocalArgs, msg::String="")
    println(msg)
    for name in fieldnames(typeof(args))
        println("-- ", rpad(name,12), " = ", getfield(args, name))
    end
end
function parse_args()
    args = LocalArgs()
    args.durmax=5
    args.use_descent=false
    args.hack_pour_eleves=false

    if length(ARGS) >= 1 && ARGS[1]=="-h"
        show_usage()
        exit(0)
    end

    # Ancienne analyse mono-argument (dir d'un projet élève)
    if length(ARGS)>=1 && isdir(ARGS[1])
        args.projdir = realpath(ARGS[1])
        popfirst!(ARGS) # sinon plantage car les ARGS sont passés à seqata qui ne reconnait pas l'action
        args.hack_pour_eleves = true
    else
        args.projdir = appdir
        args.hack_pour_eleves = false
    end

    val = "" # doit être défini pour exister dans les catch suivants
    if length(ARGS) >= 1
        try
            val = popfirst!(ARGS)
            args.durmax = parse(Int, val)
            @assert args.durmax in 1:120
        catch
            println("\n** durmax=$(val) incorrect. Doit être entier entre 1 et 120")
            show_usage()
            exit(1)
        end
    end
    if length(ARGS) >= 1
        try
            val = popfirst!(ARGS)
            args.use_descent = parse(Bool, val)
        catch
            println("\n** use_descent=$(val) incorrect. Doit être représenter un Bool")
            show_usage()
            exit(1)
        end
    end
    @show args
end
# Analyse des arguments
global args = parse_args()

##############################################################
# CHARGEMENT DU MODULE PRINCIPAL DU PROJET À TESTER
##############################################################

using Pkg
Pkg.activate(args.projdir)
Pkg.instantiate()


if args.hack_pour_eleves
    # VERRUE SPECIFIQUE POUR TESTER LES PROJETS ELEVES
    using CPLEX
    using Gurobi
    global GRB_ENV
    if !@isdefined(GRB_ENV)
        GRB_ENV = Gurobi.Env() # affiche warning pénible (academic)
    end
end


ENV["JULIA_USING_ALL"] = 1  # for loading all package
include("$(args.projdir)/src/Seqata.jl")
using .Seqata
using .Log
# @ms using Debugger # chargement du debugger (voir doc)
# @ms using Revise   # fonctionnement à vérifier avec le projet BipSolver
@ms include("$APPDIR/src/interactive.jl")
# i01()=i01
# i09()=i09
# i13()=i13
@ms using Printf


# Encapsule le résultat d'un test de la brique STP
mutable struct Result
    # Nom de l'instance
    instname::String # alp_13_a10

    # Nom comniné de l'algo STP (timing_algo_solver)
    #     [:faye, :dp, :lp, :lp2, :lp3, :lp4]
    # et du solver externe éventuel (external_lp_solver) pour les algos lpX
    #     [:cplex, :gurobi, :clp, :tulip, :glpk].
    # Soit par exemple :
    #     :faye, :dp, "dp", "lp/cplex", "lp/gurobi", ...
    fullalgo::Union{String, Symbol} # e.g. [:faye, :dp, :lp, :lp2, :lp3, :lp4]

    # Nombre d'applel à la STP par seconde (call_per_sec)
    cps::Float64
end

function to_s(res::Result)
    io = IOBuffer()
    print(io, "plane ")
    print(io, rpad(repr(res.instname), 14), " ")
    print(io, rpad(repr(res.fullalgo), 12), " ")
    print(io, lpad(repr(res.cps),      10), " ")
    String(take!(io))
end

# on veut "01" à partir de "alp_01_10" mais aussi "01" à partir de "xx_yy-zzz-01-10avions"
# shortname(inst::Instance) = match(r"^\D*?(\d+)", inst.name).captures[1]
shortname(inst::Instance) = match(r"^\D*(\d+)", inst.name).captures[1]

# get_timing_solver_perf : test d'un TimingSolver particulier.
#
# Ce teste peut se faire via un descente ou directement en manipulant une
# solution. L'utilisation de la descente peut introduire des appels parasites
# si en est mal implémentée, ou si elle effectue des tâche annexe compliqués
#
# inst : objet Instance à tester
#
# fullalgo : combinaison de timing_algo_solver et de external_lp_solver
#    Peut-être un Symbol ou un String (:lp, "lp", "lp/cplex", "lp+cplex")
#
# durmax: durée de la descente en secondes
#
# use_descent (false) : si true la descente est utilisée sinon (par défaut)
#    Une solution est crée et manipulé directement pour les tests.
#
# Exemple :
#   perf = get_stp_perf(inst, "lp/gurobi", 30)
#
function get_timing_solver_perf(inst::Instance, fullalgo; 
                durmax=5, use_descent=false, verbose=true)
    if verbose
        lg0("test $(inst.name) pour $(fullalgo) pendant  durmax=$(durmax)sec ")
        ln0("use_descent=$(use_descent) ($(ms()))... ")
    end

    # Expraction et positionnement des paramètres définiassant le timingSolver
    # e.g DynProgTimingsolver (sans solver externe), LpTiminSolver (avec cplex ou autre)
    fullalgo = String(fullalgo)
    words = split(fullalgo, r"\W")
    if length(words) == 1
        Args.set("timing_algo_solver", Symbol(words[1]))
        Args.set("external_lp_solver", nothing)
    elseif length(words) == 2
        Args.set("timing_algo_solver", Symbol(words[1]))
        Args.set("external_lp_solver", Symbol(words[2]))
    else
        error("fullalgo incorrect : #(fullalgo). Autorisé :  :dp, \"dp\", \"lp/cplex\", ...")
    end

    bestsol = nothing
    nb_calls = 0
    Log.pushlevel!(0)
    ms_start = ms()
    if use_descent
        # On teste le TimingSolver via la descente
        sv = DescentSolver(inst)
        # solve!(sv, durationmax=1.0*durmax)
        solve!(sv, durationmax=durmax)
        nb_calls = sv.bestsol.solver.nb_calls
    else
        # On teste le TimingSolver directement
        # Création de la solution triée sur le target
        bestsol = Solution(inst)
        initial_sort!(bestsol, presort = :target)
        while ( ms()-ms_start  < durmax)
            # On fait un modif minimaliste pour rester faisable
            # swap!(bestsol, 1, 2, do_update=true) # fait le même travail !
            bestsol.planes[1], bestsol.planes[2] = bestsol.planes[2], bestsol.planes[1]
            solve!(bestsol)
        end
        nb_calls = bestsol.solver.nb_calls
    end
    ms_stop = ms()
    Log.poplevel!()

    # Calcul du nombre d'appel par seconde (nb_call_per_sec)
    cps = round(nb_calls/(ms_stop-ms_start), digits=3)
    res = Result(inst.name, fullalgo, cps)

    verbose && println(to_s(res))
    return res
end

function get_stp_perfs(insts::Array{Instance}, fullalgos::Vector, durmax::Int,  use_descent)
    println("\nget_stp_perfs avec le paramètres suivants")
    println("durmax=$(durmax)")
    println("use_descent=$(use_descent)\n")
    results = Vector{Result}()
    for  inst in insts, fullalgo in fullalgos
        # res = get_stp_perf(inst, fullalgo, durmax)
        res = get_timing_solver_perf(inst, fullalgo, durmax=durmax, use_descent=use_descent)
        push!(results, res)
    end
    return results
end

function main()
    println("Test de performances de la brique STP pour tout timingsolver START ($(ms()))")

    global args
    show_fields(args, "Valeur de args :")

    start_time = ms()

    # on se débarrasse un fois pour toute du warning de Gurobi accademique
    new_lp_model(solver=:gurobi)

    # Tous les algos de timing pour Alap
    # (le terme "full" indique que le label contient aussi le mom du solver extrerne)
    # fullalgos = ["faye", "dp",
    #              "lp/cplex", "lp/gurobi", "lp/clp", "lp/cbc", "lp/tulip", "lp/glpk",
    #              "lp2/cplex", "lp2/gurobi", "lp2/cbc",
    #              "lp4/cplex", "lp4/gurobi", "lp4/cbc"]

    # Tous les algos de timing pour Seqata
    fullalgos = ["faye", "dp", "lp/cplex", "lp/gurobi", "lp/clp", "lp/cbc",
                 "lp/tulip", "lp/glpk"]

    # Quelques algos de timing pour tester Seqata
    # fullalgos = ["faye", "dp", "lp/cplex"]

    if args.hack_pour_eleves
        # VERRUE SPECIFIQUE POUR TESTER LES PROJETS ELEVES (glpk non multithreads, ...)
        fullalgos = ["lp/gurobi", "lp/cbc"]
    end

    # Les instances à tester
    insts = [i01(), i09(), i13()]

    #####################################################
    # Précompilation au cas ou sysimg non disponible
    tmp_t0 = ms()
    println("\nRécompilation de tous les solvers : BEGIN")
    for fullalgo in fullalgos
        res = get_timing_solver_perf(i01(), fullalgo, 
                durmax=0.001,   # juste pour faire une itération
                use_descent=args.use_descent, 
                verbose=false)
    end
    println("Récompilation de tous les solvers : END en $(round(ms()-tmp_t0, digits=2))sec ")
    # Précompilation fin
    #####################################################


    results = get_stp_perfs(insts, fullalgos, args.durmax, args.use_descent)

    println("\nSynthèse des résultats")
    Base.sort!(results, by=res->(res.instname,res.fullalgo), rev=false)
    println.(to_s.(results))
    println("Test de performances de la brique STP pour tout timingsolver  END ($(ms()))")

    # construction du dictionnaire des résultats
    resdict = Dict{Tuple, Float64}()
    for res in results
        resdict[(res.instname, res.fullalgo)] = res.cps
    end

    #####################################################
    # Ancien affichage du tableau horizontal
    println("\nTableau des résultats")
    print(" "^4)
    for fullalgo in fullalgos
        print(lpad(fullalgo, 10))
    end
    println()
    for inst in insts
        print(shortname(inst), ": ")  # e.g. "09: "
        for fullalgo in fullalgos
            strval = @sprintf("%10.3f", resdict[inst.name,fullalgo])
            print(strval)
            # print(lpad(repr(resdict[inst.name,fullalgo]), 10))
        end
        println()
    end
    stop_time = ms()
    println("Durée totale des tests : $(round(stop_time-start_time, digits=2))sec")

    #####################################################
    # Nouvel affichage du tableau vertical
    println("\nTableau des résultats")
    println()
    print(" "^15)
    for inst in insts
        print(lpad(shortname(inst), 10))
    end
    println()

    for fullalgo in fullalgos
        print(rpad(fullalgo, 15)) # e.g "lp3/gurobi    "
        for inst in insts
            strval = @sprintf("%10.3f", resdict[inst.name,fullalgo])
            print(strval)
        end
        println()
    end
    stop_time = ms()
    println("Durée totale des tests : $(round(stop_time-start_time, digits=2))sec")
end
main()
