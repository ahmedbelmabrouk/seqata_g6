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

# AUTRE EXEMPLE DE RÉPERTOIRE DE TEST (cf QML.jl)
#   https://github.com/barche/QML.jl/blob/master/test/runtests.jl
#

# Le realpath est nécessaire si l'un des répertoires parents est un lien symbolique :
this_appdir = dirname(dirname(realpath(@__FILE__())))
using Pkg
Pkg.activate(dirname(dirname(realpath(@__FILE__()))))
Pkg.instantiate()

ENV["JULIA_USING_ALL"] = 1  # for loading all package

pushfirst!(ARGS, "none") # On précise l'action par défaut
include("$this_appdir/src/Seqata.jl")
using .Seqata
using .Log  # pour pouvoir taper ln1() au lieu de Log.ln1()

println("BEGIN all tests at ", ms(), "sec")

using Test
using Random
using Printf
using Dates

@ms include("libtest_model_util.jl")

FILES = copy(Args.get(:files))
println("## TEST $(basename(@__FILE__)) ")

# le Log.level par défaut est 2 pour une utilisation normale de l'appli.
# En mode test, on veut qu'il soit de 1.

# On mémorise l'état initial pour le restaurer avant chaque test
TEST_LEVEL = Log.level() - 1
TEST_ARGS = copy(Args.args)
# TEST_ARGS = deepcopy(Args.args)

println("ÉTAT DE RESTAURATION DE CHAQUE TEST :")
println("ARGS_TEST :")
Args.show_args()
println("TEST_LEVEL : $TEST_LEVEL")

# Phase 1 : construction de la liste explicite du/des fichiers à tester
# - soit le/les fichier(s) passés en paramètre et aucun autre
# - soit tous les fichiers :
#   - commençant par test-*
#   - sans les fichiers de test de performance *-perf.jl (car lents !)
#   - et sans un liste de fichier explicitement à ignorer
#
global files = Vector{String}()
if !isempty(FILES)
    for file in FILES
        isfile(file) || continue
        push!(files, basename(file))
    end
else
    # Calcul de la liste des tests par défaut
    # EXEMPLE (DE ALAP) DE FICHIER A NE PAS INCLURE E.G. CAR TROP LONGGG.
    # CEUX-CI PEUVENT ÊTRE TESTÉS INDIVIDUELLEMENT EN LES PASSANT EN PARAMÈTRE.
    # Pour plus de souplesse dans les noms, ces fichiers sont définis par une
    # expression régulière.
    excluded_rpats = [
        # test-12b-dmip-cbc.jl
        r"test-.*-dmip-cbc.jl",  # long (1.5mn) (772s sur salle)
        r"test-.*-dmip-glpk.jl", # trés trés long (2647s sur salle) !!
        r"test-.*-perf.jl",      # les tests de performances sont trop long
    ]

    old_pwd = pwd()
    cd(dirname(@__FILE__))
    allfiles = readdir()
    cd(old_pwd)

    for file in allfiles
        print("\nfile=$file ?...")

        # On exclue les fichiers ne commençant pas par test-
        !startswith(file, "test-") && continue

        # On exclue le fichier s'il correspond à ceux explicitement indiqués
        if any(rpat -> occursin(rpat, file), excluded_rpats)
            print(" => (SKIPPED)")
            continue
        end

        print("  => OK")
        push!(files, file)
    end
    println()
end

println("Liste des $(length(files)) tests à effectuer :")
for file in files
    println("  $file")
end

#
# Phase 2 : exécution de chaque test dans un contexte indépendant
# - on se déplace éventuellement dans son sous-répertoire
# - on capture une erreur éventuelle pour ne pas arreter les autres tests
# - on chronomètre chaque test
#
@testset "Tests pour projet SEQATA" begin
    for file in files
        Log.level(TEST_LEVEL)
        copy!(Args.args, TEST_ARGS)

        t0 = ms()
        println("\n" * "="^80)
        println("====== Test du fichier $(file)...")
        old_pwd = pwd()
        cd(dirname(@__FILE__))
        try
            @testset "Test $(file)" begin
                # println("include(", abspath(file), ")")
                include(abspath(file))
            end
        catch err
            println()
            # println(err)
            @warn err
            # rethrow(err)
        end
        cd(dirname(old_pwd))
        dt = round(ms() - t0, digits = 3)
        println("====== Test du fichier $(file) fait en $(dt)s")
    end
    println("Fin des tests à ", ms(), "s")
end
# La suite ne serait pas exécutée en cas d'erreur
println("END all tests at ", ms(), "sec")
