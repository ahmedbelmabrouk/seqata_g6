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
# Autres options possibles pour julia
#  --depwarn=no  : pour supprimer d'éventuels warnings
#  --startup-file=no : évite de charger le fichier utilisateur
=#

#### DÉBUT DU CODE JULIA
#
# Pour lancer ce programme en mode interactif depuis un terminal :
#
# Exécution par :
#    julia -iL ./bin/run.jl
#    julia -i --color=yes -L  ./bin/run.jl
#
# En mode interactif :
#   - le main() n'est pas appelé
#   - on inclue tous les packages et les includes possibles
#   - l'analyse des arguments est faite
#   - on charge le fichier d'utilitaire dédié au mode interactif (src/interactive.jl)
#


# Le realpath est nécessaire si l'un des répertoires parents est un lien symbolique :
this_appdir = dirname(dirname(realpath(@__FILE__())))
using Pkg
Pkg.activate(this_appdir)
Pkg.instantiate()

# @show @__FILE__
# @show PROGRAM_FILE
# @show basename(PROGRAM_FILE)
# @show !empty(PROGRAM_FILE) && realpath(PROGRAM_FILE)
# @show isinteractive()

# ENV["JULIA_USING_ALL"] = 1  # for loading all package
include("$this_appdir/src/Seqata.jl")
using .Seqata

if basename(@__FILE__) == basename(PROGRAM_FILE)
    # Mode d'appel normal : on exécute le programme "bin/xxx.jl"
    if Args.get(:action) in [:none, :help]
        # Exécution normale, mais refusée car sans aucun paramètre
        print(Args.get_usage()) # new 07/12/2021
        exit(0)
    end
    # tout est prêt pour l'exécution effectif d'une action !
    main()
else
    @assert(isinteractive())

    @ms using Debugger # chargement du debugger (voir doc)
    @ms using Revise   # fonctionnement à vérifier avec le projet BipSolver
    @ms include("$APPDIR/src/interactive.jl")

    # En mode interactif, on pourrait imposer des arguments par défaut
    Args.set(:infile, "$APPDIR/data/01.alp")

    Log.lg1() && Args.show_args() # car Seqata ne réexporte pas les exports du module Log

    using .Log  # pour se passer du préfixe Log dans Log.lg1() pour la suite

    println()
    println("Début de mode interactif ($(ms()))s")
    println()
end

