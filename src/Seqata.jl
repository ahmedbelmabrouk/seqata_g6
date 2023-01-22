"""
    Seqata
Package Julia pour résoudre le problème SEQATA
(SÉQuencement d'ATterrissage d'Avions) qui est une version adaptée
du problème de la litérature ALP (Aircraft Landing Problem).
"""
module Seqata

export main, APPDIR, Args, Log

# Le realpath est nécessaire si l'un des répertoires parents est un lien
# symbolique :
global const APPDIR = dirname(dirname(realpath(@__FILE__())))


# On peut contrôler l'affichage de la macro @ms par la variable d'environnement
# JULIA_MS_LOG
# ENV["JULIA_MS_LOG"] = 1
include("time_util.jl") # CHRONOMÉTRAGE POSSIBLE LE PLUS TOT POSSIBLE

# ==================================================================
# CHARGEMENT DES FICHIERS ET SOUS-MODULES DE CETTE APPLICATION
# ==================================================================

@ms include("log_util.jl") # needed by Args module
using .Log  # for using ln1() instead of Log.ln1()

@ms include("args.jl") # define Args.get() for CLI arguments access, ...

@ms Args.parse_commandline(ARGS)
# From now we can access to programme arguments with Args.get(...)

#
# La variable USING_ALL impose le chargement de tous les packages et tous les
# fichiers julia indépendemment de l'action précisée sur la ligne de commande.
# Elle peut être positionnée par le code appelant, via la variable d'environnement
# JULIA_USING_ALL et est typiquement utilisée dans en mode interactif ou pour 
# les tests unitaires.
# Sa valeur par défaut est false.
global USING_ALL = haskey(ENV, "JULIA_USING_ALL") && 
                   occursin(r"^1|true$"i, ENV["JULIA_USING_ALL"])

# En mode interactif, on veut charger tous les packages potentiellement utiles
# USING_ALL ||= isinteractive() # ||= plante mais |= serait ok
USING_ALL = USING_ALL || isinteractive()
@show USING_ALL
MS_LOG && ln1("===== USINGS BEGIN begin($(ms()))")

@ms using Printf: @printf, @sprintf
@ms using Random # shuffle!, ...
# @ms using Dates: now, format # PLANTE !?
@ms import Dates
@ms import Glob

@ms using Crayons     # pour l'affichage en couleur ANSI
@ms using Crayons.Box # pour déclaration des constantes de couleurs RED_FG,...
@ms using Statistics  # pour mean(...) ou median(..)

@ms include("console_util.jl")
@ms include("plane.jl")
@ms include("instance.jl")
@ms include("instance_generators.jl")
@ms include("instance_read_alp.jl")
@ms include("array_util.jl")
@ms include("file_util.jl")
@ms include("solution.jl")
@ms include("solution_readsol.jl")
@ms include("mutation.jl")
    

# @ms include("$APPDIR/src/Seqata_usings.jl") # make use of Args
@ms include("Seqata_usings.jl") # make use of Args

end # module
