#!/usr/bin/env julia

# Ce script doit être spécifique à chaque application
# Il génère une image précompilée avec les modules de l'application
# (dont le module de l'application elle-même) dans le but d'accélerer 
# le lancement.
# 
# Exemple d'utilisation de ce script :
# 
# export SEQATA_SYSIMG=/tmp/julia_sysimg_seqata.so
# ./bin/build_sysimg.jl
#
# => DURÉE D'EXÉCUTION : 461.528sec
#    301585976 Dec 26 11:28 /tmp/julia_sysimg_seqata.so
#
# Utilisation de l'image précompilée
# L'exécutable `./bin/run.jl`` de ce projet utilise l'image précompilé si 
# celle-ci est spécifiés par la variable d'environnement `SEQATA_SYSIMG`
# Par contre pour l'utiliser en mode interactif, il fait la préciser sur 
# la ligne de commande par :
# julia -J$SEQATA_SYSIMG -i ./bin/run.jl \
#
# 21/01/2022 EN TEST
# De façon à permettre la précompilation d'un projet seqata d'élève ne contenant
# pas ce script, le répertoire du projet à précompiler peut être passé en paramètre
# sur la ligne de commande.
# En l'abscence de ce paramètre, c'est le courant projet séqata qui est preecompilé.
# De plus, le script contenant les appels de méthode à compiler (en plus des 
# modules eux même) doivent être dans le projet contenant ce script.

start_time = time()
@show ARGS
appdir = dirname(dirname(realpath(@__FILE__())))
if length(ARGS)==1 && isdir(ARGS[1])
  projdir = realpath(ARGS[1])
  @show ARGS
  @show typeof(ARGS)
  popfirst!(ARGS) # non plantage car les ARGS sont passés à seqata qui ne reconnait pas l'action
else
  projdir = appdir
end
println("Compilation du projet pour le répertoire suivant :")
println("projdir  $(projdir)")
# println("ON ARRET POUR TEST")
# exit(0)

# on veux "myproject" à partir de /xxx/MyProject.jl/"
projname = lowercase(splitext(basename(projdir))[1])

# On devine un nom d'image d'après le nom du répertoire à précompiler
# Seqata.jl => /tmp/julia_sysimg_seqata.so
global SYSIMG = haskey(ENV, "SEQATA_SYSIMG") ? 
                       ENV["SEQATA_SYSIMG"] : 
                       "/tmp/julia_sysimg_$projname.so"                           

println("On va re-créer :")
rm(SYSIMG, force=true)
println("SYSIMG=$SYSIMG")

using PackageCompiler
using Pkg
Pkg.activate(projdir)
Pkg.instantiate()


@show projdir

# On s'arrange pour le l'éxécutable charge le marximum de package possible
ENV["USING_ALL"]    = true
# On veut afficher la durée de chargement de chaque package
ENV["JULIA_MS_LOG"] = true


include("$projdir/src/Seqata.jl")
using .Seqata

# En principe sans spécifier de package (premier paramètre), tous les module du
# projet courant devraient être précompilés mais ce n'est pas le cas (e.g JuMP).
# Mais cela fonctionne en ajoutant le singleton :Seqata.
#  
# create_sysimage(["Plots"], sysimage_path="sys_plots.so", 
#                 precompile_execution_file="precompile_plots.jl")
#
# @time create_sysimage(:Seqata,
#                       project=projdir, # répertoire contenant le Manifest.toml
#                       sysimage_path=SYSIMG, 
#                       precompile_execution_file="$projdir/test/runtests.jl"
# )
precompile_execution_file = "$appdir/test/precompile_script.jl"
println("Fichier contenant les commandes julia à précompiler : " )
println("    precompile_execution_file = $(precompile_execution_file)")
@time create_sysimage(
    :Seqata,
    project=projdir, # répertoire contenant le Manifest.toml
    sysimage_path=SYSIMG, 
    precompile_execution_file=precompile_execution_file
)

stop_time = time()
dur = round(stop_time-start_time, digits=3)



println("FIN D'EXECUTION DU FICHIER " * @__FILE__)
println("fichier sysimg créé : ")
run(Cmd(`ls -l "$SYSIMG"`))
println("Vous pouvez utiliser l'image créée par")
println("   export SEQATA_SYSIMG=$(SYSIMG)")
println("   julia -J$(SYSIMG) ")
println("   time julia -J$(SYSIMG)  ./bin/run.jl des data/09.alp --dur 30")
println("DURÉE D'EXÉCUTION : $(dur)sec")
# test le R/09/11/2022 sur salle : compilation en 981s !!
# exit()
