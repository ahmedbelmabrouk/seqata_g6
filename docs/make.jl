#!/usr/bin/env julia --color=yes

# Utilisation du système de documentation automatique de Julia
#
# Doc du système dans Julia (présentation, principe, ...) :
# - https://docs.julialang.org/en/v1/manual/documentation/ EXCELLENT
# - https://docs.julialang.org/en/v1/stdlib/Markdown
#
# Doc du package (plus détaillé)
# - https://juliadocs.github.io/Documenter.jl/stable/
# - https://juliadocs.github.io/Documenter.jl/stable/man/examples/
# - https://juliadocs.github.io/Documenter.jl/stable/man/guide/#Package-Guide

using Documenter
using Dates
using Pkg
global const appdir = dirname(dirname(realpath(@__FILE__())))

Pkg.activate("$appdir/docs")
Pkg.instantiate()


push!(LOAD_PATH,"$appdir/src/")
@show LOAD_PATH

ENV["JULIA_MS_LOG"] = 1
ENV["JULIA_USING_ALL"] = 1

using Seqata

# makedocs(
#     root    = "<current-directory>",
#     source  = "src",
#     build   = "build",
#     clean   = true,
#     doctest = true,
#     modules = Module[],
#     repo    = "",
#     highlightsig = true,
#     sitename = "",
#     expandfirst = [],
#     page = [],
#     hide("page1.md"),
#     hide("Title" => "page2.md"),
# )


now_str = Dates.format(now(), "dd/MM/yyyy HHhmm") # NON UTIISÉ POUR L'INSTANT !

# The DOCSARGS environment variable can be used to pass additional arguments
# to make.jl.
# This is useful on CI, if you need to change the behavior of the build slightly
# but you can not change the .travis.yml or make.jl scripts any more
# (e.g. for a tag build).
if haskey(ENV, "DOCSARGS")
    for arg in split(ENV["DOCSARGS"])
        (arg in ARGS) || push!(ARGS, arg)
    end
end

# Voir liste des options de makedocs en :
#  https://juliadocs.github.io/Documenter.jl/stable/lib/public/#Documenter.makedocs
#
# On peut configurer la génération de doc pour une consultation avec
# des "prettyurl" ("xxx/mypage" au lieu de "xxx/mypage.html"
# Cependant cela ne fonctionne pas bien pour une consultation locale avec un
# navigateur.
# Du coup je choisis l'option "prettyurls = false" mais avec possibilité
# de la changer si l'on positionne  prettyurl dans l'environnement au lancement
# du programme make.jl.
#      DOCSARGS=prettyurls  julia docs/make.jl
#
@show ("prettyurls" in ARGS)  # ; exit(1)
using Dates
stamp = Dates.format(Dates.now(), "dd/mm/yyyy HH:MM")
makedocs(
    sitename = "Seqata.jl",
    format = Documenter.HTML(
        prettyurls = ("prettyurls" in ARGS),
        # prettyurls = true
        ansicolor = true,
    ),
    # format = Documenter.HTML(prettyurls = false), # pour consultation html locale
    # format = Documenter.LaTeX(), # Mais plante : à creuser...
    # format = Documenter.LaTeX(platform = "docker"),
    # format = Documenter.LaTeX(platform = "none"), # Génère tex sans le compiler
    doctest  = true, # true par défaut
    clean = true, # true par défaut
    repo     = "https://plmlab.math.cnrs.fr/diam/seqata.jl",  # github. gitlab, plmlab, ...
    modules  = [Seqata, Log, Args],
    highlightsig = true,
    expandfirst = [],
    pages    = [
    "Le prototype Seqata" => Any[
       "Index"                          => "index.md",
       "Utilisation du proto"           => "01_presentation_proto_seqata.md",
       "Organisation du code"           => "02_organisation_code.md",
       "Utilisation en mode interactif" => "03_mode_interactif.md",
       "Voisinages explicites"          => "04_voisinages_explicites.md",
    ],
    "FAQ" => Any[
        "Liens utiles"           => "faq/faq_links.md",
        "Installation"           => "faq/faq_install.md",
        "Astuces Julia"          => "faq/faq_julia.md",
        "Méthodologie"           => "faq/faq_method.md",
        "Utilisation de Seqata"  => "faq/faq_seqata.md",
    ],

    # "API" => Any[
    #    "module Seqata"      => "api/api_seqata.md",
    #    "module Seqata.Log"  => "api/api_log.md",
    #    "module Seqata.Args" => "api/api_args.md",
    # ],

    # "ADMIN/index" => Any["ADMIN/index" => "ADMIN/index.md"],
    # "ADMIN" => Any[
    #    "ADMIN/HIST"           => "ADMIN/HIST.md",
    #    "ADMIN/CHANGE"         => "ADMIN/CHANGE.md",
    #    "ADMIN/TODO"           => "ADMIN/TODO.md",
    #    "Passage de Seqata au proto" => "ADMIN/PROTO_CREATION_FROM_SEQATA.md",
    #    "Idées glouton"        => "ADMIN/IDEES_GLOUTON.md",
    #    "Idées ANT"            => "ADMIN/IDEES_ANT.md",
    #    "z_brouillon"          => "ADMIN/z_brouillon.md",
    # ],
    "."   => [], 
    "Compiled at $(stamp)"   => [], 

  ],
)
# deploydocs(
#     repo = "plmlab.math.cnrs.fr/diam/Seqata.jl.git",
#     target = "build",
#     push_preview = true,
# )
