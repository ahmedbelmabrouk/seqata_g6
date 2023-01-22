# Quelques fonctions ou macros utiles en mode interactif
# (dont des alias courts pour les types ou les méthodes)
#
# TODO :
# - faire les macros @p (pour imprimer) et @s (pour to_string générique)

export get_files_from_pats, ag, lv
export @i, @f, p, pn
export p01, p02, p03, p04, p05, p06, p07, p08, p09, p10, p11, p12, p13
export i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12, i13

using Glob
using Debugger # chargement du debugger (voir doc)
using Revise   # bon fonctionnement à vérifier ?
using Test

# Raccourcis pour la gestion des options (Args)
ag(sym::Symbol, val) = Args.set(sym, val)
ag(sym::Symbol) = Args.get(sym)
ag() = Args.show_args()

# Raccourcis pour la gestion du level (pour les messages d'erreur)
lv(L) = Log.level(L)
lv() = Log.level()

p = println
pn = print

data = dirname(dirname(@__FILE__)) * "/data"

# Les chemins d'instance (String)
p01 = "$data/01.alp"
p02 = "$data/02.alp"
p03 = "$data/03.alp"
p04 = "$data/04.alp"
p05 = "$data/05.alp"
p06 = "$data/06.alp"
p07 = "$data/07.alp"
p08 = "$data/08.alp"
p09 = "$data/09.alp"
p10 = "$data/10.alp"
p11 = "$data/11.alp"
p12 = "$data/12.alp"
p13 = "$data/13.alp"

# Définir une méthode i01() au lieu d'une variable i01 permet d'éviter de 
# rallonger de temps de chargement initial.
# i01 = Instance(p01)  # création d'une variable
i01() = Instance(p01)  # création d'une méthode
i02() = Instance(p02)
i03() = Instance(p03)
i04() = Instance(p04)
i05() = Instance(p05)
i06() = Instance(p06)
i07() = Instance(p07)
i08() = Instance(p08)
i09() = Instance(p09)
i10() = Instance(p10)
i11() = Instance(p11)
i12() = Instance(p12)
i13() = Instance(p13)


insts() = [i01(),i02(),i03(),i04(),i05(),i06(),i07(),i08(),i09(),i10(),i11(),i12(),i13()]


# Macro @i pour simplifier les "include" en interactif
# Exemple d'itilisation :
#    @i "04" "05"
# Recherche tous les fichiers correspondant à la gpatternes *04* et *05*
# puis les recharges par include (dans l'ordre)
# Recherche d'abord dans le répertoire racine, puis dans src/ puis dans test/
#
# ASTUCE :
# Si les arguments sont purement alphanumériques, alors on peut éviter des
# guillemets
#   @i inst       ok
#   => charge tous les fichiers contenant "inst" dans src/ ou dans test/
#   @i  src/inst   KO
#   @i "src/inst"  ok
#   @i  04         KO (car est équivalent à @i "4" qui couvre plus large)
#   @i "04"        ok
#
macro i(pats...)
    files = get_files_from_pats(pats...)

    if length(files) == 0
        println("Aucun fichier ne correspond.")
    else
        for absfile in files
            relfiles = replace(absfile, "$APPDIR/" => "")
            println("### include file: \"", relfiles, "\" ...")
            include(absfile)
            println("### include file: \"", relfiles, "\" FAIT")
            # println("include done.")
        end
    end
end

# macro f (find)
# Affiche les fichiers correspondants aux patternes passées (avec date de modif)
macro f(pats...)
    files = get_files_from_pats(pats...)
    wd = pwd()
    cd(APPDIR)
    # On rend le chemin des fichiers relatifs par rapport au projet pour
    # alléger le listing (quelque soit le répertoire courant)
    files = replace.(files, "$APPDIR/" => "")
    for file in files
        run(`ls -l $file`)
    end
    cd(wd)
end

function get_files_from_pats(pats...)
    # @show typeof(pats)
    # @show pats
    files = []
    for pat in pats
        append!(files, glob("$pat", "$APPDIR"))
        append!(files, glob("*$pat*", "$APPDIR"))
        append!(files, glob("*$pat*", "$APPDIR/src"))
        append!(files, glob("*$pat*", "$APPDIR/test"))
    end
    return files
end
#./
