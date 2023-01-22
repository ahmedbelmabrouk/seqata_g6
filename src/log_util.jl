"""
    Log
Module simplifiant l'écriture de logging d'une application. 
Permet d'afficher ou non un message selon le niveau de verbosité
défini par la globale Log.LEVEL
"""
module Log

using DataStructures # pour Stack (push et pop)

export lg
export lg0, lg1, lg2, lg3, lg4, lg5
export ln0, ln1, ln2, ln3, ln4, ln5
# REMARQUE : on n'exporte pas les méthodes suivantes pour obliger à qualifier
# leur appel (e.g. Log.level() )
# export level  # non exporté pour imposer de taper Log.level()
# export pushlevel!, poplevel!   # non exporté

# cette pile doit contenir au moins une valeur qui est le level courant
global LEVELS
if !@isdefined(LEVELS)
    LEVELS = Stack{Int}()
    push!(LEVELS, 2) # VALEUR PAR DÉFAUT
end

# level : modifie ou retourne le level courant (verbosité)
# - sans argument : retourne le level courant
# - avec argument : modifie le level courant et retourne sa nouvelle valeur
# 
function level(newlevel::Int)
    global LEVELS
    # LEVELS[1] = newlevel # non permis
    pop!(LEVELS)
    push!(LEVELS, newlevel)
    first(LEVELS)
end
function level()
    first(LEVELS)
end

# Mémorise la version **précédente** du level
# et positionne le level courant à newlevel
function pushlevel!(newlevel::Int)
    global LEVELS
    push!(LEVELS, newlevel)
    level()
end
# restore le level précédemment mémorisé dans la pile
function poplevel!()
    global LEVELS
    if length(LEVELS) <= 1
        # Il faut conserver au moins une valeur dans LEVELS
        error("pushlevel! should be called before poplevel!")
    end
    pop!(LEVELS)
    level()
end

# Quelques fonctions d'aide au debug
#
# La principale fonction est lg (pour log) qui affiche ou non un message selon le
# niveau de verbosité courant avec ou sans saut de ligne selon la valeur du 
# paramètre suffix.
# Dans tous les cas, cette fonction retourne un booléen indiquant si le niveau
# de verbosité est suffisant.
#
# Le message à afficher est formé par la concaténation des arguments
# vals transformés en String (par la méthode string(...)).
# 
# Exemple d'utilisation :
# 
#   ln3("coucou")
#   => affiche  "coucou" + prefix si LEVEL suffisant (ici si >= 3)
#   lg4(".")
#   => affiche un caractère "." sans saut de ligne pour illustrer une progression
#   if lg3()
#     # pré-calcul... # gros pré-calcul juste pour l'affichage qui suit
#     println("Résultat du pré-calcul")
#   end
#   lg3() n'affiche rien par lui-même, mais autorise une série d'instruction si
#         level est suffisant
#   S'écrit également :
#   lg3() && @show mavariable
#
# BUG ET DÉPENDENCES
# Pour l'instant, ces méthodes ne sont utilisables que après avoir appelé
# la méthode Args.parse_commandline() pour accéder à 
# Ceci sera à modifier (création d'un module Log paramétré de l'extérieur)
# 


# lg0(vals...) = lg(0, vals..., suffix="", doflush=true  )
lg0(vals...; kwargs...) = lg(0, vals...; kwargs...)
lg1(vals...; kwargs...) = lg(1, vals...; kwargs...)
lg2(vals...; kwargs...) = lg(2, vals...; kwargs...)
lg3(vals...; kwargs...) = lg(3, vals...; kwargs...)
lg4(vals...; kwargs...) = lg(4, vals...; kwargs...)
lg5(vals...; kwargs...) = lg(5, vals...; kwargs...)

ln0(vals...; kwargs...) = lg(0, vals...; kwargs..., suffix = "\n")
ln1(vals...; kwargs...) = lg(1, vals...; kwargs..., suffix = "\n")
ln2(vals...; kwargs...) = lg(2, vals...; kwargs..., suffix = "\n")
ln3(vals...; kwargs...) = lg(3, vals...; kwargs..., suffix = "\n")
ln4(vals...; kwargs...) = lg(4, vals...; kwargs..., suffix = "\n")
ln5(vals...; kwargs...) = lg(5, vals...; kwargs..., suffix = "\n")

function lg(
    minlevel::Int,
    vals...;
    prefix::String = "",
    suffix::String = "",
    doflush = true,
)
    if level() >= minlevel
        if !isempty(vals)
            print(prefix, join(vals, ""), suffix)
            doflush && flush(stdout)
        end
        return true
    else
        return false
    end
end

end # module

# POUR INFORMATION : Á DÉPLACER
#
# Voir aussi (capture de l'affichage d'une commande) :
#  https://stackoverflow.com/questions/54599148/in-julia-1-0-how-do-i-get-strings-using-redirect-stdout
# 
# Utilisation :
#   redirect_to_files(prefix * ".log", prefix * ".err") do
#       compute(...)
#   end
#   
# function redirect_to_files(dofunc, outfile, errfile)
#     open(outfile, "w") do out
#         open(errfile, "w") do err
#             redirect_stdout(out) do
#                 redirect_stderr(err) do
#                     dofunc()
#                 end
#             end
#         end
#     end
# end
