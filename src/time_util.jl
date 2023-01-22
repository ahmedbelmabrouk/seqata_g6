# Quelques méthodes ou macros utilitaires personnelles
# (indépendantes de tout projet julia)
# Modif le 23/04/2020

export ms, ms_reset, @ms, MS_LOG

# Log or hidde message from the @ms macro . Can be change by the caller.
# (@isdefined MS_LOG) || global MS_LOG = false  # initialize MS_LOG only if not yet defined
global MS_LOG = haskey(ENV, "JULIA_MS_LOG") && 
                occursin(r"^1|true$"i, ENV["JULIA_MS_LOG"])

"""
    ms()
démarre le chrono si nécessaire (la première fois) et retourne la durée
écoulée depuis ce démarrage (à la milli-seconde prêt)
"""
function ms()        
    global MS_START
    if !@isdefined MS_START
        MS_START = time()
    end
    ms = round(Int, 1000 * (time() - MS_START))
    return (ms / 1000)
end

"""
    ms_reset()
Réinitialise la date de référence du chrono à la date courante.
"""
function ms_reset()
    global MS_START = time()
end

"""
    :(@ms)
chronomètre l'exécution d'une commande.

- affiche la durée depuis le lancement du programme (peut-être long en
  mode interactif !)
- affiche le fichier depuis lequel cette macro est appelée
- affiche la commande à exécuter
- exécute cette commande
- affiche la durée d'exécution sur la ligne suivante

# Exemples
Exemple d'affichage en interactif :

    @ms sleep(1)
    =>
    @ms todo 9349.527s (mode repl):sleep(1) ...
    @ms DONE 9350.552s (mode repl):sleep(1) en 1.002s

Exemples depuis en script :

    @ms using JuMP
    =>
    @ms todo 1.876s Seqata_usings.jl:70: using JuMP ... 
    @ms DONE 5.864s Seqata_usings.jl:70: using JuMP en 3.988s
    """
macro ms(cmd)
    cmdstr = string(cmd)
    # info sur le fichier qui a appelé cette macro
    if isinteractive()
        # mode interactif => pas de fichier d'appel
        fileinfo = "(mode repl)"
    else
        # On récupère le nom relatif du fichier qui a appellé cette macro
        # fname = string(__source__.file)
        fname = basename(string(__source__.file))

        # On récupère le numéro de ligne de l'appel
        fline = string(__source__.line)
        fileinfo = "$(fname):$fline"
    end

    quote
        # la notation $(cmd) est réservée dans un environnement quote
        local t0 = time()
        if MS_LOG
            print("@ms todo ", round(ms(), digits = 3), "s ", $fileinfo, ": ", $(cmdstr))
            println(" ... ")
        end
        local val = $(esc(cmd)) # <= EXECUTION DU CALCUL Á CHRONOMÉTRER
        local t1 = time()
        if MS_LOG
            print("@ms DONE ", round(ms(), digits = 3), "s ", $fileinfo, ": ", $(cmdstr))
            println(" en ", round(t1 - t0, digits = 3), "s")
        end
        val
    end
end
# On démarre le chronomètre
ms()
