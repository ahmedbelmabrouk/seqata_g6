export to_sc, prt_c, prn_c, get_unicode, wait_for_key

# Return une chaine composée du caractère unicode spécifié par un symbol
# (perso et non standard)
function get_unicode(sym::Symbol)

    SYMS = Dict(
        # ligne brisée joignant le nord et l'est passant par le centre
        :NE   => "\U2517", # Nord East
        :NS   => "\U2503",
        :NW   => "\U251B",
        :ES   => "\U250F",
        :EW   => "\U2501",
        :SW   => "\U2513", # South West

        :PP   => "\U254B", # croix totale style plus "+" bold

        :ne   => "\U2514",
        :ns   => "\U2502", # barre verticale light
        :nw   => "\U2518",
        :es   => "\U250C",
        :ew   => "\U2500", # barre horizontale light
        :sw   => "\U2510",
        :pp   => "\U253C", # croix totale style plus "+" light

        :tn   => "\U2534", # Tee vers N
        :te   => "\U251C", # Tee vers E
        :ts   => "\U252C", # Tee vers S
        :tw   => "\U2524", # Tee vers W

        :cp   => "\U00B7", # point centré
        :xx   => "\U2573", # croix totale X
    )
    if haskey(SYMS, sym)
        return SYMS[sym]
    else
        error("get_unicode: symbol UNICODE inconnu : $(sym).")
    end
end

# wait_for_key : attend et retourne un seul caractère
# source: https://discourse.julialang.org/t/wait-for-a-keypress/20218/6
#
function wait_for_key(io = stdin)::Char
    setraw!(raw) = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), io.handle, raw)
    setraw!(true)
    res = read(io, Char)
    setraw!(false)
    return res
end

# GESTION DES COULEURS : NE PLUS UTILISER POUR LES COULEURS :
# voir plutôt (ou réécrire en se basant sur) Crayons.jl
# - https://github.com/KristofferC/Crayons.jl
# - http://ascii-table.com/ansi-escape-sequences.php
# - http://misc.flogisoft.com/bash/tip_colors_and_formatting
#
# Voir aussi ma commande bash : diam_color.sh
#
# Exemple d'utilisation :
#   println(to_sc("== TEST to_sc pour BLUE_STABILO ==", :BLUE_STABILO))
#
# Sur les versions récentes de julia on peut faire par exemple :
#   printstyled("Solution invalide...", color=:red)
# Cependant :
# - julia doit être lancé avec l'option --color=yes
# - cette commande ne permet pas de créer un chaine pour l'afficher plus tard.
#
function to_sc(str::AbstractString, color::Symbol=:STABILO)
    red="9";    red_l="210";  red_d="1"
    grn="2" ;   grn_l="46";   gre_d="22"
    blu="12" ;  grn_l="14" ;  blu_d="4"
    yel="226";  yel_l="228";  yel_d="3"
    ora="208";  ora_l="215";  ora_d="130"
    colors = Dict(
        :STABILO         => "\e1\e[38;5;$(red)m\e[48;5;$(yel)m",
        :RED_STABILO     => "\e1\e[38;5;$(red)m\e[48;5;$(yel)m",
        :GREEN_STABILO   => "\e1\e[38;5;$(grn)m\e[48;5;$(yel)m",
        :BLUE_STABILO    => "\e1\e[38;5;$(blu)m\e[48;5;$(yel)m",
        :ORANGE_STABILO  => "\e1\e[38;5;$(ora_d)m\e[48;5;$(yel)m",
        :RED             => "\e1\e[38;5;$(red)m",
        :GREEN           => "\e1\e[38;5;$(grn)m",
        :BLUE            => "\e1\e[38;5;$(blu)m",
        :YELLOW          => "\e1\e[38;5;$(yel)m",
        :ORANGE          => "\e1\e[38;5;$(ora_d)m",
    )
    if haskey(colors, color)
        col = colors[color]
    else
        error("to_s_color: couleur inconnue : $(color).")
    end
    def = "\e[0m"
    buf = IOBuffer()
    print(buf, "$(col)", str, "$(def)")
    # return takebuf_string(buf)
    return String(take!(buf))
end

function prt_c(str::AbstractString, color::Symbol = :STABILO)
    print(to_sc(str, color))
end
function prn_c(str::AbstractString, color::Symbol = :STABILO)
    println(to_sc(str, color))
end

