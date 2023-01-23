export AbstractMutation, Mutation, normalize, width, length
export merge_mutation
export generate_mutations_swap
export generate_mutations_shift, generate_mutations_shift2, generate_mutations_shift2g
export generate_mutations_permu
export mutate!, generate_neighborhood, generate_nbh

using Combinatorics

#noexport Main.reverse


# Représente un mouvement abstrait
# On pourra mémoriser un vecteur de de struc de ce type pour représenter un
# voisinage complet.
abstract type AbstractMutation
end


# Mutation: mutation définie par une permutation (**non mutable**)
#
# Les attributs sont les suivants
# - class:  sympole décrivant la classe (ou le type) de la mutation
#   Elle définit la sémantique des indices 
# - indices1: vecteur d'indice de départ
# - indices2: vecteur d'indice de la permuation

# Ce type est non mutable pour une raison d'efficacité (pour parallélisation 
# éventuelle).
# Le constructeur par défaut Mutation(class,indices1,indices2) est de
# bas niveau. Il n'effectuer aucune vérification ni normalisation.
# Le constructeur suivant avec paramètres pour fonctionnalité de plus haut niveau
# (dont normalisation, ...)
# 
# Idées de classes de mutations possibles (candidate à creuser) :
# - :shift : déplacement i1->i2 puis i3->i4
# - :swap : swap i1<->i2 puis i3<->i4
# - :permu : TODO 
# - :reversebranch : (TODO pour pb de routage TSP, velib, ...)
# - :indices : pour extension futures (mutations spéciales, ...)
#   La sémantique et l'utilisation des indices est à la chage de l'appelant.
# 
struct Mutation <: AbstractMutation
    class::Symbol
    indices1::Vector{Int}
    indices2::Vector{Int}
    # indices1::Tuple{Vararg{Int}} # iiiiii
    # indices2::Tuple{Vararg{Int}} # iiiiii
end

# Mutation: ce constructeur normalise la mutation crée
# e.g       P[5,4,3,2,1]->[5,2,3,4,1]
# devient : P[2,4]->[4,2]
#
# Le paramètre do_normalize permet de gagner du temps dans le cas où la
# normaisation est déjà faite par contruction d'un voisinage.
# Mais attention : cette normalisation est indispensable pour le bon 
# fonctionnement des méthodes width et length !
# 
function Mutation(indices1::Vector{Int64}, indices2::Vector{Int64};
                        class=:permu, do_normalize=true)
    if do_normalize
        # On supprime les éléments fixes des deux vecteurs
        indices = findall(!=(0), indices1 - indices2)
        reduc_indices1 = indices1[indices]
        reduc_indices2 = indices2[indices]

        # On trie les deux vecteurs selon le premier vecteur
        perm = sortperm(reduc_indices1)
        reduc_indices1 .= reduc_indices1[perm] # affectation "in place"
        reduc_indices2 .= reduc_indices2[perm] # affectation "in place"
        return Mutation(class, reduc_indices1, reduc_indices2) 
    else
        return Mutation(class, indices1, indices2) 
    end
end

# Contruction d'une Mutation à partir de deux ou quatre indices
#
# Ces indices définissent un ou deux mouvements soit de la classe swap 
# soit de la classe shift (mais pas de mixtes pour l'instant).
# 
# Principe de la conversion de deux paires d'indices en Mutation
# 
# - Exemple pour deux mouvements de classe :shift
#   - Soit mut = (i1,i2,i3,i4) = (4,6,10,8)
#     (on veut un shift : i1->i2  et i3->i4)
#     On crée la Mutation équivalente suivante :
#     => calculer lb = min(mut) et ub = max(mut)
#     => indices1 = [4, 5, 6, 7, 8, 9, 10]  e.g. par collect(lb:ub)
#     => indices2 = [5, 6, 4, 7, 10, 8, 9]
#     Les indices lb et ub de indices1 sont forcément déplacés (font partie 
#     de mut) donc la permu indices1 sera de **largeur minimale**
#   - calcul de indices1  :
#     - collect(lb:ub) => indices1 = [4, 5, 6, 7, 8, 9, 10]
#   - calcul de indices2 : 
#     - créer mut2 = [mut.i1, mut.i2, mut.i3, mut.i4] .- (lb+1)
#       => mut = (4,6,10,8) => mut2 = (1,3,7,5)
#     - appliquer mut2=(1,3,7,5) sur indice1
#       => mut2([4, 5, 6, 7, 8, 9, 10]) ==> indice2 = [5, 6, 4, 7, 10, 8, 9]
#   - un fois notre permu réduite sur le bord, on peut encore la réduire au milieu
#     en supprimant les indices du milieu non modifiés
#     => indices1 = [4, 5, 6, 7, 8, 9, 10] => [4, 5, 6, 8, 9, 10]
#     => indices2 = [5, 6, 4, 7, 10, 8, 9] => [5, 6, 4, 10, 8, 9]
#     Mais ceci sera fait par le constructeur de Mutation
#
# - Exemple pour deux mouvements de classe :swap
#     (on veut un swap : i1<->i2  et i3<->i4)
#   Soit mut = (i1,i2,i3,i4) = (4,6,10,8)
#   On peut directement constuire les tableaux :
#   - indices1 = [i1,i3]
#   - indices2 = [i3,i4]
#
# BUG : tous les indices doivent être différnts !
#
function Mutation(class::Symbol, i1::Int, i2::Int, i3::Int=1, i4::Int=1)
    # TODO TRAITER CAS PARTICULIER OÙ i1==i2 ou i3==i4
    # Par exemple : 
    #    (i1,i2,i3,i4)=(450, 453, 1, 1) 
    # créerait un tableau intermédiaire inutilement grand 
    # => on le transforme en :
    #    (i1,i2,i3,i4)=(450, 452, 450, 450)
    # ATTENTION, dans l'exemple ci-dessus le mouvement (i3,i4) doit resté nul.
    # Il faut donc maintenir l'égalité i3==i4 
    if i1 == i2
        i1 = i2 = i3 # on propage i3
    end
    if i3 == i4
        i3 = i4 = i1 # on propage i1
    end

    lb = min(i1, i2, i3, i4)
    ub = max(i1, i2, i3, i4)
    indices1 = collect(lb:ub)

    if class == :shift
        # On convertit les indices pour ce petit tableau en commencant à 1
        # On veut :  m = (4,6,10,8) => mut_local = (1,3,7,5)
        i1, i2, i3, i4 = (i1, i2, i3, i4) .- (lb - 1)
        # @show i1,i2,i3,i4
        indices2 = copy(indices1)
        shift!(indices2, i1, i2)
        shift!(indices2, i3, i4)
    elseif class == :swap
        # Ce cas est plus simple 
        # SERAIT MAUVAIS CAR DUPPLICATION POSSIBLE DES INDICES !
        #   indices1 = [i1,i2,i3,i4]
        #   indices2 = [i2,i1,i4,i3]
    
        # On convertit les indices pour ce petit tableau en commencant à 1
        # On veut :  m = (4,6,10,8) => mut_local = (1,3,7,5)
        e1, e2, e3, e4 = (i1, i2, i3, i4) .- (lb - 1)
        indices2 = copy(indices1)
        indices2[[e1,e2]] = indices2[[e2,e1]]
        indices2[[e3,e4]] = indices2[[e4,e3]]
    else
        error("\nclass :$class) not allowed in conversion from Mutation\n"*
              "should be :shift or :swap\n")
    end

    # La normalisation de la mutation est faite par le constructeur
    return Mutation(indices1, indices2, class=class)
end

# Supprime les indices inchangés dans la permutions 
# (ceux dont la valeur est identique dans indices1 et dans indices2) 
# 
# NEW 09/05/2019 : déja fait dans le constructeur de Mutation
function normalize(m::Mutation)
    return Mutation(m.indices1, m.indices2, class=m.class) # HAUT NIVEAU
end

# reverse: très rapide car utilise le construction de Mutation de bas niveau.
function Main.reverse(m::Mutation)
    return Mutation(m.class, m.indices2, m.indices1) # BAS NIVEAU
end

# Retourne la largeur du bloc de la mutation
# 
# N'exploitate pas le fait que m soit normatisée
# Sinon on pourrait faire ? :  m.indices1[end] - m.indices1[1] + 1
# 
function width(m::Mutation)
    maximum(m.indices1) - minimum(m.indices1) + 1
end

# Retourne le nombre d'éléments modifiés par la mutation
# Pour simplifier on passe par le constructeur de Mutation
# 
function Base.length(m::Mutation)
    length(m.indices1)
end

# Constructeur d'une Mutation repésentant un swap de deux éléments
# Redondant mais plus efficace que la version avec 4 indices.
# 
# function Mutation_swap(i1::Int, i2::Int)
#     i_min,i_max = minmax(i1,i2)
#     class = Symbol(string("s", i_max-i_min+1))
#     # Construction directe bas niveau normalisée
#     return Mutation(class, [i_min,i_max], [i_max,i_min])]
# end

# Fusionne un mouvement supplémentaire dans une mutation existante.
# (crée un nouvelle mutation car celle passée en paramètre est non mutable)
# 
# DEMARCHE GENERIQUE 
# 
# 1. créer in tableau fullindices1 extensif (collect(minall:maxall))
# 2. appliquer la mutation mut à fullindices1 => fullindices2
# 3. appliquer le nouveau mouvement à fullindices2 => fullindices3
# 4. créer nouvelle mutation pour fullindices1, fullindices3, class)
# 
function merge_mutation(mut::Mutation, class::Symbol, i1::Int, i2::Int)
    imin,imax = minmax(i1,i2)
    allmin = min(minimum(mut.indices1), imin)
    allmax = max(maximum(mut.indices1), imax)

    # 1. créer in tableau fullindices1 extensif (collect(minall:maxall))
    fullindices1 = collect(allmin:allmax)
    fullindices2 = copy(fullindices1)
    fullindices2[mut.indices1] = fullindices2[mut.indices2]
    # fullindices2 = apply(mut) # TODO ?
    # sol.planes[indices1] = sol.planes[indices2]

    if class == :shift
        shift!(fullindices2, i1, i2) # ordre i1 i2 important !
    elseif class == :swap
        fullindices2[[i1,i2]] = fullindices2[[i2,i1]]
    else
        error("\nclass :$class) not allowed in conversion from Mutation\n"*
              "should be :shift or :swap\n")
    end

    if mut.class != class
        # On change la classe de la nouvelle mutation
        class = :mixte
    end

    return Mutation(fullindices1, fullindices2, class=class)
end

# Fusionne deux mutations pour en créer une nouvelle.
# (crée un nouvelle mutation car celle passée en paramètre est non mutable)
# 
# DEMARCHE GENERIQUE 
# 
# 1. créer in tableau fullindices1 extensif (collect(minall:maxall))
# 2. appliquer les mutations mut1 puis mu2 à fullindices2 => fullindices2
# 3. définir la classe résultante
# 4. créer nouvelle mutation pour fullindices1, fullindices3, class)
# 
function merge_mutation(mut1::Mutation, mut2::Mutation)
    allmin = min(minimum(mut1.indices1), minimum(mut2.indices1))
    allmax = max(maximum(mut1.indices1), maximum(mut2.indices1))

    # 1. créer in tableau fullindices1 extensif (collect(minall:maxall))
    fullindices1 = collect(allmin:allmax)
    fullindices2 = copy(fullindices1)

    # 2. appliquer les mutations mut1 puis mut2 à fullindices2 => fullindices2
    fullindices2[mut1.indices1] = fullindices2[mut1.indices2]
    fullindices2[mut2.indices1] = fullindices2[mut2.indices2]

    # 3. définir la classe résultante
    if mut1.class == mut2.class
        class = mut1.class
    else
        # On change la classe de la nouvelle mutation
        class = :mixte
    end

    return Mutation(fullindices1, fullindices2, class=class)
end

# Méthode Julia pour convertir tout objet en string
function Base.show(io::IO, m::Mutation)
    if m.class == :permu
        # str = string("P(", m.indices1, "->", m.indices2, ")")
        str = string("P[", join(m.indices1,","), "]->[", join(m.indices2,","), "]")
        # str = string("P(", join(indices, ","), ")")
    elseif m.class == :shift
        str = string("T[", join(m.indices1,","), "]->[", join(m.indices2,","), "]")
    elseif m.class == :swap
        str = string("S[", join(m.indices1,","), "]->[", join(m.indices2,","), "]")
    elseif m.class == :mixte
        str = string("X[", join(m.indices1,","), "]->[", join(m.indices2,","), "]")
    else 
        str = "Mutation unknown class + $(m.class)"
    end
    Base.write(io, str)
end

# À redéfinir par pouvoir trier un vector de mutations
function Base.isless(m1::Mutation, m2::Mutation)
    if m1.indices1 < m2.indices1
        return true
    elseif m1.indices1 > m2.indices1 
        return false
    elseif m1.indices2 < m2.indices2
        return true
    elseif m1.indices2 > m2.indices2
        return false
    elseif m1.class < m2.class
        return true
    else
        return false
    end
end
# À redéfinir par pouvoir utiliser unique! sur un vector de mutations
# On ne tient pas compte de l'attribut class pour le test d'égalité car
# il ne sert que de commentaire (e.g. :shift, :swap, :permu, ...)
function Base.isequal(m1::Mutation, m2::Mutation)
    return  Base.hash(m1.indices1) == Base.hash(m2.indices1) &&
            Base.hash(m1.indices2) == Base.hash(m2.indices2)
end

# À définir par pouvoir utiliser unique! sur un vector de mutations
function Base.hash(m::Mutation)
    return Base.hash([m.indices1, m.indices2])
end

# generate_mutations_swap : génère le vecteur de toutes les mutations de type swap
# spécifiées par les paramètres.
# Cette méthode vérifie ses paramètres.
# 
# - n: taille maxi du domaine (e.g taille de l'instance)
# - idx_first, idx_last: indices limites des avions à déplacer dans la solution
# - idx_last = -1 pour indiquer la valeur maximum
# - shift_max : valeur absolue maxi d'un déplacement de i1 vers i2 (relatif)
#               (défaut -1 mour maximum)
# - shift_min=1 : écart minimum d'un swap. 
#   Une valeur shift_min=2 permet d'éviter les swaps s1 qui sont équivalents à 
#   des mutations de type shift de largeur 1 (t1)
# 
function generate_mutations_swap(n::Int;
                            idx_first::Int=1,
                            idx_last::Int=-1,
                            shift_min::Int=1,
                            shift_max::Int=-1,
                            )

    lg4() && println("L4 generate_mutations_swap BEGIN at ", ms() )
    # 
    # VÉRIFICATION DU DOMAINE DES PARAMÈTRES
    # 
    if idx_last == -1    idx_last = n end
    idx_last = min(n, idx_last)
    @assert 1 <= idx_last
    @assert 1 <= idx_first < idx_last

    # shift_max == -1 => on prend le maximum possible
    if shift_max == -1  shift_max = idx_last-idx_first  end
    @assert shift_max >= 1 "shift_max doit être strictement positif (ou nul pour max)"

    @assert 1 <= shift_min "shift_min doit être strictement positif (défaut=1)"

    # 
    # PRÉALLOCATION DU VECTEUR DE MUTATIONS
    #
    # On pré-alloue largement pour un vecteur
    # Évaluation de la taille maxi du nombre de swap i1<->i2.
    nb_i1 = idx_last-idx_first      # e.g 100-1 => 99
    nb_i2 = shift_max               # e.g 2
    N_MAX = nb_i1*nb_i2
    if (lg4())
        print("L4 generate_mutations_swap : prédimensionnement du vecteur : ") 
        @show nb_i1, nb_i2, N_MAX
    end
    muts = Vector{Mutation}()
    sizehint!(muts, N_MAX)
    
    # 
    # GÉNÉRATION DES MUTATIONS
    # 
    # On y ajoute tous les swaps possibles (i1<->i2)
    # ATTENTION : un swap de 1 est identique à un shift de 1. Il faut donc
    # prévoir la suppression de ces doublons (e.g. par shift_min==2)
    for i1 in idx_first:idx_last-1
        for i2 in (i1+shift_min):min(i1+shift_max, idx_last)
            # push!(muts, Mutation4(i1, i2, 0, 0, class=:swap))
            push!(muts, Mutation(:swap, i1, i2))
        end
    end

    # On vérifie que le prédimensionnement du tableau estsuffisant :
    if lg4()
        println("L4 generate_mutations_swap : length(muts)=", length(muts))
    end

    # # sizehint!(muts, 0)
    # # return muts
    # muts2 = Mutation.(muts)
    # return unique!(muts2)

    # size1 = length(muts)
    # unique!(muts)
    # sizehint!(muts, 0)
    # size2 = length(muts)
    # if size1!=size2
    #     ln1("ATTENTION : generate_mutations_swap doublons évités $(size1)>$(size2)")
    # end
    # return muts

    len1 = length(muts)
    muts = collect(Set(muts))
    len2 =  length(muts)

    if length(muts) > N_MAX
        ln1("ATTENTION : generate_mutations_swap: sizehint sous-dimentionné")
        ln1("            sizehint=$N_MAX mais length(muts=$(length(muts))")
    end
    if lg3()
        println("generate_mutations_swap sizehint=$N_MAX llen1=$len1 len2=$len2")
    end
    return muts

end


# generate_mutations_shift : génère le vecteur de toutes les mutations de type shift
# spécifiées par les paramètres.
# Cette méthode vérifie ses paramètres.
# 
# - sol: la solution à muter
# - idx_first, idx_last: indices limites des avions à déplacer dans la solution
# - idx_last = -1 pour indiquer la valeur maximum
# - shift_min = 1 : valeur absolue mini d'un déplacement de i1 vers i2 (relatif)
# - shift_max : valeur absolue maxi d'un déplacement de i1 vers i2 (relatif)
# 
function generate_mutations_shift(n::Int;
                            idx_first::Int=1,
                            idx_last::Int=-1,
                            shift_min::Int=1,
                            shift_max::Int=-1,
                            )

    lg4() && println("L4 generate_mutations_shift BEGIN at ", ms() )

    # 
    # VÉRIFICATION DU DOMAINE DES PARAMÈTRES
    # 
    if idx_last == -1    idx_last = n end
    idx_last = min(n, idx_last)
    @assert 1 <= idx_last
    @assert 1 <= idx_first < idx_last

    # shift_max == -1 => on prend le maximum possible
    if shift_max == -1  shift_max = idx_last-idx_first  end
    @assert shift_max >= 1 "shift_max doit être strictement positif (ou nul pour max)"
    @assert shift_min <= shift_max "shift_min doit être inférieur ou égale à shift_max"

    # 
    # PRÉALLOCATION DU VECTEUR DE MUTATIONS
    #
    # On préalloue largement pour un vecteur
    # evaluation de la taille maxi du nombre de mono_shift i1->i2 et i1->i1.
    # Il faudrait enlever des shifts non faisables pour les i1 très grands, 
    #      e.g.  si i1 == n-1 et i2 == 5
    nb_i1 = idx_last-idx_first      # e.g 100-1 => 99
    nb_i2 = shift_max               # e.g 2
    # 2->3 est possible, mais pas 3->2 car déjà fait par 2->3 ; d'où le -1
    N_MAX = nb_i1*nb_i2 + nb_i1*(nb_i2-1)
    if (lg4())
        print("L4 generate_mutations_shift : prédimensionnement du vecteur : ") 
        @show nb_i1, nb_i2, N_MAX
    end
    muts = Vector{Mutation}()
    sizehint!(muts, N_MAX)
    
    # 
    # GÉNÉRATION DES MUTATIONS
    # 
    # On y ajoute tous les mouvements uniques possibles (i1->i2)
    for i1 in idx_first:idx_last-1
        for i2 in i1+shift_min:min(i1+shift_max, idx_last)
            # @show i1,i2
            # push!(muts, Mutation4(i1, i2, 0, 0, class=:shift))
            push!(muts, Mutation(:shift, i1, i2))
            # Le recule de i2 de 1 est déjà fait par l'avance de i1 de 1
            if i2 != i1+1
                # push!(muts, Mutation4(i2, i1, 0, 0, class=:shift))
                push!(muts, Mutation(:shift, i2, i1))
            end
        end
    end

    # On vérifie que le prédimensionnement du tableau estsuffisant :
    if lg4()
        println("L4 generate_mutations_shift : length(muts)=", length(muts))
    end

    # size1 = length(muts)
    # unique!(muts)
    # sizehint!(muts, 0)
    # size2 = length(muts)
    # if size1!=size2
    #     ln1("ATTENTION : generate_mutations_shift doublons évités $(size1)>$(size2)")
    # end
    # return muts

    len1 = length(muts)
    muts = collect(Set(muts))
    len2 =  length(muts)

    if length(muts) > N_MAX
        ln1("ATTENTION : generate_mutations_shift: sizehint sous-dimentionné")
        ln1("            sizehint=$N_MAX mais length(muts=$(length(muts))")
    end
    if lg3()
        println("generate_mutations_shift sizehint=$N_MAX llen1=$len1 len2=$len2")
    end
    return muts

end

# generate_mutations_shift2g : génère le vecteur de toutes les mutations de type shift
# spécifiées par les paramètres.
# Cette méthode vérifie ses paramètres.
# 
# - sol: la solution à muter
# - idx_first, idx_last: indices limites des avions à déplacer dans la solution
# - idx_last = -1 pour indiquer la valeur maximum
# - shift_max : valeur absolue maxi d'un déplacement de i1 vers i2 (relatif)
# - gap_max : valeur absolue maxi de la différence entre  max(i1,i2) et min(i3,i4)
# 
# Principe de la génération des mutations shift2g
# -----------------------------------------------
#   - créer tous les (i2, i1, i4, i3) tels que i1 < i2 < i3 < i4
#   - pour chaque quadruplet, générer les (maxi) 4 cas possibles : i1<->i2 i3<->i4
#   le Vector des mutations (muts) est pré-alloué pour éviter les réallocations
#   dynamique (par sizehint!) puis réduit pour éviter les cases vides (toujours
#   dynamique (par sizehint!) puis réduit pour éviter les cases vides (toujours
#   par sizehint!) (ce dernier n'est peut-être pas utile pour un Set)
# 
function generate_mutations_shift2g(n::Int;
                            idx_first::Int=1,
                            idx_last::Int=-1,
                            # shift_min::Int=1, # TODO
                            shift_max::Int=-1,
                            # gap_min::Int=1, # TODO
                            gap_max::Int=-1,
                            )

    lg4() && println("L4 generate_mutations_shift2g BEGIN at ", ms() )

    # 
    # VÉRIFICATION DU DOMAINE DES PARAMÈTRES
    # 
    if idx_last == -1    idx_last = n end
    idx_last = min(n, idx_last)
    @assert 1 <= idx_last
    @assert 1 <= idx_first < idx_last

    # shift_max == -1 => on prend le maximum possible 
    if shift_max == -1  shift_max = idx_last-idx_first  end
    @assert shift_max >= 1 "shift_max doit être strictement positif (ou -1 pour max)"

    # gap_max == -1 => on prend le maximum possible 
    if gap_max == -1  gap_max = idx_last-idx_first  end
    @assert gap_max >= 1 "gap_max doit être strictement positif (ou -1 pour max)"

    # 
    # PRÉALLOCATION DU VECTEUR DE MUTATIONS
    #
    # evaluation de la taille maxi du nombre de shift2g (on prévoit large !)
    # 
    nb_i1 = idx_last-idx_first-2          # e.g 97
    nb_i2 = shift_max    # e.g 3
    nb_i3 = gap_max
    nb_i4 = shift_max
    # Pour chaque quadruplet i1..i4, on a au maximum quatre M2
    N_MAX = 4*(nb_i1 * nb_i2 * nb_i3 * nb_i4)
    if (lg4())
        print("L4 generate_mutations_shift2g : prédimensionnement du vecteur : ") 
        @show nb_i1, nb_i2, nb_i3, nb_i4, N_MAX
    end
    muts = Vector{Mutation}()
    sizehint!(muts, N_MAX)

    # 
    # GÉNÉRATION DES MUTATIONS
    # 

    # Les 4 indices de ce mouvement seront tels que : i1 < i2 < i3 < i4
    # On génère alors jusqu'à quatre mouvements possibles :  i1<->i2, i3<->i4
    for i1 in idx_first:idx_last-3
        for i2 in i1+1:min(i1+shift_max, idx_last-2)
            for i3 in i2+1:min(i2+gap_max, idx_last-1)
                for i4 in i3+1:min(i3+shift_max, idx_last)
                    # traiter les quatre cas induits potentiels
                    # do_shift (i1,i2) et (i3,i4);
                    push!(muts, Mutation(:shift, i1, i2, i3, i4))
                    if i2-i1 > 1
                        # do_shift (i2,i1) et (i3,i4);
                        push!(muts, Mutation(:shift, i2, i1, i3, i4))
                        if i4-i3 > 1
                            # do_shift (i2,i1) et (i4,i3);
                            push!(muts, Mutation(:shift, i2, i1, i4, i3))
                        end
                    end
                    if i4-i3 > 1
                        # do_shift (i1,i2) et (i4,i3)
                        push!(muts, Mutation(:shift, i1, i2, i4, i3))
                    end
                end
            end
        end
    end

    # On vérifie que le prédimensionnement du tableau estsuffisant :
    if lg4()
        println("L4 generate_mutations_shift2g : length(muts)=", length(muts))
    end

    # # sizehint!(muts, 0)
    # # return muts
    # muts2 = Mutation.(muts)
    # return unique!(muts2)

    # size1 = length(muts)
    # unique!(muts)
    # sizehint!(muts, 0)
    # size2 = length(muts)
    # if size1!=size2
    #     ln1("ATTENTION : generate_mutations_shift2g doublons évités $(size1)>$(size2)")
    # end
    # return muts

    len1 = length(muts)
    muts = collect(Set(muts))
    len2 =  length(muts)

    if length(muts) > N_MAX
        ln1("ATTENTION : generate_mutations_shift2g: sizehint sous-dimentionné")
        ln1("            sizehint=$N_MAX mais length(muts=$(length(muts))")
    end
    if lg3()
        println("generate_mutations_shift2g sizehint=$N_MAX llen1=$len1 len2=$len2")
    end
    return muts

end

# generate_mutations_shift2 : génère le vecteur de toutes les mutations de type shift
# spécifiées par les paramètres.
#
# Paramètres :
# - sol: la solution à muter
# - idx_first, idx_last: indices limites des avions à déplacer dans la solution
# - idx_last = -1 pour indiquer la valeur maximum
# - shift_max : largeur maxi des indices permutés (valeur maxi d'un déplacement)
# 
# Cacactéristiques de ce générateur :
# - Les deux shifts créés **peuvent se superposer**
# - seul le paramètre shift_max existe (gap_min et gap_max n'existe pas)
#   Ceci évite confusion par erreur avec méthode generate_mutations_shift2g
# - retourne un vecteur de Mutation canonique (minimales et uniques)
# - Cette méthode vérifie ses paramètres.
# - cette mutation (deux shift recouvrables) intègre les shift simples et 
#   les swap simples **à condition que shift_max >= 2 (ne marche pas pour 1)
#   REMARQUE : ON POURRAIT CONTOURNER CETTE RESTRICTION EN AUTORISANT LE
#   SECOND SHIFT Ã ÊTRE NUL DANS LE CAS OÙ shift_max == 1 
#   (Le cas shift_max==1 est dégénéré et rare car il correspond à un swap1)
# 
# Principe de la génération des mutations shift2
# ----------------------------------------------
# - On fait glisser une fenêtre [b1,b2] de largeur totale shift_max+1
#   (car bound b2 = b1+shift_max)
# - on prend tous les couples orientés (i1,i2) et (i3,i4) possibles
#   contenus dans [b1,b2]
# - construire les mutations Mutation4(i1,i2,i3,i4) telles que 
#    max(i1,i2,i3,i4) - min(i1,i2,i3,i4) <= shift_max 
# - éviter les redondances évidentes dès la construction des Mutation4
#   - [5 6 1 2] versus [1 2 5 6] car mutations identiques
#   - [1 3 3 1] car opération nulle
# - Attention, une même mutation4 (5->6,7->8) peut-être contenu dans plusieurs
#   intervales [b1,b2] (e.g [2,8], [3,9], ... [5,11]
# - puis supprimer toutes les redondances (dédoublonner) du tableau pécédent
#   en transformant toute mutation en sa permutation canonique équivalente.
# - Pour cela, convertir chaque Mutation4 en Mutation.
#   Exemple : soit mut::Mutation4 = (i1,i2,i3,i4) = (4,6,10,8)
#   On crée la Mutation équivalente suivante :
#   => indices1 = [4, 5, 6, 7, 8, 9, 10] => [4, 5, 6, 8, 9, 10]
#   => indices2 = [5, 6, 4, 7, 10, 8, 9] => [5, 6, 4, 10, 8, 9]
# 
#   Le Vector des mutations (muts) est pré-alloué pour éviter les réallocations
#   dynamique (par sizehint!) puis réduit pour éviter les cases vides (toujours
#   (par sizehint!)
# 
function generate_mutations_shift2(n::Int;
                            idx_first::Int=1,
                            idx_last::Int=-1,
                            shift_max::Int=-1,
                            )

    lg4() && println("L4 generate_mutations_shift2 BEGIN at ", ms() )
    lg4() && println("L4 generate_mutations_shift2 ",
                      " idx_first=$idx_first idx_last=$idx_last shift_max=$shift_max")


    # 
    # VÉRIFICATION DU DOMAINE DES PARAMÈTRES
    # 
    if idx_last == -1    idx_last = n end
    idx_last = min(n, idx_last)
    @assert 1 <= idx_last
    @assert 1 <= idx_first < idx_last

    # shift_max == -1 => on prend le maximum possible 
    if shift_max == -1  shift_max = idx_last-idx_first  end
    @assert shift_max >= 1 "shift_max doit être strictement positif (ou nul pour max)"

    # 
    # PRÉALLOCATION DU VECTEUR DE MUTATIONS
    # Exemple : 
    #   n = 10  (nb d'avions)
    #   d = 4   (shift_max d'un shift ; peut atteindre 9 si on déplace i1 en i10)
    # evaluation de la taille maxi du nombre de shift2 (on prévoit large !)
    # 
    #  hint(n,d) = nb_win * nb_i1*nb_i2 * nb_i3*nb_i4
    #  hint(n,d) = (n-d+1) * (d+1)*d * (d+1)*d
    # 
    width = idx_last-idx_first+1 # largeur de la fenêtre e.g 10-1+1 = 10
    nb_win = (width - shift_max) # nb de fenêtres possibles (10 - 4 = 6)
    nb_i1 = shift_max+1  # e.g 4+1 = 5, car cinq i1 par win
    nb_i2 = shift_max    # e.g 4 , cinq  moins le cas i1==i2
    nb_i3 = shift_max+1  
    nb_i4 = shift_max    # dont le cas i3->i4 annulant i1->i2 !
    # Nombre maxi de quadruplet possible 
    N_MAX = nb_win * nb_i1 * nb_i2 * nb_i3 * nb_i4
    if (lg4())
        print("L4 generate_mutations_shift2 : prédimensionnement du vecteur  ($(ms())) : ") 
        @show nb_win, nb_i1, nb_i2, nb_i3, nb_i4, N_MAX
    end
    muts = Vector{Mutation}()
    sizehint!(muts, N_MAX)
 
    # 
    # GÉNÉRATION DES MUTATIONS
    # 

    # b1 et b2 sont les bornes basse et haute de la fenêtre glissante
    b1_min = idx_first
    b1_max = idx_last-shift_max
    for b1 in b1_min:b1_max
        b2 = b1+shift_max
        for i1 in b1:b2, i2 in b1:b2
            if i1 == i2
                continue # shift null inutile
            end
            for i3 in b1:b2, i4 in b1:b2
                if i3 == i4 && shift_max != 1
                    continue # shift null inutile sauf dans le cas dégénéré
                end
                if i3 == i2 && i4 == i1
                    continue # shift i3->i4 inverse annulerait i1->i2
                end
                if abs(i2-i1)==1 && i1 == i3 && i2 == i4
                    continue # shift 2->3  et 2->3 !
                end
                push!(muts, Mutation(:shift, i1, i2, i3, i4))
                # if repr(muts[end]) == "T[]->[]"
                #     println("mutation nulle avec (i1,i2,i3,i4)=", (i1,i2,i3,i4))
                # end
            end
        end
    end
    if (lg4())
        print("L4 generate_mutations_shift2 : taille réelle de muts avant unique : ") 
        println(length(muts), " ($(ms()))")
    end
    unique!(muts)
    if (lg4())
        print("L4 generate_mutations_shift2 : taille réelle de muts après unique : ") 
        println(length(muts), " ($(ms()))")
    end

    # On vérifie que le prédimensionnement du tableau estsuffisant :
    if lg4()
        println("L4 generate_mutations_shift2 : length(muts)=", length(muts))
    end

    # # 24/02/2022 : il y n'a pas de doublons a priori
    # # size1 = length(muts)
    # # unique!(muts)
    # sizehint!(muts, 0)
    # # size2 = length(muts)
    # # if size1!=size2
    # #     ln1("ATTENTION : generate_mutations_shift2 doublons évités $(size1)>$(size2)")
    # # end
    # return muts

    # 24/02/2022 : il y n'a pas de doublons a priori donc inutile ?
    len1 = length(muts)
    muts = collect(Set(muts))
    len2 =  length(muts)

    if length(muts) > N_MAX
        ln1("ATTENTION : generate_mutations_permu: sizehint sous-dimentionné")
        ln1("            sizehint=$N_MAX mais length(muts=$(length(muts))")
    end
    if lg3()
        println("generate_mutations_permu sizehint=$N_MAX llen1=$len1 len2=$len2")
    end
    return muts

end

# generate_mutations_permu: Génère un ensemble de mutations de type permu
# 
# Précondition
# - idx_first est strictement positif (défaut 1)
# - idx_last est strictement positif (mais passer -1 pour valeur maxi possible)
# - idx_first < idx_last (i.e  différent)
# - permu_size vaut au moins 2 (permuter un élément n'a pas de sens)
#
# REMARQUE SUR L'EFFICACITÉ : 
# 1. Construire un Set{Mutation} [uis convertir en Vecctor est moins efficace 
# que de construire un Vector, puis de dédoublonné avec une double conversion
# Vector->Set->Vector par collect(Set(muts))
# 2. Le dédoublonnage par collect(Set(muts)) est plus efficace que par unique(muts)
# 3. le test allunique semble aussi long que le dédoublonnage.
# 
function generate_mutations_permu(n::Int;
                            idx_first::Int=1,
                            idx_last::Int=-1,
                            permu_size::Int=-1
                            )
    if idx_last == -1    idx_last = n end
    idx_last = min(n, idx_last)
    @assert 1 <= idx_last
    @assert 1 <= idx_first < idx_last
    
    if permu_size == -1    permu_size = n end
    @assert 2 <= permu_size # il faut au moins 2 éléments pour les permuter !

    muts = Vector{Mutation}()

    # On précalcule une borne sup du nombre de permus pour le sizehint
    N_MAX = (idx_last - idx_first) * factorial(permu_size)
    sizehint!(muts, N_MAX)

    # 
    # GÉNÉRATION DES MUTATIONS
    # 
    # idx1 :  premier élément (inclu) de la zone à permuter
    # idx1+permu_size-1 : dernier élément (inclu) de la zone à permuter
    for idx1 in idx_first:idx_last-permu_size+1
        range1 = idx1:idx1+permu_size-1  # de type ::UnitRange{Int}
        indices1 = collect(range1)       # de type Vector{Int}
        permu_iter = 0
        for indices2 in permutations(range1)
            permu_iter += 1
            if permu_iter==1  continue end # car identique à  l'état courant
            push!(muts, Mutation(indices1, indices2, class=:permu))
        end
    end

    len1 = length(muts)
    muts = collect(Set(muts))
    len2 =  length(muts)

    if length(muts) > N_MAX
        ln1("ATTENTION : generate_mutations_permu: sizehint sous-dimentionné")
        ln1("            sizehint=$N_MAX mais length(muts=$(length(muts))")
    end
    if lg3()
        println("generate_mutations_permu sizehint=$N_MAX llen1=$len1 len2=$len2")
    end
    return muts
end

# generate_neighborhood: génère un voisinage complet de la solution
# 
# Retourne un vecteur formé par la concaténation des vecteurs de
# chaque famille de mutation.
#
# Selon les paramètres, les mutations du type suivant sont générées
# - add_swap : si true ajoute le simple swap (i1<->i2, 0, 0)
#                 utilise shift_min et shift_max
# - add_shift : si true ajoute le simple shift (i1->i2, 0, 0)
#                 utilise shift_max
# - add_shift2g : créer les paires de mouvements **disjoints**
#                 utilise shift_max et gap_max
# - add_shift2 :  créer les paires de mouvements **chevauchables**
#                 utilise shift_max
# - add_permu : créer les paires de mouvements **chevauchables**
#                 utilise permu_size
# - ... sera complété
# 
# Les paramètres idx_first et idx_last limite la portées des mutations
# à une partie de la solution courante.
#
# IMPORTANT
# Cette méthode dédoublonne le vecteur de mutations obtenu à l'aide de la
# commande suivante :
#   muts = collect(values(Set(muts))) # 30% plus rapide que unique! (!!)
#   # unique!(all_muts) # semble plus lent que de passer pas des Set
# 
# Description des paramètres :
# ----------------------------
# 
# - TODO: shift_min  :  
#   - 1 par défaut 
#   - >1 pour interdire les mouvements proches 
# - shift_max (alias width1) : 0 pour pas de limite ; écart maxi d'un mouvement
# - TODO: gap_min 
#   - 0 pour chevaauchement total, 
#   - 1 pour contiguë => disjoint pour shift2g)
#   - >1 pour augmenter la distance entre les deux mouvements
# - gap_max (alias width2) : 0 pour pas de limite ; limt l'écart entre les deux mouvements
# 
# 
# 
function generate_neighborhood(n::Int;
                            idx_first::Int=1,
                            idx_last::Int=-1,
                            add_shift=false,
                            add_shift2=false,
                            add_shift2g=false,
                            add_swap=false,
                            add_permu=false,
                            shift_min::Int=1,
                            shift_max::Int=-1,
                            gap_min::Int=1, # TODO
                            gap_max::Int=-1,
                            permu_size=4,
                            )

    lg4() && println("L4 generate_neighborhood BEGIN at ", ms() )
    if add_swap
        # Si add_shift est aussi choisi, on doit ignorer les swap de 1 pour  
        # éviter les doublons 
        swap_shift_min = shift_min
        if add_shift && shift_min < 2 
            swap_shift_min = 2
        end
        muts_swap = generate_mutations_swap(n,
            idx_first = idx_first,
            idx_last = idx_last,
            shift_min = swap_shift_min,
            shift_max = shift_max,
        )
    else
        muts_swap = Vector{Mutation}()
    end

    if add_shift
        muts_shift = generate_mutations_shift(n,
            idx_first = idx_first,
            idx_last = idx_last,
            shift_min = shift_min, # new ajout 09/02/2022
            shift_max = shift_max,
        )
    else
        muts_shift = Vector{Mutation}()
    end

    if add_shift2
        muts_shift2 = generate_mutations_shift2(n,
            idx_first = idx_first,
            idx_last = idx_last,
            shift_max = shift_max,
        )
    else
        muts_shift2 = Vector{Mutation}()
    end

    if add_shift2g
        muts_shift2g = generate_mutations_shift2g(n,
            idx_first = idx_first,
            idx_last = idx_last,
            shift_max = shift_max,
            gap_max = gap_max,
        )
    else
        muts_shift2g = Vector{Mutation}()
    end

    if add_permu
        muts_permu = generate_mutations_permu(n,
            idx_first = idx_first,
            idx_last = idx_last,
            permu_size = permu_size,
        )
    else
        muts_permu = Vector{Mutation}()
    end

    muts = vcat(muts_swap, muts_shift, muts_shift2, muts_shift2g, muts_permu)
    if lg4()
        print("L2: muts_swap:", length(muts_swap) )
        print("  muts_shift:", length(muts_shift) )
        print("  muts_shift2:", length(muts_shift2) )
        print("  muts_shift2g:", length(muts_shift2g) )
        print("  muts_permu:", length(muts_permu) )
        println("  => length(muts)=", length(muts) , " (avant unique)")
    end

    # Si on a activé plusieurs familles de mutations, alors on impose 
    # l'unicité du tableau final
    if add_swap + add_shift + add_shift2 + add_shift2g + add_permu >= 2
        lg4() && println("L4 generate_neighborhood avant unique! ", ms() )
        # unique!(muts)
        muts = collect(values(Set(muts))) # 30% plus rapide que unique! (!!)
        ln4("L4  => length(muts)=", length(muts) , " (avant unique)")
    end

    lg4() && println("L4 generate_neighborhood end at ", ms() )
    muts
end

# mutate! : applique sur la solution une mutation de type Mutation
# 
function mutate!(sol::Solution, mut::Mutation; do_update=true)
    permu!(sol, mut.indices1, mut.indices2, do_update=false)

    # ATTENTION : ce do_update ne concerne de les costs de la solution
    do_update && solve!(sol)
    return mut
end

# generate_nbh(n::Int, label::AbstractString)
#
# Méthode de haut niveau permettant de générer des voisinages (Vector{Mutation})
# d'après leur nom (i.e. Symbol ou String définission une union de famille 
# de voisinage.
# 
# paramètees :
# n::Int: taille de vecteur à explorer (i.e. taille de l'instaance)
# label::AbstractString: nom du voisinage ("s4", "S4","D4", "p3", "p3+d4", ...)
# 
# Génère un couple (voisinage,label) en fonction du label du voisinage.
# Le paramètre "label" d'entrée est retourné en sortie pour éviter redondance 
# dans le code appelant.
# Le label définit une famille ou une réunion de famille de voisinage 
# (s4, s4+t4, D4, ...)
# 
# La construction d'un voisinage peut être la réunion de voisinages élémentaires
# même si c'est souvent un simple appel à une méthode générique de la 
# classe Mutation
#
# Liste des options disponibles de la méthode auxiliaire 
# (qui est définie dans le fichier Mutation.jl)
# 
#    muts = generate_neighborhood(n, 
#           idx_first=1,
#           idx_last=-1,
#           add_swap=false,
#           add_shift=false,
#           add_shift2g=false,
#           add_permu=false,
#           shift_min=1,
#           shift_max=-1,
#           # gap_min=1, # TODO
#           gap_max=-1,
#           permu_size=4,
#    )
#
function generate_nbh(n::Int, label::AbstractString)

    # Vérification et dispatch de type de voisinage

    # exemple de swap : s3 : swap les **extémités d'un bloc** de largeur 3
    m = match(r"^s(\d+)$", label)
    if m != nothing
        shift = parse(Int, m[1])-1
        muts = generate_neighborhood(n, 
                add_swap=true, shift_min=shift, shift_max=shift)
        return (muts, label)
    end

    # exemple de swap : S3 : tous les swaps **dans un bloc** de largeur 3
    m = match(r"^S(\d+)$", label)
    if m != nothing
        shift_max = parse(Int, m[1])-1
        muts = generate_neighborhood(n, 
                add_swap=true, shift_max=shift_max)
        return (muts, label)
    end

    # exemple de shift t3 : 
    # déplace un avion d'une **extrémité** à l'autre d'un bloc de largeur 3
    m = match(r"^t(\d+)$", label)
    if m != nothing
        shift = parse(Int, m[1])-1
        muts = generate_neighborhood(n, 
                add_shift=true, shift_min=shift, shift_max=shift)
        return (muts, label)
    end

    # exemple de shift T3 : 
    # tous les shifts d'un avion **dans un bloc** de largeur 3
    m = match(r"^T(\d+)$", label)
    if m != nothing
        shift_max = parse(Int, m[1])-1
        muts = generate_neighborhood(n, 
                add_shift=true, shift_max=shift_max)
        return (muts, label)
    end

    # exemple de shift d3 : 
    # shift ou swap d'avions d'une **extrémité** à l'autre d'un bloc de largeur 3
    m = match(r"^d(\d+)$", label)
    if m != nothing
        shift = parse(Int, m[1])-1
        muts = generate_neighborhood(n, add_swap=true, add_shift=true,
                 shift_min=shift, shift_max=shift)
        return (muts, label)
    end

    # exemple de shift D3 : 
    # tous les shifts ou swaps d'avions **dans un bloc** de largeur 3
    m = match(r"^D(\d+)$", label)
    if m != nothing
        shift_max = parse(Int, m[1])-1
        muts = generate_neighborhood(n, add_swap=true, add_shift=true,
                shift_max=shift_max)
        return (muts, label)
    end


    # exemple de permu P4 : toutes les permutation dans un bloc de largeur 4
    m = match(r"^P(\d+)$", label)
    if m != nothing
        permu_size = parse(Int, m[1])
        muts = generate_neighborhood(n, 
                add_permu=true, permu_size=permu_size)
        return (muts, label)
    end

    # exemple de shift2 2T4 :  Toutes les paires de shifts recouvrables contenus 
    # dans un bloc de largeur 4
    m = match(r"^2T(\d+)$", label)
    if m != nothing
        shift_max = parse(Int, m[1])-1
        muts = generate_neighborhood(n, 
                add_shift2=true, shift_max=shift_max)
        return (muts, label)
    end

    # exemple de shift2g 2T4g9 : Toutes les paires de shifts disjoints contenus 
    # dans deux blocs de largeur 4 distants d'au maximum 9
    # (2T4g1 correspond à des paires des shifts disjoints mais contiguës
    m = match(r"^2T(\d+)g(\d+)$", label)
    if m != nothing
        shift_max = parse(Int, m[1])-1
        gap_max = parse(Int, m[2])
        muts = generate_neighborhood(n, 
                add_shift2g=true, shift_max=shift_max, gap_max=gap_max)
        return (muts, label)
    end

    # exemple de shift et swap large : D15.6
    # tous les shift ou swap de lageur comprise entre 6 et 15 inclues
    # Motivation P5 et D15.6 sont deux voisinages non redondant
    # 
    m = match(r"^D(\d+)\.(\d+)$", label)
    if m != nothing
        shift_max = parse(Int, m[1])-1
        shift_min = parse(Int, m[2])-1
        muts = generate_neighborhood(n, 
                shift_min=shift_min, shift_max=shift_max, 
                add_shift=true, add_swap=true)
        return (muts, label)
    end

    # TODO p5: toutes la permutation dont les extrémité d'un bloc de largeur 
    # 5 sont déplacées (s5 est inclue dans p5 mais pas s4)
    #
    m = match(r"^p(\d+)$", label)
    if m != nothing
        #println("ERREUR voisinage pn non implanté"); exit(1)
        shift_max = parse(Int, m[1])-1
        muts = generate_neighborhood(n, add_swap=true, add_shift=true,
                shift_max=shift_max)
        return (muts, label)
    end

    # TODO p6.3: arrangement de 6 éléments pris 3 par 3 d'écartement 6
    # Représente tous les permutations de 3 éléments de largeur totale exactement de 6
    # (n'inclue pas s5, mais inclue s6 car on accrpte de ne déplacer que 2 éléments)
    # IMPORTANT : CONFIRMER LA SPÉCIFICATION DE LA PHRASE PRÉCÉDENTE
    m = match(r"^p(\d+)\.(\d+)$", label)
    if m != nothing
        #println("ERREUR voisinage pn non implanté"); exit(1)
        shift_max = parse(Int, m[1])-1
        muts = generate_neighborhood(n, add_swap=true, add_shift=true,
                shift_max=shift_max)
        return (muts, label)
    end


    # Vérification que la patterne est une composition de patternes élémentaires
    # S'il n'y a pas de caractère "+" c'est que cette mutation élémentaire 
    # n'existe pas (e.g erreur dans la première lettre ...)
    if match(r"^.*\+.*$", label) == nothing
        # println("ERREUR : nom de mutation non reconnue : \"$(label)\"")
        println(Crayon(foreground = :red, background = :light_yellow, bold = true))
        println("ERREUR : nom de mutation non reconnue : \"$(label)\"")
        println(Crayon(reset = true))
        exit(1)
    end

    # TODO split r"\+" et fusionner récursivement les voisinages
    # label est de la forme : "P6+T10+2T4g10"
    labs = split(label, r"\+")
    all_muts = Vector{Mutation}()
    for lab in labs
        # ln1("AVANT generate_nbh $lab $(ms())")
        muts, _ = generate_nbh(n, lab)
        # ln1("APRES generate_nbh $lab $(ms()) all_muts:$(length(all_muts)) muts:$(length(muts))")
        all_muts = vcat(all_muts, muts)
        # ln1("APRES vcat $lab $(ms()) all_muts:$(length(all_muts))")
        unique!(all_muts)
        # ln1("APRES unique! $lab $(ms())  all_muts:$(length(all_muts))")
    end
    return (all_muts, label)

    error("label de voisinage inconnu : $(label)")
end

