# Quelques méthodes de manipulation de tableau/vecteur

export shift!

# shift! : déplacement d'un élément dans un tableau
# - v : vecteur à modifier
# - idx1 : indice de l'élément à décaller
# - idx2 : indice de l'émément après déplacement
# Résultat
# - modifie le tableau passé en paramètre
# - retourne le tableau complet
#
# Exemple :
#   v = collect(1:7)
#   shift!(v, 7, 1)'
#   => 1×7 LinearAlgebra.Adjoint{Int64,Array{Int64,1}}:
#      7  1  2  3  4  5  6
#
# Voir aussi les méthodes intégrées au langage julia :
# - permutate()
# - circshift()
#
function shift!(v::Vector{T}, idx1::Int, idx2::Int) where {T}
    if idx1 < idx2
        # si idx1=2 et idx2=6
        # alors v[2,3,4,5 , 6] devient v[3,4,5,6 , 2]
        v[[idx1:idx2-1; idx2]] = v[[idx1+1:idx2; idx1]]
    else
        v[[idx2+1:idx1; idx2]] = v[[idx2:idx1-1; idx1]]
    end
    return v
end
