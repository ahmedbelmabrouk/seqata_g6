# Style de programmation Julia versus POO classique


### Différence entre approche objet de Julia et de C++

Soit la méthode `get_sep` de classe `Instance` destinée à retourner le temps de
séparation minimum entre deux avions p1 et p2.

En C++ vous créeriez une nouvelle méthode de la classe `Instance` qui 
s'utiliserait comme suit :

```
sep = inst->get_sep(p1.id, p2.id)
```
La variable `inst` étant propriétaire de l'objet, elle ne fait pas partie des deux paramètres de la méthode. Le dispatch vers le code de la méthode 
se fait selon selon le type de l'object propriétaire (simple dispatch).

En Julia, les méthodes ne font pas partie des types (i.e classe). Il n'y a pas de notion de propriétaire de méthode et le dispatch se fait selon le type de l'ensemble des paramètres de la méthode (multiple dispatch).

L'appel précédent d'écrirait comme suit :
```
sep = get_sep(inst, p1.id, p2.id)
```

L'intérêt est une banalisation de méthode qui (grâce au multiple dispacth) dispense souvent de l'héritage. Si la méthode précédente n'éstait pas prévue dans le prototype fourni, rien ne vous empêcherait de la créer sans avoir à modifier le type `Instance` du proto.



### Utiliser les objets et méthode de haut niveau

#### Dilème "collection de structures versus structure de collections"

Un dilème classique existe : collection de structures versus structure de collections.
En général, il est plus intéressant de manipuler collection de structures de haut niveau plutôt plutôt que des collections de vecteurs d'attributs séparés.



C'est pour cela que l'instance ne contient pas les 50 attributs suivants : lbs, targets, ubs, ... mais simplement le vecteur planes. Pour l'algorithmique de haut niveau il faut (généralement) raisonner en POO comme C++ ou python, et pas comme matlab ni comme en math. 

Cependant dans certains cas vous pourriez avoir besoin de faire un calcul complexe sur un seul (gros) vecteur de flottants (e.g `ubs` 
vecteur de `planes.ub`).
Il est alors toujours possible de créer ce vecteur juste avant son exploitation. Par exemple, il peut être judicieux de convertir la représentation d'un graphe sous une autre forme au début d'une fonction
rien que pour les besoins d'un algorithme particulier.

#### Éviter les méthodes de bas niveau du prototype

Par exemple, la classe Instance dispose d'une matrice `sep_mat` (qui devrait
être privée) indicées par le type d'avion (attribut `kind`).
Vous aurez besoin de l'exploiter indirectement quand vous devrez connaitre 
le temps minimum de séparation entre deux avion suxessif p1 et p2.

Une manière d'accéder à cette information pourrait-être :

```
xxx = inst.sep_mat[p1.kind, p2.kind] + xxx
```

Cependant il est bien plus simple et plus lisible d'utiliser la méthode
de plus haut niveau suivante :

```
xxx = get_sep(sv.inst, p1, p2) + xxx
```

#### Application au projet Seqata

Manipulez directement les avions (Plane) plutôt que leurs attributs (lb, ub , ...) et éviter de créer des vecteurs secondairs.

Par exemple plutot que de créer un vector `targets` pour mémoriser toutes les valeurs `target` des avions, il est bien plus lisible et généralement aussi efficace d'utiliser directement le vecteur d'avions `planes` et d'en récupérer l'attribut target (et d'autres) quand c'est nécessaire. De plus de nombreuses méthodes proposées reoivent directement un objet `p::Plane` en paramètre ou d'autre un objet `inst::Instance`.

Regarder l'exemple suivant :

Au lieu de passer par une collection intermédiaire comme :
```
# Création inutile d'un vecteur targets
targets = Dict(p => p.target for p in planes)
# puis plus tard 
xxx =  targets[planes[i]] + ..
```

C'est bien plus lisible d'utiliser directement l'attribut `ub` de l'objet
avion "
```
p1 = planes[i]
xxx =  p1.ub + ..
```

#### Illustration par un extrait de Seqata


Extrait de la classe MipDiscretTimingSolver (Pourquoi croyez vous que je 
l'ai laissée dans le proto alors qu'elle est inutile pour ce projet ? :-) 

```
for p1 in planes, p2 in planes
    if !(p1.id < p2.id)
        continue
    end # évite de traiter le cas symétrique
    # ici p1 peut être avant ou après p2
    for t_i = p1.lb:p1.ub
        t_lb_j = max(p2.lb, t_i - get_sep(sv.inst, p2, p1) + 1) 
        t_ub_j = min(p2.ub, t_i + get_sep(sv.inst, p1, p2) - 1) 
        for t_j = t_lb_j:t_ub_j
            @constraint(model, y[p1, t_i] + y[p2, t_j] <= 1)
        end
    end
end
```

L'exemple précédent est bien plus lisible que si vous aviez multiplié les tableaux intermédiaires en trainant des indices partout.



