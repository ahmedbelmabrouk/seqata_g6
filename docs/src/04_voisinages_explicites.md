# Description explicite des voisinages

### Présentation

Dans les méthaheuristiques, un voisinage est défini de manière abstraite comme 
un ensemble de voisins d'une solution particulière (e.g la solution courante).
On peut considérer un voisin particulier tiré aléatoirement (descente aléatoire), 
ou le voisinage dans son ensemble ("steepest" descente, taboo, VNS, ...).

Un opérateur de voisinage définit toute une famille de voisinages selon 
son paramétrage. Ainsi, une mutation pourra être définie simplement par la
connaissance de deux indices avec la consigne de déplacer le premier élément
vers le deuxième indice, ou encore de permuter les éléments positionnés en
chaque indices, ou encore nécessiter la connaissance de deux paires
d'indices etc...

Il est donc pratique d'associer un label à une famille de voisinages 
("2opt", swap(2), swap(3), shift(3), permu(30 ...).

En vue de son exploration, un voisinage peut-être constuit **soit implicitement**
en testant par exemple un voisin à l'intérieur d'une double boucle sur deux indices,
**soit explicitement** en construisant a priori un vecteur de mutations avant d'un 
extraire chaque voisin à tester.

La première approche (implicite) nécessite un traitement spécifique qui doit
être adapté à chaque famille de voisinage (swap, shift, permu, ...).
La seconde approche (explicipe) permet la combinaison (union, ...) de deux voisinages pour en former un plus large. Cependant elle oblige à trouver une représentation générique d'un opérateur de voisinage.


La suite de ce texte propose une méthode pour définir et gérer de manière
uniforme quelques familles de voisinage applicables au problème Seqata 
(ou ALP en général). Elle propose :
- une structure standard pour représenter une mutation quelconque
- un liste de conventions pour nommer des familles classiques de mutation.


### Une structure commune pour les mutations

Pour une permutation d'avion donnée, une mutation quelconque peut être définie
par deux vecteurs d'indices correspondant
aux anciennes et nouvelles positions des avions mutés dans la solution.
Par exemple, soit une mutation de sol1 en sol2 définie comme suit :

    sol1=[1,2,3,4,5,6,7,8,9,10]
    sol2=[1,2,8,9,5,6,7,3,4,10]

Cette mutation peut donc s'exprimer par le couple de vecteur suivant :

    [3,4,8,9]->[8,9,3,4]

Il faut noter que la mutation précédente dispose de nombreuses mutations 
équivalentes et devra donc être normalisée (sous une forme canonique) pour
des besoins éventuels de comparaison ou simplement pour éviter les doublons.
Ces deux mutations :

    [1,2,3,4,8,9]->[1,2,8,9,3,4]
    [4,3,8,9,1]->[9,8,3,4,1]

représente la même mutation que la précédente à savoir : [3,4,8,9]->[8,7,3,4].

En pratique, un voisinage est constitué de mutations appartenant à une certaine 
famille homogène (échange de deux avions, déplacement d'un avion, permutation 
d'un petit bloc d'avions) et d'une certaine largeur. Les voisinages courants 
peuvent donc être simplement caractérisés par un nom (ou label) de voisinage 
qui pourra par exemple être passé en paramètre du programme
(actuellement via l'option `--nbh` (pour neighborhood).

Enfin il peut être utile de reconnaitre à quelle famille fait partie une mutation 
particulière, ne serait ce que pour l'affichage mais aussi pour des besoins de
filtrage.  
La structure utilisée pour représenter une mutation est la suivantes (en julia).

    struct Mutation_permu <: Mutation
        class::Symbol
        indices1::Vector{Int}
        indices2::Vector{Int}
    end

La structure Mutation est un type abstrait permettant l'implantation de type dérivés
spécialisés et plus efficace pour certains voisinages
Le type (classe) Solution pourra disposer d'une méthode `apply!` lui permettant 
d'appliquer ce voisin quelque soit le type effectif de la Mutation.

    apply!(cursol, mut)
    => modifie la solution cursol selon la mutation mut, et résoud le 
       sous-problème de timing.

### Conventions pour des noms de voisinages standard

Voici ci-dessous une convention adaptée au problème d'atterrissage d'avions 
ou plus généralement à de petites perturbation d'une liste (ordonnée).

Le nom d'un voisinage simple est spécifié par :
- une lettre caractérisant la famille du voisinage (la version majuscule représente
  un ensemble plus grand),
- un premier nombre : **largeur maximale du bloc perturbé** 
  (e.g s3, S3, t3. p3),
- un second nombre éventuel dont la signification dépend de la famille de
  la mutation (par exemple  pour imposer une largeur minimale du déplacement).

Liste des principales familles de mutations standard (avec exemples) :

- t3 : déplace un avion d'une extrémité à l'autre d'un bloc de largeur 3  
  shifT ou Translate d'un avion de exactement de +/-2  
  La permutation est une rotation d'un bloc de largeur 3
- s3 : échange des extrémités d'un bloc de largeur 3 
  (Swap de deux avions espacés exactement de 2 : [4,5,6]->[6,5,4])
- t1 = s1 : bloc de 1 (pas de mutation possible)
- t2 = s2 : swap de deux avions voisins (inversion d'un bloc de 2)
- T4 = t2+t3+t4 : tous les shifts possibles dans un bloc de largeur 4 
  (donc tous les déplacements d'un avion d'au maximum de +/-3
- S4 = s2+s3+s4: tous les swaps de deux avions dans un bloc de 4
- d4 = s4+t4 : union de s4 et de t4
  (swaps et shifts des extrémités d'un bloc de 4)
- D4 = S4+T4 : tous les swaps et shifts dans un bloc de 4
- D12.7 Tous les swaps et shifts de largeurs entre 12 et 7 inclus
  D12.7 = d7+d8+d9+d10+d11+d12 = D12-D6 ou D12\D6 (was w7.14 dans alap)
- D4.1 : identique à D4 (swaps et shifts le largeur 4)
- 2T4: deux shifts simultanés **potenciellement recouvrants** dans 
  un bloc de largeur 4
- 2T4g12 : deux shifts simultanés **disjoints** de largeur maxi 4 
  avec un gap de maxi de 12 entre ces deux shifts (gap entre 1 et 12)
  TODO accepter les gap négatifs : 
  - 2T4g1: les deux blocs sont contigüs (gap entre 1 et 1)
  - 2T4g-4 : idem 2T4 (les deux blocs se superposent entièrement)
    *MAIS ATTENTION : il faudrait alors utiliser "\" au lieu de "-" 
    pour la différence ensembliste de deux voisinages*.

- P4 : toutes les permus d'un bloc de largeur 4 (moins la permu courante !)
  (intègre les mutations s2 de largeur effective 2 !)
  (Anciennement appelée p4 dans projet alap)

- p4 : toutes les permus d'un bloc de largeur *effective* de 4.  
  Autrement dit p4 ne contient aucun élément de P3  
  Autrement dit "p4" = "P4-P3" (alias "P4\P3").
  Ce voisinage peut-être intéressant pour créer une série de voisinages
  croissants et non redondants pour une méthede VNS ;-)  
  Exemple de suite VNS :  (P3,p4,p5,p6,p7,p8)

  Exemple de voisinage appartenant à p4 : [1,2,4]->[4,1,2] 
  (ne peut pas être représenté par un élément de P3).

### Composition de voisinages

Il est possible de composer des voisinages en tant qu'ensembles de mutations :

- avec + pour la réunion de voisinages élémentaires.
  Les doublons sont supprimés du résultat.
  Exemple "p6+D14.7+2T3g14" : union de trois voisinages sans les doublons,

- avec "-" ou "\" pour l'exclusion de voisinages (différence ensembliste).
  Le "-" sera à éviter si l'on souhaite pouvoir utiliser un paramètre 
  négatif dans un voisinage ; exemple "2T4G-4" définirait deux mutations
  T4 totalement recouvrable.

  *NOTA : PAS SÛR DE L'UTILITÉ DES PARAMÈTRES NÉGATIFS*,

- Sachant qu'on ne va pas gérer l'associativité, la suppression des
  ensembles par `-` se fera après les réunions des voisinages par `+`.  

  Par conséquent, V = "U1-X2+Y3-Z4" sera interprété comme "U1+Y3-X2-Z4" :  

### Idées de voisinages supplémentaires

#### `2t7` `2d7` 

De largeur effective 7 (2d7 inclue les swap contrairement à 2t7).
Exemple de mutation de type 2t7

  `s(1,4)+s(3,7)=>[1,2,3,4,5,6,7]->[4,2,7,1,5,6,3]`

#### `3d9` 

Trois shifts ou swaps simultanés de largeur totale exactement de 9

  => permet de quitter un minimum local difficile de l'instance alp13.

#### `2s10` (facile)

Motivation :
- perturbe moins le voisinage qu'un simple move T10
- maintien ensemble les avions intermédiaires
- facile à implanter
- 2s10 : deux swaps simultanés de largeur totale : exactement 10
- 2S10 : deux swaps simultanés de largeur totale : au maximum 10
