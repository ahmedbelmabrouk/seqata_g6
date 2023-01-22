# Présentation de projet Seqata

### Introduction

SEQATA (SÉQuencement d'ATterrissage d'Avions) est le nom d'un projet pour le 
cours de Recherche Opérationnelle SOD324 (ENSTA Paris).

La présentation formelle du problème SEQATA, le sujet du projet et les détails du 
travail à réaliser sont disponibles dans un fichier pdf indépendant.

Cette documentation décrit le fonctionnement et certains détails d'implantation 
du prototype de code fourni au élève pour la réalisation de ce projet.

Ce prototype est opérationnel et intègre les fonctionnalités suivantes :
- organisation multi fichiers pour simuler les conditions d'un gros projet
  (en pensant à la programmation objet et transposable en C++),
- implémentation avec choix entre plusieurs méthodes méthaheuristiques
  (dont recuit et un taboo évolué) et exactes (PLNE),
- gestion de nombreuses options dont le premier argument est une action
  (à la git) pour le choix de la méthode de résolution du problème complet,
  ou de l'algorithme du résolution du sous-problème de timing,
- choix par une option du solveur linéair externe quand un tel solveur
  est nécessaire (e.g cplex, glpk, clp, ...).

La liste des options disponibles peut-être obtenue par :

    ./bin/run.jl -h

Ou peut avoir des détails complémentaires ur les actions disponibles :

    ./bin/run.jl help

Voir en fin de fichier pour des exemples d'utilisation du code.


### Rappel du problème SEQATA 

SEQATA est une version simplifiée du problème de séquencement d'atterrissage 
d'avions (ALP pour Aircraft Landing Problem).
- prise en compte d'une seule piste d'atterrissage,
- les coûts sont linéaires par morceau en passent par 0 (en V assymétrique).

Un avion est caractérisé par :
- sa date d'atterrissage souhaitée (``T_i`` = `target`),
- ses dates d'atterrissage au plus tôt (``E_i = lb_i``) et au plus tard 
  (``L_i = ub_i``),
- sa catégorie (`kind` : plus ou moins gros, ...),
- des coûts de pénalité par unité de retard (tp = tardiness penality),
  ou d'avance (ep = earliness penality).

Une durée d'écart minimale ``S_{kl}`` doit être respectée entre un avion 
de type ``k`` et un avion de type ``l``.
En effet un gros porteur peut sans inconvénient atterrir juste après un
"delta-plane" alors que l'inverse n'est pas du tout vrai !
L'instance définit donc une matrice ``S_{kl}`` pour chaque type d'avion
possible.


L'objectif consiste à affecter la date d'atterrissage de chaque avion
en minimisant le coût de pénalité global tout en respectant les contraintes
sur les dates limites d'atterrissage.


### Exemples d'utilisation

#### Résolution du sous-problème de timing

On passe la clé "timing ou tim", une instance (option -i) et une liste d'avions
(option `--planes` ou `-p`). Le programme appelle un des solveurs de timing
disponible en fonction de l'option `--timing-algo-solver` ou `-t`
(e.g `-t lp` ou `-t earliest`)
Dans le prototype, seul le solver EarliestTimingSolver est fonctionnel.
Un fichier est générer dans le sous répertoire "_tmp/"

    ./bin/run.jl tim -t earliest -i data/01.alp  -p 3,4,5,6,7,8,9,1,10,2
    => 2830.0
    # Si le type LpTimingSolver est implanté :
    ./bin/run.jl tim -t lp -i data/01.alp  -p 3,4,5,6,7,8,9,1,10,2
    => 700.0


#### Validation d'une solution existante


On passe la clé "val", une instance (option -i) et une solution (option -i) et la
programme indique si la solution est valide ou liste les viols de contraintes.

    # test de la solution générées précédemment  avec l'algo -t earliest
    ./bin/run.jl val -t earliest -i data/01.alp -s _tmp/alp_01_p10=2830.0.sol
    => Solution correcte de coût : 2830.0
 
    # test de la solution optimale du sujet
    ./bin/run.jl val -t lp -i data/01.alp -s sols/alp_01_p10_k3=700.0.sol
    => Solution correcte de coût : 700.0

Remarque : les clés d'option -i pour l'instance et -s pour la solution sont
facultatives et seront omises dans les exemples suivants.

#### Heuristique de Monte Carlo (action carlo)

    ./bin/run.jl carlo data/01.alp  -t earliest -n 10000000 -L1
    => 2110 en 20.8s avec 250_000 calls/sec (10e6 itér)


    ./bin/run.jl carlo data/01.alp  -t lp -n 10000 -L1
    => 2080 en 36.7s avec 307 calls/sec  (10e4 itér) (l'optimum vaut 700)

#### Heuristique d'exploration (action explore)


    # test d'une exploration aléatoire (action explore) à partir d'une solution
    # aléatoire (--presort shuffle par défaut) pendant 10^7 itérations
    # On utilise ici l'algo de timing earliest (non optimal mais fourni avec
    # le proto)
    # Les solutions sont enregistrées dans le sous-répertoire _tmp.

    ./bin/run.jl explore  data/01.alp -t earliest -n 10000000 -L2
    # =>
    # ...
    ======================================================================
    Début de l'action explore
    Solution correspondant à l'ordre de l'instance
    cost=25910.0  :[p1,p2,p3,p4,p5,p6,p7,p8,p9,p10]
    Solution initiale envoyée au solver
    cost=19990.0  :[p10,p5,p1,p9,p8,p4,p3,p2,p6,p7]
    BEGIN solve!(ExploreSolver)
    iter <nb_move>=<nb_improve>+<nb_degrade> => <bestcost>
    iter 1=1+0 => cost=11150.0  :[p7,p5,p1,p9,p8,p4,p3,p2,p6,p10]
    iter 53=17+36 => cost=10180.0  :[p5,p3,p10,p9,p6,p4,p7,p8,p2,p1]
    iter 55=18+37 => cost=4540.0   :[p5,p3,p8,p4,p6,p9,p7,p10,p2,p1]
    iter 171=51+120 => cost=3870.0   :[p3,p4,p7,p5,p6,p1,p9,p8,p10,p2]
    iter 2733=796+1937 => cost=3070.0   :[p4,p5,p3,p6,p8,p9,p7,p1,p10,p2]
    iter 19724=5571+14153 => cost=3010.0   :[p4,p3,p5,p6,p7,p8,p1,p9,p10,p2]
    iter 77586=21977+55609 => cost=2950.0   :[p4,p3,p7,p6,p8,p9,p5,p1,p10,p2]
    iter 169472=48162+121310 => cost=2780.0   :[p3,p5,p4,p6,p7,p8,p1,p9,p10,p2]
    iter 242903=69022+173881 => cost=2500.0   :[p3,p5,p4,p8,p7,p6,p9,p1,p10,p2]
    iter 263272=74779+188493 => cost=2470.0   :[p4,p3,p6,p5,p7,p9,p8,p1,p10,p2]
    iter 1034240=293687+740553 => cost=2110.0   :[p4,p3,p6,p5,p7,p8,p9,p1,p10,p2]
    END solve!(ExploreSolver)

    meilleure solution trouvée :
    ======================================================================
    cost=2110.0   :[p4,p3,p6,p5,p7,p8,p9,p1,p10,p2]
    ======================================================================

Quant vous aurez implanté l'algorithme LpTimingSolver, la solution optimale
de coût 700 pourra être régulièrement atteinte en 10e7 itérations

    ./bin/run.jl explore  data/01.alp -t lp -n 10000 -L1
    # => optimum rarement atteignable (car -n10000 au lieu de -n10000000)
    => affichage (pour 10e6 iterations !)
    Début de l'action explore
    iter 2:1+/1- bestsol=cost=13850.0  :[p9,p3,p4,p1,p8,p10,p2,p5,p7,p6]
    iter 10:3+/7- bestsol=cost=9730.0   :[p8,p4,p5,p9,p6,p1,p3,p2,p7,p10]
    iter 48:13+/35- bestsol=cost=7340.0   :[p7,p5,p6,p3,p10,p9,p8,p4,p1,p2]
    iter 174:49+/125- bestsol=cost=6710.0   :[p3,p8,p7,p4,p5,p1,p6,p10,p2,p9]
    iter 175:50+/125- bestsol=cost=6470.0   :[p3,p8,p5,p4,p7,p1,p6,p10,p2,p9]
    iter 190:54+/136- bestsol=cost=6100.0   :[p6,p7,p4,p5,p9,p3,p8,p10,p2,p1]
    iter 269:77+/192- bestsol=cost=3940.0   :[p5,p3,p4,p7,p9,p8,p1,p10,p6,p2]
    iter 670:191+/479- bestsol=cost=3800.0   :[p3,p4,p1,p6,p8,p7,p5,p9,p10,p2]
    iter 1776:520+/1256- bestsol=cost=2800.0   :[p3,p5,p4,p9,p8,p7,p6,p10,p1,p2]
    iter 1915:573+/1342- bestsol=cost=2560.0   :[p4,p5,p3,p6,p9,p8,p7,p10,p1,p2]
    iter 24124:6934+/17190- bestsol=cost=1600.0   :[p3,p5,p4,p6,p7,p8,p9,p10,p1,p2]
    iter 67661:19230+/48431- bestsol=cost=1090.0   :[p3,p4,p5,p8,p6,p7,p9,p10,p1,p2]
    iter 248515:69936+/178579- bestsol=cost=700.0    :[p3,p4,p5,p7,p6,p8,p9,p1,p10,p2]
    meilleure solution trouvée :
    ======================================================================
    cost=700.0    :[p3,p4,p5,p7,p6,p8,p9,p1,p10,p2]
    ======================================================================

#### Heuristique de descente (action descent)

L'efficacité de la descente dépend beaucoup de la largeur et de la pertinence
de votre voisinage !

- avec algo de timing "earliest"
  
      ./bin/run.jl des  data/01.alp  --presort shuffle -t earliest -L1 -n 100000
      => peut trouver 2110 en moins de 200 itérations !

- avec algo de timing "lp" (à faire par les élèves !)
  
      ./bin/run.jl des  data/01.alp  --presort shuffle        -L1 -n 10000
      ./bin/run.jl des  data/01.alp  --presort shuffle  -t lp -L1 -n 10000
      => peut trouver l'opimum 700 en moins de 200 itérations

#### Méthode exacte PLNE avec coût discrétisé (action dmip)

Il y a une variable distincte pour chaque avion et chaque date d'atterrissage
possible ! Cette variable est donc directement associée à un coût d'avion
quand celui-ci atterrit à une date précise.

Avec CPLEX ou Gurobi

    ./bin/run.jl dmip data/01.alp -x cplex
    ./bin/run.jl dmip data/01.alp -x gurobi
    => 700 en 11s à 15s  (sur macbookpro mi-2015 avec précompilation faite)

Avec CBC (gratuit)

    ./bin/run.jl dmip data/01.alp -x cbc
    ./bin/run.jl dmip data/01.alp
    => 700 en... 6mn sur serveur de calcul maury (UMA) 

Avec GLPK (gratuit)

    ./bin/run.jl dmip data/01.alp -x glpk
    => 700 en... 24mn sur serveur de calcul maury (UMA)

### Tests unitaires de ce programme

Une bonne partie des fonctionnalités proposées (ou à développer) dispose de
ses propres tests unitaires. Un test unitaire consiste à appeler depuis julia
une ou plusieurs fonctions correspondant à une fonctionnalité donnée.
Ces tests sont indépendant et sont regroupés dans le répertoire test, et le
programme de lancement peut s'exécuter depuis le répertoire de travail.

    ./test/runtests.jl
    => exécute tous les tests commençant par "test/test-xxxx.jl".
    (on peut préciser des exceptions dans le fichier "runtests.jl")
    Une synthèse des tests réussis ou échoués est indiquée à la fin de l'exécution.

    ./test/runtests.jl ./test/test-07-validate.jl
    => test seulement le ou les fichiers précisés

En fonction de l'évolution de Julia, des solvers disponibles, ...
Certaines fonctionnalités de ce projet peuvent ne pas fonctionner.

L'ensembles des tests unitaires peuvent donc être exécutué par un seul
appel à Julia.

### Tests fonctionnels de l'application

Ils consistent à vérifier le bon fonctionnement du point de vue de l'utilisateur
final. Chaque test consiste donc en un appel à l'application depuis le système 
d'exploitation (e.g unix). C'est typiqument ce que fait un enseignant pour
vérifier le fonctionnement du projet, mais également le développeur pour s'assurer
que les fonctionnalités demandées sont bien opérationnelles.

Pour simplifier ces exécutions répetitives, un script Ruby générique est fourni
(`bin/funtests.jl`) ainsi qu'un fichier de configuration `Funtests.json`
spécifique au projet Seqata. Son utilisation est facultative et vous pouvez 
l'adapter à votre besoins.

Voir le contenu du script pour les prérequis d'installation.

### Précompilation du projet

Lorsque votre projet est stabilisé et que nous souhaitez effectuer une série 
d'essai (en jouant sur les options ou les instances), vous pouvez accélérer 
le lancement de l'exécutable unix en créant une image précompilée du projet.

Pour cela taper la commande suivante :

    ./bin/build_sysimg.jl

Cette commande fait appel à un script auxiliaire `test/precompile_script.jl` 
qui exécute le maximum de fonctionnalités possibles. Si celle-ci ne sont
pas encore implantées dans votre code, vous serez peut-être amené à en
commenter une partie.

Le résultat est une image julia optimisée pour votre machine et rangée 
dans un répertoire temporaire sous le nom `/tmp/julia_sysimg_seqata.so`
qui est appelé automatiquement par les scripts de lancement `bin/run.jl` 
(très utile pour les tests fonctionnels).
Mais attention la création de l'image est lente (e.g 6mn !).
### À propos des instances

Les instances d'origines (de la bibliohèque orlib) sont accessibles depuis le
site :
- <http://people.brunel.ac.uk/~mastjjb/jeb/jeb.html#aircraft>
- <http://people.brunel.ac.uk/~mastjjb/jeb/orlib/airlandinfo.html>
- <http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/>

Ces 13 instances ont été transformées dans un nouveau format plus lisible, plus dense
et mieux adapté à une évolution de la fonction de coût.
(voir sujet de projet pour sa description).




