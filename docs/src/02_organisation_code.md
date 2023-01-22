# Organisation du code julia du projet SEQATA


### Organisation globale de l'application

- l'arborescence du projet est le suivant
  - `bin/` : contient le ou les exécutables ou utilitaires de votre application 
    ```julia
      ./bin/run.jl -h  
      ./bin/run.jl help
    ```
  - `data/` : le jeu des 13 instances de références utilisées pour ce projet.
    Elle sont nommées de `01.alp` (10 avions) à `13.alp` (500 avions)
  - `docs/` : la... documentation du proto !
    (vous pourrez laisser votre rapport pdf à la racine du projet)
  - `sols/` : les meilleures solutions (validées) que vous aurez trouvées toute
    méthode confondue,
  - `src/` : l'ensemble des fichiers sources
  - `test/` : répertoire contenant des tests unitaires. Vous n'êtes pas obligés 
     de les compléter ni même de les utiliser. Mais ces tests peuvent vous faciliter 
     le maintien opérationnel ou la réparation de votre projet git après de grosse modifications.

     Son utilisation est de la forme :
    ```julia
      ./test/runtests.jl test-07-validate.jl  
      ./test/runtests.jl                    # exécute tous les tests  
    ```
  - `_tmp/` : répertoire temporaire ou seront enregistréea par défaut les solutions 
     générées lors de chaque amélioration.
     À vous de recopier celles qui méritent d'être conservées vers le répertoire `sols/`


### Description des principaux fichiers

#### À la racine

- `Rakefile` est un équivalent du Makefile pour ruby.
   Exemple d'utilisation
   - `rake help` affiche les cibles les plus utiles
   - `rake -T` affiche toutes les cibles
   - `rake dc` (distclean) vide le répertoire `_tmp`
   - `rake zip` **crée une archive datée dans le répertoire parent** :
      très utiles pour archiver ou échanger vos projets si vous n'utilisez 
      pas git !
- `Project.toml` et `Manifest.toml` propres à la gestion des packages par julia pour
  ce projet. En particulier `Project.toml` liste tous les packages requis pour 
   disposer de toutes les fonctionnalités du code (dont à la fois CPLEX, Clp et GLPK !).
- `Funtests.json` fichier décrivant une liste de **tests fonctionnels**.
   Il n'est utilisable que sur le réseau ensta avec avoir tapé `usediam ro` par la commande :
   ```
   funtests -h
   funtests list
   funtests run timing-earliest
   funtests run                   # effectue tous les tests
   ```
   Contrairement aux tests unitaires (`test/runtest.jl`) qui vérifie le fonctionnement
   en julia de méthodes individuelles, les tests fonctionnels sont plus proches de 
   l'utilisateur et permettent de valider une exécution complète de scenarios depuis 
   le terminal unix.

#### Dans le sous-répertoire bin

- le fichier `bin/run.jl` est l'exécutable principal chargé du lancement de 
  l'application (voir `LISEZ_MOI_ELEVE` pour des exemples d'utilisation).
  Son rôle est de positionner les bibliothèques et de lancer un des deux modes de fonctionnement :
  - **soit passer en mode interactif** (pratique pour le développement) en
    chargeant le maximum de bibliothèques possibles dès le départ et en affichant
    l'invite de Julia (mode REPL),
  - **soit passer en mode exécution** en lanant la fonction `main()` située dans
    le fichier `src/main.jl`

#### Dans le sous-répertoire src

- le ficher `main.jl` sous-traite le travail à une méthode  `main_xxx()`
  spécialisé en fonction de l'action détecté par l'analyse des arguments (e.g
  `carlo` => `main_carlo()`).

- les fichiers `main_xxx.jl` définissent les méthodes `main_xxx()` qui représente 
  le programme principal dédié à l'action `xxx`.
  En fonction des options, ces méthodes créent, exécutent et exploitent les
  différents solveurs.

- le fichier `instance.jl` regroupe tout ce qui concerne l'instance (aucune 
  intelligence. Il gère et manipule la collection planes d'objet de type Plane.

- le fichier `planes.jl` définit le type Plane (non modifiable) et ses méthodes
  auxiliaires .

- les fichiers de la forme `yyyy_solver.jl` sont spécialisés dans la méthode 
  de résolution `yyyy` (stupid, carlo, greedy, descent, grasp, ...).
  Les constructeurs prennent généralement une instance en paramètre (variable 
  `inst` de type `Instance`).
  Le résultat d'un appel à la méthode `solve!(...)` de ces solveurs positin
  un attribut `bestsol` de type `Solution`.

- certains solveurs sont de la forme `zzz_timing_solver.jl` car ils sont
  spécialisés pour la résolution du Sous-Problème de Timing. Leur constructeur
  prend donc en paramètre une solution et non pas une instance.
  Leur effet est de modifier directement la solution associée à ce solveur.

- le fichier `solution.jl` définit tout ce qui est nécessaire et suffisant pour
  décrire une solution complète. L'attribut principal est un vecteur d'avion
  ainsi qu'un vecteur `x` leur affectant une date d'atterrissage.
  La solution pilote elle-même son TimingSolveur par sa propre méthode `solve!()`.
  Elle propose en outre un certain nombre d'opérateurs de voisinage
  paramétrables : à vous de les exploiter judiscieusement.

- le fichier `args.jl` définit le module Args qui est responsable de l'analyse
  de la ligne de commande, de l'affectation de valeur par défaut aux paramètres
  du programme et de la vérification de ces paramètres. Il propose
  principalement les méthodes `get`, `set` pour l'accès aux paramètres et
  `show_args` pour afficher leur valeur.

D'autres fichiers existent qui sont indépendants de l'application Seqata mais en
allège le codage. Par exemple :

- `lp_model_util.jl` : fourni la méthode `new_lp_model()` qui crée et configure un
   modèle PL vide avec choix entre CPLEX, Clp ou GLPK en fonction de l'option
   `--external-lp-solver` (alias `-x`pour les fainéants ;-).
  Par défaut le solveur utilisé est `-x cplex` mais les options `-x clp` et `-x glpk`
  fonctionnent également.

- `log_util.jl` : définit le module `Log` chargé de la gestion de l'affichage pour
  déboguer (i.e. les loggins). Le niveau de verbosité dépend de l'option
  `--level` alias `-L`.
  Ce fichier fournit les méthodes `lg1`, `lg2`, ..., et leur pendant avec saut 
  de ligne `ln1`, `ln2`, ... (voir doc dans ce fichier).

- `time_util.jl` : fournit les méthode de paramétrage `ms`, `ms_reset`, et la 
  macro `@ms`.

- `console_util.jl` : quelques fonctions persos pour écrire en couleur ou 
  dessiner un graphe ou un circuit dans une console (le projet Recytom vous dit
  quelques chose ? ;-).

Enfin, les fichiers `Seqata.jl` et `Seqata_using.jl` déclarent le module principal 
Seqata et charge les fichiers nécessaires.

#### Dans le sous-répertoire docs

**À DÉTAILLER**

- présentation du système de génération de la doc automatique
- les fichiers sources de la documentation en markdown (dans `docs/src`)
- ...



### Comprendre l'application

Ce paragraphe fournit quelques conseils pour comprendre et s'approprier 
l'application complète de façon à être capable de ;a modifier la compléter... ou la refondre !

Commencer par lire le fichier `01_presentation_proto_seqata.md` et par tester
le code. Exécuter les exemples proposés pour en comprendre l'utilité.

L'idée est de comprendre le déroulement de quelques commandes :

Observer l'affichage de :

    ./bin/run.jl carlo data/01.alp  -t earliest -n 10000000 -L1


Vous pouvez alors continuer par :
- survoler le fichier `./bin/run.jl` (sans trop approfondir) qui :
  - inclue le module `src/Seqata.jl` (convention de nommage d'un package 
    pour Julia),
  - analyse les arguments de la ligne de commande,
  - finaliser l'initialisation du module Seqata,
  - appelle la fonction `main()`.
- lire le fichier `main.jl` qui :
  - détecte que l'action courante est le symbole `:carlo`,
  - appelle alors la méthode `main_carlo()` qui correspond à la fonction principale
    pour cette action.
- lire le fichier `main_carlo.jl` qui :
  - traite spécifiquement la résolution par la néthode de Monté-Carlo,
  - construit l'instance en fonction du nom du fichier passé en paramètre,
  - définit le nombre d'iétarations à réaliser,
  - initialise les objects Solution (`cursol` et `bestsol`),
  - instancie un objet `sv` du type `StupidSolver`,
  - lui demande de résoudre `itermax` itérations.
- lire le fichier `stupid_solver.jl`.
  Vous aurez alors une idée de quelques fonctionnalités des "classes" `Instance` et 
  `Solution` et de leur utilisation.
- survoler les fichiers `instance.jl` et `solution.jl` pour avoir un panorama 
  des méthodes proposées.
- regarder le fichier `plane.jl` pour vous familiariser avec les attributs et le 
  méthodes associées au type Plane.

### Créer un solveur STP

Ensuite vous aurez à coder la résolution du Sous-Problème de Timing (STP). Pour
cela vous devrez (par exemple) coder un type `LpTimingSolver` sur le modèle du
type `EarliestTimingSolver` fourni dans le prototype.

Dans l'exemple précédent utilisant le `StupidSolver`, le STP était résolu
implicipement par un appel à la méthode `solve!(cursol)`. `cursol` étant l'objet
représentatant la solution courante, c'est donc bien la classe `Solution` qui
est responsable de la résolution du STP et qui choisit le type du
`XxxTimingSolver` à utiliser au moment de l'appel à `solve!(cursol)` (grace à
l'option `--xsolver`).

Examiner le fichier `solution.jl` pour comprendre comment un `XxxTimingSolver` est
sélectionné et examiner le fichier `EarliestTimingSolver`.

Un squelette est proposé pour le `LpTimingSolver` et vous disposer d'un exemple 
d'utilisation de JuMP dans le fichier `mip_discret_solver.jl` (testable par 
l'action `dmip`).

Le prototype du projet Seqata contient pas mal de méthodes qui simplifient
l'implantation de nouveaux solveurs. Mais du coup, une difficulté que vous pourrez
rencontrer si vous n'êtes pas familier avec Julia sera de devoir deviner si la
ligne que vous lisez utilise une méthode du langage Julia ou si c'est une méthode
spécifique au projet Seqata. Un bon indice est de vérifier si un des paramètres 
de cette méthode est d'un type spécifique au projet (Instance, XxxSolver, Solution, 
Plane, ...), Mais attention cependant : une variable telle que `inst.planes` est de
type Vector et donc `sort!(inst.planes, ...)` sera une méthode standard
documentée dans la doc officielle de Julia :  
{https://docs.julialang.org/en/v1/base/sort/#Base.sort!}(https://docs.julialang.org/en/v1/base/sort/#Base.sort!)

*A contrario*, une méthode comme `sort!(sol)` qui va trie la solution selon la date
d'atterrissage des avions est propre au projet et sera à rechercher dans le fichier `solution.jl`

### Conventions de codage 

Essayez de maintenir une convention de codage relativement cohérente, garantissant
la maintenabilité du code au sein du groupe.
Dans votre cas il s'agit essentiellement de maintenir le code bien présenté et
de réfléchir pour chaque fonction créée à un nom parlant et à une signature
judicieuse.


Le **nom des types** commence par une majuscule et le nom des variables sont 
en minuscules,

Une **collection** de quelques choses se termine toujours par un `s`, que ce soit 
un tableau (Array) ou un dictionnaire (Dict). Au besoin, un liste de vecteurs d'avions
pourrait s'appeler `planess` (avec deux `s` !).

Le **nom d'une méthode** commence par un minuscule (sauf pour un constructeur).
Comme en Ruby, Julia adopte la convention d'ajouter en fin de nom un **point
d'exclamation "!"** pour indiquer que **cette méthode modifie un de ses
paramètres**.
Par exemple la commande suivant complète le vecteur planes par un avion 
supplémentaire :

    push!(planes, plane)

Contrairement à un langage objet classique, Julia supporte le **multiple-dispatch**.
Une conséquence est que les méthodes ne font pas partie d'un type. Cependant dans 
beaucoup de cas, une méthode est naturellement associée à un type particulier.
Dans ce cas le type et "ses" méthodes seront définis dans le même fichier.
C'est en particulier pour les types liés au problème lui-même : les types Instance, 
Plane et Solution.

Par exemples la méthode suivante qui retourne la durée de séparation entre les
objets avions p1 et p2, sera intégrée dans le fichier `Instance.jl` comme si
c'était une méthode de la classe Instance.

    get_sep(inst::Instance, p1::Plane, p2::Plane)

Par ailleurs, la surcharge permet d'utiliser une convention de nommage uniforme
selon la fonctionnalité de la méthode. Par exemple de nombreuses methodes
`solve!` existent dans Seqata qui sont chacune associées à un type de solveur
différent. L'objet XxxSolver est alors passé en premier paramètre
(GreedySolveur, StupidSolver, ...).

Enfin, voici quelques noms de variables locales couramment utilisées dans le
projet Seqata :

- Instance => inst
- Solution => sol, sol1, sol2, bestsol, cursol, testsol...
- XxxxSolver => solver ou sv
- Pour le grasp qui va manipuler à la fois des solveurs esclaves GreedySolver et
  DescentSolver, vour pouvez par exemple utiliser `gsv` et `dsv`.
- Plane => p, p1, p2 ou ou plane
- i1, i2 ou idx1, idx2 : indice d; avions dans une liste

### Conseils de développement

Si vous souhaitez compléter le proto fourni plutôt que de simplement vous en
inspirer, vous aurez intérêt à le comprendre suffisamment pour vous
l'approprier avant de le modifier (mais inutile de rentrer dans le détail 
d'imlantation du fichier `args.jl` par exemple !).

Vous devrez également être efficace pour valider un cycle de développement
complet "modif + test + archivage" et garantir que le code reste toujours
fonctionnel et pouvoir revenir à une version antérieure en cas de problème.

Pour finir, Julia est un langage très puissant avec de nombreuses possibilités.
Il est donc compliqué si vous voulez aborder toutes ses possibilités. Cependant 
le prototype fourni contient toutes les possibilités et syntaxes qui vous seront 
nécessaires.


./

