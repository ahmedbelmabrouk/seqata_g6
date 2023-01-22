# Utilisation du code en mode interactif

### Principe

Le lancement de l'application Seqata est plutot lent, ce qui peut s'avérer 
génant pendant le développement lorsqu'on est dans une phase avec de nombreux 
cycles "correction mineure <--> test".

Un remède classique consiste à développer l'application en mode interactif.
Cela consiste à lancer Julia en spécifiant que l'on ne veux par quitter même
quand le fichier spécifié a été exécuté.

    julia -iL ./bin/run.jl

L'application Seqata détecte qu'elle est lancée en mode interactif et charge
le maximum de packages possibles. Elle déclare également quelques alias 
pour des méthodes d'utilisation fréquentes (cf fichier src/interactive.jl)

Exemple de session : on veut mettre au point la résolution par l'action
"explore" en utilisation l'algorithme de timing "earliest" déjà fonctionnel
dans le proto.

On lance Julia en mode interactif :

    julia -iL ./bin/run.jl

On tape alors directement du julia dans la console :

### Exemples interactifs

#### Exemple 1 : exécution de l'action explore en interactif


    # on déclare le fichier d'instance à lire (Args.set(:infile))
    Args.set(:infile, "data/01.alp")
    ag(:infile, "data/01.alp") # Pour les fainéants !
    ag(:infile, p01)           # Pour les gros fainéants !!

    # On rend le solver plus verbeux pour augmenter les infos affichées
    Log.level(4)   # ou bien lv(4) pour les fainéants
    Log.level()   # retourne la valeur courante du level ; alias lv() pour les fainéants

    # On veut imposer le timing_solver :earliest au lieu de :lp 
    # (car :lp n'est pas encore fait !)
    # Zut quelle est le nom de cette option ?
    Args.show_args()  # alias ag()
    # Ah oui, c'est timing_algo_solver

    Args.set(:timing_algo_solver, :earliest)
    Args.show_args()

    # On charge la méthode main correspondant à l'action :explore
    main_explore()

    # Bon tout cela est trop long. Voyez le fichier "src/interactive.jl" pour
    # connaitre les raccourcis disponibles et testez l'exemple suivant

#### Exemple 2 : la même chose en plus rapide

    julia -iL ./bin/run.jl
    ag(:infile, p01)
    ag(:level, 4)   # ou lv(4)
    ag()
    ag(:timing_algo_solver, :earliest)
    main_explore()

    # Vous pouvez alors modifier votre code et le relancer.
    # Mais un exit() dans le code vous obligera à relancer la session Julia !

#### Exemple 3 : exécution manuelle de l'exploration et introspection

    inst = Instance(p01)
    @show inst.planes
    @show inst.planes[1] 
    dump(inst.planes[1])

    # On va utiliser l'algo de timing "earliest"
    ag(:timing_algo_solver, :earliest)
    ag()

    # On crée le solveur d'exploration pour cette instance
    sv = ExploreSolver(inst)

    @show sv.bestsol # La solution initiale aléatoire est mauvaise !

    # Une première exploration de 1000 itérations
    solve!(sv, 1000)
    @show sv.bestsol # la solution est un peu meilleure

    # On veut avoir plus de détail pendant la résolution (option :level)
    shuffle!(sv.bestsol)
    lv(3)
    solve!(sv, 1000)

    # Aller, un petit million d'itérations supplémentaire
    lv(1) # on diminue la verbosité
    shuffle!(sv.bestsol)
    solve!(sv, 2_000_000) # la solution est encore meilleure

    # Ça nous plait : on enregistre cette colution dans un fichier
    write(sv.bestsol)

### Remarques sur le mode interactif

ATTENTION, si vous modifier le fichier (e.g descent_solver.jl), 
il faut le recharger par :

    include("src/descent_solver.jl")


ou plus simplement par ma macro `@i` : 
```
@i desc  # recharge tous les fichiers du répertoire `src` contenant la chaine "desc" !
```

Il existe un package `revise.jl` qui devrait le faire automatiquement... 
mais qui n'a pas été testé pour l'instant !

D'une manière générale :

- ne pas créer de méthodes dans la console interactive de Julia, mais créer 
  ou modifier les méthodes dans un fichier du répertoire `./src/`.
- En cas de problème avec une méthode qui ne se charge pas comme prévu, 
  n'hésitez pas à relancer julia pour repartir avec un état propre !
- pour valider le code final, retester l'action qui utilise les méthodes 
  modifiées.



