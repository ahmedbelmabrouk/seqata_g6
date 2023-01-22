# Installations de Julia et des solveurs

### Installation de Julia 

#### Téléchargement de Julia

Voir la page download du site officiel : <https://julialang.org/downloads/>

#### Les packages Julia du projet


Pour la gestion des packages nécessaires à une application, Julia propose la notion de projet qui facilite la portabilité d'un code Julia en garantissant
l'utilisation des bonnes versions de package. 
L'ensemble des packages et leur version est inscrit dans le fichier `Project.toml`.

Le commandes suivantes permettent de chager les bons packages dans votre distribution Julia.

Soit en interactif 
```
julia
]           # <= on passe en mode package 
activate .  # on active le projet correspondant au répertoire courant
instantiate # on insatlle si nécessaire les package du fichier `Project.toml`
```

Soit directement dans le code de votre application (c'est le cas dans Seqata).
Dans cet exemple on y ajoute une fonctionnalité de chronométrage du processus.

```
using Dates
t1 = Dates.now()
using Pkg
Pkg.activate(appdir)
Pkg.instantiate()
println(Dates.now()-t1)
```

### Installation de CPLEX

Le packahe `CPLEX.jl` de julia permet d'appeler depuis Julia des fonctionnaltés de CPLEX supposé être préinstallé.

Il faut donc commencer par installer CPLEX (gratuit pour l'enseignement).

Soit vous vous iscriver sur le site IBM, et charger le fichier à télécharger à partir de 

<https://www.ibm.com/academic/technology/data-science>

Soit vous récupérer directement l'arche de la dernière version sur le réseau de l'ENSTA.

```ls -al 
/usr/uma/pack/cplex/LOCAL/src/download_ibm_cplex2210/
```

Vérifier l'installation de cplex en tapant `oplrun`
```
oplrun -version
=> Version IBM ILOG CPLEX Optimization Studio 22.1.0.0
```

Une fois installer CPLEX, vous devrez peut-être reconstruire le package
Julia `CPLEX.jl` pour qu'il retrouve l'installation 

En mode interactif :
```
Julia       # lance julia depuis le shell
]           # On passe en mode package
build CPLEX # recontruit le binding pour CPLEX
<backspace> # retour à julia
<Ctl-D>     # retour au shell
```
### Mise à jour du PATH

Pour accéder aux exécutables depuis la console, il faut mettre à jour la 
variable d'environnement `PATH`.

**Sous unix**
: À détailler

**Sous Macos**
: À détailler

**Sous Windows**
: À détailler
