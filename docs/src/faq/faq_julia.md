# Quelques constructions Julia (howto)

### Mélange partiel d'un tableau

    # transformation d'un intervale (implicite) en un Array (explicite)
    a = collect(1:20)
    # Modification d'une petite tranche du tableau
    a[10:15] = shuffle(a[10:15])

### Tirage aléatoire pondéré

Existe-t-il une méthode julia pour faire un choix aléatoire dans un ensemble pour lequel on donne la probabilité de tirer chaque nombre (par exemple on veut tirer un nombre entre 1 et 4 avec des probabilités [0.2,0.3,0.4,0.1]) ?

Voir fonction `sample` de `StatsBase`.

- <https://juliastats.org/StatsBase.jl/stable/sampling/>
- <https://stackoverflow.com/questions/27559958/how-do-i-select-a-random-item-from-a-weighted-array-in-julia>


```
Pkg.add("StatsBase")  # Only do this once, obviously
using StatsBase
items = ["a", 2, 5, "h", "hello", 3]
weights = [0.1, 0.1, 0.2, 0.2, 0.1, 0.3]

Echantillonnage d'un seul élèment
my_samp = sample(items, Weights(weights))

Echantillonnage de plusieurs élèments AVEC remplacement
my_samps = sample(items, Weights(weights), 10)

Echantillonnage de plusieurs élèments SANS remplacement
my_samps = sample(items, Weights(weights), 2, replace=false)
```

### Tri de vecteurs d'objets

Vous pouvez si nécessaire trier un vecteur planes en une ligne (tri direct ou reverse) sur une fonction arbitraire. Voici un exemple de tri (dans la "classe" Solution) utilisant une fonction anonyme passée au paramètre `by`  (la version de sort avec le "!" modifie le vecteur en place, sinon un nouveau vecteur est créé) :

```
Base.sort!(sol.planes, by = p -> p.lb, rev = true)
```

On peut trier sur plusieurs valeurs (cet exemple n'est pas dans Seqata !) si la fonction passés au paramètre `by` retourne un tuple

```
tmp=sort(someperfs, by=r->(r.nb_planes, r.nb_breakpoints, r.nb_calls_per_sec))
```

Dites-vous aussi que le vecteur `planes` ne prend pas plus de place qu'un vecteur `ubs` contenant les `p.ub` car il ne contient que des références 
vers des object de type Plane.


### Différence entre using et import

- `using` rend les méthodes directement disponible
- `import` rend les méthodes accessibles via leur préfixe
- `import` est nécessaire pour surcharger une méthode du module Base 
   (e.g Base.show(xxx))

```
    import PyPlot; const plt = PyPlot
    # Ceci évite de polluer l'espace de nom locaux en imposant le préfixe plt :
    function plot_costs(costs, path::String)
        println(path)
        plt.plot(costs)
        plt.xlabel("Iteration")
        plt.ylabel("Cost")
        plt.savefig("$path.pdf", bbox_inches="tight")
        plt.close()
    end
```
