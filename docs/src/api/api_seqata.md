# API du module Seqata
```@meta
CurrentModule = Seqata
```


```@docs
Seqata
```


## Méthodes générales


```@docs
Seqata.main()
Seqata.main_carlo()
Seqata.main_descent()
Seqata.main_explore()
Seqata.main_stats()
Seqata.main_test()
Seqata.main_timing()
Seqata.main_validate()
```


## Type Plane


```@docs
Plane
get_cost(p::Plane, t::Int; violcost::Float64 = 1000.0)
to_s_alp(p::Plane)
```

## Type Instance

```@docs
Seqata.Instance
```

## Type Solution

```@docs
Seqata.Solution
```

## Les solveurs STP

Voici quelques solveurs de résolution du sous-problème de timing (STP).
Il ne sont pas forcément disponible dans le prototype fourni aux élèves.

```@docs
Seqata.EarliestTimingSolver
Seqata.LpTimingSolver
```
## Les solveurs globaux

Les solveurs globaux ont pour vocation de résoudre le problème dans son 
ensemble, que se soit par une approche frontale (mathématique ou non) 
ou par une décomposition (par niveaux avec un brique STP ou par 
tranches temporelles).

```@docs
Seqata.ExploreSolver
Seqata.DescentSolver
Seqata.SteepestSolver
Seqata.MipDiscretSolver
Seqata.MipSolver
```
## Méthodes utilitaires

Ces méthode facilite de chronométrage de morceaux de code.
```@docs
Seqata.ms
Seqata.ms_reset
Seqata.@ms
```

À suivre !
