# Quelques astuces sur l'autodocumentation en Julia

### Doc de référence

- <https://juliadocs.github.io/Documenter.jl>


### Principe général

TODO...

### Quelques problèmes ou solutions

#### Insertion de html standard dans doc

On ne peut pas insérer de html natif dans ce markdown. Du coup on ne peut pas forcer un saut de ligne en insérant simplement un <br>, ou ajoute une div avec une classe personnalisée !

#### Écrire des math dans la doc de julia

Permet d'insérer des math, soit en ode inline soit en mode display. 
Quelques exemple de base :


Voici le coût de pénalité ``c_i`` de l'avion ``i`` qui atterrit à la date ``x_i`` :

```math
\begin{aligned}
  c_i(x_i) &= ep_i(T_i - x_i)^{+}  + tp_i(x_i -T_i)^{+} \\
           &  \text{avec } (X)^{+} = \max(0, X).
\end{aligned}
```


#### Créer des listes de descriptions HTML (<dl>) avec docs

Le syntaxe standard de markdown ne fonctionne pas.

```
**Sous unix**
: À détailler 1
: à détailler 2

**Sous Macos**
: À détailler

**Sous Windows**
: À détailler
```

Et comme on ne peut pas insérer de HTML natif dans ce markdown
```
<dl>
<dt>Sous unix  </dt>
<dd>À détailler</dd>
<dt>Sous Macos  </dt>
<dd>À détailler</dd>
<dt>Sous Windows  </dt>
<dd>À détailler</dd>
</dl>
````

Du coup je me comtente d'une version bidouille (proche de la version markdown 
citée en premier) qui conduit à :

**Sous unix**
: À détailler

**Sous Macos**
: À détailler

**Sous Windows**
: À détailler
