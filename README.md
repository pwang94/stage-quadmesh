# Ultimaille examples

Here you'll find a collection of examples of how to use the ultimaille library.

## Build

```sh
cmake -B build && cd build && make -j
```

## Run 

In the build directory type the following command:

```sh
examples/hello
```

## Pipeline

chart permet de générer la première partition en  taches

framefield calcule le framefield

tache_singu ajoute des sources aux singu du framefield

recentre permet de recentrer les taches, d'en ajouter ou d'en supprimer si il en manque. ça ne fait qu'une itération il faut relancer plusieur fois en changeant peut être des parametre

suppr_non_disque détecte les taches qui ne forment pas un disque topologique et rajoute des sources proches

poly_to_quad renvoie un maillage polygonal où les singus ont été découpé

annot_poly annote le maillage pour transformer le t-mesh en maillage quad.