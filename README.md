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

recentre permet de recentrer les taches, d'en ajouter ou d'en supprimer si il en manque. ça ne fait qu'une itération il faut relancer plusieur fois en ajoutant après la première fois le paramètre previous
pour ajouter ou supprimer des faces trop proche on peut ajouter le paramètre add_face ou delete_face après previous

suppr_non_disque détecte les taches qui ne forment pas un disque topologique et rajoute des sources proches

tache_to_quad renvoie un maillage polygonal où les singus ont été découpé

annot_poly annote le maillage pour transformer le t-mesh en maillage quad.
