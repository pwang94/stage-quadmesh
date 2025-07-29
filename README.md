## Pipeline

chart permet de générer la première partition en taches.
Entrée: un maillage triangulaire
Sortie: sortie_chart.geogram le maillage avec les attributs des taches

framefield calcule le framefield
Entrée: sorti un maillage triangulaire
Sortie:framefield.geogram Le maillage avec les frames pour chaque face et framefieldpl.geogram un polyline pour l'affichage du framefield

tache_singu ajoute des sources aux singu du framefield
Entrée: sortie_chart.geogram Un maillage triangulaire avec des taches et framefield.geogram le framefield
Sortie: sortie.geogram Le maillage triangulaire avec des taches en plus et les singularités sont annotés

recentre permet de recentrer les taches, d'en ajouter ou d'en supprimer si il en manque. ça ne fait qu'une itération il faut relancer plusieur fois en ajoutant après la première fois le paramètre previous
pour ajouter ou supprimer des faces trop proche on peut ajouter le paramètre add_face ou delete_face après previous
Entrée: sortie.geogram ou sortie_recentre.geogram Un maillage triangulaire avec des taches et des singu et framfield.geogram
Sortie: sortie_recentre.geogram Un maillage triangulaire avec des taches et des singu


suppr_non_disque détecte les taches qui ne forment pas un disque topologique et rajoute des sources proches
Entrée: sortie_recentre.geogram Un maillage triangulaire
Sortie: sortie_disque.geogram

tache_to_quad renvoie un maillage polygonal où les singus ont été découpé
Entrée: sortie_disque.geogram un maillage triangulaire avec taches et singus
Sortie: un maillage polygonal


annot_poly annote à chaque angle, angle plat ou angle droit pour ensuite donner le maillage à un algorithme de quantization.
Entrée: outpoly.geogram un maillage polygonal
Sortie: poly_annot.geogram le maillage avec les angles annotés

