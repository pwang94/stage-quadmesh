//on prend le centre de la source et on repropage


#include "helpers.h"
#include <ultimaille/all.h>
#include <map>
#include <set>

using namespace UM;

void find_hard_edges(Triangles& mesh, CornerAttribute<bool>& hard_edges_attr, double threshold) {

    // Iter on all mesh halfedges
    for (auto h : mesh.iter_halfedges()) {
        
        // Get opposite halfedge
        auto opposite = h.opposite();
        
        // if (!opposite.active())
        //     continue;

        // Compute normals of a face and its opposite face
        vec3 normalA = h.facet().geom<Triangle3>().normal();
        vec3 normalB = opposite.facet().geom<Triangle3>().normal();
        // Compute the dot product of normals
        double d = normalA * normalB;

        // The closer the scalar product of its normals is to zero, 
        // the closer the angle between two faces is to 90 degrees
        if (std::abs(d) <= threshold) {
            hard_edges_attr[h] = true;
        }
    }
}


int main(int argc, char** argv) {
    // --- LOAD ---

    //parametre

    
    double dist_max = 0.195; //pour supprimer des sources
    int dist_max1 = 25;  // aussi mais la distance est le nombre d'arête sur le chemin
    double dist_min = 0.05; // pour ajouter des sources
    
    double eps =1e-6;
    
    // TODO modifier les paramètres avec argv
    bool display = true;
    // bool ajoute_face = false;
    // const bool choix_supprime = false;
    const int choix_norme = 1; // 1: l_inf 0: l_2  2:l_1

    bool on_previous_out = argc >= 2 && std::string(argv[1]) == "previous";
    bool ajoute_face = argc >= 2 && std::string(argv[1]) == "add_face" || argc >= 3 && std::string(argv[2]) == "add_face" ;
    bool choix_supprime = argc >= 2 && std::string(argv[1]) == "delete_face" || argc >= 3 && std::string(argv[2]) == "delete_face" ;

    std::string path = getAssetPath();

    Triangles m; // charts
    Triangles mff; //framefield
    Triangles m_bord;
    // static const std::string entree = ; //si on prend la sortie de tache_singu
    std::string entree = on_previous_out ? "sortie_recentre.geogram" : "sortie.geogram" ; // si on veut relancer sur la nouvelle sortie
    static const std::string entree_bord = "sortie_chart.geogram";
    static const std::string sortie = "sortie_recentre.geogram";
    static const std::string ff_path = "framefield.geogram";
    SurfaceAttributes attrs = read_by_extension(entree, m);
    SurfaceAttributes attrs_framefield = read_by_extension( ff_path, mff);
    SurfaceAttributes attrs_bord = read_by_extension( entree_bord, m_bord);
    PointAttribute<int> coins_attr("coins_attr", attrs_bord, m);
    FacetAttribute<int> coins_concave("coins concave", attrs_bord, m);
    m.connect();
    m_bord.connect();

    FacetAttribute<vec3> frames0("frames0", attrs_framefield, mff);
    FacetAttribute<vec3> frames1("frames1", attrs_framefield, mff);

    FacetAttribute<bool> coins("coins", attrs, m);
    PointAttribute<int> singularite("singularite", attrs, m);
    FacetAttribute<int> source0("source", attrs, m);
    FacetAttribute<int> racine("racine", attrs, m);
    FacetAttribute<double> dist0("dist", attrs, m);

    CornerAttribute<bool> hard_edges_attr(m, false);
    find_hard_edges(m, hard_edges_attr, 0.6);
    FacetAttribute<int> source(m, -1);
    FacetAttribute<double> dist(m, 1000);
    FacetAttribute<int> new_racine(m, -1);

    //TODO si j'ai des triangles je veux les supprimer




    std::map<int, int> singu_of_source; //source -> singu
    std::map<int, int> racine_of_source; //source -> racine

    for (auto v: m.iter_vertices()) {
        if (singularite[v] != 0) {
            singu_of_source[source0[v.halfedge().facet()]] = v.halfedge().facet();
        }
    }

    for (auto f: m.iter_facets()) {
        if (racine[f] != -1) {
            racine_of_source[source0[f]] = f;
        }
    }

    int nb_source = racine_of_source.size();

    // for (auto cle: racine_of_source) {
    //     std::cout << cle.first << std::endl;
    // }
    // std::cout << "passe" <<std::endl;
    std::cout << "pase" << std::endl;

    std::vector<int> nb_cote(nb_source + 100);
    for (auto v: m.iter_vertices()) {
        // on note a chaque sources son nombre de sommets
        int compte = 0;
        std::vector<int> couleurs;
        for (auto h: v.halfedge().iter_sector_halfedges()) {
            int current_color = source0[h.facet()];
            if (std::find(couleurs.begin(), couleurs.end(), current_color) == couleurs.end()) {
                compte++;
                couleurs.push_back(current_color);
            }
        }
        if (compte >= 3) {
            for (int color: couleurs) {
                nb_cote[color]++;
            }
        }
    }
    // std::cout << "passe" <<std::endl;

    for (int i = 0 ; i< nb_source; i++) {
        if (nb_cote[i] == 3) { std::cout << i << "triangle" << std::endl;}

    } 
    // ça peut dépendre de la taille de la surface, de la face 
    // on peut enlever ou ajouter des sources si ça s'éloigne trop de l'écart moyen de la face TODO
    std::cout << "pase" << std::endl;
    

    std::vector<vec3> centres; 
    //au lieu de prendre le barycentre prendre la moyenne de toute les faces de la zone centre_tri[source] : [somme, nb point]
    std::vector<std::pair<vec3, int>> centre_tri(nb_source + 300, {vec3(0, 0, 0), 0}); //TODO nb_source

    for (auto f: m.iter_facets()) {
        if (nb_cote[source0[f]] == 3) continue;
        // if (nb_cote[i] == 4) {
        //     //si jamais ya un angle plat de plat c'est un triangle

        // }
        auto& [v, n] = centre_tri[source0[f]];
        v += f.geom<Poly3>().bary_verts();
        n++;
    }


    std::cout << "pase" << std::endl;
    // si il y a une tache qui fait un triangle on la supprime

    for (auto [v, n]: centre_tri) {
        if (n != 0) {centres.push_back(v/n);};
    }

    // si les sources sont trop loin d'un pt on en rajoute


    FacetAttribute<bool> trop_pres(m);

    if (ajoute_face) {
            std::vector<int> a_suppr, a_suppr2;

        for (auto f : m.iter_facets()) {
            if (!trop_pres[f] && dist0[f] > dist_max &&
                singu_of_source.find(source0[f]) == singu_of_source.end()) {

                nb_source++;
                centres.push_back(f.geom<Poly3>().bary_verts());
                trop_pres[f] = true;
                a_suppr = {f};
                for (int i = 0; i < dist_max1; ++i) {
                    a_suppr2.clear();
                    for (auto fid: a_suppr) { // f à distance i
                        Surface::Facet f(m, fid);
                        for (auto h: f.iter_halfedges()) {
                            if (hard_edges_attr[h]) continue; //suit les features
                            auto f2 = h.opposite().facet();
                            if (!trop_pres[f2]) {
                                a_suppr2.push_back(f2);
                                trop_pres[f2] = true;
                            }
                        }
                    }
                    a_suppr = a_suppr2;
                }
            }
        }

    }



    BVHTriangles bvh(m);

    std::vector<int> a_traite;
    std::vector<int> a_traite2;
    int compte_source = 0;



    for (vec3& v: centres) {
        auto new_f = bvh.nearest_point(v).f;
        // Surface::Facet f(m, racine_of_source[source0[new_f]]); //une seule source entre les 2
        // new_f = bvh.nearest_point(new_f.geom<Poly3>().bary_verts() + f.geom<Poly3>().bary_verts()).f;
        if (singu_of_source.find(source0[new_f]) != singu_of_source.end()) new_f = singu_of_source[source0[new_f]]; //on bouge pas la source si c'est une singu
        //Surface::Facet facet_source(m, racine_of_source[source0[new_f]]);
        // if (coins[facet_source]) new_f = racine_of_source[source0[new_f]]; //on bouge pas la source si c'est un coin
        if (dist[new_f] == 0) continue;
        a_traite.push_back(new_f);
        new_racine[new_f] = compte_source;
        dist[new_f] = 0;
        source[new_f] = compte_source;
        // if (racine_of_source[source0[new_f]] != new_f) {compte_source++;}
        racine_of_source[compte_source]= new_f;
        // source[racine_of_source[source0[new_f]]] = compte_source;
        compte_source++;
    }



    FacetAttribute<bool> trop_proche(m);
    std::vector<bool> a_garder(nb_source);
    std::vector<bool> a_jeter(nb_source);
    //chart_gen
    do {
        for (auto fid: a_traite) { // f à distance i
            Surface::Facet f(m, fid);
            Surface::Facet f_source(m, racine_of_source[source[f]]); 
            for (auto h: f.iter_halfedges()) {
                if (hard_edges_attr[h]) {continue;} //suit les features
                
                auto f2 = h.opposite().facet();

                //on calcule la distance linf suivant le framefield de f2 à la source inf TODO que si la source est pas une singu
                vec3 p1 = f_source.geom<Poly3>().bary_verts(), p2 = f2.geom<Poly3>().bary_verts();
                vec3 v =  p1 - p2;
                //si c'est une singu on prend juste toutes les frames qui pointent vers elle
                double distance;
                // augmenter un peu la distance pour les trucs qui sont loins mais dont la proj est proche  
                double correction = abs(f_source.geom<Poly3>().normal() * v);
                switch (choix_norme) {
                    case 1:
                        distance = std::max(abs(v*frames0[f_source]), abs(v*frames1[f_source])) +correction * 1 /*+  0.01 *v.norm() */;
                        break;
                    case 2:
                        distance = abs(v*frames0[f_source]) + abs(v*frames1[f_source]) + correction * 1;
                        break;
                    default:
                        distance = v.norm();
                        break;
                }

                if (singu_of_source.find(source0[f_source]) != singu_of_source.end()) {
                    distance = 0.7 * v.norm();
                }
                // on peut vérifier ici si ca touche une autre face trop proche
                if (choix_supprime) {
                    if (source[f2] != - 1 && source[f2] != source[f] && distance < dist_min) {
                        //on garde la plus petite des deux
                        int f_min = std::min(f, f2), f_max = std::max(f, f2);;
                        trop_proche[f_max] = true;

                        if (!a_garder[source[f_max]] && !a_jeter[source[f_min]] && singu_of_source.find(source0[racine_of_source[source[f_min]]]) == singu_of_source.end() && !a_jeter[source[f_max]]) {
                            a_jeter[source[f_max]] = true;
                            std::cout << "passage0" << std::endl;
                            //on supprime source[f2] et on met tout couleur source[f1]
                            a_garder[source[f]] = true;
                            int to_kill = source[f_max];
                            std::cout << source[f] <<std::endl;
                            for (auto f0: m.iter_facets()) {
                                if (source[f0] == to_kill) {source[f0] = source[f_min];}
                                
                            }
                            new_racine[racine_of_source[source[f_max]]] = -1;
                        }
                    }
                }
                if (dist[f2] > distance) {
                    dist[f2] = distance;
                    a_traite2.push_back(f2);
                    source[f2] = source[f];
                }
            }
        }
        a_traite = a_traite2;
        a_traite2.clear();
    } while (a_traite.size() != 0);






    write_by_extension(sortie, m, {{{"singularite", singularite.ptr}, {"coins_attr", coins_attr.ptr}}, {{"source", source.ptr}, {"racine", new_racine.ptr}, {"dist", dist.ptr}, {"porche", trop_proche.ptr}, {"coins concave",coins_concave.ptr}}, {}});

    if (display) int result = system((getGraphitePath() + " " + sortie).c_str());
    // --- END ---

    return 0;
}
