
#include "helpers.h"
#include <ultimaille/all.h>
#include <map>
#include <set>

using namespace UM;

bool deuxcommun(std::vector<int> v1, std::vector<int> v2) {
    int compte = 0;
    for (int c : v2) {
        if (std::find(v1.begin(), v1.end(), c) != v1.end()) {compte++;}
    }
    return (compte == 2);
}

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

    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Triangles m;
    Triangles m_bord;
    static const std::string entree = "framefield.geogram";
    static const std::string entree_bord = "sortie_chart.geogram";
    static const std::string sortie = "sortie.geogram";
    SurfaceAttributes attrs = read_by_extension(entree, m);
    PointAttribute<int> singularite("singularite", attrs, m);
    SurfaceAttributes attrs_bord = read_by_extension(entree_bord, m_bord);
    FacetAttribute<bool> seed_bord("seed", attrs_bord, m);

    m.connect();

    // les hard edges
    CornerAttribute<bool> hard_edges_attr(m, false);
    find_hard_edges(m, hard_edges_attr, 0.6);
    std::vector<int> a_traite;
    std::vector<int> a_traite2;
    int compte_face = 0;
    FacetAttribute<int> zone(m, -1);
    int zone_num = 0;
    std::vector<std::vector<int>> de_zone; // a chaque zone donne les faces dedans
    // obtenir les grosses zones engendrées par les hards edges

    while (compte_face < m.nfacets()) { 
        int frand;
        de_zone.push_back({});
        do  {
            frand = rand() % m.nfacets();
        } while (zone[frand] != -1);
        // on propage sur toute la zone
        a_traite.push_back(frand);
        while (a_traite.size() != 0) {
            for (int fid: a_traite) {
                if (zone[fid] != -1) continue;
                compte_face++;
                zone[fid] = zone_num;
                de_zone[zone_num].push_back(fid);
                Surface::Facet f(m, fid);
                for (auto h : f.iter_halfedges()) {
                    int f_opp = h.opposite().facet();
                    if (!hard_edges_attr[h] && zone[f_opp] == -1) {
                        a_traite2.push_back(f_opp);
                    }
                }
            }
            a_traite = a_traite2;
            a_traite2.clear();
        }
        zone_num++;
    }

    a_traite.clear();
    a_traite2.clear();




    // on stock les coins de chaques faces 
    
    double tolerance = 0.3;
    PointAttribute coins_attr(m, 0);
    FacetAttribute<bool> coins_attr_facets(m, 0);
    FacetAttribute<int> source(m, -1);
    FacetAttribute<int> racine(m, -1);
    
    FacetAttribute<double> dist(m, 1000);
    int compte_source = 0;

    for (auto v: m.iter_vertices()) {
        //si il y a 3 facets où c'est 2 à 2 de proche de 90 deg
        for (auto h1: v.iter_halfedges()) {
            for (auto h2: v.iter_halfedges()) {
                for (auto h3: v.iter_halfedges()) {
                    auto f1 = h1.facet();
                    auto f2 = h2.facet();
                    auto f3 = h3.facet();
                    vec3 normalA = f1.geom<Triangle3>().normal();
                    vec3 normalB = f2.geom<Triangle3>().normal();
                    vec3 normalC = f3.geom<Triangle3>().normal();
                    if (normalA * normalB < tolerance && normalB * normalC < tolerance && normalC * normalA < tolerance) {
                        //coins[zone[f1]].push_back(v);
                        //on met les coins dans les source
                        //on veut éviter les coins concaves
                        if (coins_attr[v] != 1) {
                            source[f1] = compte_source;
                            racine[f1] = compte_source;
                            dist[f1] = 0;
                            compte_source++;
                            a_traite.push_back(f1);
                            source[f2] = compte_source;
                            racine[f2] = compte_source;
                            dist[f2] = 0;
                            compte_source++;
                            a_traite.push_back(f2);
                            source[f3] = compte_source;
                            racine[f3] = compte_source;
                            dist[f3] = 0;
                            compte_source++;
                            a_traite.push_back(f3);



                            coins_attr_facets[f3] = true;
                            coins_attr_facets[f2] = true;
                            coins_attr_facets[f1] = true;
                            coins_attr[v] = 1;
                        }

                        
                    }
                }
            }
        }
    }


    //seed

    for (auto v: m.iter_vertices()) {
        if (singularite[v] != 0) {
            a_traite.push_back(v.halfedge().facet());
            source[v.halfedge().facet()] = compte_source;
            racine[v.halfedge().facet()] = compte_source;
            dist[v.halfedge().facet()] = 0;
            compte_source++;
        }
    }



    // peut etre sampler toute la surface en utilisant les 1ere tache pour faire des trucs

    //on autorise les nouvelles source loins des anciennes:

    FacetAttribute<int> dist_interdit(m, -1);
    std::vector<int> a_lire = a_traite;
    std::vector<int> a_lire2;
    const int dist_max = 25;

    do {
        for (auto fid: a_lire) { // f à distance i
            Surface::Facet f(m, fid);
            for (auto h: f.iter_halfedges()) {
                if (hard_edges_attr[h]) {continue;} //suit les features
                auto f2 = h.opposite().facet();
                if (dist_interdit[f2] == -1 && dist_interdit[f] < dist_max) {
                    dist_interdit[f2] = dist_interdit[f] + 1;
                    a_lire2.push_back(f2);
                }
            }
        }
        a_lire = a_lire2;
        a_lire2.clear();
    } while (a_lire.size() != 0);
    
    // on lit les nouvelles sources de bord.cpp

    for (auto f: m_bord.iter_facets()) {
        if (seed_bord[f] && dist_interdit[f] == -1) {

            a_traite.push_back(f);
            source[f] = compte_source;
            racine[f] = compte_source;
            dist[f] = 0;
            compte_source++;
        }

    }
    //charts gen





    do {
        for (auto fid: a_traite) { // f à distance i
            Surface::Facet f(m, fid);
            for (auto h: f.iter_halfedges()) {
                if (hard_edges_attr[h]) {continue;} //suit les features
                auto f2 = h.opposite().facet();
                if (dist[f2] > dist[f] + 1) {
                    dist[f2] = dist[f] + 1;
                    a_traite2.push_back(f2);
                    source[f2] = source[f];
                }
            }
        }
        a_traite = a_traite2;
        a_traite2.clear();
    } while (a_traite.size() != 0);

    





    //si faces vides
    std::map<int, int> num;
    for (auto f: m.iter_facets()) {
        if (source[f] == -1) {
            if (!(num.find(zone[f]) !=num.end())) {num[zone[f]] = compte_source; compte_source++;}
            source[f] = num[zone[f]];
        }
    }


    write_by_extension(sortie, m, {{{"singularite", singularite.ptr}}, {{"source", source.ptr}, {"dist", dist.ptr}, {"racine", racine.ptr}, {"coins", coins_attr_facets.ptr}}, {{"hard_edge", hard_edges_attr.ptr}}});

    int result = system((getGraphitePath() + " " + sortie).c_str());
    // --- END ---

    return 0;
}
