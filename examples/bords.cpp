
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
    // on met argv[0] en entrée TODO et faire un script shell qui fait tout bien
    static const std::string entree = "joint_sub_sub.geogram";
    static const std::string sortie = "sortie_bord.geogram";
    read_by_extension(path + entree, m);
    m.connect();

    // les hard edges
    CornerAttribute<bool> hard_edges_attr(m, false);
    find_hard_edges(m, hard_edges_attr, 0.6);
    FacetAttribute<int> dist2(m, -1);

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
                if (zone[fid] != -1) {continue;}
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


    //compute distance from bord

    for (auto f : m.iter_facets()) {
        if (hard_edges_attr[f.halfedge(0)] ||hard_edges_attr[f.halfedge(1)] ||hard_edges_attr[f.halfedge(2)]) {
            a_traite.push_back(f);
            dist2[f] = 0;
        }
    }

    int current_dist = 0;
    FacetAttribute coupe(m, 0);
    FacetAttribute ligne(m, -1);
    int espace = 25 ; //TODO essayer
    std::vector<int> maybe_seed;
    while (!a_traite.empty()) {
        a_traite2.clear();
        for (int fid: a_traite) {
            Surface::Facet f(m, fid);
            for (auto h: f.iter_halfedges()) {
                auto newf = h.opposite().facet();
                if (dist2[newf] != -1) {continue;}
                a_traite2.push_back(newf);
                dist2[newf] = current_dist;
                // if (current_dist % espace == 0 && current_dist != espace ) {coupe[newf] = 1; maybe_seed.push_back(newf); ligne[newf] = current_dist/espace;}
                if (current_dist == 0 ) {coupe[newf] = 1; maybe_seed.push_back(newf); ligne[newf] = current_dist/espace;}
            }
        }    
        a_traite = a_traite2;
        current_dist++;
    }


    // on stock les coins de chaques faces 
    
    //std::vector<std::vector<int>> coins(zone_num, {}); // stock les coins de chaque zone
    //std::vector<int> toutcoins;
    double tolerance = 0.3;
    PointAttribute coins_attr(m, 0);
    FacetAttribute coins_facet(m, 0);
    // FacetAttribute source(m, -1);
    
    //FacetAttribute dist(m, 1000);
    // int compte_source = 0;

    for (auto v: m.iter_vertices()) {
        //si il y a 3 facets ou c'est 2 à 2 de proche de 90 deg
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
                        // on regarde si c'est concave.
                        double angle1 = 0, angle2 = 0, angle3 = 0;
                        int zone1 = zone[f1], zone2 = zone[f2], zone3 = zone[f3];
                        vec3 p_prec;
                        for (auto h: v.iter_halfedges()) {
                            if (zone[h.facet()] == zone1) {
                                angle1 += 1; // vsy on suppose que si ya 4 triangle c"est concave moins robuste
                            }
                            if (zone[h.facet()] == zone2) {
                                angle2 += 1; 
                            }
                            if (zone[h.facet()] == zone3) {
                                angle3 += 1; 
                            }
                        }
                        if (angle1 >= 4) {coins_facet[f1] = true; coins_attr[v] = 1;}
                        if (angle2 >= 4) {coins_facet[f2] = true; coins_attr[v] = 1;}
                        if (angle3 >= 4) {coins_facet[f3] = true; coins_attr[v] = 1;}
                        

                    }
                }
            }
        }

    }

    FacetAttribute source2(m, -1);
    FacetAttribute dist3(m, 1000);
    FacetAttribute too_close(m, -1);
    int current_source = 0;
    a_traite.clear();
    a_traite2.clear();
    int dist_min = 30; //TODO essayer
    std::vector<int> a_voir;
    std::vector<int> a_voir2;
    FacetAttribute<bool> seed(m, false);
    


    for (int f: maybe_seed) {
        //on voudrait que ceux d'une ligne soit alignés avec ceux de la ligne d'après
        if (too_close[f] != ligne[f] && (ligne[f] == 0 || ligne[f] != -1)) {  //on veut éviter que ceux sur la même ligne
            too_close[f] = ligne[f];    
            source2[f] = current_source;
            dist3[f] = 0;
            current_source++;
            a_traite.push_back(f);
            seed[f] = true;
            a_voir = {f};
            for (int compte = 0; compte < dist_min; compte++) {
                for (auto fid: a_voir) { // f à distance i
                    Surface::Facet f1(m, fid);
                    for (auto h: f1.iter_halfedges()) {
                        if (hard_edges_attr[h]) {continue;}
                        auto f2 = h.opposite().facet();
                        if (too_close[f2] != ligne[f]) {
                            too_close[f2] = ligne[f],
                            a_voir2.push_back(f2);
                        }
                    }
                }
                a_voir = a_voir2;
                a_voir2.clear();
            }
        }
    }
    //charts gen
    do {
        for (auto fid: a_traite) { // f à distance i
            Surface::Facet f(m, fid);
            for (auto h: f.iter_halfedges()) {
                if (hard_edges_attr[h]) {break;}
                auto f2 = h.opposite().facet();
                if (dist3[f2] > dist3[f] + 1) {
                    dist3[f2] = dist3[f] + 1;
                    a_traite2.push_back(f2);
                    source2[f2] = source2[f];
                }
            }
        }
        a_traite = a_traite2;
        a_traite2.clear();
    } while (a_traite.size() != 0);



    
    write_by_extension(sortie, m, {{{"coins_attr", coins_attr.ptr}}, {{"coupe", coupe.ptr}, {"source", source2.ptr}, {"ligne", ligne.ptr}, {"seed", seed.ptr}, {"coins concave", coins_facet.ptr}/* {{"dist", dist.ptr},  */}, {{"hard_edge", hard_edges_attr.ptr}}});

    int result = system((getGraphitePath() + " " + sortie).c_str());
    // --- END ---

    return 0;
}
