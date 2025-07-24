
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

void print(auto s) {
    std::cout<< s << std::endl;
}

int main(int argc, char** argv) {
    // --- LOAD ---

    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Triangles m;
    static const std::string entree = "joint_sub_sub.geogram";
    static const std::string sortie = "sortie.geogram";
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

    PointAttribute bord(m, false);
    //compute distance from bord
    for (auto h: m.iter_halfedges()) {
        if (hard_edges_attr[h]) {bord[h.from()] = true;}
    }

    for (auto f : m.iter_facets()) {
        if (hard_edges_attr[f.halfedge(0)] ||hard_edges_attr[f.halfedge(1)] ||hard_edges_attr[f.halfedge(2)]) {
            a_traite.push_back(f);
            dist2[f] = 0;
        }
    }
    int nb_sourcebord = a_traite.size();

    // on sample le bord d'une face

    FacetAttribute interdit(m, false);
    int dist_min = 13;
    std::vector<int> a_voir;
    std::vector<int> a_voir2;
    std::vector<std::vector<int>> sample;
    FacetAttribute source(m, -1);
    FacetAttribute dist(m, 1000);
    int current_source = 0;
    std::vector<int> face_from_source;

    for (int fid: a_traite) {
        //on propage de distance_min interdit
        if (interdit[fid]) {continue;}
        a_traite2.push_back(fid);
        source[fid] = current_source;
        face_from_source.push_back(fid);
        current_source++;
        dist[fid] = 0;
        a_voir = {fid};
        for (int compte = 0; compte < dist_min; compte++) {
            for (auto fid: a_voir) { // f à distance i
                Surface::Facet f(m, fid);
                for (auto h: f.iter_halfedges()) {
                    if (hard_edges_attr[h]) {continue;}
                    auto f2 = h.opposite().facet();
                    if (interdit[f2] == false) {
                        interdit[f2] = true;
                        a_voir2.push_back(f2);
                    }
                }
            }
            a_voir = a_voir2;
            a_voir2.clear();
        }
    }

    a_traite = a_traite2;
    a_traite2.clear();
    //on fait grossir les source au bord 
    do {
        for (auto fid: a_traite) { // f à distance i
            Surface::Facet f(m, fid);
            for (auto h: f.iter_halfedges()) {
                if (hard_edges_attr[h]) {continue;}
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
        std::reverse(a_traite.begin(), a_traite.end());
    } while (a_traite.size() != 0);

    FacetAttribute<bool> axe_median(m, false);
    std::vector<int> maybe_seed(nb_sourcebord, false);

    //on note toutes les paires de sources voisines
    std::vector<std::vector<int>> voisins(nb_sourcebord, std::vector<int>(nb_sourcebord, false));
    for (auto v: m.iter_vertices()) {
        if (bord[v]) {
            for (auto h: v.iter_halfedges()) {
                    voisins[source[h.facet()]][source[h.opposite().facet()]] = true;
                    voisins[source[h.opposite().facet()]][source[h.facet()]] = true;
            }
        }
    }
    

    // axe médian
    for (auto h: m.iter_halfedges()) {
        int s1 = source[h.facet()];
        int s2 = source[h.opposite().facet()];
        if (s1 != s2) {
            Surface::Facet f1(m, face_from_source[s1]); 
            Surface::Facet f2(m, face_from_source[s2]); 
            //on veut verifier que les sources sont pas voisines
            if (voisins[s1][s2] == false) {

                axe_median[h.facet()]= true;
                maybe_seed.push_back(h.facet());
            }
        }
    }

    //seed gen
    FacetAttribute source2(m, -1);
    FacetAttribute dist3(m, 1000);
    FacetAttribute too_close(m, false);
    current_source = 0;
    a_traite.clear();
    a_traite2.clear();
    dist_min = 5;

    for (int f: maybe_seed) {
        if (too_close[f] == false) {  
            too_close[f] = true;    
            source2[f] = current_source;
            dist3[f] = 0;
            current_source++;
            a_traite.push_back(f);
            a_voir = {f};
            for (int compte = 0; compte < dist_min; compte++) {
                for (auto fid: a_voir) { // f à distance i
                    Surface::Facet f1(m, fid);
                    for (auto h: f1.iter_halfedges()) {
                        if (hard_edges_attr[h]) {continue;}
                        auto f2 = h.opposite().facet();
                        if (too_close[f2] == false) {
                            too_close[f2] = true;
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



    write_by_extension(sortie, m, {{}, {{"dist", dist3.ptr},  {"source", source2.ptr},{"too_close", too_close.ptr} ,{"source0", source.ptr},{"axe_median", axe_median.ptr}}, {{"hard_edge", hard_edges_attr.ptr}}});

    int result = system((getGraphitePath() + " " + sortie).c_str());
    // --- END ---
    

    return 0;
}
