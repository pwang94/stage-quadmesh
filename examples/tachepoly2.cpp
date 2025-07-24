
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
    static const std::string entree = "B72sub.geogram";
    static const std::string sortie_chart = "charts.geogram";
    static const std::string sortie_tri = "charts_to_polygon.geogram";
    read_by_extension(path + entree, m);
    m.connect();

    // les hard edges
    CornerAttribute<bool> hard_edges_attr(m, false);
    find_hard_edges(m, hard_edges_attr, 0.6);

    FacetAttribute<int> source(m, -1);
    FacetAttribute<int> dist(m, 1000);
    FacetAttribute<int> hard_face(m);

    for (auto f : m.iter_facets()) {
        if (hard_edges_attr[f.halfedge(0)] ||hard_edges_attr[f.halfedge(1)] ||hard_edges_attr[f.halfedge(2)]) {
            hard_face[f] = true;
        }
    }

    std::vector<int> a_traite;
    std::vector<int> a_traite2;
    int compte_face = 0;
    FacetAttribute<int> zone(m, -1);
    int zone_num = 0;
    std::vector<std::vector<int>> de_zone; // a chaque zone donne les faces dedans

    // obtenir les grosses zone engendrées par les hards edges

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

    //generate seeds
    const int source_add = 75;
    const int source_mini = 5;
    //on met source_mini source dans chaque zone
    int current_source = 0;
    for (int c = 0; c < zone_num; c++) {
        for (int i = 0; i< source_mini; i++) {
            int frand = de_zone[c][rand() % de_zone[c].size()];
            if (std::find(a_traite.begin(), a_traite.end(), frand) != a_traite.end() || hard_face[frand]){i--;}
            else {        
                source[frand] = current_source;
                dist[frand] = 0;
                current_source++;
                a_traite.push_back(frand);
            }
        }
    }

    int nb_source = source_mini * zone_num + source_add;
    FacetAttribute<int> interdit(m, false);
    int dist_mini = 30;

    for (int i = source_mini* zone_num; i < nb_source; i++) {
        int frand = rand() % m.nfacets();
        if (std::find(a_traite.begin(), a_traite.end(), frand) != a_traite.end() || hard_face[frand] || interdit[frand]){i--;}
        // on enleve les sources sur les bords
        // on veux qu'ils soit loins des bords? et loins les uns des autres
        else {
            source[frand] = current_source;
            current_source++;
            dist[frand] = 0;
            a_traite.push_back(frand);
            //on regarde la distance aux autres: on peut propager de 20 cases chaque fois qu'on rencontre une nouvelle source
            std::vector<int> a_voir = {frand};
            std::vector<int> a_voir2;
            for (int compte = 0; compte < dist_mini; compte++) {
                for (auto fid: a_voir) { // f à distance i
                    Surface::Facet f(m, fid);
                    for (auto h: f.iter_halfedges()) {
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
    }

    compte_face = 0;

    //generate charts
    do {
        for (auto fid: a_traite) { // f à distance i
            compte_face++;
            Surface::Facet f(m, fid);
            for (auto h: f.iter_halfedges()) {
                if (hard_edges_attr[h]) {break;}
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
    } while (compte_face < m.nfacets());

    //generation de la nouvelle surface

    Polygons m2;
    std::vector<std::vector<int>> sommetsface(nb_source);
    PointAttribute<std::vector<int>> charts(m2.points);
    std::set<int> charts_unordered;
    int num_vertice = 0;

    //on fait un tableau qui à chaque futur sommet donne toute ses sources voisines

    for (auto v: m.iter_vertices()) {
        //check si c'est une intersection de 3 charts
        charts_unordered.clear();
        for (auto h: v.iter_halfedges()) {
            charts_unordered.insert(source[h.facet()]);
        } //donne pas le bon ordre

        if (charts_unordered.size() >= 3) { 
            //on rajoute le point
            m2.points.create_points(1);
            m2.points[num_vertice] = v.pos();
            for (auto f: charts_unordered) {
                sommetsface[f].push_back(num_vertice);
                charts[num_vertice].push_back(f);
            }
        
            num_vertice++;
            //pour chacun des sommets on les lie si ils ont 2 sources en commun
        }
    }

    // creer facet

    for (int f = 0; f< nb_source; f++) {
        int s0 = sommetsface[f][0];
        int s_curr = s0;
        int s_passe = s0;
        m2.create_facets(1, sommetsface[f].size());
        for (int j = 0; j < sommetsface[f].size(); j++) {
            //on passe au sommet suivant, il doit avoir 2 couleur en commun avec le cote précédent 
            for (int s : sommetsface[f]) {
                if (s != s_curr &&  s != s_passe && deuxcommun(charts[s], charts[s_curr])) {    
                    m2.vert(f, j) = s_curr;
                    s_passe = s_curr;
                    s_curr = s;
                    break;
                }
            }
        }
    }
    
    //projeter

    /*BVHTriangles bvh(m);

    for (auto v : m2.iter_vertices()) {
        m2.points[v] = vec3(bvh.nearest_point(v.pos()));
    } */


    write_by_extension(sortie_chart, m, {{}, {{"source", source.ptr}, {"dist", dist.ptr}, {"hard", zone.ptr}, {"interdit", interdit.ptr}}, {{"hard_edge", hard_edges_attr.ptr}} });
    write_by_extension(sortie_tri, m2, {{}, {}, {} });

    int result = system((getGraphitePath() + " " + sortie_tri).c_str());
    // --- END ---

    return 0;
}
