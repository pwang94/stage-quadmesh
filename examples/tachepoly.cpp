
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

    //TODO choisir le placement des sources





    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Triangles m;
    static const std::string entree = "joint_sub_sub.geogram";
    static const std::string sortie_chart = "charts.geogram";
    static const std::string sortie_tri = "charts_to_polygon.geogram";
    read_by_extension(path + entree, m);
    m.connect();

    // les hard edges
    CornerAttribute<bool> hard_edges_attr(m);
    find_hard_edges(m, hard_edges_attr, 0.6);

    FacetAttribute<int> source(m, -1);
    FacetAttribute<int> dist(m, 1000);
    FacetAttribute<int> hard_face(m);

    for (auto f : m.iter_facets()) {
        if (hard_edges_attr[f.halfedge(0)] || hard_edges_attr[f.halfedge(1)] || hard_edges_attr[f.halfedge(2)]) {
            hard_face[f] = true;
        }
    }

    std::vector<int> a_traite;
    std::vector<int> sommets;
    //generate seeds
    const int nb_source = 120;

    for (int i = 0; i< nb_source; i++) {
        int frand = rand() % m.nfacets();
        if (std::find(a_traite.begin(), a_traite.end(), frand) != a_traite.end() || hard_face[frand]){i--;}
        // on enleve les sources sur les bords
        // il faut mettre assez de sources dans chaque zone meme si les zones sont petites
        else {
        source[frand] = i;
        dist[frand] = 0;
        a_traite.push_back(frand);
        }
    }
    
    

    sommets = a_traite;
    std::vector<int> a_traite2;
    int compte_face = 0;


    
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
    } while (compte_face < m.nfacets() );



    write_by_extension(sortie_chart, m, {{}, {{"source", source.ptr}, {"dist", dist.ptr}}, {{"hard_edge", hard_edges_attr.ptr}} });


    //generation de la nouvelle sphere

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
    compte_face = 0;

    for (int f = 0; f < nb_source; f++) {
        //si nb paire de sommets on peut découper en quad

        int s0 = sommetsface[f][0];
        int s_curr = s0;
        int s_passe = s0;
        std::vector<int> sommetsordonne;
        sommetsordonne.clear();
        for (int j = 0; j < sommetsface[f].size(); j++) {
            //on passe au sommet suivant, il doit avoir 2 couleur en commun avec le cote précédent 
            for (int s : sommetsface[f]) {
                if (s != s_curr &&  s != s_passe && deuxcommun(charts[s], charts[s_curr])) {    
                    sommetsordonne.push_back(s_curr);
                    s_passe = s_curr;
                    s_curr = s;
                    break;
                }
            }
        }
        
        if (sommetsface[f].size() % 2 == 0) {
            m2.create_facets(1, 4);
            m2.vert(compte_face, 0) = sommetsordonne[0];
            m2.vert(compte_face, 1) = sommetsordonne[1];
            for (int i = 1; i < sommetsface[f].size() / 2 - 1; i++) {
                m2.vert(compte_face, 2) = sommetsordonne[1+i];
                m2.vert(compte_face, 3) = sommetsordonne[sommetsordonne.size() - i];
                compte_face++;
                m2.create_facets(1, 4);
                m2.vert(compte_face, 0) = sommetsordonne[sommetsordonne.size() - i];
                m2.vert(compte_face, 1) = sommetsordonne[1+i];
            }
                m2.vert(compte_face, 2) = sommetsordonne[sommetsface[f].size() / 2];
                std::cout << sommetsordonne.size() << " " << sommetsface[f].size() << " " << sommetsordonne.size() - (sommetsface[f].size() / 2 - 1) << std::endl;
                m2.vert(compte_face, 3) = sommetsordonne[sommetsordonne.size() - (sommetsface[f].size() / 2 - 1)];
            compte_face++;
        }
        else {
            m2.create_facets(1, sommetsface[f].size());
            for (int i= 0; i < sommetsface[f].size(); i++) {
                m2.vert(compte_face, i) = sommetsordonne[i];
            }
            compte_face++;
        }
    }
    


    //projeter

    /* BVHTriangles bvh(m); */

    /*  for (auto v : m2.iter_vertices()) {
        m2.points[v] = vec3(bvh.nearest_point(v.pos()));
    } */


    write_by_extension(sortie_chart, m, {{}, {{"source", source.ptr}, {"dist", dist.ptr}}, {{"hard_edge", hard_edges_attr.ptr}} });
    write_by_extension(sortie_tri, m2, {{}, {}, {} });

    int result = system((getGraphitePath() + " " + sortie_tri).c_str());
    // --- END ---

    return 0;
}
