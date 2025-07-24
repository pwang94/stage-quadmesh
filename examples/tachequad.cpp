
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

    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Triangles m;
    static const std::string entree = "cube_tri.geogram";
    static const std::string sortie_chart = "sphere2.geogram";
    static const std::string sortie_tri = "a_transformer.geogram";
    read_by_extension(path + entree, m);
    m.connect();

    // les hard edges
    CornerAttribute<bool> hard_edges_attr(m);
    find_hard_edges(m, hard_edges_attr, 0.1);



    FacetAttribute<int> source(m, -1);
    FacetAttribute<int> dist(m, 1000);
    FacetAttribute<int> hard_face(m);
    
    for (auto f: m.iter_facets()) {
        for (int i= 0; i<3 ; i++) {
            if (hard_edges_attr[f.halfedge(i)]) {
                hard_face[f] = true;
            }
        }
    }

    std::vector<int> a_traite;
    std::vector<int> sommets;
    //generate seeds
    const int nb_source = 60;
    /*    int compteur = 0;
    
     for (auto f: m.iter_facets())
        {
            if (f % (m.nfacets()/nb_source) == 3) {
                source[f] = compteur;
                dist[f] = 0;
                compteur++;
                a_traite.push_back(f);
                //std::cout << source[f] << std::endl;
            }
            if (compteur >= nb_source) {break;}
        }
    */
    
    for (int i = 0; i<nb_source; i++) {
        int frand = rand() % m.nfacets();
        if (std::find(a_traite.begin(), a_traite.end(), frand) != a_traite.end() || !(hard_face[frand])){i--;}
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

    //generation de la nouvelle sphere

    Triangles m2;
    m2.points.create_points(nb_source);

    //generate vertices
    for (int i: sommets) {
        auto f = Triangles::Facet(m, i);
        auto tri_geo = f.geom<Triangle3>();
        auto b = tri_geo.bary_verts();
        m2.points[source[i]] = b;
    }
    
    int current_facet = 0;
    //generate facets

    std::vector<int> charts;
    std::set<int> charts_unordered;

    for (auto v: m.iter_vertices()) {
        //check si c'est une intersection de 3 charts
        charts.clear();
        charts_unordered.clear();
        for (auto h: v.iter_halfedges()) {
            charts_unordered.insert(source[h.facet()]);
        } //donne pas le bon ordre


        if (charts_unordered.size() >= 3) { 
            //on mets les sources dans le bon ordre
            Surface::Halfedge h0 = *v.iter_halfedges().begin() ;
            Surface::Halfedge h_curr = h0 ;

            int c0 = source[h_curr.facet()];
            int curr_chart = c0;
            charts.push_back(c0);
            Surface::Facet f = h_curr.facet();

            do {
                h_curr = h_curr.opposite();
                f = h_curr.facet();

                if (source[f] != curr_chart && source[f]!= c0) {
                    curr_chart = source[f];
                    charts.push_back(curr_chart);
                }
                //on cherche le 3eme sommet pour avoir la face d'après;
                int s;
                for (int i = 0; i < 3; i ++)
                    if (f.vertex(i) != v && f.vertex(i) != h_curr.from()) {
                        s = f.vertex(i);
                        //on cherche le halfedge v->s
                        for (auto h : v.iter_halfedges()) {
                            if (h.to() == s) {
                                h_curr = h;
                                 break;}
                        }
                        break;
                    }
                std::cout << charts.size() << std::endl;
            } while (charts.size() < charts_unordered.size());


            std::cout << charts.size() << std::endl;
            int base = charts[0];
            int suivant = charts[1];
            int compte = 2;
            while (compte < charts.size())
            { 
                m2.create_facets(1);
                m2.vert(current_facet, 0) = base;
                m2.vert(current_facet, 1) = suivant;
                suivant = charts[compte];
                m2.vert(current_facet, 2) = suivant;
                current_facet++;
                compte++;
            }
        }
    }

    

    write_by_extension(sortie_chart, m, {{}, {{"source", source.ptr}, {"dist", dist.ptr}, {"hard_face", hard_face.ptr}}, {{"hard_edge", hard_edges_attr.ptr}} });
    //int result = system((getGraphitePath() + " sphere2.geogram").c_str());
    write_by_extension(sortie_tri, m2, {{}, {}, {} });

    int result = system((getGraphitePath() + " " + sortie_tri).c_str());
    // --- END ---

    return 0;
}

