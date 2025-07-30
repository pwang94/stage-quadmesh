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

    std::string path = getAssetPath();

    Triangles m; // charts
    Triangles mff; //framefield
    static const std::string entree = "sortie_recentre.geogram";
    static const std::string entree_poly = "outpoly.geogram";
    static const std::string sortie = "sortie_disque.geogram";
    static const std::string ff_path = "framefield.geogram";
    SurfaceAttributes attrs = read_by_extension( entree, m);
    SurfaceAttributes attrs_framefield = read_by_extension( ff_path, mff);
    m.connect();


    PointAttribute<int> coins_attr("coins_attr", attrs, m);
    FacetAttribute<int> coins_concave("coins concave", attrs, m);

    PointAttribute<int> singularite("singularite", attrs, m);
    FacetAttribute<int> source("source", attrs, m);
    FacetAttribute<int> racine("racine", attrs, m);
    FacetAttribute<double> dist("dist", attrs, m);

    CornerAttribute<bool> hard_edges_attr(m, false);
    find_hard_edges(m, hard_edges_attr, 0.6);
    FacetAttribute<int> new_racine(m, -1);

    std::map<int, int> singu_of_source; //source -> singu
    std::map<int, int> racine_of_source; //source -> racine

    for (auto v: m.iter_vertices()) {
        if (singularite[v] != 0) {
            singu_of_source[source[v.halfedge().facet()]] = v.halfedge().facet();
        }
    }

    for (auto f: m.iter_facets()) {
        if (racine[f] != -1) {
            racine_of_source[source[f]] = f;
        }
    }

    int nb_source = racine_of_source.size();

    //on check pour chaque charts si c'est un non disque topologique avecc la formule d'euler
    // on calcule sommets arêtes faces

    
    std::vector<int> sommets(nb_source);
    std::vector<double> aretes(nb_source);
    std::vector<int> faces(nb_source, 1);

    //faces
    for (auto f: m.iter_facets()) {
        faces[source[f]]++;
    }

    //sommets
    std::vector<int> deja_vu;
    for (auto v: m.iter_vertices()) {
        deja_vu.clear();
        for (auto h: v.halfedge().iter_sector_halfedges()) {
            if (std::find(deja_vu.begin(), deja_vu.end(), source[h.facet()]) == deja_vu.end()) {
                deja_vu.push_back(source[h.facet()]);
                sommets[source[h.facet()]]++;
            }
        }
    }

    //aretes
    for (auto h: m.iter_halfedges()) {
        aretes[source[h.facet()]]+= 1/2.;
        if (source[h.facet()] != source[h.opposite().facet()]) {
            aretes[source[h.facet()]]+= 1/2.;
        }
    }

    int nb_source_init = nb_source;

    for (int i = 0; i < nb_source_init; i++) {
        if (sommets[i] - aretes[i] + faces[i] != 2) { // par la formule d'euler c'est un disque que dans ce cas
            //on ajoute une source dans ce charts
            for (auto f: m.iter_facets()) {
                if (source[f] == i) {
                    racine[f] = nb_source;
                    source[f] = nb_source;
                    nb_source++;
                    break;
                }
            }
        }
    }




    write_by_extension(sortie, m, {{{"singularite", singularite.ptr}, {"coins_attr", coins_attr.ptr}}, {{"source", source.ptr}, {"racine", racine.ptr}, {"dist", dist.ptr}, {"coins concave",coins_concave.ptr}}, {}});

    int result = system((getGraphitePath() + " " + sortie).c_str());
    // --- END ---

    return 0;
}
