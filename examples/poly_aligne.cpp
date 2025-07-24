#include "helpers.h"
#include <ultimaille/all.h>
#include <map>

using namespace UM;

void find_hard_edges(Polygons& mesh, CornerAttribute<bool>& hard_edges_attr, PointAttribute<bool>& hard_point ,double threshold) {

    // Iter on all mesh halfedges
    for (auto h : mesh.iter_halfedges()) {
        
        // Get opposite halfedge
        auto opposite = h.opposite();

        if (!opposite.active()) {std::cout << (int)h << std::endl; continue;}
        
        // Compute normals of a face and its opposite face
        vec3 normalA = h.facet().geom<Poly3>().normal();
        vec3 normalB = opposite.facet().geom<Poly3>().normal();
        // Compute the dot product of normals
        double d = normalA * normalB;

        // The closer the scalar product of its normals is to zero, 
        // the closer the angle between two faces is to 90 degrees
        if (std::abs(d) <= threshold) {
            hard_edges_attr[h] = true;
            hard_point[h.from()] = true;
            hard_point[h.to()] = true;
        }
    }
}


int main(int argc, char** argv) {

    // --- LOAD ---

    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Polygons m;
    Triangles m0;
    read_by_extension(path + "joint_sub_sub.geogram", m0);
    
    read_by_extension(path + "\\..\\Debug\\outpoly.geogram", m);
    m.connect();
    // on veut aligner les pts presque aligné
    //peut être les suppr ?
    FacetAttribute<int> nbcote(m); //on ne compte pas les côtés ou les sommets sont presque aligné 
    PointAttribute<std::vector<vec3>> proche(m); 
    PointAttribute<int> proche2(m); 
    double marge = 0.06;
    for (auto f: m.iter_facets()) {
        nbcote[f] = m.facet_size(f);
        if (nbcote[f] <= 4) {continue;}
        for (int i = 0; i < nbcote[f] && nbcote[f] > 4; i++) {
            vec3 v1 = f.vertex(i).pos();
            vec3 v2 = f.vertex((i + 1) % nbcote[f]).pos();
            vec3 v3 = f.vertex((i + 2) %nbcote[f]).pos();
            Segment3 s({v1, v3});
            vec3 vapprox = s.nearest_point(v2);

            //ou

            if ((v2 - vapprox).norm()/((v1 - v3).norm()) < marge) {
                proche[f.vertex((i + 1) % nbcote[f])] = {vapprox};
                proche2[f.vertex((i + 1) % nbcote[f])] = 1;
                nbcote[f]--;
            }
        }
    }

    write_by_extension("newpolya.geogram", m, {{{"proche", proche2.ptr}}, {}, {}});
    int result = system((getGraphitePath() + " newpolya.geogram").c_str()); 
    // --- END ---
    
    return 0;
}