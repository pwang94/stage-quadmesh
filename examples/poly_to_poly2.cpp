#include "helpers.h"
#include <ultimaille/all.h>
#include <map>

using namespace UM;

void find_hard_edges(Polygons& mesh, CornerAttribute<bool>& hard_edges_attr, PointAttribute<bool>& hard_point ,double threshold) {

    // Iter on all mesh halfedges
    for (auto h : mesh.iter_halfedges()) {
        
        // Get opposite halfedge
        auto opposite = h.opposite();

        if (!opposite.active()) {std::cout << "inactive" << (int)h << std::endl; continue;}
        
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
    
    SurfaceAttributes attrs = read_by_extension(path + "\\..\\Debug\\outpoly.geogram", m);
    m.connect();
    
    FacetAttribute<bool> singu("singu", attrs, m);
    PointAttribute<bool> hard_point(m, false);
    CornerAttribute<bool> hard_edges_attr(m, false);
    find_hard_edges(m, hard_edges_attr, hard_point, 0.6);
    // on peut supprimer des petites arêtes
    PointAttribute<bool> to_change(m, 0);

    std::vector<int> a_change;

    for (auto v: m.iter_vertices()) {
        // on check si c'est valence 4 et tous ces voisins sont valences 3
        int valence_v = 0;
        bool rajoute = true;
        if (!v.halfedge().active()) continue;
        for (auto h: v.halfedge().iter_sector_halfedges()) {
            if (h.active()) continue;
            valence_v++;
            int valence_voisins = 0;
            for (auto h2: h.to().halfedge().iter_sector_halfedges()) {
                valence_voisins++;
            }
            if (valence_voisins != 3) rajoute = false;
        } 
        if (valence_v != 4) rajoute = false;
        if (rajoute) {
            to_change[v] = true;
            a_change.push_back(v);
            std::cout << "test" <<std::endl;}
    }

    write_by_extension("polyatt2.geogram", m, {{{"to_change", to_change.ptr}}, {}, {}});

    // on ajoute une face et on modifie les 4autres
    //m.conn->create_facet(new_face.data(), m.facet_size(f1) - 1);
    std::vector<int> new_face;
    std::vector<std::pair<std::vector<int>, int>> new_faces;

    for (int v: a_change) {
        Surface::Vertex sommet(m, v);
        for (auto h: sommet.halfedge().iter_sector_halfedges()) {
            new_face.clear();
            for (auto h: h.facet().iter_halfedges()) {
                //on enlève le sommet v à la face
                if (h.to() != v) new_face.push_back(h.to());
            }
            //TODO faire tout a la fin
            new_faces.push_back({new_face, h.facet()});
            m.conn->create_facet(new_face.data(), new_face.size());
            m.conn->active[h.facet()] = false; 
        }
        // on ajoute la dernière face
        new_face.clear();
        for (auto h: sommet.halfedge().iter_sector_halfedges()) {
            new_face.push_back(h.to());
        }
        m.conn->create_facet(new_face.data(), new_face.size());
    }

    int result = system((getGraphitePath() + " polyatt2.geogram").c_str()); 
    // --- END ---
    
    return 0;
}