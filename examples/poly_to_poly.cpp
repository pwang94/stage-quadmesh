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
    // on veut plus de quad !
    FacetAttribute<bool> singu("singu", attrs, m);
    PointAttribute<bool> hard_point(m, false);
    CornerAttribute<bool> hard_edges_attr(m, false);
    find_hard_edges(m, hard_edges_attr, hard_point, 0.6);
    // on peut supprimer des petites arêtes
    CornerAttribute<int> to_kill(m, 0);
    PointAttribute<int> to_kill_ext(m, 0);
    const double dist_max = 0.12;
    std::vector<int> small_he;
    CornerAttribute<int> bizarre(m, 0);
    double marge = 0.01;
    PointAttribute<std::vector<vec3>> proche(m); 
    PointAttribute<int> proche2(m); 
    FacetAttribute<int> nbcote(m);

    for (auto f: m.iter_facets()) {
        nbcote[f] = f.size();
        for (int i= 0; i < f.size(); i ++) {
            vec3 v1 = f.vertex(i).pos();
            vec3 v2 = f.vertex((i + 1) % f.size()).pos();
            vec3 v3 = f.vertex((i + 2) % f.size()).pos();
            Segment3 s({v1, v3});
            vec3 vapprox = s.nearest_point(v2);
            if ((v2 - vapprox).norm()/((v1 - v3).norm()) < marge) {
                proche[f.vertex((i + 1) % nbcote[f])] = {vapprox};
                proche2[f.vertex((i + 1) % nbcote[f])] = 1;
                nbcote[f]--;
            }
        }
    }

    for (auto h : m.iter_halfedges()) {
        if (h.opposite() == - 1 || !h.active()) {bizarre[h] = 1; std::cout << "oups" <<std::endl;}
        else {
            int f1 = h.facet();
            int f2 = h.opposite().facet();
            int n1 = nbcote[f1]; 
            int n2 = nbcote[f2];
            int s1 = h.from(); 
            int s2 = h.to();

            if (n1 > 4 && n2 > 4 && 
                // TODO: mieux définir le nb de sommet si y'a 3 pts presque alignés on veut pas compter celui du milieux
                // !hard_edges_attr[h] &&
                // !hard_point[h.from()]&&
                // !hard_point[h.to()] &&
                h.active() && 
                !(hard_point[s1] && !hard_point[s2]) && //????? ok
                !(hard_point[s2] && !hard_point[s1]) &&
                (h.from().pos() - h.to().pos()).norm() < dist_max) {
                to_kill[h] = 1;
                to_kill[h.opposite()] = 1;
                to_kill_ext[h.from()] = 1;
                to_kill_ext[h.to()] = 1;
                small_he.push_back(h);
            }
        }
    }

    std:: cout << "passage" << std::endl;
    write_by_extension("polyatt.geogram", m, {{{"to_killpt", to_kill_ext.ptr},  {"hard_point", hard_point.ptr}}, {}, {{"to_kill", to_kill.ptr}, {"bizarre", bizarre.ptr}} });
    //TODO on transforme du plus petit écart au plus grand si c'est encore possible

    std::sort(small_he.begin(), small_he.end(), 
        [&m](auto h1, auto h2) {
            Surface::Halfedge he1(m, h1);
            Surface::Halfedge he2(m, h2);
            return (he1.from().pos() - he1.to().pos()).norm() < (he2.from().pos() - he2.to().pos()).norm();
        } );

    for (int hid : small_he) {
        Surface::Halfedge h(m, hid);
        if (!h.active()) continue;
        auto s1 = h.from();
        auto s2 = h.to();
        if (!h.active() || !s1.halfedge().active() || !s2.halfedge().active()) {continue;}
        //on check que toutes les arêtes sont pas inactives
        auto new_point_pos = (s1.pos() + s2.pos())/2;
        if (hard_point[s1] && !hard_point[s2]) {new_point_pos = s1.pos();}
        if (hard_point[s2] && !hard_point[s1]) {new_point_pos = s2.pos();}
        
        m.points.create_points(1);
        int new_point = m.points.size()-1;
        m.points[new_point] = new_point_pos;
        //pour les 2 faces qui ont les 2 pts on remplace
        Surface::Facet f1 = h.facet();
        Surface::Facet f2 = h.opposite().facet();
        m.conn->active[f1] = false;
        std::vector<int> new_face;
        for (int n = 0; n < m.facet_size(f1); n++) {
            int v = f1.vertex(n);
            if (v != s1 && v != s2) {new_face.push_back(v);}
            else if (v == s1) {new_face.push_back(new_point);}
        }
        m.conn->create_facet(new_face.data(), m.facet_size(f1) - 1);
        //face 2
    std:: cout << "passage" << std::endl;
        
        m.conn->active[f2] = false;
        new_face.clear();
        for (int n = 0; n < m.facet_size(f2); n++ ) {
            int v = f2.vertex(n);
            if (v != s1 && v != s2) {new_face.push_back(v);}
            else if (v == s1) {new_face.push_back(new_point);}
        }
        m.conn->create_facet(new_face.data(), m.facet_size(f2) - 1);
        //s1
        for (auto he: s1.iter_halfedges()) { 
            new_face.clear();
            Surface::Facet f = he.facet();
            if (int(he) != int(h)) {
                for (int n = 0; n < m.facet_size(f); n++ ) {
                    int v = f.vertex(n);
                    if (v != s1) {new_face.push_back(v);}
                    else {new_face.push_back(new_point);}
                }
            }
            m.conn->active[f] = false;
            m.conn->create_facet(new_face.data(), m.facet_size(f));
        }
        //s2
        for (auto he: s2.iter_halfedges()) { 
            new_face.clear();
            Surface::Facet f = he.facet();
            if (int(he) != int(h)) {
                for (int n = 0; n < m.facet_size(f); n++ ) {
                    int v = f.vertex(n);
                    if (v != s2) {new_face.push_back(v);}
                    else {new_face.push_back(new_point);}
                }
            }
            m.conn->active[f] = false;
            m.conn->create_facet(new_face.data(), m.facet_size(f));
    std:: cout << "passage" << "2" << std::endl;
            
        }
    }
    m.compact();

    // on récupère toutes les faces qui touchent un des 2 point et on remplace par newpoint

    //TODO si une tache touche le bord en 1 point on agrandit

    write_by_extension("newpoly.geogram", m, {{{"to_killpt", to_kill_ext.ptr}, {"hard_point", hard_point.ptr}}, {}, {{"to_kill", to_kill.ptr}} });
    int result = system((getGraphitePath() + " newpoly.geogram").c_str()); 
    // --- END ---
    
    return 0;
}