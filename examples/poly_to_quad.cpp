
#include "helpers.h"
#include <ultimaille/all.h>
#include <map>
#include <set>

using namespace UM;

int main(int argc, char** argv) {
    // --- LOAD ---

    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Polygons m;
    Polygons m2;
    Triangles m_r;
    static const std::string entree = "\\..\\Debug\\newpoly.geogram";
    static const std::string sortie = "sortie_quad.geogram";
    read_by_extension(path + entree, m);

    m.connect();
    m2.connect();
    
    //TODO on regarde si la face est issue d'une singu: si le projeté du bary dans le modele tri à une bonne couleur de tache
    SurfaceAttributes attrs = read_by_extension(path + "\\..\\Debug\\sortie_recentre.geogram", m_r);
    PointAttribute<int> singularite("singularite", attrs, m_r);
    FacetAttribute<int> source0("source", attrs, m_r);
    BVHTriangles bvh(m_r);
    m_r.connect();
    std::map<int, int> singu_of_source; //source ->singu

    for (auto v: m_r.iter_vertices()) {
        if (singularite[v] != 0) {
            std::cout << singularite[v] <<"singu?" <<std::endl;
            singu_of_source[source0[v.halfedge().facet()]] = v.halfedge().facet();
        }
    }

    for (auto v: m.iter_vertices()) {
        m2.points.create_points(1);
        m2.points[v] = v.pos();
    }


    for (auto f: m.iter_facets()) {

        std::vector<int> new_face;
        if (f.size() % 2 == 0 /* et pas issue d'une singu */ &&
            singu_of_source.find(source0[bvh.nearest_point (f.geom<Poly3>().bary_verts()).f]) == singu_of_source.end()
            )
             {
            std::cout << f.size() <<std::endl;

            new_face.push_back(f.vertex(0));
            new_face.push_back(f.vertex(1));

            //TODO bien choisir qui mettre avec qui, que les valences soit bien

            for (int i = 1; i < f.size() / 2 - 1; i++) {
                new_face.push_back(f.vertex(i+1));
                new_face.push_back(f.vertex(-i+f.size()));
                
                m2.conn->create_facet(new_face.data(), 4); 

                new_face.clear();

                new_face.push_back(f.vertex(-i+f.size()));

                new_face.push_back(f.vertex(i+1)); 
            }
            new_face.push_back(f.vertex(f.size() / 2));
            new_face.push_back(f.vertex(f.size() - (f.size()/2 - 1)));
            m2.conn->create_facet(new_face.data(), 4); 
        }
        else {
            //on garde la face >_<
            for (int i = 0; i < f.size(); i++) {new_face.push_back(f.vertex(i));}
            m2.conn->create_facet(new_face.data(), new_face.size()); 
            new_face.clear();
        }
    }

    write_by_extension(sortie, m2, { {}, {}, {} });
    int result = system((getGraphitePath() + " " + sortie).c_str()); 

    return 0;
    
}
