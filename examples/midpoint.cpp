/**
 * This example shows how to iterate primitives over a mesh without connectivity
*/
#include "helpers.h"
#include <ultimaille/all.h>
#include <map>

using namespace UM;

int main(int argc, char** argv) {

    // --- LOAD ---

    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Triangles m;
    // Loading catorus.geogram into m
    read_by_extension(path + "a_transformer.geogram", m);

    //init quad
    Quads q;
    q.points.create_points(m.nfacets()*7);
    q.create_facets(m.nfacets()*3);

    /**/
    // Iterate over facets
    int compteur = 0;
    int compteur2 = 0;
    for (auto f : m.iter_facets()) {
        q.points[7*compteur] = (f.vertex(0).pos() + f.vertex(1).pos() + f.vertex(2).pos()) / 3;
        q.points[7*compteur+1] = f.vertex(0).pos();
        q.points[7*compteur+2] = f.vertex(1).pos();
        q.points[7*compteur+3] = f.vertex(2).pos();
        q.points[7*compteur+4] = (f.vertex(0).pos() + f.vertex(1).pos())/2;
        q.points[7*compteur+5] = (f.vertex(1).pos() + f.vertex(2).pos())/2;
        q.points[7*compteur+6] = (f.vertex(0).pos() + f.vertex(2).pos())/2;
        q.vert(compteur2*3, 0) = 7*compteur+1;
        q.vert(compteur2*3, 1) = 7*compteur+4;
        q.vert(compteur2*3, 2) = 7*compteur;
        q.vert(compteur2*3, 3) = 7*compteur+6;
        q.vert(compteur2*3+1, 0) = 7*compteur+2;
        q.vert(compteur2*3+1, 1) = 7*compteur+5;
        q.vert(compteur2*3+1, 2) = 7*compteur;
        q.vert(compteur2*3+1, 3) = 7*compteur+4;
        q.vert(compteur2*3+2, 0) = 7*compteur+3;
        q.vert(compteur2*3+2, 1) = 7*compteur+6;
        q.vert(compteur2*3+2, 2) = 7*compteur;
        q.vert(compteur2*3+2, 3) = 7*compteur+5;
        compteur++;
        compteur2++;
    }

    std::vector<vec3> pt_list;
    std::vector<int> oldtonew;
    //std::vector<bool> to_kill(taille, false);

    //build pt_list
    for (auto v : q.iter_vertices()) {
        pt_list.push_back(v.pos());
    } 

    colocate(pt_list, oldtonew, 1e-20);

    
    //BVHTriangles bvh(m);

    std::map<int, int> taken;

    q.connect();

    for (auto f : q.iter_facets()) {
        for (int i = 0; i <= 3;i++) {
            if (taken.find(oldtonew[(int)f.vertex(i)]) == taken.end()) {
                taken[oldtonew[(int)f.vertex(i)]] = (int)f.vertex(i);
                //to_kill[(int)f.vertex(i)] = false; //
            }
            else {
                if (taken.find(oldtonew[(int)f.vertex(i)])-> second != (int)f.vertex(i)) {
                    //to_kill[(int)f.vertex(i)] = true;
                    q.conn->active[f] = false;
                }
            }
        }
        if (q.conn->active[f] == false) {
            q.conn->create_facet({taken.find(oldtonew[f.vertex(0)])->second, taken.find(oldtonew[f.vertex(1)])-> second, taken.find(oldtonew[f.vertex(2)])-> second, taken.find(oldtonew[f.vertex(3)])-> second});
            //nouvelle face avec que des sommets gardés
        }
    }
    
    //q.delete_vertices(to_kill);

    q.compact();
    q.disconnect();
    /* FacetAttribute<int> chart(q, -1);
    chart[3] = 12; */
    
    write_by_extension("sphere_quad.geogram", q, {{}, {}, {} });
    int result = system((getGraphitePath() + " sphere_quad.geogram").c_str());
    //save
    // --- END ---

    return 0;
}