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
    Polygons m;
    Triangles m0;
    read_by_extension(path + "joint_sub_sub.geogram", m0);
    // Loading catorus.geogram into m
    std::cout << path << std::endl;

    read_by_extension(path + "\\..\\Debug\\outpoly.geogram", m);

    //init quad
    Quads q;

    // Iterate over facets
    int compteur = 0;
    int compteur_face = 0;
    m.connect();
    for (auto f : m.iter_facets()) {
        //recuperer nb sommets
        auto f_geom = f.geom<Poly3>();
        int n = 0;
        for (auto _ : f.iter_halfedges()) {n++;}
        q.create_facets(n);
        q.points.create_points(2*n+1);
        q.points[compteur + 2*n] = f_geom.bary_verts(); //TODO faire mieux
        q.points[compteur] = f.vertex(0).pos();
        for (int i = 0; i < n; i++) {
            q.points[compteur + i ] = f.vertex(i).pos();
            q.points[compteur + n + i ] = (f.vertex(i).pos() + f.vertex((i+1) % n).pos())/2; //TODO faire mieux
        }
        for (int i = 0; i < n; i++) {
            q.vert(compteur_face, 0) = compteur + i ;
            q.vert(compteur_face, 1) = compteur + n + i;
            q.vert(compteur_face, 2) = compteur + 2 * n;
            q.vert(compteur_face, 3) = compteur + n  +  ((n + i - 1)% n);
            compteur_face++;
        }
        compteur += 2*n + 1;

    }

    std::vector<vec3> pt_list;
    std::vector<int> oldtonew;

    //build pt_list
    for (auto v : q.iter_vertices()) {
        pt_list.push_back(v.pos());
    } 

    colocate(pt_list, oldtonew, 1e-20);

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
    

    q.compact();
    q.disconnect();

    BVHTriangles bvh(m0);

    for (auto v : q.iter_vertices()) {
        q.points[v] = vec3(bvh.nearest_point(v.pos()));
    }
    
    write_by_extension("sphere_quad.geogram", q, {{}, {}, {} });
    int result = system((getGraphitePath() + " sphere_quad.geogram").c_str());
    // --- END ---

    return 0;
}