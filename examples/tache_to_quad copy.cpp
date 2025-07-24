

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
    return (compte >= 2);
}


int encommun(std::vector<int> v1, std::vector<int> v2, int f) {
    for (int c : v2) {
        if (std::find(v1.begin(), v1.end(), c) != v1.end() && c != f) {return c;};
    }
    assert(false);
    return 0;
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

int mid(std::map<std::pair<int, int>, int> &milieu, std::map<std::pair<int, int>, vec3> &milieu1, int c1, int c2, int &compte, Polygons& m
     /*, PointAttribute<std::pair<int,int>> &milieu2 */) {
    if (milieu.find({c1, c2}) == milieu.end()) {
        m.points.create_points(1);
        m.points[compte] = milieu1[{c1, c2}];
        std::cout << "nouveau pt "<< compte << std::endl;
        std::cout << "c1, c2 " << c1 << " " << c2 <<std::endl;

        std::cout << "pos " << milieu1[{c1, c2}]<<std::endl;
        milieu[{c1, c2}] = compte;
        //milieu2[compte] = {s1, s2};
        compte++;
    }
    return milieu[{c1, c2}];
}


int main(int argc, char** argv) {
    // --- LOAD ---

    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Triangles m;
    static const std::string entree = "\\..\\Debug\\sortie.geogram";
    static const std::string sortie = "outpoly.geogram";
    read_by_extension(path + entree, m);
    m.connect();

    //il faut avoir source, nb_source
    SurfaceAttributes attrs = read_by_extension(path + entree, m);


    FacetAttribute<int> source("source", attrs, m);
    FacetAttribute<int> racine("racine", attrs, m); 

    int nb_source = 0;

    for (auto f:m.iter_facets()) {
        if (source[f] + 1>nb_source) {nb_source = source[f] + 1;}
    }
    

    //generation de la nouvelle surface
    Polygons m2;
    std::vector<std::vector<int>> sommetsface(nb_source);
    std::map<int, int> oldtonew;
    PointAttribute<std::vector<int>> charts(m2.points);
    std::set<int> charts_unordered; // ancien -> nouveau
    int num_vertice = 0;

    //on fait un tableau qui à chaque futur sommet donne toute ses sources voisines
    std::cout << "test0" << std::endl;

    for (auto v: m.iter_vertices()) {
        //check si c'est une intersection de 3 charts
        charts_unordered.clear();
        for (auto h: v.iter_halfedges()) {
            charts_unordered.insert(source[h.facet()]);
        } 
        if (charts_unordered.size() >= 3) { 
            //on rajoute le point
            m2.points.create_points(1);
            m2.points[num_vertice] = v.pos();
            oldtonew[v] = num_vertice;

             for (auto f: charts_unordered) {
                 sommetsface[f].push_back(v);
                 charts[num_vertice].push_back(f);
                }
            num_vertice++;
            //pour chacun des sommets on les lie si ils ont 2 sources en commun
        }
    }
    
    std::vector<int> racine_of_source(nb_source);//on fait un tableau qui a chaque source renvoie la racine    
    
    //racines
    for (auto f: m.iter_facets()) {
        if (racine[f] != -1) {
            //on rajoute le sommet a m2
            m2.points.create_points(1);
            m2.points[num_vertice] = f.geom<Triangle3>().bary_verts();//bary
            racine_of_source[source[f]] = num_vertice;
            num_vertice++;
        }
    }
    
    
    // creer facet

    m2.connect();
    
    std::vector<int> new_face;
    std::vector<int> side_color;
    std::map<std::pair<int, int>, int> milieu;
    //PointAttribute<std::pair<int, int>> milieu2(m2);

    for (int color = 0; color < nb_source; color++) {
        new_face.clear();
        side_color.clear();
        std::cout << sommetsface[color].size() <<std::endl;
        Surface::Vertex s0(m, sommetsface[color][0]);
        int current_vertice = s0;
        int last_vertice = s0;
        int compte = 0;
        double dist_min = 100;
        Segment3 current_droite;
        std::map<std::pair<int, int>, vec3> milieu1;
        while (compte < sommetsface[color].size() + 1){

            Surface::Vertex current_vertice_v(m, current_vertice);
            //on cherche un halfedge qui suit le bord de la zone
            //pour chaque coté on cherche le point sur la frontiere le plus proche de la droite source1 : source2
            for (auto h: (current_vertice_v.halfedge()).iter_sector_halfedges()) {
                if (source[h.opposite().facet()] == color && source[h.facet()] != color && h.to() != last_vertice) {
                    if (std::find(sommetsface[color].begin(), sommetsface[color].end(), current_vertice) != sommetsface[color].end()) {
                        new_face.push_back(oldtonew[current_vertice]);
                        side_color.push_back(source[h.facet()]);
                        Segment3 new_droite(m2.points[racine_of_source[color]], m2.points[racine_of_source[source[h.facet()]]]);
                        current_droite = new_droite;
                        compte++;
                        dist_min = 100;
                        if (color == 6) { std::cout <<"colorside" << source[h.facet()]<<std::endl; }
                    }
                    last_vertice = current_vertice;
                    current_vertice = h.to();
                    
                    //on met a jour le pt le plus proche
                    if (milieu1.find({color, source[h.facet()]}) == milieu1.end() ||
                        (current_droite.nearest_point(h.to().pos()) - h.to().pos()).norm()
                        < (current_droite.nearest_point(milieu1[{color,source[h.facet()]}])- milieu1[{color,source[h.facet()]}]).norm()
                    ) 
                    {
                        
                        milieu1[{color,source[h.facet()]}] = h.to().pos();
                    }
                    if (source[h.facet()] == 7 && color == 6) {std::cout << "entre six et sept " << milieu1[{color,source[h.facet()]}]<< std::endl; }
                    break;
                }
            }
        }
        
        // m2.conn->create_facet(new_face.data(), new_face.size()); //sans mp
        std::vector<int> new_quad;
        for (int i= 0; i < sommetsface[color].size(); i++) {
            std::cout << "nouveau quad"<< std::endl;
            std::cout << "ccouleur " <<color <<std::endl;
            std::cout << "voisin1 " << side_color[(sommetsface[color].size() + i - 1) % sommetsface[color].size()] <<std::endl;
            std::cout << "voisin2 " << side_color[i] <<std::endl;
             
            new_quad.clear();
            new_quad.push_back(racine_of_source[color]);
            //on cherche ce point sur la frontière
            new_quad.push_back(mid(milieu, milieu1, color, side_color[(sommetsface[color].size() + i - 1) % sommetsface[color].size()], num_vertice, m2));
            new_quad.push_back(new_face[i]);
            new_quad.push_back(mid(milieu, milieu1, color, side_color[i], num_vertice, m2));
            for (auto v: new_quad) {
                std::cout << m2.points[v] <<std::endl;
            }
            m2.conn->create_facet(new_quad.data(), 4);
        }
    }
    
    m2.disconnect();

    write_by_extension(sortie, m2, {{/* {"milieu2", milieu2.ptr} */}, {}, {}});
    int result = system((getGraphitePath() + " " + sortie).c_str());
    return 0;
}

