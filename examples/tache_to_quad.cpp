#include "helpers.h"
#include <ultimaille/all.h>
#include <map>
#include <set>

using namespace UM;

// bool deuxcommun(std::vector<int> v1, std::vector<int> v2) {
//     int compte = 0;
//     for (int c : v2) {
//         if (std::find(v1.begin(), v1.end(), c) != v1.end()) {compte++;}
//     }
//     return (compte >= 2);
// }


// int encommun(std::vector<int> v1, std::vector<int> v2, int f) {
//     for (int c : v2) {
//         if (std::find(v1.begin(), v1.end(), c) != v1.end() && c != f) {return c;};
//     }
//     assert(false);
//     return 0;
// }

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

// int mid(std::map<std::pair<int, int>, int> &milieu, int s1, int s2, int &compte, Polygons& m) {
//     if (milieu.find({s1, s2}) == milieu.end()) {
//         m.points.create_points(1);
//         m.points[compte] = (m.points[s1] + m.points[s2]) / 2;
//         milieu[{s1, s2}] = compte;
//         compte++;
//         std::cout << "nouveau pt" << std::endl;
//     }
//     return milieu[{s1, s2}];
// }


int main(int argc, char** argv) {
    // --- LOAD ---

    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Triangles m;
    static const std::string entree = "\\..\\Debug\\sortie_disque.geogram";
    static const std::string ff_path = "\\..\\Debug\\framefield.geogram";
    static const std::string sortie = "outpoly.geogram";

    //il faut avoir source, nb_source
    SurfaceAttributes attrs = read_by_extension(path + entree, m);
    m.connect();
    
    
    Triangles mff; //framefield
    SurfaceAttributes attrs_framefield = read_by_extension(path + ff_path, mff);
    FacetAttribute<vec3> frames0("frames0", attrs_framefield, mff);
    FacetAttribute<vec3> frames1("frames1", attrs_framefield, mff);

    FacetAttribute<int> source("source", attrs, m);
    FacetAttribute<int> racine("racine", attrs, m); 

    PointAttribute<int> coins_attr("coins_attr", attrs, m);
    FacetAttribute<int> coins_concave("coins concave", attrs, m);

    int nb_source = 0;

    for (auto f:m.iter_facets()) {
        if (source[f] + 1>nb_source) {nb_source = source[f] + 1;}
    }
    

    //generation de la nouvelle surface
    Polygons m2;
    
    std::vector<std::vector<int>> sommetsface(nb_source + 5);
    FacetAttribute<std::vector<int>> sommets_face(m2);
    std::map<int, int> oldtonew;
    PointAttribute<std::vector<int>> charts(m2.points);
    std::set<int> charts_unordered; // ancien -> nouveau
    int num_vertice = 0;

    FacetAttribute<bool> facet_with_coins(m2);
    PointAttribute<bool> coins_new(m2);



    FacetAttribute<bool> singu(m2);
    FacetAttribute<int>  index_singu(m);
    PointAttribute<int> singularite("singularite", attrs, m);

    std::map<int, int> coin_of_source; //source ->singu
    std::map<int, int> singu_of_source; //source ->singu
    std::map<int, int> singu_vertice_of_source; //source ->singu

    PointAttribute<int> pt_singu(m2, 4);

    for (auto v: m.iter_vertices()) {
        if (singularite[v] != 0) {
            index_singu[v.halfedge().facet()] = singularite[v];
            singu_of_source[source[v.halfedge().facet()]] = v.halfedge().facet();
            singu_vertice_of_source[source[v.halfedge().facet()]] = v;
        }
    }
    std::cout << "test130" << std::endl; 


    for (auto f: m.iter_facets()) {
        if (coins_concave[f]) {
            coin_of_source[source[f]] = f; 
        }
    }

    //on fait un tableau qui à chaque futur sommet donne toute ses sources voisines

    for (auto v: m.iter_vertices()) {
        //check si c'est une intersection de 3 charts
        charts_unordered.clear();
        for (auto h: v.iter_halfedges()) {
            charts_unordered.insert(source[h.facet()]);
        } 

        if (charts_unordered.size() >= 3) { 
            //on rajoute le point
            //on vérif si c'est un coin
            m2.points.create_points(1);
            m2.points[num_vertice] = v.pos();
            oldtonew[v] = num_vertice;
            
            if (coins_attr[v]) {
                coins_new[num_vertice] = true;
            }

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
    // for (auto f: m.iter_facets()) {
    //     if (racine[f] != -1) {
    //         //on rajoute le sommet a m2
    //         m2.points.create_points(1);
    //         m2.points[num_vertice] = f.geom<Triangle3>().bary_verts();//bary
    //         racine_of_source[source[f]] = num_vertice;
    //         num_vertice++;
    //     }
    // }
    
    
    

    m2.connect();
    
    std::vector<int> new_face;
    std::map<std::pair<int, int>, int> milieu;
    int taille = 0;
    // const int taille_min = 4;
    //TODO ici: supprimer les petites arêtes
    std::map<std::pair<int, int>, int> mid_se;

    for (int color = 0; color < nb_source; color++) {
        new_face.clear();
        if (sommetsface[color].size() == 0) {continue;}
        Surface::Vertex s0(m, sommetsface[color][0]);
        int current_vertice = s0;
        int last_vertice = s0;
        int compte = 0;
        taille = 0;
        int current_color;
        while (compte < sommetsface[color].size()){
            Surface::Vertex current_vertice_v(m, current_vertice);
            //on cherche un halfedge qui suit le bord de la zone
            for (auto h: (current_vertice_v.halfedge()).iter_sector_halfedges()) {
                if (source[h.opposite().facet()] == color && source[h.facet()] != color && h.to() != last_vertice) {
                    if (std::find(sommetsface[color].begin(), sommetsface[color].end(), current_vertice) != sommetsface[color].end()) {
            
                        // if (taille < taille_min && taille != 0) {
                        //     //on change le dernier sommet par la moyenne avec le nouveau
                        //     //on a un dico qui donne le point qu'on garde
                        //     //après on remplace toute les faces
                        //     if (mid_se.find({color, current_color}) == mid_se.end()) {
                        //         m2.points.create_points(1);
                        //         m2.points[num_vertice] = (current_vertice_v.pos() + m2.points[new_face[new_face.size() - 1]])/2;
                        //         std::cout << "test2" <<std::endl;

                        //         mid_se[{color, current_color}] = num_vertice;
                        //         num_vertice++;
                        //     }
                        //     new_face.pop_back();
                        //     new_face.push_back(mid_se[{color, current_color}]);
                        // }
                        // else {new_face.push_back(oldtonew[current_vertice]);}
                        // if (color == 68) {std::cout << "sommets " <<oldtonew[current_vertice] << std::endl;}

                        new_face.push_back(oldtonew[current_vertice]);
                        compte++;
                        taille = 0;
                        current_color = source[h.facet()];
                    }
                    taille++;
                    last_vertice = current_vertice;
                    current_vertice = h.to();
                    break;
                }
            }
        }
        //on remplit sommet_face pour recentre.cpp
        m2.conn->create_facet(new_face.data(), new_face.size()); //sans mp
        if (singu_of_source.find(color) != singu_of_source.end()) {
            singu[m2.nfacets()-1] = true;}
        if (coin_of_source.find(color) != coin_of_source.end()) facet_with_coins[m2.nfacets()-1] = true;
        // std::vector<int> new_quad;
        // for (int i= 0; i < new_face.size(); i++) {
        //     std::cout << "nouveau quad"<< std::endl;
            
        //     new_quad.clear();
        //     new_quad.push_back(racine_of_source[color]);
        //     new_quad.push_back(mid(milieu, new_face[(new_face.size() + i - 1) % new_face.size()], new_face[i], num_vertice, m2));
        //     new_quad.push_back(new_face[i]);
        //     new_quad.push_back(mid(milieu, new_face[i], new_face[(i+1)%new_face.size()], num_vertice, m2));
        //     m2.conn->create_facet(new_quad.data(), 4);
        // }
    }











    //faire les subdivisions des charts singus
    
    PolyLine test;
    int num_vert = 0;
    std::vector<std::vector<vec3>> lignes_de_source(nb_source);
    for (auto v: m.iter_vertices()) {
        if (singularite[v] != 0) {
            // on regarde tous les frames et prends les cinq directions les plus différentes
            std::vector<vec3> directions;
            for (auto h: v.halfedge().iter_sector_halfedges()) {
                // parmi les deux frames on prend celle qui s'aligne le plus avec v - v_bary_facet
                Surface::Facet current_facet = h.facet();
                vec3 direction = - v.pos() + current_facet.geom<Poly3>().bary_verts();
                if (abs(frames0[current_facet] * direction) > abs(frames1[current_facet] * direction)) {
                    if (frames0[current_facet] * direction > 0) directions.push_back(frames0[current_facet]);
                    else directions.push_back(-frames0[current_facet]);
                }
                else {
                    if (frames1[current_facet] * direction > 0) directions.push_back(frames1[current_facet]);
                    else directions.push_back(-frames1[current_facet]);
                }
            }

            while (directions.size() != singularite[v]) {
                // on supprime le plus proche d'un autre ou on prend la moyenne des 2 plus proche ?
                // pas les plus proche mais les plus proches apres rotation 
                double max = -1;
                std::pair<int,int> couple_min;
                for (int i = 0;  i < directions.size(); i++) {
                    for (int j = 0;  j < directions.size(); j++) {
                        if (i >= j) {continue;}
                        if (directions[i] * directions[j] > max) {max = abs(directions[i] * directions[j]);couple_min = {i, j};}
                    } 
                }
                directions.push_back((directions[couple_min.first] + directions[couple_min.second])/2);
                directions.erase(directions.begin() + couple_min.second);
                directions.erase(directions.begin() + couple_min.first);
            }
            test.points.create_points(1);
            int centre = test.points.size() - 1;
            test.points[centre] = v.pos();
            lignes_de_source[source[v.halfedge().facet()]] = directions;
            for (auto dir: directions) {
                test.points.create_points(1);
                test.points[test.points.size() - 1] = v.pos() + dir*0.08;
                test.create_edges(1);
                test.vert(num_vert, 0) = centre;
                test.vert(num_vert, 1) = test.points.size() -1;
                num_vert++;
            }
        }
    }


    

    // on subdivise les poly singu selon ces lignes
    std::vector<int> a_supprimer_au_plus_vite;
    std::vector<std::vector<int>> a_creer_au_plus_vite;
    std::vector<std::pair<int,Surface::Halfedge>> a_recoller;
    std::map<int, int> new_edge;

    m2.disconnect();
    m2.connect();
    PointAttribute<int> test_ordre(m2);
    int facet_depart = m2.nfacets();





    for (int n = 0;  n < facet_depart;n++) {
        Surface::Facet f(m2, n);
        
        if (singu_of_source.find(f) != singu_of_source.end()) {
            //il faut trouver les cotes sur lesquelles se coupent les directions ligne_de_source[f]
            
            Surface::Vertex v(m, singu_vertice_of_source[f]);
            
            m2.points.create_points(1);
            m2.points[m2.points.size()-1] = v.pos();
            int v_num = m2.points.size()-1;
            pt_singu[v_num] = singularite[v];


            if (singularite[v] >= f.size()) {} //oups plus d'une subdivision par faceon 
            // on récupère le nouveaux pt de la subdivision 
            //en prenant le pt le plus proche d'une arête de v.pos() + bcp * direction
            std::map<int, bool> deja_pris;
            std::map<int, bool> deja_pris_vraiment;
            for (auto h: f.iter_halfedges()) {
                deja_pris[h] = false;
                deja_pris_vraiment[h] = false;
            }
            
            std::vector<std::pair<int, int>> new_points; // num_pt, he
            for (vec3 dir: lignes_de_source[f]) {
                vec3 pt;
                vec3 sub_ici;
                Surface::Halfedge h_choisit = f.halfedge();
                
                // on cherche le points sur la frontière qui coup la direction !
                
                double dist_min = 1000;
                for (auto h: f.iter_halfedges()) {
                    Segment3 s({h.to(), h.from()});
                    for (int x= 1; x< 50; x++) {
                        pt = v.pos() + 0.01*x*dir;
                        vec3 v_subd = s.nearest_point(pt);
                        if ((v_subd - pt).norm() <dist_min) {
                            dist_min = (v_subd - pt).norm();
                            sub_ici = v_subd;
                            h_choisit = h;
                        }
                    }
                }

                // std::cout << "new_dir" << std::endl;
                // //on veut eviter les petites arêtes
                if (deja_pris[h_choisit] /* || (h_choisit.to().pos() - h_choisit.from().pos()).norm() <0.08 */) {
                    //on cherche une autre arête qui convient
                    dist_min = 1000;
                    for (auto h: f.iter_halfedges()) {
                        if (!(deja_pris[h] /* || (h.to().pos() - h.from().pos()).norm() <0.08 */))  {
                            Segment3 s({h.to(), h.from()});
                            for (int x= 1; x<50; x++) {
                                pt = v.pos() + 0.01*x*dir;
                                vec3 v_subd = s.nearest_point(pt);
                                if ((v_subd - pt).norm() <dist_min) {
                                    dist_min = (v_subd - pt).norm();
                                    sub_ici = v_subd;
                                    h_choisit = h;
                                }
                            }
                        }
                    }
                }
                if (!deja_pris_vraiment[h_choisit] || true) {
                    // si on a un truc droit on le met en déjà pris
                    deja_pris[h_choisit] = true;
                    deja_pris_vraiment[h_choisit] = true;
                    for (auto h: f.iter_halfedges()) {
                        vec3 v1 = h.to().pos() - h.from().pos();
                        v1 = v1/v1.norm();
                        vec3 v2 = h_choisit.to().pos() - h_choisit.from().pos();
                        v2 = v2/v2.norm();
                        if ((h.from() == h_choisit.to() || h.to() == h_choisit.from()) && abs(v1*v2) > 0.8) {
                            deja_pris[h] = true;
                        }
                    }
                    m2.points.create_points(1);
                    int pts = m2.points.size()-1;
                    m2.points[pts] = (sub_ici + (h_choisit.to().pos() + h_choisit.from().pos())/4)/(3/2.);
                    new_points.push_back({pts, h_choisit});
                    // il faut connecter avec la face opposer h_choisit.opposite
                    a_recoller.push_back({pts, h_choisit.opposite()});
                    new_edge[h_choisit.opposite()] = h_choisit.opposite();
                }

                
            }

            //on crée toute les facets de la face
            //on ordonne new_points :
            //on parcours les arêtes de la face et on 
            Surface::Halfedge current_h(m2, new_points[0].second);
            int compte = 0;
            new_face.clear();
            new_face.push_back(v_num);
            new_face.push_back(new_points[0].first);
            while (compte < new_points.size()) {
                new_face.push_back(current_h.to());
                //on obtiens l'arête suivante

                for (auto h: f.iter_halfedges()) {
                    if (h.from() == current_h.to()) {
                        current_h = h;
                        break;
                    }
                }
                // on cherche l'arete
                for (auto& [pt, h]: new_points) {
                    if (h == current_h) {
                        compte++;
                        test_ordre[pt] = compte;
                        //on cree la nouvelle face
                        new_face.push_back(pt);
                        m2.conn->create_facet(new_face.data(), new_face.size());
                        new_face.clear();
                        new_face.push_back(v_num);
                        new_face.push_back(pt);
                    }
                }
            }
            // new_face.push_back(current_h.to());
            // new_face.push_back(new_points[0].first);
            m2.conn->create_facet(new_face.data(), new_face.size());
            if (facet_with_coins[f]) facet_with_coins[m2.nfacets()-1] = true;
            
            m2.conn->active[f] = false; 
        }
    }


    // on veut découper les points concave avec les segment  -3 et +3

    for (auto f: m2.iter_facets()) {
        if (facet_with_coins[f]) {
            //on découpe en 3 quad
            int i0;
            for (int i = 0; i < f.size(); i++) {
                if (coins_new[f.vertex(i)]) {
                    i0 = i;
                }
            }
            // on veut prendre le segment +3 et -3 sans compter les segments plats
            int i_curr = i0;
            new_face.clear();
            vec3 v1;
            vec3 v2;
            vec3 v1_n;
            vec3 v2_n;
            // 1er cote
            int compte = 0;
            std :: cout << f << " f" << std::endl;
            do {
                new_face.push_back(f.vertex(i_curr% f.size()));
                std::cout << "sommet " << f.vertex(i_curr%f.size()) << std::endl;
                std::cout << i_curr % f.size() << std::endl;
                std::cout << f.size() << std::endl;
                std::cout << f.halfedge(i_curr%f.size()).from() << std::endl;
                std::cout << f.halfedge(i_curr%f.size()).to() << std::endl;
                v1 = f.halfedge(i_curr % f.size()).to().pos() - f.halfedge(i_curr % f.size()).from().pos();
                // std::cout << v1 << std::endl;
                v1_n = v1/ v1.norm();
                v2 = f.halfedge((i_curr + 1)% f.size()).to().pos() - f.halfedge((i_curr + 1) % f.size()).from().pos();
                // std::cout << v2 << std::endl;
                v2_n = v2/ v2.norm();
                std::cout << v1 * v2 << std::endl;
                if (v1_n * v2_n < 0.9) compte++;
                i_curr++;
            } while (compte < 2);

            new_face.push_back(f.vertex(i_curr % f.size()));
            // milieu du 3eme coté il faudra ensuite fixer l'autre bord parce qu'on crée un tmesh
            m2.points.create_points(1);
            m2.points[m2.points.size() - 1] = (f.halfedge((i_curr)% f.size()).to().pos() + f.halfedge((i_curr) % f.size()).from().pos())/2;
            new_face.push_back(m2.points.size() - 1);
            m2.conn->create_facet(new_face.data(), new_face.size());
            // on fait la deuxieme face:
            new_face.clear();
            new_face.push_back(f.vertex(i0 % f.size()));
            new_face.push_back(m2.points.size()-1);
            
            a_recoller.push_back({m2.points.size() - 1, f.halfedge((i_curr)% f.size()).opposite()});
            i_curr++;
            std::cout << "deuxieme face \n" <<std::endl;
            do {
                new_face.push_back(f.vertex(i_curr % f.size()));
                std::cout << "sommet " << f.vertex(i_curr%f.size()) << std::endl;
                // std::cout << i_curr % f.size() << std::endl;
                // std::cout << f.size() << std::endl;
                // std::cout << f.halfedge(i_curr%f.size()).from() << std::endl;
                // std::cout << f.halfedge(i_curr%f.size()).to() << std::endl;
                v1 = f.halfedge(i_curr % f.size()).to().pos() - f.halfedge(i_curr % f.size()).from().pos();
                // std::cout << v1 << std::endl;
                v1_n = v1/ v1.norm();
                v2 = f.halfedge((i_curr - 1 + f.size())% f.size()).to().pos() - f.halfedge((i_curr - 1 +f.size()) % f.size()).from().pos();
                // std::cout << v2 << std::endl;
                v2_n = v2/ v2.norm();
                // std::cout << v1 * v2 << std::endl;
                // std::cout << v1_n * v2_n << std::endl;
                i_curr++;
            } while (v1_n * v2_n > 0.9);
            i_curr--;
            m2.points.create_points(1);
            m2.points[m2.points.size() - 1] = (f.halfedge((i_curr)% f.size()).to().pos() + f.halfedge((i_curr) % f.size()).from().pos())/2;
            new_face.push_back(m2.points.size() - 1);
            std::cout << "sommet " << m2.points.size() - 1 << std::endl;
            std::cout << "nb sommet " << new_face.size() << std::endl;
            std::cout << m2.nverts() << std::endl;

            m2.conn->create_facet(new_face.data(), new_face.size());

            //on fait la 3eme face
            new_face.clear();
            new_face.push_back(f.vertex(i0 % f.size()));
            new_face.push_back(m2.points.size()-1);
            a_recoller.push_back({m2.points.size() - 1, f.halfedge((i_curr)% f.size()).opposite()});

            i_curr++;
            do {
                new_face.push_back(f.vertex(i_curr % f.size()));
                i_curr++;
            } while (i_curr % f.size() != i0);
            m2.conn->create_facet(new_face.data(), new_face.size());
            m2.conn->active[f] = false; 

        }
    }


    //on fait une map qui donne toute les arêtes à recoller et leur nouvelle arête asso si elle a été suppr
    //on recolle les morceaux
    std::vector<int> sommet_toujours_mauvais;
    PointAttribute<int> sommet_toujours_mauvais_attr(m2.points);
    CornerAttribute<int> tienstiens(m2);
    for (auto& [pts, h]: a_recoller){
        // si jamais !h.active(),  c'est qu'on a supprimé la face soit parce qu'elle est singu soit pendant le recollage
        // il faut donc donner la nouvelle arête
        
        //on supprime l'arête de new_edge, plus besoin:
        new_face.clear();
        Surface::Halfedge new_h(m2, new_edge[h]);
        if (!new_h.active()) {
            sommet_toujours_mauvais.push_back(pts);
            sommet_toujours_mauvais_attr[pts] = true;
            continue;
        }
        Surface::Facet bord = new_h.facet();
        
        new_edge.erase(h);


        int nv_pt;
        for (int i = 0; i < bord.size(); i++) {
            new_face.push_back(bord.vertex(i));
            if (bord.vertex(i) == h.from()) {new_face.push_back(pts); nv_pt = i;}
        }
        // a_creer_au_plus_vite.push_back(new_face);
        m2.conn->create_facet(new_face.data(), new_face.size());
        Surface::Facet new_f(m2, m2.nfacets()-1);
        for (int i = 0; i < bord.size(); i++) {
            if (new_edge.find(bord.halfedge(i)) != new_edge.end()) {
                new_edge[bord.halfedge(i)] = new_f.halfedge(i) + (nv_pt < i);
            }
   
        }
        if (singu_of_source.find(bord) == singu_of_source.end())  {
            // a_supprimer_au_plus_vite.push_back(bord);
            if (facet_with_coins[bord]) facet_with_coins[m2.nfacets()-1] = true;
            m2.conn->active[bord] = false; 

            
        }
    }
    m2.compact();






    //il reste des sommets mal connecté
    
    for (int s: sommet_toujours_mauvais) {

        // on cherche l'arête qui ne contient pas s mais qui doit contenir s en vrai
        for (auto h: m2.iter_halfedges()) {
            vec3 pt = m2.points[s];
            if (h.to() != s && h.from() !=s) {
                Segment3 seg({h.to(), h.from()});
                vec3 proche = seg.nearest_point(pt);
                if ((proche - pt).norm() < 1e-14) {
                    // on a trouver le halfedge il faut mtn le couper en deux supprimer et ajouter la nouvelle face
                    new_face.clear();
                    Surface::Facet bord = h.facet();
                    int num;
                    for (int i = 0; i < bord.size(); i++) {
                        new_face.push_back(bord.halfedge(i).from());
                        if (bord.halfedge(i) == h) {new_face.push_back(s); num = i;}
                    }
                    // a_creer_au_plus_vite.push_back(new_face);
                    m2.conn->create_facet(new_face.data(), new_face.size());
                    if (facet_with_coins[bord]) facet_with_coins[m2.nfacets()-1] = true;

                    m2.conn->active[bord] = false; 
                    m2.compact();
                    break;
                }
            }
        }
    }
    
    for (auto f: m2.iter_facets()) {
        if (f.size() < 3) {
            // std::cout << "ah bon"<< std::endl;
            m2.conn->active[f] = false; 
            
        }
    }
    m2.compact();




    //on voudra mettre un 3 sur les coins concaves 








    // write_by_extension("test_subdivise.geogram", test, {{}, {}});


    for (auto v: m2.iter_vertices()) {
        if (coins_new[v]) {pt_singu[v]++;}
    }

    write_by_extension(sortie, m2, {{{"pt_singu", pt_singu.ptr}, {"coins_new",coins_new.ptr} }, {{"singu", singu.ptr}, {"facet_with_coins", facet_with_coins.ptr}}, {{"tienstiens", tienstiens.ptr}}});
    int result = system((getGraphitePath() + " " + sortie).c_str());
    return 0;
}
