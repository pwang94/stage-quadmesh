#include "helpers.h"
#include <ultimaille/all.h>
#include <map>
#include <fstream>

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

    SurfaceAttributes attrs = read_by_extension(path + "\\..\\Debug\\outpoly.geogram", m);
    static const std::string sortie = "poly_annot.geogram";

    m.connect();
    PointAttribute<int> pt_singu("pt_singu", attrs, m);
    PointAttribute<bool> coins("coins_new", attrs, m);
    
    FacetAttribute<int> facet_obj_sommet(m, 4);
    // FacetAttribute<bool> facet_with_coins("facet_with_coins", attrs, m);
    PointAttribute<bool> hard_point(m, false);
    CornerAttribute<bool> hard_edges_attr(m, false);
    find_hard_edges(m, hard_edges_attr, hard_point, 0.6);

    // for (auto f: m.iter_facets()) {
    //     if (facet_with_coins[f]) facet_obj_sommet[f]++;
    // }
    //on note 2 si c'est angle pi 1 si c'est pi/2

    //on veut chaque sommet d'indice total 4

    // on fait par face on donne 1 au 4 meilleurs 2 au reste
    // si il y ambiguité on le fait plus tard, ambiguité à partir de 3pi/4
    // plus tard: on prend les meilleurs selon qqch qui prend en compte que indice 4 c'est mieux


    PointAttribute<int> compte_sommet(m);

    FacetAttribute<std::vector<std::pair<double, int>>> angles(m);
    FacetAttribute<std::vector<int>> annotation(m);
    CornerAttribute<int> annotations(m);
    FacetAttribute<int> valide(m, 0);



    for (auto f: m.iter_facets()) {
        std::vector<std::pair<double, int>> angle;
        if (f.size() <= 3) continue; 
        if (f == 190) std::cout << "190 eme FACE" << std::endl;

        annotation[f] = std::vector<int>(f.size(), -1);
        // if (facet_with_coins[f]) {
        //     //on met à 3 le coin concave
        //     for (auto h: f.iter_halfedges()) {
        //         if (coins[h.from()]) {
        //             annotation[f][h.id_in_facet()] = 3;
        //             annotations[h] = 3;
        //         }
        //     }
        // }

        for (int i = 0; i <f.size(); i++) {
            //on calcule l'angle
            vec3 v1 = f.vertex(i).pos() - f.vertex((i + 1) % f.size()).pos();
            v1 = v1/v1.norm();
            vec3 v2 = f.vertex(i).pos() - f.vertex((i - 1 + f.size()) % f.size()).pos();
            v2 = v2/v2.norm();
            angle.push_back({v1 * v2, i});
        }

        // on regarde si les 4 meilleurs sont bien meilleurs que le reste

        std::sort(angle.begin(), angle.end(), 
        [](auto p1, auto p2) {

            return abs(p1.first) < abs(p2.first);
        } );
        //heuristique comme ça 
        std::cout << "test0"<< std::endl;
        std::cout << f.size()<< std::endl;
        if ((f.size() == 4 || (abs(angle[3].first) < 0.3 )|| (abs(angle[3].first) < 0.50 && abs(angle[4].first) > 0.7))
            /* && !facet_with_coins[f] */) {
            std::cout << "test2"<< std::endl;
            valide[f] = 4;
            for (int i= 0; i < 4; i ++) {
                compte_sommet[f.vertex(angle[i].second)]++;
                //annotation[f][i] = 1;
                //annotations[f.halfedge(i)] = 1;
                annotation[f][angle[i].second] = 1;
                annotations[f.halfedge(angle[i].second)] = 1;
            }
            for (int i= 4; i < f.size(); i ++) {
                compte_sommet[f.vertex(angle[i].second)] += 2;
                annotation[f][angle[i].second] = 2;
                annotations[f.halfedge(angle[i].second)] = 2;

            }
        }
        std::cout << "test3"<< std::endl;
        angles[f] = angle;
    }
    // on regarde les faces qui n'ont pas été validée

    //on annote directe tous les angles ou c'est vraiment proche de 90 ou 180
    for (auto f: m.iter_facets()) {
        if (f.size() <= 3 || valide[f] == facet_obj_sommet[f]) continue;
        if (f == 190) std::cout << "190 PASSAGE 2" << std::endl;
        for (auto& [alpha, num]: angles[f]) {
            if (f == 190) std::cout << "alpha" << alpha << std::endl;

            if (abs(alpha) < 0.3 && valide[f] < 4 && annotation[f][num] == -1) {
                annotation[f][num] = 1;
                annotations[f.halfedge(num)] = 1;
                valide[f]++;
                compte_sommet[f.vertex(num)]++;
            }
            if (abs(alpha) > 0.8 && annotation[f][num] == -1 && compte_sommet[f.vertex(num)] < pt_singu[f.vertex(num)] - 1) {
                annotation[f][num] = 2;
                annotations[f.halfedge(num)] = 2;

                compte_sommet[f.vertex(num)] +=2;
            }
        }
        if (valide[f] == facet_obj_sommet[f]) {
            for (int i = 0; i < f.size(); i++) {
                if (annotation[f][i] == -1) {
                    annotation[f][i] = 2;
                    annotations[f.halfedge(i)] = 2;
                    
                    compte_sommet[f.vertex(f.vertex(i))]+=2;
                }
            }
        }
    }

    for (auto f: m.iter_facets()) {
        if (f.size() <= 3 || valide[f] >= facet_obj_sommet[f]) continue;
        // on prend dans l'ordre si c'est légal
        auto candidats = angles[f];
        int pris_nb = 0;
        for (int i = 0; i< f.size(); i++) {  
            auto [angle, vert_num] =candidats[i];
            //on doit prendre les sommets qui ont une valence tq +2 ça fait dépasser sinon :<

            if (compte_sommet[f.vertex(vert_num)] + 2 > pt_singu[f.vertex(vert_num)] && annotation[f][vert_num] == -1) {
                // il est pris 
                pris_nb += 1;
                annotation[f][vert_num] = 1; 
                annotations[f.halfedge(vert_num)] = 1;

                valide[f]++;
                compte_sommet[f.vertex(vert_num)]++;
            }
            if (valide[f] >= 4) break;
        }
        //le compte n'est pas encore bon mais on n'a pas tout essayé on prend le reste dans l'ordre
        int i = 0;

        while (valide[f] <  facet_obj_sommet[f] && i < f.size()) {
            if (annotation[f][candidats[i].second] == -1) {
                valide[f]++;
                annotation[f][candidats[i].second] = 1; 
                annotations[f.halfedge(candidats[i].second)] = 1; 
                compte_sommet[f.vertex(candidats[i].second)]++;
            }
            i++;
        }
        for (int i = 0; i < annotation[f].size(); i++) {
            if (annotation[f][i] == -1) {
                annotation[f][i] = 2;
                annotations[f.halfedge(i)] = 2; 
                if (f == 124) {std::cout << i << " = 1" << std::endl;};

                compte_sommet[f.vertex(i)]+=2;
            }
        }
    }

    // si il y a un sommet à 3 et un autre à 5 et celui à trois est 90 et celui à 5 180 on peut les échanger

    for (auto f: m.iter_facets()) {
        int a3 = -1;
        int a5 = -1;
        for (int i = 0; i< f.size(); i++) {
            if (compte_sommet[f.vertex(i)] < pt_singu[f.vertex(i)] && annotation[f][i] == 1) {
                a3 = i;
            }
            if (compte_sommet[f.vertex(i)] > pt_singu[f.vertex(i)] && annotation[f][i] == 2) {
                a5 = i;
            }
        }
        if (a5 != -1 && a3 != -1) {
            annotation[f][a3] = 2;
            annotation[f][a5] = 1;
            annotations[f.halfedge(a3)] = 2; 
            annotations[f.halfedge(a5)] = 1; 

            compte_sommet[f.vertex(a5)]--;
            compte_sommet[f.vertex(a3)]++;
        }
    }

    //si il y a un sommet avec trop d'angle et il a un angle 2 dans une face qui à pas assez d'angle 1 on met un angle 1

    for (auto f: m.iter_facets()) {
        if (valide[f] >= 4) continue;
        for (int i = 0; i < f.size(); i++) {
            if (f == 190) {
                
            std:: cout << "annot"<< annotation[f][i] << std::endl;
            std:: cout << "somm"<< compte_sommet[f.vertex(i)] << std::endl;
            std:: cout << "sing"<< pt_singu[f.vertex(i)] << std::endl;}
            if (annotation[f][i] == 2 && compte_sommet[f.vertex(i)] > pt_singu[f.vertex(i)]) {

                annotation[f][i] = 1;
                compte_sommet[f.vertex(i)]--;
                valide[f]++;
                if (valide[f] == 4) break;
            }
        }
    }

    std::ofstream file("file_annotation.txt");
    // for (auto f: m.iter_facets()) {
    //     file << f << "\n" ;
    //     for (int i = 0; i< annotation[f].size(); i++) {
    //         file << annotation[f][i] << " ";
    //     }
        
    //     file << "\n";
    // }
    // file.close();

    CornerAttribute<vec2> mapquad(m);

    for (auto f: m.iter_facets()) {
        int i = 0;
        int compte = 0;
        file << f << "\n" ;

        while (i < 4 && compte < f.size()) {
            if (annotation[f][compte] == 1) {
                if (i == 0) mapquad[f.halfedge(compte)] = vec2(0, 1);  
                if (i == 1) mapquad[f.halfedge(compte)] = vec2(1, 1);  
                if (i == 2) mapquad[f.halfedge(compte)] = vec2(1, 0);  
                if (i == 3) mapquad[f.halfedge(compte)] = vec2(0, 0);  
                i++;
                file << " coins: " << f.halfedge(compte).from();
                
            }
            compte++;
        }
        file << "\n";

    }
    file.close();

    // write_by_extension("poly_annot.geogram", m, {{}, {}, {{"annotations", annotations.ptr}} });
    write_by_extension(sortie, m, {{{"compte", compte_sommet.ptr}, }, {{"compte_face", valide.ptr}}, {{"annotations", annotations.ptr}, {"map", mapquad.ptr}} });


    int result = system((getGraphitePath() + " " + sortie).c_str()); 
    // --- END ---
    
    return 0;
}