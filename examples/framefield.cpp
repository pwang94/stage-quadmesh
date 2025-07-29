#include "helpers.h"
#include <ultimaille/all.h>
#include <map>
#include <cmath>
// example: construct a quadratic program from given iterators
// the QP below is the first quadratic program example in the user manual

using namespace UM;
using namespace Linear;

void find_hard_edges(Triangles& mesh, CornerAttribute<bool>& hard_edges_attr, PointAttribute<bool>& hard_point ,double threshold) {

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


std::vector<vec3> repere_ref(Surface::Halfedge h1, Surface::Halfedge h2) {
    vec3 v1 = h1.to().pos() - h1.from().pos();
    v1 = v1/(v1.norm());
    vec3 v2 = h2.to().pos() - h2.from().pos();
    v2 -= v1*(v1*v2);
    v2 = v2/(v2.norm());
    return std::vector{v1, v2};
}

int main(int argc, char** argv) {

    // --- LOAD ---
    const double pi = 3.141592653589793;
    // Get path of current executable
    std::string path = getAssetPath();

    // Declare a mesh with triangle surface
    Triangles m;
    static const std::string entree = "joint_sub_sub.geogram";
    static const std::string sortie = "framefield.geogram";
    static const std::string sortiepl = "framefieldpl.geogram"; // pour l'affichage avec un polyline
    read_by_extension(entree, m);
    m.connect();
    //calcul des faces
    CornerAttribute<bool> hard_edges_attr(m, false);
    PointAttribute<bool> hard_point(m, false);
    find_hard_edges(m, hard_edges_attr, hard_point, 0.6);
    FacetAttribute<bool> hard_face(m, false);

    for (auto f : m.iter_facets()) {
        if (hard_edges_attr[f.halfedge(0)]||hard_edges_attr[f.halfedge(1)]||hard_edges_attr[f.halfedge(2)]) {
            hard_face[f] = true;
        }
    }

    std::vector<int> a_traite;
    std::vector<int> a_traite2;
    FacetAttribute<int> zone(m, -1);
    int zone_num = 0;
    std::vector<std::vector<int>> de_zone;
    int compte_face = 0;
    // a chaque zone donne les faces dedans
    // obtenir les grosses zones engendrées par les hards edges

    while (compte_face < m.nfacets()) { 
        int frand;
        de_zone.push_back({});
        do  {
            frand = rand() % m.nfacets();
        } while (zone[frand] != -1);
        // on propage sur toute la zone
        a_traite.push_back(frand);
        while (a_traite.size() != 0) {
            for (int fid: a_traite) {
                if (zone[fid] != -1) continue;
                zone[fid] = zone_num;
                compte_face++;
                de_zone[zone_num].push_back(fid);
                Surface::Facet f(m, fid);
                for (auto h : f.iter_halfedges()) {
                    int f_opp = h.opposite().facet();
                    if (!hard_edges_attr[h] && zone[f_opp] == -1) {
                        a_traite2.push_back(f_opp);
                    }
                }
            }
            a_traite = a_traite2;
            a_traite2.clear();
        }
        zone_num++;
    }

    PolyLine test;

    LeastSquares ls(2*m.nfacets());
    //on fixe les contraintes au bord
    for (auto f: m.iter_facets()) {
        if (hard_face[f]) { // TODO: traiter le cas plusieurs halfedge
            auto repere = repere_ref(f.halfedge(0), f.halfedge(1));
            for (auto h : f.iter_halfedges()) {
                if (hard_edges_attr[h]) {
                    vec3 e = h.to().pos() - h.from().pos();
                    e = e/(e.norm());
                    // il faut avoir 4x l'angle
                    // double theta = fmod(acos(co) * 4, 2 * pi);
                    double theta = atan2(repere[1] * e,repere[0] * e)* 4;
                    ls.fix(2*f, cos(theta));
                    ls.fix(2*f + 1, sin(theta));
                    break;
                } 
            }
        }
    }



    //on calcule les angle pour Rij
    
    for (auto h: m.iter_halfedges()) {
        if (hard_edges_attr[h]) continue;
        Surface::Facet f1 = h.facet();
        Surface::Facet f2 = h.opposite().facet();
        auto repere1 = repere_ref(f1.halfedge(0), f1.halfedge(1));
        auto repere2 = repere_ref(f2.halfedge(0), f2.halfedge(1));

        //TODO tourner le repere1 pour qu'il soit dans le meme plan quel le 2 selon l'arête qui les relie
        // besoin de l'axe et de l'angle
        vec3 axe = h.to().pos() - h.from().pos();
        axe = axe/axe.norm();
        vec3 normal1 = f1.geom<Triangle3>().normal();
        vec3 normal2 = f2.geom<Triangle3>().normal();
        vec3 ve_normal2_axe = {normal2[1] * axe[2] - normal2[2] * axe[1], normal2[2] * axe[0] -normal2[0] * axe[2], normal2[0] * axe[1] -normal2[1] * axe[0]};
        vec3 ve_normal1_axe = {normal1[1] * axe[2] - normal1[2] * axe[1], normal1[2] * axe[0] -normal1[0] * axe[2], normal1[0] * axe[1] -normal1[1] * axe[0]};
        double alpha = atan2(ve_normal2_axe * normal1,ve_normal2_axe*ve_normal1_axe);

        //on applique la rotation à repere2 pour les mettre sur le mm plan
        //rotation:
        double a = cos(alpha/2);
        double b = axe[0]*sin(alpha/2);
        double c = axe[1]*sin(alpha/2);
        double d = axe[2]*sin(alpha/2);
        // double m[3][3];
        // m[0][0] = a*a + b*b - c*c - d*d;
        // m[0][1] = 2 * (b*c - a*d);
        // m[0][2] = 2 * (b*d + a*c);
        // m[1][0] = 2 * (b*c + a*d);
        // m[1][1] = a*a + c*c - b*b - d*d;
        // m[1][2] = 2 * (c*d - a*b);
        // m[2][0] = 2 * (b*d - a*c);
        // m[2][1] = 2 * (c*d + a*b);
        // m[2][2] = a*a + d*d - b*b - c*c;
        // std::vector<vec3> new_repere2;
        // double rx = m[0][0] * repere2[0][0] + m[0][1] * repere2[0][1] + m[0][2] * repere2[0][2];
        // double ry = m[1][0] * repere2[0][0] + m[1][1] * repere2[0][1] + m[1][2] * repere2[0][2];
        // double rz = m[2][0] * repere2[0][0] + m[2][1] * repere2[0][1] + m[2][2] * repere2[0][2];
        // new_repere2.push_back(vec3(rx, ry, rz));
        // rx = m[0][0] * repere2[1][0] + m[0][1] * repere2[1][1] + m[0][2] * repere2[1][2];
        // ry = m[1][0] * repere2[1][0] + m[1][1] * repere2[1][1] + m[1][2] * repere2[1][2];
        // rz = m[2][0] * repere2[1][0] + m[2][1] * repere2[1][1] + m[2][2] * repere2[1][2];
        // new_repere2.push_back(vec3(rx, ry, rz));

        // double theta = acos(abs(repere1[0] * repere2[0])); //cos truc il faut avoir l'angle orienté
        // if (abs(repere1[0] * repere2[0]) < 0) theta = pi - theta;
        // if (repere2[0] * repere1[1] < 0) theta = -theta;
        double theta = atan2(repere2[0]*repere1[1], repere2[0] *repere1[0])*4;

        ls.add_to_energy(-cos(theta)*X(2*f2) + sin(theta)*X(2*f2 + 1) + X(2*f1));
        ls.add_to_energy(sin(theta)*X(2*f2) + cos(theta)*X(2*f2 + 1) - X(2*f1+1));
        
        //test tourner repere1 de theta degré

        if (true) {
            test.points.create_points(3);
            test.points[test.points.size() - 3] = (f1.vertex(0).pos() + f1.vertex(1).pos() + f1.vertex(2).pos()) /3.;
            test.points[test.points.size() - 2] = test.points[test.points.size() - 3] + (repere1[0]*cos(theta) + repere1[1]*sin(theta))*0.08;//repere2[0]
            test.points[test.points.size() - 1] = test.points[test.points.size() - 3] + (-repere1[0]*sin(theta) + repere1[1]*cos(theta))*0.08;//repere2[1]
            test.create_edges(2);
            test.vert(2*test.points.size()/3 -1, 0) = test.points.size() - 3;
            test.vert(2*test.points.size()/3 -1, 1) = test.points.size() - 1;
            test.vert(2*test.points.size()/3 -2, 0) = test.points.size() - 3;
            test.vert(2*test.points.size()/3 -2, 1) = test.points.size() - 2;
        }

    }

    write_by_extension("test_frame.geogram", test, {{}, {}});

    ls.solve();
    
    // normaliser et faire alpha = alpha / 4
    


    FacetAttribute<std::vector<vec3>> frames(m);
    FacetAttribute<vec3> frames0(m);
    FacetAttribute<vec3> frames1(m);

    PolyLine p;
    // display frame
    const double scale = 0.2e-2;
    for (auto f: m.iter_facets()) {
        auto repere = repere_ref(f.halfedge(0), f.halfedge(1));
        p.points.create_points(5);
        p.points[5*f] =  (f.vertex(0).pos() + f.vertex(1).pos() + f.vertex(2).pos()) /3.;
        vec3 frame = repere[0] * ls.value(2*f) + repere[1] * ls.value(2*f + 1);
        vec3 frame2 = repere[0] * ls.value(2*f + 1)  - repere[1] * ls.value(2*f) ;
        //normaliser
        if (std::pow((std::pow(ls.value(2*f),2) +std::pow(ls.value(2*f +1),2)), 1./2) > 1e-20) {
        frame = frame/std::pow((std::pow(ls.value(2*f),2) +std::pow(ls.value(2*f +1),2)), 1./2);
        frame2 = frame2/std::pow((std::pow(ls.value(2*f),2) +std::pow(ls.value(2*f +1),2)), 1./2);
        }
        // faire alpha/4
        double theta = atan2(frame * repere[1], frame * repere[0]);
        theta = fmod(theta/4., 2* pi) ;
        frame = repere[0] * cos(theta) + repere[1] * sin(theta);
        frame2 = repere[0] * sin(theta) - repere[1] * cos(theta);
        frames[f].push_back(frame);
        frames0[f] = frame;
        frames[f].push_back(frame2);
        frames1[f] = frame2;
        frames[f].push_back(-frame);
        frames[f].push_back(-frame2);
        p.points[5*f+1] =  p.points[5*f] + frame * scale;
        p.points[5*f+2] =  p.points[5*f] + frame2 * scale;
        p.points[5*f+3] =  p.points[5*f] - frame * scale;
        p.points[5*f+4] =  p.points[5*f] - frame2 * scale;
        // p.points[3*f+1] =  p.points[3*f] + repere[0] * scale;
        // p.points[3*f+2] =  p.points[3*f] + repere[1] * scale;
        // p.points[3*f+3] =  p.points[3*f] + repere[0] * scale;
        // p.points[3*f+4] =  p.points[3*f] + repere[1] * scale;
        p.create_edges(4);
        p.vert(4*f, 0) = 5*f;
        p.vert(4*f, 1) = 5*f + 1;
        p.vert(4*f + 1, 0) = 5*f;
        p.vert(4*f + 1, 1) = 5*f+2;
        p.vert(4*f + 2, 0) = 5*f;
        p.vert(4*f + 2, 1) = 5*f+3;
        p.vert(4*f + 3, 0) = 5*f;
        p.vert(4*f + 3, 1) = 5*f+4;

    }

    //détection des singularite
    PointAttribute<int> singularite(m, 0);

    for (auto v: m.iter_vertices()) {
        double total = 0;
        Surface::Halfedge h0 = v.halfedge();
        Surface::Facet current_face = h0.facet();
        vec3 current_fr = frames[current_face][0];
        for (auto h : h0.iter_sector_halfedges()) {
            //on prend le meilleur vecteur parmis les 4 et on regarde si c'est le mm à la fin
            double max = 0;
            vec3 newframe;
            for (auto fr: frames[h.facet()]) {
                //on applique la rotation sur fr
                // vec3 axe = h.to().pos() - h.from().pos();
                // axe = axe/axe.norm();
                // vec3 normal1 = current_face.geom<Triangle3>().normal();
                // vec3 normal2 = h.facet().geom<Triangle3>().normal();
                // vec3 ve_normal2_axe = {normal2[1] * axe[2] - normal2[2] * axe[1], normal2[2] * axe[0] -normal2[0] * axe[2], normal2[0] * axe[1] -normal2[1] * axe[0]};
                // vec3 ve_normal1_axe = {normal1[1] * axe[2] - normal1[2] * axe[1], normal1[2] * axe[0] -normal1[0] * axe[2], normal1[0] * axe[1] -normal1[1] * axe[0]};
                // double alpha = atan2(ve_normal2_axe * normal1,ve_normal2_axe*ve_normal1_axe);

                // //on applique la rotation à repere2 pour les mettre sur le mm plan
                // //rotation:
                // double a = cos(alpha/2);
                // double b = axe[0]*sin(alpha/2);
                // double c = axe[1]*sin(alpha/2);
                // double d = axe[2]*sin(alpha/2);
                // float m[3][3];
                // m[0][0] = a*a + b*b - c*c - d*d;
                // m[0][1] = 2 * (b*c - a*d);
                // m[0][2] = 2 * (b*d + a*c);
                // m[1][0] = 2 * (b*c + a*d);
                // m[1][1] = a*a + c*c - b*b - d*d;
                // m[1][2] = 2 * (c*d - a*b);
                // m[2][0] = 2 * (b*d - a*c);
                // m[2][1] = 2 * (c*d + a*b);
                // m[2][2] = a*a + d*d - b*b - c*c;
                // std::vector<vec3> new_repere2;
                // double rx = m[0][0] * fr[0] + m[0][1] * fr[1] + m[0][2] * fr[2];
                // double ry = m[1][0] * fr[0] + m[1][1] * fr[1] + m[1][2] * fr[2];
                // double rz = m[2][0] * fr[0] + m[2][1] * fr[1] + m[2][2] * fr[2];
                // vec3 newfr = vec3(rx, ry, rz);
                if (fr*current_fr > max) {newframe = fr; max = fr*current_fr;};
            }
            current_face = h.facet();
            current_fr = newframe;
        }
        // on regarde si on a la bonne frame à la fin
        for (int i = 1; i<4; i ++) {
            if (current_fr * frames[h0.facet()][i] > current_fr * frames[h0.facet()][0]) {
                // std::cout << current_fr * frames[h0.facet()][i] << " " << current_fr * frames[h0.facet()][0] << std::endl;
                // on note si la frame la plus proche est 1 ou 3 
                if (!hard_point[v]) {
                    singularite[v] = i == 1? 5: 3;
                }
            }
        }
    }

    write_by_extension(sortiepl, p, {{}, {}});
    write_by_extension(sortie, m, {{{"singularite", singularite.ptr}}, {{"frames0", frames0.ptr}, {"frames1", frames1.ptr}}, {}});
    // int result = system((getGraphitePath() + " " + sortie).c_str()); 

}