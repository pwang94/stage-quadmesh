
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
    read_by_extension(path + "spherep.geogram", m);
    m.connect();

    FacetAttribute<int> source(m, -1);
    FacetAttribute<int> dist(m, 1000);

    std::vector<int> a_traite; 

    //generate seeds
    const int nb_source = 7;
    /*    int compteur = 0;
    
     for (auto f: m.iter_facets())
        {
            if (f % (m.nfacets()/nb_source) == 3) {
                source[f] = compteur;
                dist[f] = 0;
                compteur++;
                a_traite.push_back(f);
                //std::cout << source[f] << std::endl;
            }
            if (compteur >= nb_source) {break;}
        }
    */
    
    for (int i = 0; i<nb_source; i++) {
        int frand = rand() % m.nfacets();
        if (std::find(a_traite.begin(), a_traite.end(), frand) != a_traite.end()){i--;}
        else {
        source[frand] = i;
        dist[frand] = 0;
        a_traite.push_back(frand);
        }
    }

    std::vector<int> a_traite2;

    int compte_face = 0;

    do {
        for (auto fid: a_traite) { // f à distance i
            compte_face++;
            Surface::Facet f(m, fid);
            for (auto h: f.iter_halfedges()) {
                auto f2 = h.opposite().facet();
                if (dist[f2] > dist[f] + 1) {
                    dist[f2] = dist[f] + 1;
                    a_traite2.push_back(f2);
                    source[f2] = source[f];
                }
            }
        }
        a_traite = a_traite2;
        a_traite2.clear();
    } while (compte_face < m.nfacets());

    std::cout << "test" << std::endl;

    write_by_extension("sphere2.geogram", m, {{}, {{"source", source.ptr}, {"dist", dist.ptr}}, {} });
    int result = system((getGraphitePath() + " sphere2.geogram").c_str());
    // --- END ---

    return 0;
}

