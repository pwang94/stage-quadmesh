
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
    read_by_extension(path + "charts_to_polygon", m);

    write_by_extension("sphere2.geogram", m, {{}, {{"source", source.ptr}, {"dist", dist.ptr}}, {} });
    int result = system((getGraphitePath() + " sphere2.geogram").c_str());
    // --- END ---

    return 0;
}

