#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"

enum EMatid {ENone,EDomain,EBC};

TPZGeoMesh* CreateGMesh(int ndiv);

int main() {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    
    // Create gmesh
    int ndiv = 2;
    TPZGeoMesh* gmesh = CreateGMesh(ndiv);
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    // Create compmeshes
    
    // Analysis
    
    // Post Process
    
    // deleting stuff
    delete gmesh;
}

TPZGeoMesh* CreateGMesh(int ndiv) {
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    const TPZManVector<REAL, 3> x0 = {0., 0., 0.};
    const TPZManVector<REAL, 3> x1 = {1., 1., 1.};
    const TPZManVector<int, 3> ndivvec = {ndiv, ndiv, ndiv};
    TPZGenGrid2D gen(ndivvec, x0, x1);
    gen.Read(gmesh);

    gen.SetBC(gmesh, 4, EBC); // bot
    gen.SetBC(gmesh, 5, EBC); // right
    gen.SetBC(gmesh, 6, EBC); // top
    gen.SetBC(gmesh, 7, EBC); // left
    
    return gmesh;
}
