#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include "TPZMixedModelProblem.h"
#include <TPZNullMaterial.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <TPZMultiPhysicsCompMesh.h>
#include <pzstepsolver.h>
#include <TPZLinearAnalysis.h>


enum EMatid {ENone,EDomain,EBC};

TPZGeoMesh* CreateGMesh(int ndiv);
TPZCompMesh* CreateFluxCMesh(TPZGeoMesh* gmesh);
TPZCompMesh* CreatePressureCMesh(TPZGeoMesh* gmesh);
TPZCompMesh* CreateMPMesh(TPZManVector<TPZCompMesh*,2> &meshvec);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);


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
    TPZCompMesh* cmeshflux = CreateFluxCMesh(gmesh);
    TPZCompMesh* cmeshp = CreatePressureCMesh(gmesh);
    TPZManVector<TPZCompMesh*,2> cmeshvec = {cmeshflux,cmeshp};
    TPZCompMesh* mpcmesh = CreateMPMesh(cmeshvec);
    
    // Analysis
    //Solve Multiphysics
    TPZLinearAnalysis an(mpcmesh,true);
    SolveProblemDirect(an,mpcmesh);
    
    // Post Process
    std::cout << "--------- PostProcess ---------" << std::endl;
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmeshvec, mpcmesh);
    PrintResultsMultiphysic(cmeshvec,an,mpcmesh);

    
    // deleting stuff
    delete mpcmesh;
    delete cmeshflux;
    delete cmeshp;
    delete gmesh;
    
    
    std::cout << "--------- Simulation finished ---------" << std::endl;
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

TPZCompMesh* CreateFluxCMesh(TPZGeoMesh* gmesh) {
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(1);
    cmesh->SetAllCreateFunctionsHDiv();
    
    // Domain null mat
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EDomain);
    cmesh->InsertMaterialObject(mat);
    mat->SetDimension(dim);
    
    // BC null mat
    const STATE imposedpressure = 2.;
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,imposedpressure);
    
    const int diri = 0;
    auto* BCCond0 = mat->CreateBC(mat, EBC, diri, val1, val2);
    cmesh->InsertMaterialObject(BCCond0);
    
    // Constructs mesh
    cmesh->AutoBuild();
    
    return cmesh;
}

TPZCompMesh* CreatePressureCMesh(TPZGeoMesh* gmesh) {
    gmesh->ResetReference();
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(1);
    
    // Setting L2 spaces
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    // Create domain mat
    TPZNullMaterial<>* mat = new TPZNullMaterial<>(EDomain);
    cmesh->InsertMaterialObject(mat);
    
    // Create the compels
    cmesh->AutoBuild();
    
    // Set the lagrange multiplier of this mesh
    const int64_t ncon = cmesh->NConnects();
    const int lag = 1;
    for(int64_t i = 0 ; i < ncon ; i++) {
        TPZConnect& con = cmesh->ConnectVec()[i];
        con.SetLagrangeMultiplier(lag);
    }
    
    return cmesh;
}

TPZCompMesh* CreateMPMesh(TPZManVector<TPZCompMesh*,2> &meshvec) {
    TPZGeoMesh* gmesh = meshvec[0]->Reference();
    const int dim = gmesh->Dimension();
    TPZMultiphysicsCompMesh* mpcmesh = new TPZMultiphysicsCompMesh(gmesh);
    mpcmesh->SetDefaultOrder(1);
    mpcmesh->SetDimModel(dim);
    
    // Create domain mat
//    TPZMixedDarcyFlow* mat = new TPZMixedDarcyFlow(EDomain, dim);
    TPZMixedModelProblem* mat = new TPZMixedModelProblem(EDomain, dim);
    mat->SetConstantPermeability(1.);
    mpcmesh->InsertMaterialObject(mat);
    
    // BCs
    const int diri = 0;
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,2.);
    auto* BCCond = mat->CreateBC(mat, EBC, diri, val1, val2);
    mpcmesh->InsertMaterialObject(BCCond);
    
    // Setting flux and p meshes
    TPZManVector<int,2> activevec(2,1);
    mpcmesh->SetAllCreateFunctionsMultiphysicElem();
    mpcmesh->BuildMultiphysicsSpace(activevec, meshvec);
    mpcmesh->LoadReferences();
    
    return mpcmesh;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{0};
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    an.Assemble();
    
    ///solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    an.Solve();

    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResultsMultiphysic(TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    
    TPZManVector<std::string,10> scalnames(1), vecnames(1);
    const int dim = cmesh->Dimension();
    
    scalnames[0] = "Pressure";
//    scalnames[1] = "ExactPressure";
    vecnames[0]= "Flux";
//    vecnames[1]= "ExactFlux";
    
    int div = 0;
    std::string plotfile = "solutionHdiv.vtk";
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    
    return;
}
