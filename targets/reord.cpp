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
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <TPZSimpleTimer.h>
#include "pzvisualmatrix.h"
#include "TPZSYSMPMatrix.h"

enum EMatid {ENone,EDomain,EBC};

TPZGeoMesh* CreateGMesh(int ndiv);
TPZCompMesh* CreateFluxCMesh(TPZGeoMesh* gmesh, const int pord);
TPZCompMesh* CreatePressureCMesh(TPZGeoMesh* gmesh, const int pord);
TPZCompMesh* CreateMPMesh(TPZManVector<TPZCompMesh*,2> &meshvec);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);

auto exactSol = [](const TPZVec<REAL> &loc,
                   TPZVec<STATE>&u,
                   TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];
    u[0] = 1.;
    gradU(0,0) = 0.;
    gradU(1,0) = 0.;
};

auto rhsfunc = [](const TPZVec<REAL> &loc,
                   TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];
    u[0] = 10000.;
};


int main() {
#ifdef PZ_LOG
//    TPZLogger::InitializePZLOG("log4cxx.cfg");
    TPZLogger::InitializePZLOG();
#endif
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    
    // Create gmesh
    int ndiv = 60;
    const int pord = 4;
    TPZGeoMesh* gmesh = CreateGMesh(ndiv);
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    // Create compmeshes
    TPZCompMesh* cmeshflux = CreateFluxCMesh(gmesh,pord);
    TPZCompMesh* cmeshp = CreatePressureCMesh(gmesh,pord);
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

TPZCompMesh* CreateFluxCMesh(TPZGeoMesh* gmesh, const int pord) {
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
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

TPZCompMesh* CreatePressureCMesh(TPZGeoMesh* gmesh, const int pord) {
    gmesh->ResetReference();
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    
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
//    mpcmesh->SetDefaultOrder(1);
    mpcmesh->SetDimModel(dim);
    
    // Create domain mat
//    TPZMixedDarcyFlow* mat = new TPZMixedDarcyFlow(EDomain, dim);
    TPZMixedModelProblem* mat = new TPZMixedModelProblem(EDomain, dim);
    mat->SetConstantPermeability(1.);
    mat->SetForcingFunction(rhsfunc, 4);
    mpcmesh->InsertMaterialObject(mat);
    
    
    // BCs
    const int diri = 0;
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    auto* BCCond = mat->CreateBC(mat, EBC, diri, val1, val2);
//    BCCond->SetForcingFunctionBC(rhsfunc, 4);
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
    constexpr int nThreads{16};
//    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    TPZSimpleTimer time_ass;
    an.Assemble();
    std::cout << "Total time = " << time_ass.ReturnTimeDouble()/1000. << " s" << std::endl;

    
    const bool isUsePermVecFromFile = false;
    TPZMatrix<STATE>* mat = an.MatrixSolver<STATE>().Matrix().operator->();
    TPZSYsmpMatrix<STATE>* pardisomat = dynamic_cast<TPZSYsmpMatrix<STATE>*>(mat);
    if(isUsePermVecFromFile) {
        std::ifstream in("permvec.txt");
        int64_t npos = -1;
        in >> npos;
        TPZVec<int64_t> perm(npos);
        for (int i = 0; i < npos; i++) {
            int64_t pos = -1;
            in >> pos;
            perm[i] = pos;
        }
        pardisomat->Permute(perm);
    }
    const int64_t resolution = 300;
    TPZFMatrix<REAL> fillin;
    pardisomat->ComputeFillIn(resolution, fillin);
//    cmesh->ComputeFillIn(resolution, fillin);
    std::string outvisual("visual_mat.vtk");
    VisualMatrix(fillin,outvisual);
    
    ///solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble()/1000. << " s" << std::endl;
    
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
