/* interface to TrilinosKokkos

This interface can in principal be used both in parallel and serially,
but currently it is only called in SParIterSolver.
It is compiled into Elmer by adding -DHAVE_TRILINOS
and linking with the TrilinosKokkos libraries. We currently
allow creating an iterative solver (from the Belos library)
and a preconditioner (ifpack or ML). To use a direct solver,
set "Ifpack Method" to "Amesos" and "Iterative Solver" to "None".

Older versions (below 10.10) of TrilinosKokkos demand a different call 
of the XML input file. This can be force by -DOLD_TRILINOS in the
CXXFLAGS 

To activate these linear solvers, set
'Linear System Use TrilinosKokkos' = Logical True
'TrilinosKokkos Input File' = String <xml filename>

see elmerfem/fem/examples/trilinos for an example.

*/

#include "config.h"


#ifdef HAVE_TRILINOS

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef CHECK_ZERO
#undef CHECK_ZERO
#endif

// enable this to store matrices, generate debugging output etc
//#define DEBUG_TRILINOS_INTERFACE
//#define DUMP_IN_TRILINOS_INTERFACE

#define CHECK_ZERO(funcall) {ierr = funcall;\
if (ierr) {std::cout<<"TrilinosKokkos Error "<<ierr<<" returned from call "<<#funcall<<std::endl; return;}}

#define FCHECK_ZERO(funcall) {ierr = funcall;\
if (ierr) {std::cout<<"TrilinosKokkos Error "<<ierr<<" returned from call "<<#funcall<<std::endl; return Teuchos::null;}}

#define ERROR(msg,file,line)  {std::cerr << "Error in file "<<file<<", line "<<line<<":"<<std::endl; \
  std::cerr << msg << std::endl; \
  return;}

#define FERROR(msg,file,line)  {std::cerr << "Error in file "<<file<<", line "<<line<<":"<<std::endl; \
  std::cerr << msg << std::endl; \
  return Teuchos::null;}

#define WARNING(msg,file,line)  {std::cerr << "Warning in file "<<file<<", line "<<line<<":"<<std::endl; \
  std::cerr << msg << std::endl;}

#define PRINTLEVEL 0

#define OUT(s) if (am_printer) std::cout << s << std::endl;

#ifdef HAVE_MPI
#define HAD_MPI
#undef HAVE_MPI
#endif

#ifdef HAVE_HYPRE
#define HAD_HYPRE
#undef HAVE_HYPRE
#endif

#include "Tpetra_Map.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Operator.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Utils.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultSerialComm.hpp"

#include "BelosTpetraAdapter.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "Ifpack2_Factory.hpp"
#include "Ifpack2_AdditiveSchwarz.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Kokkos_SerialNode.hpp"

#ifdef HAVE_CUDA
#include "Kokkos_ThrustGPUNode.hpp"
#endif

//#ifdef DEBUG_TRILINOS_INTERFACE
// only for debugging
//#include "TpetraExt_RowMatrixOut.h"
//#include "TpetraExt_MultiVectorOut.h"
//#include "TpetraExt_BlockMapOut.h"
//#endif

#ifdef HAD_MPI
#define HAVE_MPI
#endif
#ifdef HAD_HYPRE
#define HAVE_HYPRE
#endif


typedef float ST;
typedef float CMST;
typedef int OT;
typedef OT GOT;
typedef OT LOT;

typedef Kokkos::SerialNode serial_node_type;
#ifdef HAVE_CUDA
typedef Kokkos::ThrustGPUNode gpu_node_type;
#endif

#if defined(HAVE_CUDA) && defined(USE_GPU_SOLVER)
typedef gpu_node_type node_type;
#else
typedef serial_node_type node_type;
#endif

typedef Tpetra::Map<LOT,GOT,node_type> map_type;
typedef Tpetra::Export<OT,OT,node_type> export_type;
typedef Tpetra::CrsMatrix<CMST,LOT,GOT,node_type> crsmatrix_type;
typedef Teuchos::Comm<OT> comm_type;
typedef Tpetra::Vector<ST,LOT,GOT,node_type> vector_type;
typedef Tpetra::MultiVector<ST,LOT,GOT,node_type> multivector_type;
typedef Tpetra::Operator<ST,LOT,GOT,node_type> operator_type;

typedef Teuchos::ScalarTraits<ST> SCT;
typedef SCT::magnitudeType MT;
typedef Belos::MultiVecTraits<ST,multivector_type> MVT;
typedef Belos::OperatorTraits<ST,multivector_type,operator_type> OPT;

typedef Ifpack2::ILUT< crsmatrix_type > LocalInverseType;

#ifdef DEBUG_TRILINOS_INTERFACE
class MemTest 
  {
  public:
  
  MemTest() {std::cerr <<  "Elmer TrilinosKokkos container constructed"<<std::endl;}
  ~MemTest() {std::cerr << "Elmer TrilinosKokkos container destroyed"<<std::endl;}
  };
#endif

typedef struct ElmerTrilinosKokkosContainer {

Teuchos::RCP< comm_type > comm_;
Teuchos::RCP< const map_type > assemblyMap_; // map with 'overlap' of nodes
Teuchos::RCP< const map_type > solveMap_; // map with each node on one proc
Teuchos::RCP< export_type > exporter_;
Teuchos::RCP< crsmatrix_type > matrix_;
Teuchos::RCP< vector_type > rhs_;
Teuchos::RCP< vector_type > sol_;
Teuchos::RCP< multivector_type > coords_; // node coordinates can be used by ML
Teuchos::RCP< Teuchos::Array<OT> > GDOFs_;

Teuchos::RCP< node_type > node_;

double scaleFactor_; // scale entire system by a scalar constant, scaleFactor*Ax = scaleFactor*b
                     // enabled by setting "Scale Factor" in xml file

Teuchos::RCP<Teuchos::ParameterList> params_;

Teuchos::RCP< operator_type > prec_;
Teuchos::RCP<Belos::SolverManager<ST,multivector_type,operator_type> > solver_;

Teuchos::RCP<struct ElmerTrilinosKokkosContainer> previous_;
Teuchos::RCP<struct ElmerTrilinosKokkosContainer> next_;

#ifdef DEBUG_TRILINOS_INTERFACE
Teuchos::RCP<MemTest> memtest_;
#endif
} ElmerTrilinosKokkosContainer;

// to avoid warnings from RCP's in debug mode when 
// temporarily returning the control to Fortran,   
// we store an extern pointer to a (doubly) linked 
// list of all 'Container' objects here
static Teuchos::RCP<ElmerTrilinosKokkosContainer> containerListHead=Teuchos::null;



// some auxiliary functions implemented below:

template<class S1,class S2>
void ElmerTrilinosKokkos_convertS1ToS2(S1* s1arr, S2* s2arr, size_t length);

// creates the map without overlap
Teuchos::RCP< const map_type > ElmerTrilinosKokkos_createSolveMap
        (Teuchos::RCP< comm_type > comm,
        int n, const Teuchos::ArrayView< const OT >& GID, int* owner);

// creates the matrix with overlap (e.g. shared nodes)
Teuchos::RCP< crsmatrix_type > ElmerTrilinosKokkos_createMatrix
        (Teuchos::RCP< const map_type > assemblyMap,
         Teuchos::RCP< Teuchos::ParameterList > params,
        int *rows, int *cols, double* vals);

Teuchos::RCP< operator_type > ElmerTrilinosKokkos_createPreconditioner(
        Teuchos::RCP< crsmatrix_type > A,
        Teuchos::ParameterList& params,
        Teuchos::RCP< multivector_type > coords);

Teuchos::RCP< operator_type > ElmerTrilinosKokkos_createIfpackPreconditioner(
        Teuchos::RCP< const crsmatrix_type > A, Teuchos::ParameterList& params);

Teuchos::RCP<Belos::SolverManager<ST,multivector_type,operator_type> > ElmerTrilinosKokkos_createSolver
        (Teuchos::RCP<operator_type> A, Teuchos::RCP<operator_type> P,
        Teuchos::RCP<multivector_type> x, Teuchos::RCP<multivector_type> b,
        Teuchos::ParameterList& params);

static bool am_printer = true; // for output


// need to declare this extern so that the function names are
// consistent with Fortran
extern "C" {


// the calling order of SolveTrilinosKokkos1..4 is the same as for SolveHYPRE1..4

// construct matrix, solver and preconditioner
// nrows - number of local rows
// ncols - number of column indices on local partition
void SolveTrilinosKokkos1
 (
  int *n, int *nnz,
  int *rows, int *cols, double *vals,
  int *globaldof, int *owner, char* xmlfile, int* verbosityPtr, int** ContainerPtr,
  int *num_nodes, 
  double* xcoords, double* ycoords, double* zcoords,
  int *returnCode)
{

   bool serial=false;
   
   int verbose=*verbosityPtr;
   
   int& ierr=*returnCode;
   ierr=0;
   bool success=true;

   OUT("starting TrilinosKokkos setup");
   
   // Elmer has the unpleasant habit of finalizing MPI
   // if only one process is involved, this has to be 
   // changed in the source code if we want to use 
   // TrilinosKokkos for sequential runs as well. Check to make sure:
   int mpi_state;
   MPI_Initialized(&mpi_state);
   if (mpi_state==0)
     {
     OUT("MPI_Init has not been called, using serial comm");
     serial=true;
     }
   else
     {
     MPI_Finalized(&mpi_state);
     if (mpi_state!=0)
       {
       OUT("MPI_Finalize has already been called, using serial comm");
       serial=true;
       }
     }

   Teuchos::RCP< comm_type > comm;
   
   //Teuchos::RCP< Kokkos::ThrustGPUNode > node = Teuchos::rcp(new Kokkos::ThrustGPUNode(params)) //RMV
   
   // create communicator
   if (serial)
     {
     comm = Teuchos::rcp(new Teuchos::SerialComm<OT>());
     }
   else
     {
     comm = Teuchos::rcp(new Teuchos::MpiComm<OT>(MPI_COMM_WORLD));
     }
   
   am_printer = (comm->getRank()==0)&&(verbose>PRINTLEVEL);

   Teuchos::RCP< Teuchos::Array<OT> > GDOFs;
   GDOFs = Teuchos::rcp(new Teuchos::Array<OT>(*n));
   for ( Teuchos::Array<OT>::size_type k = 0 ; k < *n ; k++ )
       (*GDOFs)[k] = globaldof[k];
   
   // sanity check, we expect standard CRS format
   if (*nnz != rows[*n])
     {
     ERROR("number of nonzeros incorrect",__FILE__,__LINE__);
     }

   Teuchos::RCP< const map_type > assemblyMap;
   Teuchos::RCP< crsmatrix_type > A;

   // first create a map which has nodes shared between partitions
   // (e.g. an 'Elmer'-map)
   try {
     assemblyMap = Teuchos::rcp
        (new map_type (-1, *GDOFs, 1, comm));
   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)

   if (!success)
     {
     WARNING("Failed to construct Elmer map",__FILE__,__LINE__);
     ierr = -1;
     return;
     }


#ifdef DUMP_IN_TRILINOS_INTERFACE 
CHECK_ZERO(TpetraExt::BlockMapToMatrixMarketFile("ElmerMap.txt",*assemblyMap));
#endif

//////////////////////////////////////////////////////////////
// read parameters from input file 
//////////////////////////////////////////////////////////////

   // read solver parameters from a file
   Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList());
   Teuchos::RCP< Teuchos::ParameterList > fillparams = Teuchos::sublist(params, "Local Sparse Ops");
   fillparams->set("Prepare Solve", true);

   std::string filename(xmlfile);
   if (filename=="none")
     {
     WARNING("no parameter file specified, using default settings in TrilinosKokkos",__FILE__,__LINE__);
     }
   else
     {
     OUT("reading parameters from '"+filename+"'");
     try {
#ifdef OLD_TRILINOS
       Teuchos::updateParametersFromXmlFile(filename,params.get());
#else
// for TrilinosKokkos 10.10 and later     
       Teuchos::updateParametersFromXmlFile(filename,params.ptr());
#endif       
       } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
       if (!success)
         {
         WARNING("failed to read your TrilinosKokkos input file, using default values",
                __FILE__,__LINE__);
         ierr = 1;
         }
     }
     
   // based on this map, create a sparse matrix. 
   Teuchos::RCP< crsmatrix_type > A_elmer;
   try {
   A_elmer = 
        ElmerTrilinosKokkos_createMatrix(assemblyMap, params, rows, cols, vals);
   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)

   if (!success || A_elmer==Teuchos::null)
     {
     WARNING("Failed to construct Elmer matrix",__FILE__,__LINE__);
     ierr = -1;
     return;
     }

#ifdef DUMP_TRILINOS_INTERFACE
   CHECK_ZERO(TpetraExt::RowMatrixToMatrixMarketFile("ElmerMatrix.txt",*A_elmer));
#endif
     
   // now construct a map with each node owned by one partition (a solver-map)
   Teuchos::RCP< const map_type > solveMap;
   try {   
    solveMap = ElmerTrilinosKokkos_createSolveMap(comm, *n, *GDOFs, owner);
   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)

   if (!success || solveMap==Teuchos::null)
     {
     WARNING("Failed to construct map",__FILE__,__LINE__);
     ierr = -1;
     return;
     }

#ifdef DUMP_IN_TRILINOS_INTERFACE 
CHECK_ZERO(TpetraExt::BlockMapToMatrixMarketFile("TrilinosKokkosMap.txt",*solveMap));
#endif

   Teuchos::RCP< export_type > exporter;

   // construct an exporter to transfer data between the two object types
   try {
   exporter = Teuchos::rcp(new export_type
        (assemblyMap, solveMap));
   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)

   if (!success || exporter==Teuchos::null)
     {
     WARNING("Failed to construct exporter",__FILE__,__LINE__);
     ierr = -1;
     return;
     }


try {   
   // build the non-overlapping matrix
   A = Tpetra::createCrsMatrix<CMST,LOT,GOT,node_type>(solveMap, A_elmer->getNodeMaxNumRowEntries());
   // create the matrix from the overlapping input matrix
   A->doExport(*A_elmer, *exporter, Tpetra::ADD); //TODO:PTW:CHECK_ZERO
   Teuchos::RCP< Teuchos::ParameterList > fillparams = Teuchos::sublist(params, "Local Sparse Ops");
   fillparams->set("Prepare Solve", true);
   A->fillComplete(params); //TODO:PTW:CHECK_ZERO

   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)
   
   if (!success)
     {
     WARNING("Failed to construct matrix",__FILE__,__LINE__);
     ierr = -1;
     return;
     }
   
   // Create the rhs and solution
   Teuchos::RCP< vector_type > sol=Teuchos::rcp(new vector_type(A->getRowMap()));
   Teuchos::RCP< vector_type > rhs=Teuchos::rcp(new vector_type(A->getRowMap()));

   OUT("matrix constructed")


   bool print_matrix = params->get("Dump Matrix",false);

   if (print_matrix)
   {
#if 1
//#ifdef DEBUG_TRILINOS_INTERFACE
//   std::string filename = params->get("Filename Base","TrilinosKokkos")+"Matrix.mtx";
//   CHECK_ZERO(TpetraExt::RowMatrixToMatrixMarketFile(filename.c_str(),*A));
   WARNING("you have specified 'Dump Matrix', but DEBUG_TRILINOS_INTERFACE is not working for Tpetra",
        __FILE__,__LINE__);
#else
   WARNING("you have specified 'Dump Matrix', but DEBUG_TRILINOS_INTERFACE is not defined",
        __FILE__,__LINE__);
#endif
   }

   double scaleFactor = params->get("Scale Factor",1.0);
   if (scaleFactor!=1.0) A->scale(scaleFactor); //TODO:PTW:CHECK_ZERO

   Teuchos::RCP< multivector_type > coords0=Teuchos::rcp(new multivector_type(assemblyMap,3));
   Teuchos::RCP< multivector_type > coords=Teuchos::rcp(new multivector_type(solveMap,3));
   std::cout << (OT)coords->getLocalLength() << "@@" << (OT)coords->getGlobalLength() << std::endl;
   std::cout << (OT)coords0->getLocalLength() << "@@" << (OT)coords0->getGlobalLength() << std::endl;
   
   int k = *num_nodes;
   int dof = (int)(assemblyMap->getNodeNumElements()/k);
   if (dof*k != assemblyMap->getNodeNumElements())
     {
     ERROR("size mismatch of coord arrays",__FILE__,__LINE__);
     }

   Teuchos::ArrayRCP< Teuchos::ArrayRCP<ST> > coordsview = coords->get2dViewNonConst();

   //TODO:PTW:What about parallel considerations? RMV
   Teuchos::ArrayRCP< Teuchos::ArrayRCP<ST> > coords0view = coords0->get2dViewNonConst();
   Teuchos::ArrayRCP<ST> xcoords0view = coords0view[0];
   Teuchos::ArrayRCP<ST> ycoords0view = coords0view[1];
   Teuchos::ArrayRCP<ST> zcoords0view = coords0view[2];

   for (int i=0;i<k; i++)
     {
     for (int j=0;j<dof;j++)
       {
       xcoords0view[dof*i+j] = xcoords[i];
       ycoords0view[dof*i+j] = ycoords[i];
       zcoords0view[dof*i+j] = zcoords[i];
       }
     }

   coords->doExport(*coords0, *exporter, Tpetra::INSERT); //TODO:PTW:CHECK_ZERO && Tpetra::ZERO //RMV
   
   //TODO:PTW: check if updating k is the correct response?
   k = (int)(solveMap->getNodeNumElements() / dof);

   Teuchos::RCP< const map_type > auxMap = Tpetra::createContigMapWithNode<LOT,GOT,node_type>(-1, k, comm, Teuchos::rcp<node_type>(new node_type()));
   Teuchos::RCP< multivector_type > auxVec = Teuchos::rcp(new multivector_type(auxMap,3));
   //When leaving this scope, the view should update the vector (auxVec)
   {
     Teuchos::ArrayRCP< Teuchos::ArrayRCP<ST> > auxVecview = auxVec->get2dViewNonConst();
     for (int i=0;i<k;i++)
       { 
       for (int j=0;j<3;j++)
         {
         ST l = coordsview[j][dof*i];
         auxVecview[j][i] =l;
         }
       }
   }
   
   coords = auxVec;
   
  ///////////////////////////////////////////////////////////////////
  // create/setup preconditioner                                   //
  ///////////////////////////////////////////////////////////////////
  Teuchos::RCP< operator_type > prec = ElmerTrilinosKokkos_createPreconditioner(A, *params, coords);

  //////////////////////////////////////////////////////////////
  // Krylov subspace method setup (Belos)                     //
  //////////////////////////////////////////////////////////////
  Teuchos::RCP<Belos::SolverManager<ST,multivector_type,operator_type> > solver = ElmerTrilinosKokkos_createSolver(A,prec,
        sol, rhs, *params);

  // print parameter list so we can see default values
  // and if some of our input is unused.
  if (am_printer && verbose>=5)
    {
    std::cout << "TrilinosKokkos parameter list: "<<std::endl;
    std::cout << *params<<std::endl;
    }


   //////////////////////////////////////////////////////////////
   // construct a container to return to Fortran               //
   //////////////////////////////////////////////////////////////
   if (*ContainerPtr != 0)
     {
     WARNING("pointer passed into SolveTrilinosKokkos1 not NULL, possible memory leak.",
        __FILE__,__LINE__);
     }

   Teuchos::RCP<ElmerTrilinosKokkosContainer> Container = 
        Teuchos::rcp(new ElmerTrilinosKokkosContainer());

   *ContainerPtr=(int*)(Container.get());
   
   // store a pointer to the container in a static 
   // linked list to avoid confusing the garbage   
   // collection (Teuchos RCPs won't understand that
   // we have stored a raw pointer in Elmer somewhere)
   if (containerListHead!=Teuchos::null)
      {
      containerListHead->previous_=Container;
      }
    Container->next_=containerListHead;
    containerListHead = Container;
   
   // put pointers to all the important objects 
   // in the container.
   Container->comm_=comm;
   Container->params_=params;
   Container->assemblyMap_=assemblyMap;
   Container->solveMap_=solveMap;
   Container->exporter_=exporter;
   Container->matrix_=A;
   Container->coords_=coords;
   Container->scaleFactor_=scaleFactor;
   Container->rhs_=rhs;
   Container->sol_=sol;
   Container->GDOFs_ =GDOFs;
  Container->prec_=prec;
  Container->solver_=solver;
  
#ifdef DEBUG_TRILINOS_INTERFACE
  Container->memtest_=Teuchos::rcp(new MemTest());
#endif
  
  return;

  }
////////////////////////////////////////////////////////////

void SolveTrilinosKokkos2
 (
  int *n, double *xvec, double *rhsvec, int *Rounds, ST *TOL,
  int *verbosityPtr, int** ContainerPtr,
  int* returnCode
 )
  {
  int verbose = *verbosityPtr;
   int& ierr=*returnCode;   
   ierr=0;
   bool success=true;

ElmerTrilinosKokkosContainer* Container = (ElmerTrilinosKokkosContainer*)(*ContainerPtr);
if (Container==NULL) ERROR("invalid pointer passed to SolveTrilinosKokkos2",__FILE__,__LINE__);

   am_printer = (Container->comm_->getRank()==0)&&(verbose>PRINTLEVEL);
   
   // get the data structures built in SolveTrilinosKokkos1():
   Teuchos::RCP< const map_type > assemblyMap = Container->assemblyMap_;
   Teuchos::RCP< const map_type > solveMap = Container->solveMap_;
   Teuchos::RCP< export_type > exporter = Container->exporter_;
   Teuchos::RCP< crsmatrix_type > A = Container->matrix_;
   Teuchos::RCP< vector_type > x = Container->sol_;
   Teuchos::RCP< vector_type > b = Container->rhs_;
   Teuchos::RCP< operator_type > prec = Container->prec_;
   Teuchos::RCP<Belos::SolverManager<ST,multivector_type,operator_type> > solver = Container->solver_;
   Teuchos::RCP<Teuchos::ParameterList> params = Container->params_;
   size_t num_elts = assemblyMap->getGlobalNumElements();
   ST scaleFactor = Container->scaleFactor_;
   float rhsvecf[num_elts];
   float xvecf[num_elts];

   // import the vectors
   //TODO: For view approach: Teuchos::ArrayRCP<const ST> bview(rhsvec, 0, assemblyMap->getGlobalNumElements(), false);
   ElmerTrilinosKokkos_convertS1ToS2<double,float>(rhsvec, rhsvecf, num_elts);
   Teuchos::ArrayView<ST> bview(rhsvecf, num_elts);
   //TODO:Maybe find some way of restoring this for GPU, or at least making it optional
   //Teuchos::RCP< vector_type > bvec = Tpetra::createVectorFromView<ST,LOT,GOT,node_type>(assemblyMap, bview);
   Teuchos::RCP< vector_type > bvec = Teuchos::rcp(new vector_type(assemblyMap, bview));
   b->doExport(*bvec, *exporter, Tpetra::ADD); //TODO:PTW:CHECK_ZERO

   // import the vectors
   //Teuchos::ArrayRCP<ST> xview(xvec, 0, assemblyMap->getGlobalNumElements(), false);
   ElmerTrilinosKokkos_convertS1ToS2<double,float>(xvec, xvecf, num_elts);
   Teuchos::ArrayView<ST> xview(xvecf, num_elts);
   //TODO:Maybe find some way of restoring this for GPU, or at least making it optional
   //Teuchos::RCP< vector_type > xv = Tpetra::createVectorFromView<ST,LOT,GOT,node_type>(assemblyMap, xview);
   Teuchos::RCP< vector_type > xv = Tpetra::rcp(new vector_type(assemblyMap, xview));
   x->doExport(*xv, *exporter, Tpetra::INSERT); //TODO:PTW:CHECK_ZERO

   if (scaleFactor!=1.0)
     { 
     b->scale(scaleFactor); //TODO:PTW:CHECK_ZERO
     }

   std::cout << "A norm : " << A->getFrobeniusNorm() << std::endl;
   std::cout << "b norm : " << b->normInf() << std::endl;
   std::cout << "bvec norm : " << bvec->normInf() << std::endl;
   std::cout << "xv norm : " << xv->normInf() << std::endl;

   // override the settings for tconv tol and num iter using Elmer inut data:
   if (*TOL>=0.0) params->sublist("Belos").set("Convergence Tolerance",*TOL);
   if (*Rounds>0) params->sublist("Belos").set("Maximum Iterations",*Rounds);
   
  // check initial residual - do not start solver if already converged.
  
  int numrhs=1;
  std::vector<ST> actual_resids( numrhs );
  std::vector<ST> rhs_norm( numrhs );
  multivector_type resid(A->getRangeMap(), numrhs);

#if 0
  OPT::Apply( *A, *x, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *b, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *b, rhs_norm );  
  bool converged = true;
  
  for ( int i=0; i<numrhs; i++)
    {
    double actRes = actual_resids[i]/rhs_norm[i];
    if (actRes >= *TOL) converged = false;
    } 
  
  if (converged)
    {
    OUT("initial solution passed to TrilinosKokkos is already good enough - returning to Elmer.");
    return;
    }
#else
// always start with a 0 initiaial guess because Belos
// does a relative residual check only
//x->PutScalar(0.0);
#endif   

/////////////////////////////////////////////////////////////
// Perform solve                                           //
/////////////////////////////////////////////////////////////

  if (solver!=Teuchos::null)
    {
    OUT("start iterative solver");
    solver->setParameters(Teuchos::rcp(&(params->sublist("Belos")),false));
    solver->reset(Belos::Problem);
    Belos::ReturnType ret;
    try {
    ret = solver->solve();
    } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

    if (!success) 
      {
      WARNING("TrilinosKokkos solve failed!",__FILE__,__LINE__);
      ierr=-1; // caught an exception -> error
      return;
      }

    // check for loss of accuracy
    bool loa = solver->isLOADetected();

    if (loa)
      {
      WARNING("loss of accuracy in Belos solve!",__FILE__,__LINE__);
      ierr=1;
      }

    if (ret!=Belos::Converged)
      {
      WARNING("Belos did not converge!",__FILE__,__LINE__);
      ierr=2;
      }
    ierr=0;//RMV
    }
  else if (prec!=Teuchos::null)
    {
    OUT("can't apply operator inverse");
    //prec->applyInverse(*b,*x); //TODO:PTW:CHECK_ZERO && REINSTATE!!! //RMV
    }
  else
    {
   WARNING("no solver or preconditioner available",__FILE__,__LINE__);
    ierr=3;
    *x=*b;
    }

   bool print_vectors = params->get("Dump Vectors",false);
   if (print_vectors)
   {
#ifdef DEBUG_TRILINOS_INTERFACE   
   string filebase = params->get("Filename Base","TrilinosKokkos");
   TpetraExt::MultiVectorToMatrixMarketFile((filebase+"Rhs.txt").c_str(),*b);
   TpetraExt::MultiVectorToMatrixMarketFile((filebase+"Sol.txt").c_str(),*x);
#else
   WARNING("you have specified 'Dump Vectors', but DEBUG_TRILINOS_INTERFACE is not defined",
        __FILE__,__LINE__);
#endif
   }

   

#ifdef DEBUG_TRILINOS_INTERFACE   
    
  //
  // Compute actual residual
  //
  bool badRes = false;

  OPT::Apply( *A, *x, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *b, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *b, rhs_norm );  

  if (verbose>=3)
    {
    std::cout << "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    }
  for ( int i=0; i<numrhs; i++) 
    {
    ST actRes = actual_resids[i]/rhs_norm[i];
    if (verbose>=3)
      {
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      }
    if (actRes > *TOL) badRes = true;
    } 

  if (badRes) 
    {
    WARNING("bad actual residual found!",__FILE__,__LINE__);
    std::cerr << "(required accuracy was "<<*TOL<<")"<<std::endl;
    ierr=4;
    }

#endif

   // import the vectors
   xv->doImport(*x, *exporter, Tpetra::INSERT); //TODO:PTW:CHECK_ZERO & Work out Zero
   std::cout << "xv norm : " << xv->normInf() << std::endl;
   xv->get1dCopy(xview, assemblyMap->getGlobalNumElements());
   ElmerTrilinosKokkos_convertS1ToS2<float,double>(xvecf, xvec, num_elts);

// TrilinosKokkos cleans up itself (because of RCP's)
return;
}

////////////////////////////////////////////////////////////
// destructor                                             //
////////////////////////////////////////////////////////////

void SolveTrilinosKokkos4(int** ContainerPtr)
  {
  ElmerTrilinosKokkosContainer* Container = 
        (ElmerTrilinosKokkosContainer*)(*ContainerPtr);

#ifdef DEBUG_TRILINOS_INTERFACE  
std::cerr << "PID "<<Container->comm_->getRank()<<": destroy TrilinosKokkos object "<<std::endl;
#endif

  // remove this container from the list
  if (Container->next_!=Teuchos::null)
    {
    Container->next_->previous_=Container->previous_;
    }
  if (Container->previous_!=Teuchos::null)
    {
    Container->previous_->next_=Container->next_;
    }

  // nullify the pointer in fortran
  if (Container == containerListHead.get())
    {
    containerListHead=Container->next_;
    }

  *ContainerPtr=0;
  Container = NULL;
  
  
  
  // do some global barriers after destruction
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  }                                                                                                                                         

// this function deletes ALL TrilinosKokkos objects created by Elmer that     
// have not been destroyed properly, yet (by SolveTrilinosKokkos4).           
// It should only be called at the very end of an Elmer run.            
void TrilinosKokkosCleanup(void)
  {
#ifdef DEBUG_TRILINOS_INTERFACE  
  std::cout << "Destroying all remaining TrilinosKokkos objects..."<<std::endl;
#endif
  Teuchos::RCP<struct ElmerTrilinosKokkosContainer> c = containerListHead;
  while (c!=Teuchos::null)
    {
    c=c->next_;
    if (c!=Teuchos::null) c->previous_ = Teuchos::null;
    }
  containerListHead = Teuchos::null;
  }

}//extern "C"

// creates the map without overlap
Teuchos::RCP< const map_type > ElmerTrilinosKokkos_createSolveMap
        (Teuchos::RCP< Tpetra::Comm<OT> > comm,
        int n, const Teuchos::ArrayView< const OT >& GID, int* owner)
  {
  int ierr;
  // first figure out how many nodes are 'owned' by this process
  int nloc = 0;
  for (int i=0;i<n;i++) nloc+=owner[i];
  
  // construct an array with only owned global node id's
  Teuchos::RCP< Teuchos::Array<OT> > MyGIDs;
  MyGIDs = Teuchos::rcp(new Teuchos::Array<OT>(nloc));
  int k=0;
  for (int i=0;i<n;i++)
    {
    if (owner[i])
      {
      (*MyGIDs)[k++] = GID[i];
      }
    }
  
  Teuchos::RCP< const map_type > map = Teuchos::rcp(new map_type(-1, *MyGIDs, 1, comm, Teuchos::rcp<node_type>(new node_type())));
  
  //delete [] MyGIDs; //This should be removed now by the reference counter...
  return map;
  }

// creates the matrix with overlap (e.g. shared nodes)
Teuchos::RCP< crsmatrix_type > ElmerTrilinosKokkos_createMatrix
        (Teuchos::RCP< const map_type > assemblyMap,
         Teuchos::RCP< Teuchos::ParameterList > params,
        int *rows, int *cols, double* vals)
  {
  int ierr=0;
  bool success=true;
  Teuchos::RCP< crsmatrix_type > A;
  
  if (assemblyMap==Teuchos::null)
    {
    FERROR("map passed to ElmerTrilinosKokkos_createMatrix is null",__FILE__,__LINE__);
    }
  
  int nrows = assemblyMap->getNodeNumElements();
  
  int *row_size = new int[nrows];
  
  int max_row_size=0;
  
  for (int i=0;i<nrows;i++)
    {
    row_size[i]=rows[i+1]-rows[i];
    if (row_size[i]>max_row_size) max_row_size=row_size[i];
    }
  
  int *gcols=new int[max_row_size];

  A = Teuchos::rcp(new crsmatrix_type(assemblyMap, max_row_size, Tpetra::StaticProfile));
  
  // Now go through my local rows and set the matrix entries.
  try {
  for (int i = 0; i < nrows; i++)
    {
    for (int j=0;j<row_size[i];j++) gcols[j]=assemblyMap->getGlobalElement(cols[rows[i]-1+j]-1);
    Tpetra::ArrayView<OT> gcolsview(gcols, row_size[i]);
    float valsrowf[row_size[i]];
    ElmerTrilinosKokkos_convertS1ToS2<double,float>(&(vals[rows[i]-1]), valsrowf, row_size[i]);
    Tpetra::ArrayView<const ST> valsview(valsrowf, row_size[i]);
    A->insertGlobalValues(assemblyMap->getGlobalElement(i),gcolsview,valsview); //TODO:PTW:CHECK_ZERO
    }
  } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)
  delete [] row_size;
  delete [] gcols;
  // Assemble after setting the coefficients
  A->fillComplete(params); //TODO:PTW:CHECK_ZERO
  return A;
  }

//TODO:Ensure this doesn't run if S1==S2
template<class S1,class S2>
void ElmerTrilinosKokkos_convertS1ToS2(S1* s1arr, S2* s2arr, size_t length)
{
    std::copy(s1arr, s1arr+length, s2arr); //TODO: improve end arg
}

Teuchos::RCP< operator_type > ElmerTrilinosKokkos_createPreconditioner(
        Teuchos::RCP< crsmatrix_type > A, 
        Teuchos::ParameterList& params,
        Teuchos::RCP< multivector_type > coords)
        {
        std::string type="None";
        type=params.get("Preconditioner",type);
        if (type=="Ifpack") return ElmerTrilinosKokkos_createIfpackPreconditioner(A,params);
        if (type=="None") return Teuchos::null;
        if (type=="ML") WARNING("ML not available for Tpetra (waiting on release of MueLu), returning null",__FILE__,__LINE__);
        WARNING("invalid 'Preconditioner', returning null",__FILE__,__LINE__);
        return Teuchos::null;
        }


Teuchos::RCP< operator_type > ElmerTrilinosKokkos_createIfpackPreconditioner(
        Teuchos::RCP< const crsmatrix_type > A,
        Teuchos::ParameterList& params)
  {
  int ierr;
  Teuchos::RCP< operator_type > prec = Teuchos::null;
  string type = "None";
  type=params.get("Preconditioner",type);
  
  if (type=="Ifpack")
    {
    Teuchos::RCP< Ifpack2::AdditiveSchwarz< crsmatrix_type, LocalInverseType > > ifPrec;
    Teuchos::ParameterList& ifpackList = params.sublist("Ifpack");
    string ifpackType = "RELAXATION";
    ifpackType=params.get("Ifpack Preconditioner",ifpackType);
    OUT(ifpackType);
    int OverlapLevel = params.get("Ifpack Overlap Level",0); 
    
    OUT("construct ifpack preconditioner object");
    //Ifpack2::Factory factory;
    
    //ifPrec = factory.create(ifpackType, A);
    ifPrec = Teuchos::rcp(new Ifpack2::AdditiveSchwarz< crsmatrix_type, LocalInverseType >(A, OverlapLevel));
    ifpackList.set("fact: drop-tolerance", 1e-9);
    ifpackList.set("fact: level-of-fill", 1);
    //ifPrec = Teuchos::rcp(factory.Create(ifpackType,
    //           A.get(), OverlapLevel) );

    OUT("set parameters");
    ifPrec->setParameters(ifpackList);
    OUT("initialize");
    ifPrec->initialize(); //TODO:PTW:CHECK_ZERO
 
    OUT("compute");
    ifPrec->compute(); //TODO:PTW:CHECK_ZERO
    
    prec = ifPrec;
    }
  else
    {
    FERROR(" 'Preconditioner' must be 'Ifpack' for this function.",
        __FILE__,__LINE__);
    }
  return prec;
  }

 
// create an iterative solver for the linear system Ax=b,
// preconditioned by P.
Teuchos::RCP<Belos::SolverManager<ST,multivector_type,operator_type> > ElmerTrilinosKokkos_createSolver
        (Teuchos::RCP<operator_type> A, Teuchos::RCP<operator_type> P,
         Teuchos::RCP<multivector_type> x, Teuchos::RCP<multivector_type> b,
         Teuchos::ParameterList& params)
  {
  Teuchos::RCP<Belos::SolverManager<ST,multivector_type,operator_type> > belosSolverPtr;


  // retrieve User's Belos list, add more things
  Teuchos::ParameterList& belosList=params.sublist("Belos");

  //! Belos linear problem interface
  Teuchos::RCP<Belos::LinearProblem<ST,multivector_type,operator_type> > belosProblemPtr;

  string linearSolver=params.get("Iterative Solver","GMRES");
  if (linearSolver=="None")
    {
    return Teuchos::null;
    }
  bool verbose = true;
#ifdef DEBUG_TRILINOS_INTERFACE
  bool debug = true;
#else
  bool debug = false;
#endif  

  int verbosity = Belos::Errors + Belos::Warnings;
  if (verbose)
    { //TODO: where to put which option? how do we get readable output?
    verbosity+=Belos::TimingDetails+Belos::IterationDetails;
    verbosity+=Belos::StatusTestDetails+Belos::OrthoDetails+Belos::FinalSummary;
    }
  if (debug) verbosity+=Belos::Debug;
  // User is allowed to override these settings
  if (belosList.isParameter("Verbosity")==false)
    {
    belosList.set("Verbosity",verbosity);
    }

  if (belosList.isParameter("Output Stream")==false)
    {
    belosList.set("Output Stream",Teuchos::rcp(&std::cout, false));
    }

  belosList.set("Output Style",(int)Belos::Brief);
  
  // by default, Belos checks only ||r_imp||_2/||r_0||_2, but if we're
  // solving a sequence of systems and therefore have a very good initial
  // guess, it is better to check ||r||_2/||b||_2 (actually a more fancy 
  // convergence test like the one used in Elmer would be better, but I  
  // haven't figured out how to tell Belos to do that)
  if (linearSolver=="GMRES" && belosList.isParameter("Implicit Residual Scaling")==false)
    {
    belosList.set("Implicit Residual Scaling","Norm of RHS");
    }

  // create Belos interface to preconditioner.
  // This is simply an operator_type with 'Apply' and 'ApplyInverse' switched.
  //if (!Teuchos::is_null(P))
  //  {
  //  belosPrecPtr = Teuchos::rcp(new Belos::TpetraPrecOp(P));
  //  }
  // create Belos problem interface
  belosProblemPtr = Teuchos::rcp(new Belos::LinearProblem<ST,multivector_type,operator_type>(A,x,b));

  if (P!=Teuchos::null)
    {
    // set preconditioner
    OUT("set preconditioner...");
    belosProblemPtr->setLeftPrec(P);
    }
  bool set = belosProblemPtr->setProblem();
  if (set == false) {
    FERROR("Belos::LinearProblem failed to set up correctly!",__FILE__,__LINE__);
    }

  // create the solver
  Teuchos::RCP<Teuchos::ParameterList> belosListPtr=rcp(&belosList,false);
  if (linearSolver=="CG")
    {
    belosSolverPtr = Teuchos::rcp(new Belos::BlockCGSolMgr<ST,multivector_type,operator_type>(belosProblemPtr,belosListPtr));
    }
  else if (linearSolver=="GMRES")
    {
    Teuchos::RCP<Teuchos::ParameterList> belosListPtr=Teuchos::rcp(&belosList,false);
    belosSolverPtr = Teuchos::rcp(new Belos::BlockGmresSolMgr<ST,multivector_type,operator_type>(belosProblemPtr,belosListPtr));
    }
  else if (linearSolver=="None")
    {
    belosSolverPtr = Teuchos::null;
    }
  else
    {
    FERROR("Currently 'CG', 'GMRES' and 'None' \n"
        " are supported as 'Iterative Solver'",__FILE__,__LINE__);
    }
  return belosSolverPtr;
  }

void TrilinosKokkosSetNodeCoords(int* num_nodes, double* x, double* y, double* z)
  {
  }

#endif
