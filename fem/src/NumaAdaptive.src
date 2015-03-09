! *****************************************************************************/
! *  Adaptive Meshing routines
! *****************************************************************************/
! *  This file is based on the file "Adaptive.src" of Elmer (Finite Element 
! *  Software for Multiphysical Problems)
! *  Copyright 1st April 1995 - , CSC - Scientific Computing Ltd., Finland
! *****************************************************************************/
! *****************************************************************************/
! * Version 1b (02/09/09)
! * Main changes with Elmer file Adaptive.src  (02/09/09)
! *****************************************************************************/
! *  - Computation of indicator
! *  - Refinement only by splitting
! *****************************************************************************/
! *****************************************************************************/
! * List of subroutines:
! *****************************************************************************/
! *  - NumaRefineMesh
! *  - SplitMesh
! *  - SplitOneLevel
! *  - RGBRefinement
! *  - SetParents
! *  - ComputeError
! *****************************************************************************/
! *****************************************************************************/
! *  --/--/--:
! *  - ...
! *****************************************************************************/ 

MODULE NumaAdaptive

	USE GeneralUtils
	USE SolverUtils
	USE NumaSolverUtils
	USE ModelDescription

	IMPLICIT NONE

CONTAINS

!******************************************************************************
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
	SUBROUTINE NumaRefineMesh( Model,Solver,Quant,Perm, &
		InsideResidual, EdgeResidual, BoundaryResidual, NDOFs )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Adaptive mesh refinement methods.
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!	REAL(KIND=dp) :: Quant(:)
!  	INPUT: Variable computed on the mesh on which the error indicator is based
!
!	INTEGER :: Perm(:)
!		INPUT: Element local numbering
!
!	INTEGER :: NDOFs
!		INPUT: Degrees of freedom of the variable Quant
!
!******************************************************************************
		IMPLICIT NONE

		TYPE(Solver_t), TARGET :: Solver
		INTEGER :: Perm(:), NDOFs
		REAL(KIND=dp) :: Quant(:)
		TYPE( Model_t ) :: Model
!------------------------------------------------------------------------------
		INTERFACE
			FUNCTION BoundaryResidual( Model,Edge,Mesh,Quant,Perm,Gnorm,NDOFs) RESULT(Indicator)
				USE Types
				TYPE(Element_t), POINTER :: Edge
				TYPE(Model_t) :: Model
				TYPE(Mesh_t), POINTER :: Mesh
				REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
				INTEGER :: Perm(:)
				INTEGER :: NDOFs
			END FUNCTION BoundaryResidual


			FUNCTION EdgeResidual( Model,Edge,Mesh,Quant,Perm,NDOFs) RESULT(Indicator)
				USE Types
				TYPE(Element_t), POINTER :: Edge
				TYPE(Model_t) :: Model
				TYPE(Mesh_t), POINTER :: Mesh
				REAL(KIND=dp) :: Quant(:), Indicator(2)
				INTEGER :: Perm(:)
				INTEGER :: NDOFs
			END FUNCTION EdgeResidual


			FUNCTION InsideResidual( Model,Element,Mesh,Quant,Perm,Fnorm,NDOFs) RESULT(Indicator)
				USE Types
				TYPE(Element_t), POINTER :: Element
				TYPE(Model_t) :: Model
				TYPE(Mesh_t), POINTER :: Mesh
				REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
				INTEGER :: Perm(:)
				INTEGER :: NDOFs
			END FUNCTION InsideResidual
		END INTERFACE
!------------------------------------------------------------------------------
		TYPE(Mesh_t), POINTER   :: RefMesh,NewMesh, Mesh
		TYPE( Nodes_t ) :: Nodes
		TYPE( Matrix_t ), POINTER :: NewMatrix
		INTEGER, POINTER :: Permutation(:)
		LOGICAL, POINTER       :: EdgeSplitted(:)
		INTEGER, POINTER       :: Referenced(:)
		TYPE( Element_t ), POINTER :: RefElement
		INTEGER :: i,j,k,n,nn,MarkedElements
		TYPE( Variable_t ), POINTER :: Var, Var1, NewVar
		REAL(KIND=dp) :: MaxError, ErrorLimit, minH, maxH, MaxChangeFactor, &
			LocalIndicator,ErrorEstimate,t,TotalTime,CPUTime,RealTime,RemeshTime,s
		LOGICAL :: BandwidthOptimize, Found, Coarsening, GlobalBubbles
		INTEGER :: MaxDepth, NLen
		CHARACTER(LEN=1024) :: Path
		CHARACTER(LEN=MAX_NAME_LEN) :: VarName
		REAL(KIND=dp), POINTER  :: Time(:), NodalError(:), PrevValues(:), &
			Hvalue(:),PrevNodalError(:), PrevHValue(:), hConvergence(:), ptr(:), tt(:)
		REAL(KIND=dp), POINTER  :: ErrorIndicator(:), eRef(:), hRef(:), Work(:)
!---------------------------------------------------------------------------------
!   	Initialize:
!------------------------------------------------------------------------------
		CALL Info( 'NumaRefineMesh', ' ', Level=5 )
		CALL Info( 'NumaRefineMesh', &
			'-----------N U M A   M E S H   R E F I N E M E N T --------------', Level=5 )
		TotalTime = CPUTime()
		RemeshTime = 0.0d0

		RefMesh => Solver % Mesh

		MaxDepth = ListGetInteger( Solver % Values, 'Adaptive Max Depth', Found )
		IF ( Found .AND. Refmesh % AdaptiveDepth > MaxDepth ) THEN
			WRITE( Message, * ) 'Max adaptive depth reached.'
			CALL Info( 'NumaRefineMesh', Message, Level = 6 )
			GOTO 20
		END IF

		DO i=1,RefMesh % NumberOfBulkElements
			RefMesh % Elements(i) % Splitted = 0
		END DO
		WRITE( Message, * ) 'RefMesh % NumberOfBulkElements= ',RefMesh % NumberOfBulkElements
		CALL Info( 'NumaRefineMesh', Message, Level = 6 )
!------------------------------------------------------------------------------
!     Compute the local error indicators:
!------------------------------------------------------------------------------
		t = CPUTime()
		CALL AllocateVector( ErrorIndicator, RefMesh % NumberOfBulkElements )

		MaxError = ComputeError( Model, ErrorIndicator, RefMesh, &
			Quant, Perm, InsideResidual, EdgeResidual, BoundaryResidual, NDOFs )
		WRITE( Message, * ) 'Error computation time (cpu-secs):               ',CPUTime()-t
		CALL Info( 'NumaRefineMesh', Message, Level = 6 )
!------------------------------------------------------------------------------
!   	Global error estimate:
!------------------------------------------------------------------------------
		ErrorEstimate =  SQRT( SUM( ErrorIndicator**2  ) )

		WRITE( Message, * ) 'Max error      =                                 ',MaxError
		CALL Info( 'NumaRefineMesh', Message, Level = 6 )
		WRITE( Message, * ) 'Error estimate =                                 ',ErrorEstimate
		CALL Info( 'NumaRefineMesh', Message, Level = 6 )
		WRITE(12,*) RefMesh % NumberOfBulkElements,ErrorEstimate,MaxError
!------------------------------------------------------------------------------
!   	Add nodal average of the h-value to the mesh variable list:
!------------------------------------------------------------------------------
		NN = RefMesh % NumberOfNodes

		Var => VariableGet( RefMesh % Variables, 'Hvalue', ThisOnly=.TRUE. )

		IF ( ASSOCIATED( Var ) ) THEN
			Hvalue      => Var % Values
			Var % PrimaryMesh => RefMesh
		ELSE
			CALL AllocateVector( Hvalue, nn )

			CALL VariableAdd( RefMesh % Variables, RefMesh, Solver, &
				'Hvalue', 1, Hvalue, Output=.FALSE. )
			Hvalue = 0.0d0
		END IF

		CALL AllocateVector( PrevHvalue, nn )
		IF ( RefMesh % AdaptiveDepth > 0 ) THEN
			PrevHvalue(1:nn) = Hvalue(1:nn)
		ELSE
			PrevHvalue(1:nn)= 0.0d0
		END IF

		CALL AllocateVector( Referenced, nn )

		Hvalue = 0.0d0
		Referenced = 0
		CALL AllocateVector( Nodes % x, RefMesh % MaxElementNodes )
		CALL AllocateVector( Nodes % y, RefMesh % MaxElementNodes )
		CALL AllocateVector( Nodes % z, RefMesh % MaxElementNodes )

		DO i=1,RefMesh % NumberOfBulkElements
			RefElement => RefMesh % Elements(i)
			n = RefElement % Type % NumberOfNodes

			Nodes % x(1:n) = RefMesh % Nodes % x(RefElement % NodeIndexes)
			Nodes % y(1:n) = RefMesh % Nodes % y(RefElement % NodeIndexes)
			Nodes % z(1:n) = RefMesh % Nodes % z(RefElement % NodeIndexes)
			s = ElementDiameter( RefElement, Nodes )
			DO j=1,n
				k = RefMesh % Elements(i) % NodeIndexes(j)
				Hvalue(k) = Hvalue(k) + s
				Referenced(k) = Referenced(k) + 1
			END DO
		END DO

		DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

		WHERE( Referenced(1:nn) > 0 )
			Hvalue(1:nn) = Hvalue(1:nn) / Referenced(1:nn)
		END WHERE
!------------------------------------------------------------------------------
!		Add estimate of the convergence with respecto to h:
!------------------------------------------------------------------------------
		Var => VariableGet( RefMesh % Variables, 'hConvergence', ThisOnly=.TRUE. )

		IF ( ASSOCIATED( Var ) ) THEN
			hConvergence => Var % Values
			Var % PrimaryMesh => RefMesh
		ELSE
			CALL AllocateVector( hConvergence, nn )
			hConvergence = 1.0d0

			CALL VariableAdd( RefMesh % Variables, RefMesh, Solver, &
				'hConvergence', 1, hConvergence, Output=.FALSE. )
		END IF
!------------------------------------------------------------------------------
!		Add nodal average of the computed estimate to the
!		solution error to the mesh variable list:
!------------------------------------------------------------------------------
		NLen = Solver % Variable % NameLen
		VarName = Solver % Variable % Name(1:NLen)
		Var => VariableGet( RefMesh % Variables, &
			VarName(1:NLen)  // '.error', ThisOnly=.TRUE. )

		IF ( ASSOCIATED( Var ) ) THEN
			NodalError  => Var % Values
			Var % PrimaryMesh => RefMesh
		ELSE
			CALL AllocateVector( NodalError, nn )

			CALL VariableAdd( RefMesh % Variables, RefMesh, Solver, &
				VarName(1:NLen) // '.error', 1, NodalError, Solver % Variable % Perm )
		END IF

		Var => VariableGet( RefMesh % Variables, &
			VarName(1:NLen) // '.perror', ThisOnly=.TRUE. )

		IF ( ASSOCIATED( Var ) ) THEN
			PrevNodalError  => Var % Values
			Var % PrimaryMesh => RefMesh
		ELSE
			CALL AllocateVector( PrevNodalError, RefMesh % NumberOfNodes )
			PrevNodalError = 0.0d0

			CALL VariableAdd( RefMesh % Variables, RefMesh, Solver, &
				VarName(1:NLen) // '.perror', 1, PrevNodalError, Output=.FALSE. )
		END IF

		NodalError = 0.0d0
		Referenced = 0
		DO i = 1, RefMesh % NumberOfBulkElements
			DO j=1,RefMesh % Elements(i) % Type % NumberOfNodes
				k = RefMesh % Elements(i) % NodeIndexes(j)
				Referenced(k) = Referenced(k) + 1
				NodalError(k) = NodalError(k) + ErrorIndicator(i)
			END DO
		END DO

		WHERE( Referenced(1:nn) > 0 )
			NodalError(1:nn) = NodalError(1:nn) / Referenced(1:nn)
		END WHERE
!------------------------------------------------------------------------------
!		Smooth error, if requested:
!------------------------------------------------------------------------------
		k = ListGetInteger( Solver % Values, 'Adaptive Pre Smoothing', Found )
		IF ( Found .AND. k > 0 ) THEN 
			CALL AllocateVector( eRef, nn )
			DO j=1,k
				eRef(1:nn) = NodalError(1:nn)
				Referenced = 0
				NodalError = 0
				DO i=1,RefMesh % NumberOfBulkElements
					n = RefMesh % Elements(i) % Type % NumberOfNodes
					NodalError(RefMesh % Elements(i) % NodeIndexes) = &
						NodalError(RefMesh % Elements(i) % NodeIndexes) + &
						SUM( eRef(RefMesh % Elements(i) % NodeIndexes) ) / n
					Referenced( RefMesh % Elements(i) % NodeIndexes ) = &
						Referenced( RefMesh % Elements(i) % NodeIndexes ) + 1
				END DO
				WHERE( Referenced(1:nn) > 1 )
					NodalError(1:nn) = NodalError(1:nn) / Referenced(1:nn)
				END WHERE
			END DO
			DEALLOCATE( eRef )
		END IF

		DEALLOCATE( Referenced )
!------------------------------------------------------------------------------
!		Add reference error to variable list:
!------------------------------------------------------------------------------
		Var => VariableGet( RefMesh % Variables, &
			VarName(1:NLen) // '.eRef', ThisOnly=.TRUE. )

		IF ( ASSOCIATED( Var ) ) THEN
			eRef => Var % Values
			Var % PrimaryMesh => RefMesh
		ELSE
			CALL AllocateVector( eRef, nn )
			eRef(1:nn) = NodalError(1:nn)

		CALL VariableAdd( RefMesh % Variables, RefMesh, Solver, &
			VarName(1:NLen) // '.eRef',1,eRef, Output=.FALSE. )
		END IF
!------------------------------------------------------------------------------
!		Mesh projection may alter the values somewhat!
!------------------------------------------------------------------------------
		eRef = MAX( eRef, 1.0d-12 )
!------------------------------------------------------------------------------
!		Add reference h to variable list:
!------------------------------------------------------------------------------
		Var => VariableGet( RefMesh % Variables, 'hRef', ThisOnly=.TRUE. )

		IF ( ASSOCIATED( Var ) ) THEN
			hRef => Var % Values
			Var % PrimaryMesh => RefMesh
		ELSE
			CALL AllocateVector( hRef, nn )
			hRef(1:nn) = Hvalue(1:nn)

			CALL VariableAdd( RefMesh % Variables, RefMesh, Solver, &
				'hRef', 1, hRef, Output=.FALSE. )
		END IF
!------------------------------------------------------------------------------
!		Mesh projection may alter the values somewhat!
!------------------------------------------------------------------------------
		hRef = MAX( hRef, 1.0d-12 )
!------------------------------------------------------------------------------
!		Check for convergence:
!------------------------------------------------------------------------------
		ErrorLimit = ListGetConstReal( Solver % Values, &
			'Adaptive Error Limit', Found )

		IF ( .NOT.Found ) ErrorLimit = 0.5d0

		IF ( MaxError < ErrorLimit ) THEN ! ErrorEstimate < ErrorLimit ) THEN
			CALL Info( 'NumaRefineMesh', &
				'Mesh convergence limit reached. I will do nothing.', Level=6 )
			GOTO 10
		END IF
!------------------------------------------------------------------------------
!   	Get additional parameters:
!------------------------------------------------------------------------------
		minH = ListGetConstReal( Solver % Values, 'Adaptive Min H', Found )
		maxH = ListGetConstReal( Solver % Values, 'Adaptive Max H', Found )

		MaxChangeFactor = ListGetConstReal( Solver % Values, &
			'Adaptive Max Change', Found )
		IF ( .NOT.Found .OR. MaxChangeFactor <= AEPS ) MaxChangeFactor = 3.0d0

		Coarsening = ListGetLogical( Solver % Values, 'Adaptive Coarsening', Found )
		if( .not.Found ) Coarsening = .TRUE.
!------------------------------------------------------------------------------
!		Compute local convergence of the solution with respect to h:
!------------------------------------------------------------------------------
		WHERE( eRef(1:nn) > 0 )
			PrevNodalError(1:nn) = PrevNodalError(1:nn) + &
				LOG( HValue(1:nn) / hRef(1:nn) ) * LOG( NodalError(1:nn) / eRef(1:nn) )
		END WHERE

		PrevHvalue(1:nn) = PrevHvalue(1:nn) + LOG( HValue(1:nn) / hRef(1:nn) )**2

		IF ( RefMesh % AdaptiveDepth > 0 ) THEN
			WHERE( PrevHValue(1:nn) > 0 )
				hConvergence(1:nn)  = MAX( PrevNodalError(1:nn) / PrevHValue(1:nn), 0.25d0 )
			ELSEWHERE
				hConvergence(1:nn)  = 0.25d0
			END WHERE
		END IF
!------------------------------------------------------------------------------
!   	Generate the new mesh:
!------------------------------------------------------------------------------
		NewMesh => SplitMesh( RefMesh, ErrorIndicator, ErrorLimit, &
		hConvergence, MaxChangeFactor )

		Hvalue(1:nn) = PrevHValue(1:nn)
		!NodalError = PrevNodalError

		IF ( .NOT.ASSOCIATED( NewMesh ) ) THEN
			CALL Info( 'NumaRefineMesh', &
				'Current mesh seems fine. I will do nothing.', Level=6 )
			GOTO 10
		END IF
		CALL Info( 'NumaRefineMesh', '-------------------------------------------', Level=5 )
		CALL Info( 'NumaRefineMesh', 'The new mesh consists of: ', Level=5 )
		WRITE( Message, * ) NewMesh % NumberOfNodes,' nodal points'
		CALL Info( 'NumaRefineMesh', Message, Level = 5 )
		WRITE( Message, * ) NewMesh % NumberOfBulkElements,' bulk elements'
		CALL Info( 'NumaRefineMesh', Message, Level = 5 )
		WRITE( Message, * ) NewMesh % NumberOfBoundaryElements,' boundary elements'
		CALL Info( 'NumaRefineMesh', Message, Level = 5 )
		CALL Info( 'NumaRefineMesh', '-------------------------------------------', Level=5 )
!-------------------------------------------------------------------
!		All the mesh geometry related tables are ready now,
!		next we update model and solver related tables:
!------------------------------------------------------------------------------
		t = CPUTime()
!------------------------------------------------------------------------------
!		Add the new mesh to the global list of meshes:
!------------------------------------------------------------------------------
		NewMesh % Next   => Model % Meshes 
		Model % Meshes   => NewMesh
		RefMesh % Child  => NewMesh
		NewMesh % Parent => RefMesh
		NULLIFY( NewMesh % Child )

		NewMesh % MaxBDOFs = RefMesh % MaxBDOFs

		NewMesh % Name = ListGetString( Solver % Values, &
			'Adaptive Mesh Name', Found )
		IF ( .NOT. Found ) NewMesh % Name = 'RefinedMesh'

		NewMesh % AdaptiveDepth = RefMesh % AdaptiveDepth + 1

		i = NewMesh % AdaptiveDepth
		n = FLOOR(LOG10(REAL(i))) + 1.5d0 
		Nlen = LEN_TRIM(NewMesh % Name)
		DO j=n,1,-1
			k = i / 10**(j-1)
			NewMesh % Name = NewMesh % Name(1:NLen) // CHAR(k+ICHAR('0'))
			i = i - k*10**(j-1)
		END DO

		Nlen = LEN_TRIM(OutputPath)
		IF ( Nlen > 0 ) THEN
			Path = OutputPath(1:Nlen) // '/' // TRIM(NewMesh % Name)
		ELSE
			Path = TRIM(NewMesh % Name)
		END IF
		CALL MakeDirectory( TRIM(path) // CHAR(0) )

		IF ( ListGetLogical( Solver % Values, 'Adaptive Save Mesh', Found ) ) &
			CALL WriteMeshToDisk( NewMesh, Path )
!------------------------------------------------------------------------------
!		Initialize local variables for the new mesh:
!------------------------------------------------------------------------------
		NULLIFY( NewMesh % Variables )

		CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
			'Coordinate 1', 1, NewMesh % Nodes % x )

		CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
			'Coordinate 2', 1, NewMesh % Nodes % y )

		CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
			'Coordinate 3', 1, NewMesh % Nodes % z )
!------------------------------------------------------------------------------
!		Time must always be there:
!------------------------------------------------------------------------------
		Var => VariableGet( RefMesh % Variables,'Time',ThisOnly=.TRUE. )
		CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
			'Time', 1, Var % Values )

		Var => VariableGet( RefMesh % Variables,'Timestep',ThisOnly=.TRUE. )
		CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
			'Timestep', 1, Var % Values )

		Var => VariableGet( RefMesh % Variables,'Timestep size',ThisOnly=.TRUE. )
		CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
			'Timestep size', 1, Var % Values )

		Var => VariableGet( RefMesh % Variables,'Timestep interval',ThisOnly=.TRUE. )
		CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
			'Timestep interval', 1, Var % Values )
!------------------------------------------------------------------------------
! 		Initialize the field variables for the new mesh. These are
! 		interpolated from the old meshes variables. Vector variables
! 		are in the variable lists in two ways: as vectors and as
! 		vector components. We MUST update the vectors (i.e. DOFs>1) first!!!!!
!------------------------------------------------------------------------------
		CALL SetCurrentMesh( Model, NewMesh )
		Var => RefMesh % Variables
		DO WHILE( ASSOCIATED( Var ) )
			IF ( Var % DOFs > 1 ) THEN
				NewVar => VariableGet( NewMesh % Variables,Var % Name,.FALSE. )
				k = SIZE( NewVar % Values )
				IF ( ASSOCIATED( NewVar % Perm ) ) THEN
					k = COUNT( NewVar % Perm > 0 )
				END IF
				IF ( NewVar % Name == 'flow solution' ) THEN
					NewVar % Norm = 0.0d0
					DO i=1,NewMesh % NumberOfNodes
						DO j=1,NewVar % DOFs-1
							NewVar % Norm = NewVar % Norm + &
								NewVar % Values( NewVar % DOFs*(i-1)+j )**2
						END DO
					END DO
					NewVar % Norm = SQRT( NewVar % Norm / k )
				ELSE
					NewVar % Norm = SQRT( SUM(NewVar % Values**2) / k )
				END IF
			END IF
			Var => Var % Next
		END DO
!------------------------------------------------------------------------------
!		Second time around, update scalar variables and vector components:
!------------------------------------------------------------------------------
		Var => RefMesh % Variables
		DO WHILE( ASSOCIATED( Var ) )
			SELECT CASE( Var % Name )
				CASE( 'coordinate 1', 'coordinate 2', 'coordinate 3', 'time' , &
					'timestep', 'timestep size', 'timestep interval' )
				CASE DEFAULT
					IF ( Var % DOFs == 1 ) THEN
						Found = .FALSE.
						Found = Found .OR. INDEX( Var % Name, '.error'  ) > 0
						Found = Found .OR. INDEX( Var % Name, '.eref'   ) > 0
						Found = Found .OR. INDEX( Var % Name, '.perror' ) > 0
						IF ( Found ) THEN
							k = Solver % Variable % NameLen
							IF ( Var % Name(1:k) /= Solver % Variable % Name(1:k) ) THEN
								Var => Var % Next
							CYCLE
							END IF
						END IF

						NewVar => VariableGet( NewMesh % Variables, Var % Name, .FALSE. )
						k = SIZE( NewVar % Values )
						IF ( ASSOCIATED( NewVar % Perm ) ) THEN
							k = COUNT( NewVar % Perm > 0 )
						END IF
						NewVar % Norm = SQRT( SUM(NewVar % Values**2) / k )
					END IF
			END SELECT
			Var => Var % Next
		END DO
!------------------------------------------------------------------------------
		WRITE( Message, * )  &
			'Mesh variable update time (cpu-secs):            ',CPUTime()-t
		CALL Info( 'NumaRefineMesh', Message, Level = 6 )
!------------------------------------------------------------------------------
!		Update Solver structure to use the new mesh:
!------------------------------------------------------------------------------ 
		Solver % Mesh => NewMesh
		CALL MeshStabParams( NewMesh )
!------------------------------------------------------------------------------
!		Nothing computed on this mesh yet:
!------------------------------------------------------------------------------
		NewMesh % SavesDone    = 0  ! start new output file
		NewMesh % OutputActive = .FALSE.
		NewMesh % Changed   = .TRUE.
!------------------------------------------------------------------------------
!		Update the solvers variable pointer:
!------------------------------------------------------------------------------
		Solver % Variable => VariableGet( Solver % Mesh % Variables, &
			Solver % Variable % Name, ThisOnly=.TRUE. )
		Solver % Variable % PrimaryMesh => NewMesh
!------------------------------------------------------------------------------
!		Create matrix structures for the new mesh:
!------------------------------------------------------------------------------
		t = CPUTime()
!------------------------------------------------------------------------------
!		Try to account for the reordering of DOFs
!		due to bandwidth optimization:
!------------------------------------------------------------------------------
		GlobalBubbles = ListGetLogical( Solver % Values, &
			'Bubbles in Global System', Found )
		IF ( .NOT. Found ) GlobalBubbles = .TRUE.

		BandwidthOptimize = ListGetLogical( Solver % Values, &
			'Optimize Bandwidth', Found )
		IF ( .NOT. Found ) BandwidthOptimize = .TRUE.

		IF ( BandwidthOptimize ) THEN
			n = NewMesh % NumberOfNodes
			IF ( GlobalBubbles ) &
				n = n + NewMesh % MaxBDOFs*NewMesh % NumberOFBulkElements
				CALL AllocateVector( Permutation,  n )
			ELSE
				Permutation => Solver % Variable % Perm
		END IF
!------------------------------------------------------------------------------
!		Create the CRS format matrix tables for solving the
!		current equation on the new mesh. Also do bandwidth
!		optimization, if requested:
!------------------------------------------------------------------------------
		NewMatrix => CreateMatrix( Model, Solver, Solver % Mesh,  &
			Permutation, Solver % Variable % DOFs, MATRIX_CRS, &
			BandwidthOptimize, ListGetString( Solver % Values, 'Equation' ), &
			GlobalBubbles=GlobalBubbles )

		IF ( ASSOCIATED( Solver % Matrix ) ) THEN
			CALL FreeMatrix( Solver % Matrix )
			NULLIFY( Solver % Matrix )
		END IF

		!Solver % Matrix % Child  => NewMatrix
		NewMatrix % Parent => Solver % Matrix
		NULLIFY( NewMatrix % Child )
		Solver % Matrix => NewMatrix
!------------------------------------------------------------------------------
!		Reorder the primary variable for bandwidth optimization:
!------------------------------------------------------------------------------
		IF ( BandwidthOptimize ) THEN
			n = Solver % Variable % DOFs
			ALLOCATE( Work(SIZE(Permutation)) )
			DO i=0,n-1
#if 0
					WHERE( Permutation > 0 )
						Solver % Variable % Values(n*Permutation-i) = &
							Solver % Variable % Values(n*Solver % Variable % Perm-i)
					END WHERE

					IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
						DO j=1,SIZE( Solver % Variable % PrevValues,2)
							WHERE( Permutation > 0 )
								Solver % Variable % PrevValues(n*Permutation-i,j) = &
									Solver % Variable % PrevValues(n*Solver % Variable % Perm-i,j)
							END WHERE
						END DO
					END IF
#else
					Work = Solver % Variable % Values(i+1::n)
					DO j=1,SIZE(Permutation)
						IF ( Permutation(j) > 0 ) THEN
							Solver % Variable % Values(n*Permutation(j)-i) = &
								Work(Solver % Variable % Perm(j))
						END IF
					END DO
					IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
						DO j=1,SIZE( Solver % Variable % PrevValues,2)
							Work = Solver % Variable % PrevValues(i+1::n,j)
							DO k=1,SIZE(Permutation)
								IF ( Permutation(k) > 0 ) THEN
									Solver % Variable % PrevValues(n*Permutation(k)-i,j) = &
										Work(Solver % Variable % Perm(k))
								END IF
							END DO
						END DO
					END IF
#endif
			END DO
			k = SIZE(Permutation)
			Solver % Variable % Perm(1:k) = Permutation
			DEALLOCATE( Permutation, Work )
		END IF
!------------------------------------------------------------------------------
!		TODO: CreateMatrix should do these
!------------------------------------------------------------------------------
		Solver % Matrix % Lumped = ListGetLogical( Solver % Values, &
			'Lumped Mass Matrix', Found )

		Solver % Matrix % Symmetric = ListGetLogical( Solver % Values, &
			'Linear System Symmetric', Found )

		CALL AllocateVector( Solver % Matrix % RHS, SIZE(Solver % Variable % Values) )
		Solver % Matrix % RHS = 0.0d0
!------------------------------------------------------------------------------
!		Transient case additional allocations:
!------------------------------------------------------------------------------
		IF (ListGetString(Model % Simulation,'Simulation Type')=='transient') THEN
			n = SIZE( Solver % Variable % Values )

			CALL AllocateArray( Solver % Matrix % Force, n, Solver % TimeOrder+1 )
			Solver % Matrix % Force = 0.0d0
		END IF
!------------------------------------------------------------------------------
!		Eigen analysis case additional allocations:
!------------------------------------------------------------------------------
		IF ( ListGetLogical(Solver % Values, 'Eigen Analysis', Found ) ) THEN
			n = SIZE( Solver % Variable % Values )
			CALL AllocateArray( Solver % Variable % EigenVectors, &
				Solver % NOFEigenValues, n )
			CALL AllocateVector( Solver % Variable % EigenValues, &
				Solver % NOFEigenValues )
		END IF

		CALL ParallelInitMatrix( Solver, Solver % Matrix )

		WRITE( Message, * ) 'Matrix structures update time (cpu-secs):        ',CPUTime()-t
		CALL Info( 'NumaRefineMesh', Message, Level=6 )
!------------------------------------------------------------------------------
!		Release previous meshes. Keep only the original mesh, and
!		the last two meshes:
!------------------------------------------------------------------------------
		Mesh => RefMesh % Parent
		DO WHILE( ASSOCIATED( Mesh ) )
			IF ( Mesh % AdaptiveDepth /= 0 ) THEN
				IF ( ASSOCIATED( Mesh % Parent ) ) THEN
					Mesh % Parent % Child => Mesh % Child
				END IF

				IF ( ASSOCIATED( Mesh % Child ) ) THEN
					Mesh % Child % Parent => Mesh % Parent
				END IF

				CALL ReleaseMesh( Mesh )
			END IF
			Mesh => Mesh % Parent
		END DO
!------------------------------------------------------------------------------
10		CONTINUE
!------------------------------------------------------------------------------
!		Comment the next calls, if you want to keep the edge tables:
!------------------------------------------------------------------------------
		CALL ReleaseMeshEdgeTables( RefMesh )
		CALL ReleaseMeshFaceTables( RefMesh )

		CALL SetCurrentMesh( Model, RefMesh )
		DEALLOCATE( ErrorIndicator, PrevHvalue )
!------------------------------------------------------------------------------
20		CONTINUE
!------------------------------------------------------------------------------
		WRITE( Message, * ) 'Mesh refine took in total (cpu-secs):           ', &
			CPUTIme() - TotalTime 
		CALL Info( 'NumaRefineMesh', Message, Level=6 )
		IF ( RemeshTime > 0 ) THEN
			WRITE( Message, * ) 'Remeshing took in total (real-secs):            ', &
				RemeshTime
			CALL Info( 'NumaRefineMesh', Message, Level=6 )
		END IF
		CALL Info( 'NumaRefineMesh', &
		'----------- E N D   M E S H   R E F I N E M E N T --------------', Level=5 )
!------------------------------------------------------------------------------
	CONTAINS
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
	FUNCTION SplitMesh( RefMesh,ErrorIndicator,ErrorLimit,  &
		hConvergence, MaxChange ) RESULT(NewMesh)
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Check which elements have to be refined and call the splitting method
!
!  ARGUMENTS:
!
!  TYPE(Mesh_t), POINTER :: RefMesh
!		INPUT: Current mesh to be refined
!
!	REAL(KIND=dp) :: ErrorIndicator(:)
!		INPUT: Refinement indicator computed on the current mesh
!
!	REAL(KIND=dp) :: ErrorLimit
!		INPUT: Limit value required for the refinement indicator
!
!	REAL(KIND=dp) :: hConvergence(:)
!		INPUT: ?
!
!	REAL(KIND=dp) :: MaxChange
!		INPUT: Maximum density change factor
!
!  TYPE(Mesh_t), POINTER :: NewMesh
!		OUTPUT: Refined mesh
!
!******************************************************************************
		REAL(KIND=dp) :: hConvergence(:), MaxChange
		TYPE(Mesh_t), POINTER :: NewMesh, RefMesh
		REAL(KIND=dp) :: ErrorIndicator(:),ErrorLimit
		!------------------------------------------------------------------------------
		TYPE(Mesh_t), POINTER :: NewMesh1
		REAL(KIND=dp) :: Lambda, EhConvergence
		INTEGER :: i,j,k,n,MarkedElements
		TYPE(Element_t), POINTER :: RefElement
!------------------------------------------------------------------------------
		NULLIFY( NewMesh )
!------------------------------------------------------------------------------
!		Determine the marked elements:
!------------------------------------------------------------------------------
		MarkedElements = 0
!------------------------------------------------------------------------------
!		Go through the bulk elements:
!------------------------------------------------------------------------------
		DO i = 1,RefMesh % NumberOfBulkElements
			RefElement => RefMesh % Elements(i)
!------------------------------------------------------------------------------
!			Check if the element is a linear triangle:
!------------------------------------------------------------------------------
			IF ( RefElement % Type % ElementCode /= 303 ) THEN
				CALL Fatal( 'SplitMesh', 'Internal splitting implemented only for linear triangles.' )
			END IF
!------------------------------------------------------------------------------
!			Get element description:
!------------------------------------------------------------------------------
			n = RefElement % Type % NumberOfNodes

			IF( RefMesh % AdaptiveDepth < 1 ) THEN
				EhConvergence = 1.0d0 ! First round: Assume full convergence speed
			ELSE
				EhConvergence = SUM( hConvergence( RefElement % Nodeindexes(1:n) ) ) / n
			END IF
!------------------------------------------------------------------------------
!			Check if the element has to be refined according to MaxChange:
!------------------------------------------------------------------------------
			RefElement % Splitted = 0
			IF( ErrorIndicator(i) > 100*AEPS ) THEN
				Lambda = ( ErrorLimit / ErrorIndicator(i) ) ** ( 1.0d0 / EhConvergence )
				RefElement % Splitted = MIN( MaxChange, 1.0d0/Lambda )
			END IF

			IF ( RefElement % Splitted > 0 ) MarkedElements = MarkedElements  + 1
!------------------------------------------------------------------------------
		END DO
!------------------------------------------------------------------------------
!		If no refinement required, quit:
!------------------------------------------------------------------------------
		IF ( MarkedElements == 0 ) THEN
			RefMesh % Changed = .FALSE.
			RETURN
		END IF
!------------------------------------------------------------------------------
!		Refine until all elements splitted specified times:
!------------------------------------------------------------------------------
		NewMesh => SplitOneLevel( RefMesh )
!------------------------------------------------------------------------------		
		DO WHILE( .TRUE. )
!------------------------------------------------------------------------------
			MarkedElements = 0
			DO i=1,NewMesh % NumberOfBulkElements
				IF ( NewMesh % Elements(i) % Splitted > 0 ) THEN
					MarkedElements = MarkedElements + 1
				END IF
			END DO

			IF ( MarkedElements == 0 ) EXIT

			NewMesh1 => SplitOneLevel( NewMesh )
			CALL ReleaseMesh( NewMesh )
			DEALLOCATE( NewMesh )

			NewMesh => NewMesh1
!------------------------------------------------------------------------------
		END DO
!------------------------------------------------------------------------------
	END FUNCTION SplitMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
	FUNCTION SplitOneLevel( RefMesh ) RESULT( NewMesh )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Compute the refined mesh from the current one by a splitting method
!
!  ARGUMENTS:
!
!  TYPE(Mesh_t), POINTER :: RefMesh
!		INPUT: Current mesh to be refined
!  
!  TYPE(Mesh_t), POINTER :: NewMesh
!		OUTPUT: Refined mesh 
!
!******************************************************************************
		IMPLICIT NONE

		TYPE( Mesh_t ), POINTER :: RefMesh, NewMesh
		!------------------------------------------------------------------------------
		REAL(KIND=dp) :: CPUTime,t

		INTEGER :: EdgeNumber,LongestEdge,Node1,Node2
		INTEGER :: i,j,k,l,n,NewElCnt,NewNodeCnt,MarkedEdges

		TYPE(Element_t), POINTER :: RefElement,Parent,Child,Edge

		LOGICAL, POINTER :: EdgeSplitted(:)
		INTEGER, POINTER :: MarkedOrder(:), Children(:,:)

		TYPE(PElementDefs_t), POINTER :: PD
		REAL(KIND=dp) :: x1, x2, y1, y2, EdgeLength, MaxLength
!------------------------------------------------------------------------------
!		Localization of the mesh edges
!------------------------------------------------------------------------------
		t = CPUTime()
		CALL Info( 'SplitOneLevel', '--------------------------------------------', Level=6 )
		CALL FindMeshEdges( RefMesh )
		WRITE( Message, * ) 'Find mesh edges time (cpu-secs):                 ',CPUTime()-t
		CALL Info( 'SplitOneLevel', Message, Level=6 )
!------------------------------------------------------------------------------
!		RGB Refinement:
!------------------------------------------------------------------------------
		t = CPUTime()
		CALL AllocateVector( EdgeSplitted, RefMesh % NumberOfEdges )
		MarkedEdges = RGBRefinement( EdgeSplitted,RefMesh )
		WRITE( Message, * ) 'RGB Refinement time (cpu-secs):                  ',CPUTime()-t
		CALL Info( 'SplitOneLevel', Message, Level=6 )
!------------------------------------------------------------------------------
!		Initialize the new mesh:
!------------------------------------------------------------------------------
		NewMesh => AllocateMesh()
		NewMesh % MaxElementNodes = 3
		NewMesh % MaxElementDOFs  = 3
!------------------------------------------------------------------------------
!		Create node tables for the new mesh:
!------------------------------------------------------------------------------ 
		t = CPUTime()
		NewMesh % NumberOfNodes = RefMesh % NumberOfNodes + MarkedEdges
		CALL AllocateVector( NewMesh % Nodes % x, NewMesh % NumberOfNodes )
		CALL AllocateVector( NewMesh % Nodes % y, NewMesh % NumberOfNodes )
		CALL AllocateVector( NewMesh % Nodes % z, NewMesh % NumberOfNodes )
!------------------------------------------------------------------------------
!		Add old nodes to the new mesh:
!------------------------------------------------------------------------------   
		NewMesh % Nodes % x(1:RefMesh % NumberOfNodes) = &
			RefMesh % Nodes % x(1:RefMesh % NumberOfNodes)
		NewMesh % Nodes % y(1:RefMesh % NumberOfNodes) = &
			RefMesh % Nodes % y(1:RefMesh % NumberOfNodes)
		NewMesh % Nodes % z(1:RefMesh % NumberOfNodes) = &
			RefMesh % Nodes % z(1:RefMesh % NumberOfNodes)
!------------------------------------------------------------------------------
!		Add new nodes to the new mesh:
!------------------------------------------------------------------------------
		NewNodeCnt = RefMesh % NumberOfNodes
		DO i = 1,RefMesh % NumberOfEdges
			IF ( EdgeSplitted(i) ) THEN
				Node1 = RefMesh % Edges(i) % NodeIndexes(1)
				Node2 = RefMesh % Edges(i) % NodeIndexes(2)
				x1 = RefMesh % Nodes % x(Node1)
				x2 = RefMesh % Nodes % x(Node2)
				y1 = RefMesh % Nodes % y(Node1)
				y2 = RefMesh % Nodes % y(Node2)
				NewNodeCnt = NewNodeCnt + 1
				NewMesh % Nodes % x(NewNodeCnt) = (x1+x2) / 2.0d0
				NewMesh % Nodes % y(NewNodeCnt) = (y1+y2) / 2.0d0
				NewMesh % Nodes % z(NewNodeCnt) = 0.0d0
			END IF
		END DO
		WRITE( Message, * ) 'Node tables generation time (cpu-secs):          ',CPUTime()-t
		CALL Info( 'SplitOneLevel', Message, Level=6 )
!------------------------------------------------------------------------------
!		Count the new number of bulk elements:
!------------------------------------------------------------------------------
		CALL AllocateVector( MarkedOrder, RefMesh % NumberOfEdges )
		MarkedOrder = 0

		k = 0
		NewElCnt = 0
		DO i = 1,RefMesh % NumberOfBulkElements
			MarkedEdges = 0
			DO j = 1,3
				EdgeNumber = RefMesh % Elements(i) % EdgeIndexes(j)
				IF( EdgeSplitted(EdgeNumber) ) THEN
					MarkedEdges = MarkedEdges + 1
					IF ( MarkedOrder(EdgeNumber) == 0 ) THEN
						k = k + 1
						MarkedOrder(EdgeNumber) = k + RefMesh % NumberOfNodes
					END IF
				END IF
			END DO
			NewElCnt = NewElCnt + MarkedEdges + 1
		END DO
		NewMesh % NumberOfBulkElements = NewElCnt
!------------------------------------------------------------------------------
!		Count the new number of boundary elements:
!------------------------------------------------------------------------------
		NewElCnt = 0
		DO i = RefMesh % NumberOfBulkElements+1,RefMesh % NumberOfBulkElements+&
			RefMesh % NumberOfBoundaryElements

			RefElement => RefMesh % Elements(i) % BoundaryInfo % Left
			IF ( .NOT.ASSOCIATED( RefElement) ) &
				RefElement => RefMesh % Elements(i) % BoundaryInfo % Right

			IF ( ASSOCIATED( RefElement ) ) THEN
				NULLIFY( Edge )

				DO j=1,3
					Edge => RefMesh % Edges(RefElement % EdgeIndexes(j))

					IF ( Edge % NodeIndexes(1) == RefMesh % Elements(i) % NodeIndexes(1) .AND. &
						Edge % NodeIndexes(2) == RefMesh % Elements(i) % NodeIndexes(2) .OR.  &
						Edge % NodeIndexes(2) == RefMesh % Elements(i) % NodeIndexes(1) .AND. &
						Edge % NodeIndexes(1) == RefMesh % Elements(i) % NodeIndexes(2) ) EXIT
				END DO

				IF ( EdgeSplitted( RefElement % EdgeIndexes(j) ) ) THEN
					NewElCnt = NewElCnt + 2
				ELSE
					NewElCnt = NewElCnt + 1
				END IF
			ELSE
				NewElCnt = NewElCnt + 1
			END IF
		END DO

		NewMesh % NumberOfBoundaryElements = NewElCnt
!------------------------------------------------------------------------------
!		Allocate element tables:
!------------------------------------------------------------------------------
		t = CPUTime()
		CALL AllocateVector( NewMesh % Elements, NewMesh % NumberOfBulkElements + &
			NewMesh % NumberOfBoundaryElements )

		CALL AllocateArray( Children, RefMesh % NumberOfBulkElements + &
			RefMesh % NumberOfBoundaryElements, 4 )
		Children = 0
!------------------------------------------------------------------------------
!		Find the new bulk elements:
!------------------------------------------------------------------------------
		NewElCnt    = 0
!------------------------------------------------------------------------------
		DO i = 1,RefMesh % NumberOfBulkElements
!------------------------------------------------------------------------------
			RefElement => RefMesh % Elements(i)
			n = RefElement % Type % NumberOfNodes

			MarkedEdges = 0
			DO j = 1,3
				EdgeNumber = RefElement % EdgeIndexes(j)
				IF ( EdgeSplitted(EdgeNumber) ) THEN
					MarkedEdges = MarkedEdges + 1
				END IF
			END DO
!------------------------------------------------------------------------------
!			Make elements for the new mesh:
!------------------------------------------------------------------------------
			SELECT CASE(MarkedEdges)
!------------------------------------------------------------------------------
				CASE(0)
!------------------------------------------------------------------------------
!					Just copy of the old one:
!------------------------------------------------------------------------------
					NewElCnt = NewElCnt + 1
					NewMesh % Elements(NewElCnt) = RefElement
					NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
					CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
					NewMesh % Elements(NewElCnt) % NodeIndexes(1:n) = &
						RefElement % NodeIndexes(1:n)

					Children(i,1) = NewElCnt
!------------------------------------------------------------------------------
				CASE(1)
!------------------------------------------------------------------------------				
!					Bisect the longest edge to give two triangles:
!------------------------------------------------------------------------------
					DO j = 1,3
						EdgeNumber = RefElement % EdgeIndexes(j)
						IF ( EdgeSplitted( EdgeNumber ) ) EXIT
					END DO
!------------------------------------------------------------------------------            
!					Find node (k) opposite to the splitted edge:
!------------------------------------------------------------------------------
					DO k = 1,3
						IF ( RefElement % NodeIndexes(k) /= &
							RefMesh % Edges(EdgeNumber) % NodeIndexes(1) .AND. &
							RefElement % NodeIndexes(k) /= &
							RefMesh % Edges(EdgeNumber) % NodeIndexes(2) ) EXIT
					END DO
!------------------------------------------------------------------------------
!					New element 1
!------------------------------------------------------------------------------
					NewElCnt = NewElCnt + 1
					NewMesh % Elements(NewElCnt) = RefElement
					NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
					CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

					NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
						RefElement % NodeIndexes(k)
					NewMesh % Elements(NewElCnt) % NodeIndexes(2) = &
						RefMesh % Edges(EdgeNumber) % NodeIndexes(1)

					NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
						MarkedOrder(RefElement % EdgeIndexes(j))

					Children(i,1) = NewElCnt
!------------------------------------------------------------------------------
!					New element 2
!------------------------------------------------------------------------------
					NewElCnt = NewElCnt + 1
					NewMesh % Elements(NewElCnt) = RefElement
					NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
					CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

					NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
						RefElement % NodeIndexes(k)

					NewMesh % Elements(NewElCnt) % NodeIndexes(2) = &
						MarkedOrder(RefElement % EdgeIndexes(j))

					NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
						RefMesh % Edges(EdgeNumber) % NodeIndexes(2)

					Children(i,2) = NewElCnt
!------------------------------------------------------------------------------
			CASE(2)
!------------------------------------------------------------------------------
!				Bisect two of the edges to give three new elements:
!------------------------------------------------------------------------------
!				Find the edge NOT splitted:
!------------------------------------------------------------------------------
				DO j = 1,3
					EdgeNumber = RefElement % EdgeIndexes(j)
					IF ( .NOT.EdgeSplitted( EdgeNumber ) ) EXIT             
				END DO
!------------------------------------------------------------------------------
!				Find node (k) opposite to the edge NOT splitted:
!------------------------------------------------------------------------------
				DO k = 1,3
					IF (RefElement % NodeIndexes(k) /= &
						RefMesh % Edges(EdgeNumber) % NodeIndexes(1) .AND. &
						RefElement % NodeIndexes(k) /= &
						RefMesh % Edges(EdgeNumber) % NodeIndexes(2) ) EXIT
				END DO
!------------------------------------------------------------------------------
!				New element 1
!------------------------------------------------------------------------------
				NewElCnt = NewElCnt + 1
				NewMesh % Elements(NewElCnt) = RefElement
				NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
				CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

				NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
					RefElement % NodeIndexes(k)

				l = 1
				DO k = 1,3
					IF ( k /= j ) THEN
						l = l + 1
						NewMesh % Elements(NewElCnt) % NodeIndexes(l) = &
							MarkedOrder(RefElement % EdgeIndexes(k))
					END IF
				END DO

				Children(i,1) = NewElCnt
!------------------------------------------------------------------------------
!				New element 2
!------------------------------------------------------------------------------
				NewElCnt = NewElCnt + 1
				NewMesh % Elements(NewElCnt) = RefElement
				NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
				CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

				l = 0
				DO k = 1,3
					IF ( k /= j ) THEN
						l = l + 1
						NewMesh % Elements(NewElCnt) % NodeIndexes(l) = &
							MarkedOrder(RefElement % EdgeIndexes(k))
					END IF
				END DO

				MaxLength = 0.0d0
				DO k = 1,3
					IF ( k /= j ) THEN
						EdgeNumber = RefElement % EdgeIndexes(k)
						Node1 = RefMesh % Edges( EdgeNumber ) % NodeIndexes(1)
						Node2 = RefMesh % Edges( EdgeNumber ) % NodeIndexes(2)
						x1 = RefMesh % Nodes % x( Node1 )
						x2 = RefMesh % Nodes % x( Node2 )
						y1 = RefMesh % Nodes % y( Node1 )
						y2 = RefMesh % Nodes % y( Node2 )
						EdgeLength = ((x2-x1)**2+(y2-y1)**2)**0.5
						IF (EdgeLength >= MaxLength) THEN
							MaxLength = EdgeLength
							LongestEdge = k
						END IF
					END IF
				END DO
				k = LongestEdge
				if ( k <= 0 .or. k > 3 ) print*,k

				IF ( RefMesh % Edges(RefElement % EdgeIndexes(j)) % NodeIndexes(1) ==  &
					RefMesh % Edges(RefElement % EdgeIndexes(k)) % NodeIndexes(1) .OR.&
					RefMesh % Edges(RefElement % EdgeIndexes(j)) % NodeIndexes(1) ==  &
					RefMesh % Edges(RefElement % EdgeIndexes(k)) % NodeIndexes(2) ) THEN
					NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
					RefMesh % Edges(RefElement % EdgeIndexes(j)) % NodeIndexes(2)
				ELSE
					NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
					RefMesh % Edges(RefElement % EdgeIndexes(j)) % NodeIndexes(1)
				END IF

				Children(i,2) = NewElCnt
!------------------------------------------------------------------------------
!				New element 3
!------------------------------------------------------------------------------
				NewElCnt = NewElCnt + 1
				NewMesh % Elements(NewElCnt) = RefElement
				NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
				CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

				DO j = 1,3
					EdgeNumber = RefElement % EdgeIndexes(j)
					IF ( .NOT.EdgeSplitted( EdgeNumber ) ) EXIT             
				END DO

				DO k = 1,2
					NewMesh % Elements(NewElCnt) % NodeIndexes(k) = &
						RefMesh % Edges(EdgeNumber) % NodeIndexes(k)
				END DO

				NewMesh % Elements(NewElCnt) % NodeIndexes(3) = &
					MarkedOrder(RefElement % EdgeIndexes(LongestEdge))

				Children(i,3) = NewElCnt
!------------------------------------------------------------------------------
			CASE(3)
!------------------------------------------------------------------------------		
!				Bisect all the edges to give four new elements:
!------------------------------------------------------------------------------
!				New element 1
!------------------------------------------------------------------------------
				NewElCnt = NewElCnt + 1
				NewMesh % Elements(NewElCnt) = RefElement
				NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
				CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

				NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
					RefElement % NodeIndexes(1)

				j = RefElement % EdgeIndexes(1)
				NewMesh % Elements(NewElCnt) % NodeIndexes(2) = MarkedOrder(j)

				j = RefElement % EdgeIndexes(3)
				NewMesh % Elements(NewElCnt) % NodeIndexes(3) = MarkedOrder(j)

				Children(i,1) = NewElCnt
!------------------------------------------------------------------------------
!				New element 2
!------------------------------------------------------------------------------
				NewElCnt = NewElCnt + 1
				NewMesh % Elements(NewElCnt) = RefElement
				NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
				CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

				NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
					RefElement % NodeIndexes(2)

				j = RefElement % EdgeIndexes(2)
				NewMesh % Elements(NewElCnt) % NodeIndexes(2) = MarkedOrder(j)

				j = RefElement % EdgeIndexes(1)
				NewMesh % Elements(NewElCnt) % NodeIndexes(3) = MarkedOrder(j)

				Children(i,2) = NewElCnt
!------------------------------------------------------------------------------
!				New element 3
!------------------------------------------------------------------------------
				NewElCnt = NewElCnt + 1
				NewMesh % Elements(NewElCnt) = RefElement
				NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
				CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

				NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
					RefElement % NodeIndexes(3)

				j = RefElement % EdgeIndexes(3)
				NewMesh % Elements(NewElCnt) % NodeIndexes(2) = MarkedOrder(j)

				j = RefElement % EdgeIndexes(2)
				NewMesh % Elements(NewElCnt) % NodeIndexes(3) = MarkedOrder(j)

				Children(i,3) = NewElCnt
!------------------------------------------------------------------------------
!				New element 4
!------------------------------------------------------------------------------
				NewElCnt = NewElCnt + 1
				NewMesh % Elements(NewElCnt) = RefElement
				NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
				CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )

				DO j=1,n
					NewMesh % Elements(NewElCnt) % NodeIndexes(j) = &
						MarkedOrder( RefElement % EdgeIndexes(j) )
				END DO

				Children(i,4) = NewElCnt
!------------------------------------------------------------------------------
			END SELECT
!------------------------------------------------------------------------------
			DO j=1,4
				k = Children(i,j)
				IF ( k > 0 ) THEN
					NewMesh % Elements(k) % Splitted = RefElement % Splitted-1
				END IF
			END DO
!------------------------------------------------------------------------------
		END DO ! Bulk elements
!------------------------------------------------------------------------------
		WRITE( Message, * ) 'Bulk element tables generation time (cpu-secs):  ',CPUTime()-t
		CALL Info( 'SplitOneLevel', Message, Level=6 )
!------------------------------------------------------------------------------
!		Update boundary elements:
!------------------------------------------------------------------------------
		t = CPUTime()
		NewElCnt = NewMesh % NumberOfBulkElements
!------------------------------------------------------------------------------
		DO j = RefMesh % NumberOfBulkElements + 1, &
			RefMesh % NumberOfBulkElements + RefMesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
			RefElement => RefMesh % Elements(j) % BoundaryInfo % Left
			IF ( .NOT.ASSOCIATED( RefElement) ) &
				RefElement => RefMesh % Elements(j) % BoundaryInfo % Right
!------------------------------------------------------------------------------
			IF ( ASSOCIATED( RefElement ) ) THEN
!------------------------------------------------------------------------------
				NULLIFY( Edge )
				DO i=1,3
					Edge => RefMesh % Edges(RefElement % EdgeIndexes(i))
					IF ( Edge % NodeIndexes(1) == RefMesh % Elements(j) % NodeIndexes(1) .AND. &
						Edge % NodeIndexes(2) == RefMesh % Elements(j) % NodeIndexes(2) .OR.  &
						Edge % NodeIndexes(2) == RefMesh % Elements(j) % NodeIndexes(1) .AND. &
						Edge % NodeIndexes(1) == RefMesh % Elements(j) % NodeIndexes(2) ) EXIT
				END DO
				EdgeNumber = RefElement % EdgeIndexes(i)

				RefElement => RefMesh % Elements(j)
				n = RefElement % Type % NumberOfNodes
!------------------------------------------------------------------------------
				IF ( EdgeSplitted(EdgeNumber) ) THEN
!------------------------------------------------------------------------------
!					New element 1:
!------------------------------------------------------------------------------
					NewElCnt = NewElCnt + 1
					NewMesh % Elements(NewElCnt) = RefElement
					NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
					CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
					NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
						RefElement % NodeIndexes(1)
					NewMesh % Elements(NewElCnt) % NodeIndexes(2) = &
						MarkedOrder(EdgeNumber)

					ALLOCATE( NewMesh % Elements(NewElCnt) % BoundaryInfo )
					NewMesh % Elements(NewElCnt) % BoundaryInfo = &
						RefElement % BoundaryInfo

					CALL SetParents( NewMesh % Elements(NewElCnt), &
						NewMesh, Children, Edge )

					Children(j,1) = NewElCnt             
!------------------------------------------------------------------------------
!					New element 2:
!------------------------------------------------------------------------------
					NewElCnt = NewElCnt + 1
					NewMesh % Elements(NewElCnt) = RefElement
					NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
					CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
					NewMesh % Elements(NewElCnt) % NodeIndexes(1) = &
						MarkedOrder(EdgeNumber)
					NewMesh % Elements(NewElCnt) % NodeIndexes(2) = &
						RefElement % NodeIndexes(2)

					ALLOCATE( NewMesh % Elements(NewElCnt) % BoundaryInfo )
					NewMesh % Elements(NewElCnt) % BoundaryInfo = &
						RefElement % BoundaryInfo

					CALL SetParents( NewMesh % Elements(NewElCnt), &
						NewMesh, Children, Edge )

					Children(j,2) = NewElCnt
!------------------------------------------------------------------------------					
				ELSE
!------------------------------------------------------------------------------
!					New element 1:
!------------------------------------------------------------------------------
					NewElCnt = NewElCnt + 1
					NewMesh % Elements(NewElCnt) = RefElement
					NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
					CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
					NewMesh % Elements(NewElCnt) % NodeIndexes = &
						RefElement % NodeIndexes

					ALLOCATE( NewMesh % Elements(NewElCnt) % BoundaryInfo )

					NewMesh % Elements(NewElCnt) % BoundaryInfo = &
						RefElement % BoundaryInfo

					CALL SetParents( NewMesh % Elements(NewElCnt), &
						NewMesh, Children, Edge )

					Children(j,1) = NewElCnt
!------------------------------------------------------------------------------
				END IF
!------------------------------------------------------------------------------
			ELSE
!------------------------------------------------------------------------------
!				New element 1, this is point element:
!------------------------------------------------------------------------------
				NewElCnt = NewElCnt + 1
				RefElement => RefMesh % Elements(j)
				n = RefElement % Type % NumberOfNodes

				NewMesh % Elements(NewElCnt) = RefElement
				NewMesh % Elements(NewElCnt) % ElementIndex = NewElCnt
				CALL AllocateVector( NewMesh % Elements(NewElCnt) % NodeIndexes,n )
				NewMesh % Elements(NewElCnt) % NodeIndexes = &
					RefElement % NodeIndexes

				ALLOCATE( NewMesh % Elements(NewElCnt) % BoundaryInfo )

				NewMesh % Elements(NewElCnt) % BoundaryInfo = &
					RefElement % BoundaryInfo

				NULLIFY( NewMesh % Elements(NewElCnt) % BoundaryInfo % Left )
				NULLIFY( NewMesh % Elements(NewElCnt) % BoundaryInfo % Right )

				Children(j,1) = NewElCnt
!------------------------------------------------------------------------------
			END IF
!------------------------------------------------------------------------------
		END DO
!------------------------------------------------------------------------------
		NewMesh % MaxBDOFs = RefMesh % MaxBDOFs
		DO i = 1,NewMesh % NumberOfBulkElements+NewMesh % NumberOfBoundaryElements
			RefElement => NewMesh % Elements(i)
			NULLIFY( RefElement % PDefs )
			NULLIFY( RefElement % DGIndexes )
			NULLIFY( RefElement % EdgeIndexes )
			NULLIFY( RefElement % FaceIndexes )
			NULLIFY( RefElement % BubbleIndexes )
			IF ( RefElement % BDOFs > 0 ) THEN
				ALLOCATE( RefElement % BubbleIndexes(RefElement % BDOFs) )
				DO j=1,RefElement % BDOFs
					RefElement % BubbleIndexes(j) = NewMesh % MaxBDOFs*(i-1)+j
				END DO
			END IF

			IF ( ASSOCIATED(RefElement % PDefs) ) THEN
				PD => RefElement % PDefs
				CALL AllocatePDefinitions( RefElement )
				RefElement % PDefs = PD
			END IF
		END DO
!------------------------------------------------------------------------------
		WRITE( Message, * ) 'Bndry element tables generation time (cpu-secs): ',CPUTime()-t
		CALL Info( 'SplitOneLevel', Message, Level=6 )

		DEALLOCATE( EdgeSplitted, MarkedOrder, Children )
!------------------------------------------------------------------------------
	END FUNCTION SplitOneLevel
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
	FUNCTION RGBRefinement(  EdgeSplitted,RefMesh ) RESULT(MarkedEdges)
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Decide which edges have to be splitted in the current mesh, according to 
!	RGB-refinement method
!
!  ARGUMENTS:
!
!	LOGICAL :: EdgeSplitted(:)
!		OUTPUT: Indicate which edge has to be splitted in SplitOneLevel(...) 
!
!  TYPE(Mesh_t), POINTER :: RefMesh
!		INPUT: Current mesh to be refined
!
!	INTEGER :: MarkedEdges
!		OUTPUT: Number of edges to be splitted in SplitOneLevel(...) 
!		(= number of new nodes in the mesh)
!		
!******************************************************************************
		IMPLICIT NONE

		LOGICAL :: EdgeSplitted(:)
		INTEGER :: MarkedEdges
		TYPE(Mesh_t), POINTER :: RefMesh
!------------------------------------------------------------------------------
		LOGICAL :: MarkedEdgesFound
		INTEGER :: i,j,EdgeNumber,HangingNodes,RGBIterations,Node1,Node2,&
			LongestEdge
		REAL(KIND=dp) :: x1,y1,x2,y2,EdgeLength,MaxLength
!------------------------------------------------------------------------------
		EdgeSplitted = .FALSE.
!------------------------------------------------------------------------------
!		Mark all three edges of the marked elements (RED refinement):
!------------------------------------------------------------------------------
!     DO i = 1,RefMesh % NumberOfBulkElements
!        IF ( RefMesh % Elements(i) % Splitted > 0 ) THEN
!           DO j = 1,3
!              EdgeNumber = RefMesh % Elements(i) % EdgeIndexes(j)
!              EdgeSplitted( EdgeNumber ) = .TRUE.
!           END DO
!        END IF
!     END DO
!------------------------------------------------------------------------------
!		Mark the longest edges of the marked elements (GREEN refinement):
!------------------------------------------------------------------------------
		DO i = 1,RefMesh % NumberOfBulkElements
!------------------------------------------------------------------------------
			IF ( RefMesh % Elements(i) % Splitted > 0 ) THEN
				MaxLength   = 0.0D0
				LongestEdge = 0
				DO j = 1,3
					EdgeNumber = RefMesh % Elements(i) % EdgeIndexes(j)
					Node1 = RefMesh % Edges( EdgeNumber ) % NodeIndexes(1)
					Node2 = RefMesh % Edges( EdgeNumber ) % NodeIndexes(2)
					x1 = RefMesh % Nodes % x( Node1 )
					x2 = RefMesh % Nodes % x( Node2 )
					y1 = RefMesh % Nodes % y( Node1 )
					y2 = RefMesh % Nodes % y( Node2 )
					EdgeLength = SQRT((x2-x1)**2+(y2-y1)**2)
					IF (EdgeLength >= MaxLength) THEN
						MaxLength = EdgeLength
						LongestEdge = EdgeNumber
					END IF
				END DO
				EdgeSplitted( LongestEdge ) = .TRUE.
			END IF
!------------------------------------------------------------------------------
		END DO
!------------------------------------------------------------------------------
		MarkedEdges = 0
		DO i = 1,RefMesh % NumberOfEdges
			IF ( EdgeSplitted(i) ) THEN
				MarkedEdges = MarkedEdges + 1
			END IF
		END DO
!------------------------------------------------------------------------------
!		Mark longest edges until we have a RGB-refinement:
!------------------------------------------------------------------------------
		RGBiterations = 0 
!------------------------------------------------------------------------------
		DO WHILE( .TRUE. )
!------------------------------------------------------------------------------
			HangingNodes = 0
			RGBiterations = RGBiterations+1
!------------------------------------------------------------------------------
			DO i = 1,RefMesh % NumberOfBulkElements
!------------------------------------------------------------------------------
!				Check for marked edges and find the longest edge:
!------------------------------------------------------------------------------
				MarkedEdgesFound = .FALSE.
				LongestEdge      = 0
				MaxLength        = 0.0d0
				DO j = 1,3
					EdgeNumber = RefMesh % Elements(i) % EdgeIndexes(j)
					MarkedEdgesFound = MarkedEdgesFound.OR.EdgeSplitted(EdgeNumber)
					Node1 = RefMesh % Edges(EdgeNumber) % NodeIndexes(1)
					Node2 = RefMesh % Edges(EdgeNumber) % NodeIndexes(2)
					x1 = RefMesh % Nodes % x( Node1 )
					x2 = RefMesh % Nodes % x( Node2 )
					y1 = RefMesh % Nodes % y( Node1 )
					y2 = RefMesh % Nodes % y( Node2 )
					EdgeLength = ((x2-x1)**2+(y2-y1)**2)**0.5
					IF (EdgeLength >= MaxLength) THEN
						MaxLength = EdgeLength
						LongestEdge = EdgeNumber
					END IF
				END DO
!------------------------------------------------------------------------------
!				If there are marked edges, the longest edge must be one of them:
!------------------------------------------------------------------------------
				IF ( MarkedEdgesFound.AND.(.NOT.EdgeSplitted(LongestEdge)) ) THEN
					HangingNodes = HangingNodes + 1
					EdgeSplitted( LongestEdge ) = .TRUE.
				END IF
!------------------------------------------------------------------------------
			END DO
!------------------------------------------------------------------------------
			IF( HangingNodes > 0) THEN
				WRITE( Message, * ) 'RGB ',RGBiterations,' : ',HangingNodes,' new nodes'
				CALL Info( 'RGBRefinement', Message, Level=6 )
				MarkedEdges = MarkedEdges + HangingNodes
			ELSE
				EXIT
			END IF
!------------------------------------------------------------------------------
		END DO
!------------------------------------------------------------------------------
	END FUNCTION RGBRefinement
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
	SUBROUTINE SetParents( Element, Mesh, Children, Edge )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Find the parent elements to the splitted boundary element among the children 
!	of the original parent element:
!
!  ARGUMENTS:
!
!	TYPE(Element_t) :: Element
!		INPUT/OUTPUT: 
!
!	TYPE(Mesh_t), POINTER :: Mesh
!		INPUT: New mesh
!
!	INTEGER :: Children(:,:)
!  	INPUT: Children of the original parent element
!
!	TYPE(Element_t), POINTER :: Edge
!		INPUT: 
!		
!
!******************************************************************************
		TYPE(Element_t) :: Element
		TYPE(Element_t), POINTER :: Edge

		INTEGER :: Children(:,:)
		TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
		INTEGER j,k,l,n,i0,j0,k0

		TYPE(Element_t), POINTER :: Child
!------------------------------------------------------------------------------
!		Get elememt description:
!------------------------------------------------------------------------------
		n = Element % Type % NumberOfNodes
!------------------------------------------------------------------------------
!		Get edge description:
!------------------------------------------------------------------------------
		k = Edge % BoundaryInfo % Left % ElementIndex
		
		NULLIFY( Child )
		
		DO l=1,4
			IF ( Children(k,l)>0 ) THEN
				Child => Mesh % Elements( Children(k,l) )
				i0 = 0
				DO j0=1,n
					DO k0=1,Child % Type % NumberOfNodes
						IF ( Child % NodeIndexes(k0) == Element % NodeIndexes(j0) ) THEN
							i0 = i0 + 1 
							EXIT
						END IF
					END DO
				END DO
				IF ( i0 == n ) EXIT
			END IF
		END DO

		IF ( l > 4 ) STOP 'Adaptive: parent 1 not found'
        
		Element % BoundaryInfo % Left  => Child
		NULLIFY( Element % BoundaryInfo % Right )
!------------------------------------------------------------------------------
		NULLIFY( Child )
		IF ( ASSOCIATED(Edge % BoundaryInfo % Right) ) THEN
			k = Edge % BoundaryInfo % Right % ElementIndex
			DO l=1,4
				IF ( Children(k,l)>0 ) THEN
					Child => Mesh % Elements( Children(k,l) )
					i0 = 0
					DO j0=1,n
						DO k0=1,Child % Type % NumberOfNodes
							IF ( Child % NodeIndexes(k0) == Element % NodeIndexes(j0) ) THEN
								i0 = i0 + 1 
								EXIT
							END IF
						END DO
					END DO
					IF ( i0 == n ) EXIT
				END IF
			END DO

			Element % BoundaryInfo % Right => Child
		END IF
!------------------------------------------------------------------------------
	END SUBROUTINE SetParents
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
 END SUBROUTINE NumaRefineMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
	FUNCTION ComputeError( Model, ErrorIndicator, RefMesh,  &
		Quant, Perm, InsideResidual, EdgeResidual, BoundaryResidual, NDOFs ) RESULT(MaxError)
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Compute a refinement indicator from the residuals computed by InsideResidual, 
!	EdgeResidual and BoundaryResidual.
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!	REAL(KIND=dp) :: ErrorIndicator(:)
!		OUTPUT: Refinement Indicator for the splitting method
!
! 	TYPE(Mesh_t), POINTER :: RefMesh
!		INPUT: Current mesh to be refined
!
!	REAL(KIND=dp) :: Quant(:)
!  	INPUT: Variable computed on the mesh on which the error indicator is based
!
!	INTEGER :: Perm(:)
!		INPUT: Element local numbering
!
!	INTEGER :: NDOFs
!		INPUT: Degrees of freedom of the variable Quant
!
!******************************************************************************
		USE CRSMATRIX
		IMPLICIT NONE
!------------------------------------------------------------------------------
		TYPE(Mesh_t), POINTER :: RefMesh
		TYPE(Model_t) :: Model
		INTEGER :: Perm(:),NDOFs
		REAL(KIND=dp) :: ErrorIndicator(:), Quant(:), MaxError
!------------------------------------------------------------------------------
		INTERFACE
			FUNCTION BoundaryResidual( Model,Edge,Mesh,Quant,Perm,Gnorm,NDOFs ) RESULT(Indicator)
				USE Types
				TYPE(Element_t), POINTER :: Edge
				TYPE(Model_t) :: Model
				TYPE(Mesh_t), POINTER :: Mesh
				REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
				INTEGER :: Perm(:), NDOFs
			END FUNCTION BoundaryResidual

			FUNCTION EdgeResidual( Model,Edge,Mesh,Quant,Perm,NDOFs ) RESULT(Indicator)
				USE Types
				TYPE(Element_t), POINTER :: Edge
				TYPE(Model_t) :: Model
				TYPE(Mesh_t), POINTER :: Mesh
				REAL(KIND=dp) :: Quant(:), Indicator(2)
				INTEGER :: Perm(:), NDOFs
			END FUNCTION EdgeResidual

			FUNCTION InsideResidual( Model,Element,Mesh,Quant,Perm,Fnorm,NDOFs ) RESULT(Indicator)
				USE Types
				TYPE(Element_t), POINTER :: Element
				TYPE(Model_t) :: Model
				TYPE(Mesh_t), POINTER :: Mesh
				REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
				INTEGER :: Perm(:), NDOFs
			END FUNCTION InsideResidual
		END INTERFACE
!------------------------------------------------------------------------------
		TYPE(Element_t), POINTER :: Edge, Element
		INTEGER :: i, j, k, Parent
		REAL(KIND=dp), POINTER :: TempIndicator(:,:)
		REAL(KIND=dp) :: LocalIndicator(2), Fnorm, LocalFnorm,s
!------------------------------------------------------------------------------
!		Localize the mesh edges and do some initializations:
!------------------------------------------------------------------------------
		CALL FindMeshEdges( RefMesh )

		Fnorm = 0.0d0
		ErrorIndicator = 0.0d0

		CALL AllocateArray( TempIndicator, 2,SIZE(ErrorIndicator) )
		TempIndicator = 0.0d0
!------------------------------------------------------------------------------
!		Bulk equation residuals:
!------------------------------------------------------------------------------
		DO i=1,RefMesh % NumberOfBulkElements
			Element => RefMesh % Elements(i)
			CurrentModel % CurrentElement => Element

			LocalIndicator = InsideResidual( Model, Element, &
				RefMesh, Quant, Perm, LocalFnorm, NDOFs )

			Fnorm = Fnorm + LocalFnorm
			TempIndicator(:,i) = TempIndicator(:,i) + LocalIndicator
			
			!PRINT *, 'Ind= ',TempIndicator(:,i)
		END DO
!------------------------------------------------------------------------------
!		Boundary condition residuals:
!------------------------------------------------------------------------------
		DO i = RefMesh % NumberOfBulkElements + 1,  &
			RefMesh % NumberOfBulkElements + RefMesh % NumberOfBoundaryElements
			Edge => RefMesh % Elements(i)
			CurrentModel % CurrentElement => Edge

			IF ( Edge % Type % ElementCode == 101 ) CYCLE

			LocalIndicator = BoundaryResidual( Model, Edge, &
				RefMesh, Quant, Perm, LocalFnorm, NDOFs )

			Fnorm = Fnorm + LocalFnorm

			IF ( ASSOCIATED( Edge % BoundaryInfo % Left) ) THEN
				Parent = Edge % BoundaryInfo % Left % ElementIndex
				IF ( Parent > 0 ) TempIndicator( :,Parent ) = &
					TempIndicator( :,Parent ) + LocalIndicator
			END IF

			IF ( ASSOCIATED( Edge % BoundaryInfo % RIght) ) THEN
				Parent = Edge % BoundaryInfo % Right % ElementIndex
				IF ( Parent > 0 ) TempIndicator( :,Parent ) = &
					TempIndicator( :,Parent ) + LocalIndicator
			END IF
		END DO
!------------------------------------------------------------------------------
		s = SQRT( SUM(TempIndicator(2,:)) ) / SQRT( SUM(TempIndicator(1,:)) )
		!PRINT *,'s= ',s
!------------------------------------------------------------------------------
		ErrorIndicator = SQRT( TempIndicator(1,:)/(2*s) + s*TempIndicator(2,:)/2 )
		!PRINT *,'ErrorIndicator= ', ErrorIndicator
!------------------------------------------------------------------------------
		IF ( Fnorm > AEPS ) THEN
		ErrorIndicator = ErrorIndicator / SQRT( Fnorm )
		END IF
!------------------------------------------------------------------------------
		MaxError = MAXVAL( ErrorIndicator )
		DEALLOCATE( TempIndicator )
!------------------------------------------------------------------------------
	END FUNCTION ComputeError
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END MODULE NumaAdaptive
!-----------------------------------------------------------------------------
