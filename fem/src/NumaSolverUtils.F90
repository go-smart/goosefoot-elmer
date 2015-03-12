
MODULE NumaSolverUtils

   USE DirectSolve
   USE Multigrid
   USE IterSolve
   USE ElementUtils
   USE TimeIntegrate
   USE ModelDescription
   USE MeshUtils
   USE SolverUtils
   USE ParallelUtils
   USE ParallelEigenSolve

   IMPLICIT NONE

   CHARACTER(LEN=MAX_NAME_LEN), PRIVATE :: NormalTangentialName
   INTEGER, PRIVATE :: NormalTangentialNOFNodes
   INTEGER, POINTER, PRIVATE :: NTelement(:,:)
   LOGICAL, POINTER, PRIVATE :: NTzeroing_done(:,:)
   INTEGER, POINTER, PRIVATE :: BoundaryReorder(:)
   REAL(KIND=dp), POINTER, PRIVATE :: BoundaryNormals(:,:),  &
                                      BoundaryTangent1(:,:), &
                                      BoundaryTangent2(:,:)

   SAVE BoundaryReorder, NormalTangentialNOFNodes, BoundaryNormals, &
              BoundaryTangent1, BoundaryTangent2, NormalTangentialName

CONTAINS


!------------------------------------------------------------------------------
   SUBROUTINE NumaSetPeriodicBoundariesPass1( Model, StiffMatrix, ForceVector, &
                      Name, DOF, NDOFs, Perm, This, Done )
!------------------------------------------------------------------------------
!******************************************************************************
!
! Set dirichlet boundary condition for given dof
!
! TYPE(Model_t) :: Model
!   INPUT: the current model structure
!
! TYPE(Matrix_t), POINTER :: StiffMatrix
!   INOUT: The global matrix
!
! REAL(KIND=dp) :: ForceVector(:)
!   INOUT: The global RHS vector
! 
! CHARACTER(LEN=*) :: Name
!   INPUT: name of the dof to be set
!
! INTEGER :: DOF, NDOFs
!   INPUT: The order number of the dof and the total number of DOFs for
!          this equation
!
! INTEGER :: Perm(:)
!   INPUT: The node reordering info, this has been generated at the
!          beginning of the simulation for bandwidth optimization
!******************************************************************************
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Matrix_t), POINTER :: StiffMatrix

    REAL(KIND=dp) :: ForceVector(:)

    CHARACTER(LEN=*) :: Name
    LOGICAL :: Done(:)
    INTEGER :: This, DOF, NDOFs, Perm(:)
!------------------------------------------------------------------------------

    INTEGER :: i,j,k,l,m,n,nn,ii,nlen
    LOGICAL :: GotIt
    REAL(KIND=dp) :: Scale
    TYPE(Matrix_t), POINTER :: Projector, Projector1
!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)

    Scale = -1.0d0
    IF ( .NOT. ListGetLogical( Model % BCs(This) % Values, &
       'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
       IF ( .NOT. ListGetLogical( Model % BCs(This) % Values, &
          'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) RETURN
       Scale = 1.0d0
    END IF

    Projector => Model % BCs(This) % PMatrix
    IF ( .NOT. ASSOCIATED(Projector) ) RETURN

!
!   Do the assembly of the projector:
!   ---------------------------------
!
    DO i=1,Projector % NumberOfRows
       ii = Projector % InvPerm(i)
       k = Perm(ii)
       IF ( .NOT. Done(ii) .AND. k > 0 ) THEN
          k = NDOFs * (k-1) + DOF
          DO l=Projector % Rows(i),Projector % Rows(i+1)-1
            IF ( Projector % Cols(l) <= 0 .OR. Projector % Values(l)==0.0d0 ) CYCLE

            m = Perm(Projector % Cols(l))
            IF ( m > 0 ) THEN
              m = NDOFs*(m-1) + DOF
              DO nn=StiffMatrix % Rows(k),StiffMatrix % Rows(k+1)-1
                 CALL AddToMatrixElement( StiffMatrix, m, StiffMatrix % Cols(nn), &
                        Projector % Values(l) * StiffMatrix % Values(nn) )
              END DO
              ForceVector(m) = ForceVector(m)+Projector % Values(l)*ForceVector(k)
            END IF
          END DO
       END IF
       Done(ii) = .TRUE.
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE NumaSetPeriodicBoundariesPass1
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE NumaSetPeriodicBoundariesPass2( Model, StiffMatrix, ForceVector, &
                      Name, DOF, NDOFs, Perm, This, Done )
!------------------------------------------------------------------------------
!******************************************************************************
!
! Set dirichlet boundary condition for given dof
!
! TYPE(Model_t) :: Model
!   INPUT: the current model structure
!
! TYPE(Matrix_t), POINTER :: StiffMatrix
!   INOUT: The global matrix
!
! REAL(KIND=dp) :: ForceVector(:)
!   INOUT: The global RHS vector
! 
! CHARACTER(LEN=*) :: Name
!   INPUT: name of the dof to be set
!
! INTEGER :: DOF, NDOFs
!   INPUT: The order number of the dof and the total number of DOFs for
!          this equation
!
! INTEGER :: Perm(:)
!   INPUT: The node reordering info, this has been generated at the
!          beginning of the simulation for bandwidth optimization
!******************************************************************************
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Matrix_t), POINTER :: StiffMatrix

    REAL(KIND=dp) :: ForceVector(:)

    CHARACTER(LEN=*) :: Name
    LOGICAL :: Done(:)
    INTEGER :: This, DOF, NDOFs, Perm(:)
!------------------------------------------------------------------------------

    INTEGER :: i,j,k,l,m,n,nn,ii,nlen
    LOGICAL :: GotIt
    REAL(KIND=dp) :: Scale,s
    TYPE(Matrix_t), POINTER :: Projector
!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)

    Scale = -1.0d0
    IF ( .NOT. ListGetLogical( Model % BCs(This) % Values, &
       'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
       IF ( .NOT. ListGetLogical( Model % BCs(This) % Values, &
          'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) RETURN
       Scale = 1.0d0
    END IF

    Projector => Model % BCs(This) % PMatrix
    IF ( .NOT. ASSOCIATED(Projector) ) RETURN
!
!   Do the assembly of the projector:
!   ---------------------------------
    DO i=1,Projector % NumberOfRows
       ii = Projector % InvPerm(i)
       k = Perm( ii )
       IF ( .NOT. Done(ii) .AND. k > 0 ) THEN
          k = NDOFs * (k-1) + DOF
          CALL ZeroRow( StiffMatrix,k )
          DO l = Projector % Rows(i), Projector % Rows(i+1)-1
             IF ( Projector % Cols(l) <= 0 ) CYCLE
             m = Perm( Projector % Cols(l) )
             IF ( m > 0 ) THEN
               m = NDOFs * (m-1) + DOF
               CALL AddToMatrixElement( StiffMatrix,k,m,Projector % Values(l) )
             END IF
          END DO
          ForceVector(k) = 0.0d0
          CALL AddToMatrixElement( StiffMatrix, k, k, scale )
       END IF
       Done(ii) = .TRUE.
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE NumaSetPeriodicBoundariesPass2
!------------------------------------------------------------------------------

END MODULE NumaSolverUtils
