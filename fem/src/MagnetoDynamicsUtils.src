!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/
    
!------------------------------------------------------------------------------
!>  Solve Maxwell equations in vector potential formulation (or the A-V
!>  formulation) and (relatively)low frequency approximation using lowest
!>  order Withney 1-forms (edge elements).
!> \ingroup Solvers
!-------------------------------------------------------------------------------
MODULE MagnetoDynamicsUtils

   USE DefUtils

   COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)

   INTERFACE SetDOFtoValue
     MODULE PROCEDURE SetDOFtoValueR, SetDOFtoValueC
   END INTERFACE

   INTERFACE GetReluctivity
     MODULE PROCEDURE GetReluctivityR, GetReluctivityC
   END INTERFACE

CONTAINS

!------------------------------------------------------------------------------
  FUNCTION GetBoundaryEdgeIndex(Boundary,nedge) RESULT(n)
!------------------------------------------------------------------------------
    INTEGER :: n,nedge
    TYPE(Element_t) :: Boundary
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,jb1,jb2,je1,je2
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Parent, Edge, Face
!------------------------------------------------------------------------------
    Mesh => GetMesh()
    n = 0
    SELECT CASE(GetElementFamily(Boundary))
    CASE(1)
      RETURN
    CASE(2)
      IF ( nedge==1 ) THEN
        Parent => Boundary % BoundaryInfo % Left
        IF ( .NOT. ASSOCIATED(Parent) ) &
            Parent => Boundary % BoundaryInfo % Right
 
        jb1 = Boundary % NodeIndexes(1)
        jb2 = Boundary % NodeIndexes(2)
        DO i=1,Parent % TYPE % NumberOfEdges
          Edge => Mesh % Edges(Parent % EdgeIndexes(i))
          je1 = Edge % NodeIndexes(1)
          je2 = Edge % NodeIndexes(2)
          IF ( jb1==je1.AND.jb2==je2 .OR. jb1==je2.AND.jb2==je1) EXIT
        END DO
        n = Parent % EdgeIndexes(i)
      END IF
    CASE(3,4)
      j = GetBoundaryFaceIndex(Boundary)
      Face => Mesh % Faces(j)
      IF ( nedge>0.AND.nedge<=Face % TYPE % NumberOfEdges ) &
        n = Face % EdgeIndexes(nedge) 
    END SELECT
!------------------------------------------------------------------------------
  END FUNCTION GetBoundaryEdgeIndex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetBoundaryFaceIndex(Boundary) RESULT(n)
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t) :: Boundary
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,m
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Parent, Face
!------------------------------------------------------------------------------
    Mesh => GetMesh()
    Parent => Boundary % BoundaryInfo % Left
    IF ( .NOT. ASSOCIATED(Parent) ) &
       Parent => Boundary % BoundaryInfo % Right

    DO i=1,Parent % TYPE % NumberOfFaces
      Face => Mesh % Faces(Parent % FaceIndexes(i))
      m = 0
      DO j=1,Face % TYPE % NumberOfNodes
        DO k=1,Boundary % TYPE % NumberOfNodes
          IF ( Face % NodeIndexes(j)==Boundary % NodeIndexes(k)) m=m+1
        END DO
      END DO
      IF ( m==Boundary % TYPE % NumberOfNodes) EXIT
    END DO
    n = Parent % FaceIndexes(i)
!------------------------------------------------------------------------------
  END FUNCTION GetBoundaryFaceIndex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SetDOFToValueR(Solver,k,VALUE)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: VALUE,v
    TYPE(Solver_t) :: Solver
    INTEGER :: n,k
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh

    Mesh => GetMesh(Solver)
    n = Solver % Variable % Perm(k+Mesh % NumberOfNodes)
    A => GetMatrix()
    CALL CRS_SetSymmDirichlet(A,A % RHS,n,VALUE)
!------------------------------------------------------------------------------
  END SUBROUTINE SetDOFToValueR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SetDOFToValueC(Solver,k,VALUE)
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: VALUE
    TYPE(Solver_t) :: Solver
    INTEGER :: n,k
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh

    Mesh => GetMesh(Solver)
    n = Solver % Variable % Perm(k+Mesh % NumberOfNodes)
    A => GetMatrix()
    CALL CRS_SetSymmDirichlet(A,A % RHS,2*(n-1)+1,REAL(VALUE))
    CALL CRS_SetSymmDirichlet(A,A % RHS,2*(n-1)+2,AIMAG(VALUE))
!------------------------------------------------------------------------------
  END SUBROUTINE SetDOFToValueC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivityR(Material,Acoef,n)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(:)
!------------------------------------------------------------------------------
    LOGICAL :: Found, FirstTime = .TRUE.
    REAL(KIND=dp) :: Avacuum

    SAVE Avacuum 

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF(.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF
  
    Acoef(1:n) = GetReal( Material, 'Relative Permeability', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Avacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permeability', Found )
    END IF
    IF ( Found ) THEN
      Acoef(1:n) = 1._dp / Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityR
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivityC(Material,Acoef,n)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    COMPLEX(KIND=dp) :: Acoef(:)
!------------------------------------------------------------------------------
    LOGICAL :: Found, FirstTime = .TRUE.
    REAL(KIND=dp) :: Avacuum

    SAVE Avacuum 

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF(.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF

    Acoef(1:n) = GetReal( Material, 'Relative Permeability', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Avacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permeability', Found )
    END IF
    IF ( Found ) THEN
      Acoef(1:n) = 1._dp / Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found )
      Acoef(1:n) = CMPLX( REAL(Acoef(1:n)), &
         GetReal( Material, 'Reluctivity im', Found ), KIND=dp )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetPermittivity(Material,Acoef,n)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(:)
!------------------------------------------------------------------------------
    LOGICAL :: Found, FirstTime = .TRUE.
    REAL(KIND=dp) :: Pvacuum = 0._dp

    IF ( FirstTime ) THEN
      Pvacuum = GetConstReal( CurrentModel % Constants, &
              'Permittivity of Vacuum', Found )
      FirstTime = .FALSE.
    END IF
    

    Acoef(1:n) = GetReal( Material, 'Relative Permittivity', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Pvacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permittivity', Found )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetPermittivity
!------------------------------------------------------------------------------

END MODULE MagnetoDynamicsUtils

