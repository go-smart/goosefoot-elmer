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
! *  Module containing a solver for heat equation
! *
! ******************************************************************************
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
!> Subroutine for solving the energy a.k.a. heat equation in various coordinate systems.
!> \ingroup Solvers
!------------------------------------------------------------------------------

   SUBROUTINE HeatSolver_init( Model,Solver,Timestep,TransientSimulation )

          USE DefUtils

          IMPLICIT NONE
          TYPE(Solver_t) :: Solver  
          TYPE(Model_t) :: Model    
          REAL(KIND=dp) :: Timestep
          LOGICAL :: TransientSimulation 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
          TYPE(ValueList_t),POINTER :: SolverParams
          LOGICAL :: Found

          SolverParams => GetSolverParams()
   END SUBROUTINE HeatSolver_init

   SUBROUTINE HeatSolver( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
    USE DiffuseConvective
    USE DiffuseConvectiveGeneral
    USE Differentials
    USE Radiation
    USE MaterialModels
    USE NumaAdaptive
    USE DefUtils

!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: PHASE_SPATIAL_1 = 1
    INTEGER, PARAMETER :: PHASE_SPATIAL_2 = 2
    INTEGER, PARAMETER :: PHASE_TEMPORAL  = 3
    
    TYPE(Model_t)  :: Model
    TYPE(Solver_t), TARGET :: Solver

    LOGICAL :: TransientSimulation
    REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: StiffMatrix

    INTEGER :: i,j,k,l,m,n,nd,t,tt,iter,k1,k2,body_id,eq_id,istat,LocalNodes,bf_id, &
        ierr, PowerControl_Integ_Length

    TYPE(Nodes_t)   :: ElementNodes
    TYPE(Element_t),POINTER :: Element,Parent,RadiationElement
    TYPE(Variable_t), POINTER :: TimeVar
    REAL(KIND=dp) :: RelativeChange, &
           Norm,PrevNorm,Text,S,C,Emissivity,StefanBoltzmann, &
           ReferencePressure=0.0d0, SpecificHeatRatio, &
           PowerControl_Kp, PowerControl_Ki, PowerControl_Kd, &
           Time, PrevTime, TargetTemperature, tmpPower, &
           TotalDistanceToTip, x, y, z, xProjTip, yProjTip, zProjTip

    CHARACTER(LEN=MAX_NAME_LEN) :: RadiationFlag,ConvectionFlag

    INTEGER :: PhaseChangeModel
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseModel, StabilizeFlag, VarName

    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Stabilize = .TRUE., Bubbles = .TRUE., UseBubbles,NewtonLinearization = .FALSE., &
         Found, GotIt, HeatFluxBC, HeatGapBC, GotMeltPoint, IsRadiation, InfBC, &
         CellsDeath = .FALSE.,ControlMaxTemperature = .FALSE., VisualizePerfusion = .FALSE., &
         ControlTotalPower = .FALSE., InterpolatedElectricPower = .FALSE., UseElectricPower = .FALSE., &
         TemperatureControlledPower = .FALSE., VisualizeHeatSource = .FALSE.
    INTEGER :: MaxCouplingIterDone,CouplingIter, ADOFs, JHDOFs, ios, CellStateModel, NbElectricTips, MaxNbElectricPoints
    INTEGER, ALLOCATABLE :: nbElectricpoints(:)

    CHARACTER(LEN=MAX_NAME_LEN) :: char_MyPe, MaxTemperatureFilename,iterationsFilename, &
        TestName, char_ElectricTip, ElectricTipFile,TotalPowerFilename, line, str1, str2

! Which compressibility model is used
    CHARACTER(LEN=MAX_NAME_LEN) :: CompressibilityFlag, ConvectionField
    INTEGER :: CompressibilityModel

    LOGICAL :: AllocationsDone = .FALSE.,PhaseSpatial=.FALSE., &
        PhaseChange=.FALSE., CheckLatentHeatRelease=.FALSE., FirstTime, &
        ElectricPowerCutOff = .FALSE., SmartHeaterControl, IntegralHeaterControl, HeaterControlLocal, SmartTolReached=.FALSE., &
        TransientHeaterControl, SmartHeaterAverage, ConstantBulk, SaveBulk, &
	TransientAssembly, Converged
    LOGICAL, POINTER :: SmartHeaters(:), IntegralHeaters(:)

    TYPE(Variable_t), POINTER :: TempSol,FlowSol,HeatSol,CurrentSol, MeshSol, DensitySol, &
        CellStateSol
    TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants

    INTEGER, POINTER :: TempPerm(:),FlowPerm(:),CurrentPerm(:),MeshPerm(:), &
        CellStatePerm(:)

    INTEGER :: NSDOFs,NewtonIter,NonlinearIter,MDOFs, &
         SmartHeaterBC, SmartHeaterNode, DoneTime=0
    REAL(KIND=dp) :: NonlinearTol,NewtonTol,SmartTol,Relax, &
            SaveRelax,dt,dt0,CumulativeTime, VisibleFraction, PowerScaling=1.0, PrevPowerScaling=1.0, &
            allocrealtime,alloctime,srealt,arealt, TotalElectricTipsMesure, ElectricPowerEpsilon, DistanceToTip, &
            BodyTemperature, DeadThreshold, TotalPower, &
            PowerRelax, PowerTimeScale, PowerSensitivity, xave, yave, Normal(3), &
            dist, mindist, ControlPoint(3), &
            DeathCapacity, MaxTemperature

    REAL(KIND=dp), POINTER :: Temperature(:),PrevTemperature(:),FlowSolution(:), &
        ElectricCurrent(:), PhaseChangeIntervals(:,:),ForceVector(:), &
        PrevSolution(:), HC(:), Hwrk(:,:,:),MeshVelocity(:), XX(:), YY(:),ForceHeater(:),&
        RealWork(:,:), &
        HistoryErrorToTargetTemperatureVariable(:), HistoryErrorToTargetTemperature(:,:), &
        IntegErrorToTargetTemperature(:),DerivErrorToTargetTemperature(:), &
        HeatSource(:), ElectricPowerVisualization(:), ElectricPower(:), &
        CellState(:), C0(:), C1(:)

    REAL(KIND=dp), ALLOCATABLE :: vals(:)
    REAL(KIND=dp) :: Jx,Jy,Jz,JAbs, MeltPoint, IntHeatSource

    INTEGER, ALLOCATABLE, SAVE :: Indexes(:), SaveIndexes(:)

    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
        STIFF(:,:), LOAD(:), HeatConductivity(:,:,:), &
        FORCE(:), U(:), V(:), W(:), MU(:,:),TimeForce(:), &
        Density(:), LatentHeat(:), HeatTransferCoeff(:), &
        HeatCapacity(:), Enthalpy(:), Viscosity(:), LocalTemperature(:), &
        NodalEmissivity(:), ElectricConductivity(:), Permeability(:), Work(:), &
        Pressure(:), dPressuredt(:), GasConstant(:),AText(:), HeaterArea(:), &
        HeaterTarget(:), HeaterScaling(:), HeaterDensity(:), HeaterSource(:), &
        HeatExpansionCoeff(:), ReferenceTemperature(:), PressureCoeff(:), &
        PhaseVelocity(:,:), HeatConductivityIso(:), &
        PerfusionRate(:), PerfusionDensity(:), PerfusionHeatCapacity(:), PerfusionRefTemperature(:), &
        Integ_Force(:), power(:), &
        x_ElectricTip(:,:),y_ElectricTip(:,:),z_ElectricTip(:,:),ElectricTipMesure(:), &
        LocalHeatSource(:), &
        PreviousPower(:),CurrentPower(:), LocalCellState(:)

    SAVE U, V, W, MU, MASS, STIFF, LOAD, PressureCoeff, &
        FORCE, ElementNodes, HeatConductivity, HeatCapacity, HeatTransferCoeff, &
        Enthalpy, Density, LatentHeat, PhaseVelocity, AllocationsDone, Viscosity, TimeForce, &
        LocalNodes, LocalTemperature, Work, ElectricConductivity, &
        NodalEmissivity, Permeability, C0, C1, dPressuredt, Pressure, &
        GasConstant,AText,Hwrk, XX, YY, ForceHeater, HeaterArea, HeaterTarget, &
        HeaterScaling, HeaterDensity, HeaterSource, SmartHeaters, IntegralHeaters, SmartTolReached,    &
        ReferenceTemperature, HeatExpansionCoeff, PrevPowerScaling, PowerScaling, &
        MeltPoint, DoneTime, SmartHeaterNode, SmartHeaterBC, SmartHeaterAverage, &
        HeatConductivityIso, TransientAssembly, &
        PerfusionRate, PerfusionDensity, PerfusionHeatCapacity, PerfusionRefTemperature, &
        DeathCapacity, Integ_Force, power, &
        x_ElectricTip,y_ElectricTip,z_ElectricTip,ElectricTipMesure, &
        LocalHeatSource,TotalPower, &
        IntegErrorToTargetTemperature,DerivErrorToTargetTemperature,HistoryErrorToTargetTemperature, &
        PreviousPower,CurrentPower, LocalCellState,ElectricPower,HeatSource, &

        MaxTemperature, &
        Time,PrevTime,MaxCouplingIterDone,CouplingIter, &
        TestName,CellStateModel, &
        ElectricPowerCutOff,BodyTemperature, &
        DeadThreshold,Stabilize,UseBubbles,NonlinearIter,NonlinearTol,Relax, &
        CellsDeath,ControlMaxTemperature, UseElectricPower, &
        InterpolatedElectricPower,VisualizeHeatSource,ControlTotalPower,&
        VisualizePerfusion, TemperatureControlledPower, &
        PowerControl_Kp, PowerControl_Ki, PowerControl_Kd,TargetTemperature,&
        PowerControl_Integ_Length


    INTERFACE
        FUNCTION HeatBoundaryResidual( Model,Edge,Mesh,Quant,Perm,Gnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
          INTEGER :: Perm(:)
        END FUNCTION HeatBoundaryResidual

        FUNCTION HeatEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2)
          INTEGER :: Perm(:)
        END FUNCTION HeatEdgeResidual

        FUNCTION HeatInsideResidual( Model,Element,Mesh,Quant,Perm, Fnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Element
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
          INTEGER :: Perm(:)
        END FUNCTION HeatInsideResidual
    END INTERFACE

    REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime


!------------------------------------------------------------------------------
!    The View and Gebhardt factors may change. If this is necessary, this is 
!    done within this subroutine. The routine is called in the
!    start as it may affect the matrix topology.
!    Newton lineariarization option is needed only when there is radiation.
!------------------------------------------------------------------------------
    IsRadiation = ListCheckPresentAnyBC( Model,'Radiation')
    IF( IsRadiation ) THEN
        CALL RadiationFactors( Solver, .FALSE.)
    END IF

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

    IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

    StiffMatrix => GetMatrix()
    ForceVector => Solver % Matrix % RHS

    TempSol => Solver % Variable
    TempPerm    => TempSol % Perm
    Temperature => TempSol % Values
    VarName = GetVarName( TempSol ) 
  
    CellStateSol => VariableGet( Solver % Mesh % Variables, 'CellState' )
    IF ( ASSOCIATED( CellStateSol ) ) THEN
         CellStatePerm => CellStateSol % Perm
         CellState => CellStateSol % Values
         ADOFs =  CellStateSol % DOFs
    END IF

    LocalNodes = COUNT( TempPerm > 0 )
    IF ( LocalNodes <= 0 ) RETURN

    SolverParams => GetSolverParams()
    ConvectionField = GetString( SolverParams, 'Temperature Convection Field', Found )

    IF ( Found ) THEN
        FlowSol => VariableGet( Solver % Mesh % Variables, ConvectionField )
    ELSE
        FlowSol => VariableGet( Solver % Mesh % Variables, 'Flow Solution' )
    END IF

    IF ( ASSOCIATED( FlowSol ) ) THEN
        FlowPerm     => FlowSol % Perm
        NSDOFs       =  FlowSol % DOFs
        FlowSolution => FlowSol % Values
    END IF

    DensitySol => VariableGet( Solver % Mesh % Variables, 'Density' )

    ! Check whether we have some heater controls. This will affect initialization stuff. 
    SmartHeaterControl = ListCheckPresentAnyBodyForce( Model,'Smart Heater Control')
    IntegralHeaterControl = ListCheckPresentAnyBodyForce( Model,'Integral Heat Source')
   
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
    IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
        N = Solver % Mesh % MaxElementDOFs

        IF ( AllocationsDone ) THEN
          DEALLOCATE(  &
                 U, V, W, MU,           &
                 Pressure,              &
                 dPressureDt,           &
                 PressureCoeff,        &
                 ElementNodes % x,      &
                 ElementNodes % y,      &
                 ElementNodes % z,      &
                 Density,Work,          &
                 LatentHeat,            &
                 PhaseVelocity,         &
                 ElectricConductivity,  &
                 Permeability,          &
                 Viscosity,C0,C1,          &
                 HeatTransferCoeff,     &
                 HeatExpansionCoeff,    &
                 ReferenceTemperature,  &
                 MASS,       &
                 LocalTemperature,      &
                 HeatCapacity,Enthalpy, &
                 NodalEmissivity,       &
                 GasConstant, AText,    &
                 HeatConductivity,      &
                 STIFF,LOAD,            &
                 Indexes, SaveIndexes,  &
                 FORCE, TimeForce,      &
                 HeatConductivityIso,   &
                 PerfusionRate,         &
                 PerfusionDensity,      &
                 PerfusionHeatCapacity, &
                 PerfusionRefTemperature, &
                 power,                                          &
                 LocalCellState,        &
                 Integ_Force,                                &
                 ElectricPower,                          &
                 HeatSource,                  &
                 LocalHeatSource,                        &
                 PreviousPower,                          &
                 CurrentPower                                    )
        END IF

        ALLOCATE( &
                 Indexes(N), SaveIndexes(N),           &
                 U( N ),   V( N ),  W( N ),            &
                 MU( 3,N ),                            &
                 Pressure( N ),                        &
                 dPressureDt( N ),                     &
                 PressureCoeff( N ),                   &
                 ElementNodes % x( N ),                &
                 ElementNodes % y( N ),                &
                 ElementNodes % z( N ),                &
                 Density( N ),Work( N ),               &
                 LatentHeat( N ),                      &
                 PhaseVelocity(3, N),                  &
                 ElectricConductivity(N),              &
                 Permeability(N),                      &
                 Viscosity(N),C0(N),C1(N),                   &
                 HeatTransferCoeff( N ),               &
                 HeatExpansionCoeff( N ),              &
                 ReferenceTemperature( N ),            &
                 MASS(  2*N,2*N ),                     &
                 LocalTemperature( N ),                &
                 HeatCapacity( N ),Enthalpy( N ),      &
                 NodalEmissivity( N ),                 &
                 GasConstant( N ),AText( N ),          &
                 HeatConductivity( 3,3,N ),            &
                 HeatConductivityIso( N ),             &
                 STIFF( 2*N,2*N ),LOAD( N ), &
                 FORCE( 2*N ), TimeForce(2*N), &
                 PerfusionRate( N ),     &
                 PerfusionDensity( N ),      &
                 PerfusionHeatCapacity( N ), &
                 PerfusionRefTemperature( N ), &
                 power( N ),                     &  
                 LocalCellState(N),              &
                 ElectricPower(LocalNodes),      &
                 HeatSource(LocalNodes),         &
                 LocalHeatSource(N),             &
                 PreviousPower(LocalNodes),          &
                 CurrentPower(LocalNodes),         &
                 STAT=istat )

        IF ( istat /= 0 ) THEN
         CALL Fatal( 'HeatSolve', 'Memory allocation error' )
        END IF


        IF ( SmartHeaterControl .OR. IntegralHeaterControl) THEN          
          n = Model % NumberOfBodyForces
          IF ( AllocationsDone ) DEALLOCATE( HeaterArea, HeaterDensity, HeaterSource, &
               HeaterScaling, HeaterTarget, SmartHeaters, IntegralHeaters )
          ALLOCATE( HeaterArea(n), HeaterDensity(n), HeaterSource(n), &
               HeaterScaling(n), HeaterTarget(n), SmartHeaters(n), &
               IntegralHeaters(n) )
          IF ( istat /= 0 ) THEN
             CALL Fatal( 'HeatSolve', 'Memory allocation error' )
          END IF
          SmartHeaters = .FALSE.
          IntegralHeaters = .FALSE.
        END IF
       
        IF( SmartHeaterControl ) THEN
          IF ( AllocationsDone ) DEALLOCATE( XX, YY, ForceHeater  )
          n = SIZE( Temperature )
          ALLOCATE( XX( n ), YY(n), ForceHeater( n ), STAT=istat )
          IF ( istat /= 0 ) THEN
             CALL Fatal( 'HeatSolve', 'Memory allocation error' )
          END IF
          XX = 0.0d0 
          YY = 0.0d0
          ForceHeater = 0.0d0
        END IF
       
        NULLIFY( Hwrk )
        AllocationsDone = .TRUE.

    !------------------------------------------------------------------------------ 
    !   Check if Cells Death Modelling specified in input file 
    !------------------------------------------------------------------------------ 
        CellsDeath = GetLogical( SolverParams,'Cells Death Modelling',Found )
        IF ( .NOT.Found ) CellsDeath = .TRUE.
    !------------------------------------------------------------------------------
    !       Physical time of the current and previous coupled system iterations: 
    !------------------------------------------------------------------------------             
        Time = 0.0
        PrevTime = 0.0
    !------------------------------------------------------------------------------
    !       If specified in the input file, compute the maximum temperature over the model: 
    !------------------------------------------------------------------------------ 
        ControlMaxTemperature = GetLogical( SolverParams,'Control Max Temperature',Found )      
        IF ( .NOT.Found ) ControlMaxTemperature = .FALSE.
        IF (ControlMaxTemperature) THEN
            MaxTemperature = 0.0D0
            DO i=1,LocalNodes
                    IF (MaxTemperature < Temperature(i)) THEN
                        MaxTemperature = Temperature(i)
                    END IF
            END DO
    !------------------------------------------------------------------------------ 
    !       Write the header of the max-temperature control file
    !------------------------------------------------------------------------------     
            IF(ParEnv % PEs>1) THEN
                WRITE(char_MyPe,*) ParEnv % MyPe
                char_MyPe = ADJUSTL(char_MyPe)
                MaxTemperatureFilename = 'maxtemp'//'_'//TRIM(TestName)//'_'//TRIM(char_MyPe)//'.dat'       
            ELSE 
                MaxTemperatureFilename = 'maxtemp'//'_'//TRIM(TestName)//'.dat'         
            END IF
            
            OPEN(UNIT=1,FILE=MaxTemperatureFilename)
            WRITE(UNIT=1,FMT=*) 'Time    ', '    Blood Temperature    ', '    Tissue Temperature'
            CLOSE(1)
    !------------------------------------------------------------------------------
        END IF ! ControlMaxTemperature
    !------------------------------------------------------------------------------ 
    !       If saving the total electric power in a control file, write header:
    !------------------------------------------------------------------------------ 
        TotalPower = 0.0
        ControlTotalPower = GetLogical( SolverParams,'Control Total Power',Found )      
        IF ( .NOT.Found ) ControlTotalPower = .TRUE.
        IF ((ControlTotalPower) .AND. (ParEnv % MyPe==0)) THEN  
            TotalPowerFilename = 'totalpower'//'_'//TRIM(TestName)//'.csv'          
            OPEN(UNIT=1,FILE=TotalPowerFilename)
            WRITE(UNIT=1,FMT=*) 'Time,', 'Blood-Total-Power,', 'Tissue-Total-Power'
            CLOSE(1)
        END IF ! ControlTotalPower      
    !------------------------------------------------------------------------------ 
    !       Get electric power geometry if interpolated (multi)line source:
    !------------------------------------------------------------------------------ 
        InterpolatedElectricPower = GetLogical( Model % Simulation, &
            'Interpolated Electric Power',Found )
        IF (.NOT. Found) InterpolatedElectricPower = .FALSE.

        UseElectricPower = GetLogical( SolverParams, &
            'Use Electric Power',Found )
        IF (.NOT. Found) UseElectricPower = InterpolatedElectricPower

    !------------------------------------------------------------------------------ 
        IF (InterpolatedElectricPower) THEN
    !------------------------------------------------------------------------------
    !           Get the number of tips:
    !------------------------------------------------------------------------------ 
            NbElectricTips = GetInteger( Model % Simulation, &
                'Electric Tips Number', Found )
            IF ( .NOT.Found ) NbElectricTips = 10
    !------------------------------------------------------------------------------
    !           Get the (numerical) width of the tips:
    !------------------------------------------------------------------------------ 
            ElectricPowerEpsilon = GetConstReal( Model % Simulation, &
                'Electric Power Epsilon', Found )
            IF ( .NOT.Found ) ElectricPowerEpsilon = 1.5
    !------------------------------------------------------------------------------
    !           Read the coordinates of the points defining the tips in text files
    !------------------------------------------------------------------------------
            MaxNbElectricPoints = 10        
            ALLOCATE(x_ElectricTip(NbElectricTips,MaxNbElectricPoints), &
                y_ElectricTip(NbElectricTips,MaxNbElectricPoints), &
                z_ElectricTip(NbElectricTips,MaxNbElectricPoints), &
                ElectricTipMesure(NbElectricTips), &
                nbElectricpoints(NbElectricTips)) 

            x_ElectricTip = 0.0D0
            y_ElectricTip = 0.0D0
            z_ElectricTip = 0.0D0
            ElectricTipMesure = 0.0D0
            TotalElectricTipsMesure = 0.0D0
            NbElectricPoints = 0
            !------------------------------------------------------------------------------
            !   Go through the tips
            !------------------------------------------------------------------------------     
            DO j=1,NbElectricTips
            !------------------------------------------------------------------------------     
            !   Open the Tips file:
            !------------------------------------------------------------------------------  
                ElectricTipFile = GetString( Model % Simulation,'Electric Tips Filename Root',Found )   
                !------------------------------------------------------------------------------ 
                IF (Found) THEN
                !------------------------------------------------------------------------------ 
                    WRITE(char_ElectricTip,*) j
                    char_ElectricTip = ADJUSTL(char_ElectricTip)
                    ElectricTipFile = TRIM(ElectricTipFile) // "_" // TRIM(char_ElectricTip) // ".txt"
                    OPEN(UNIT=10,FILE=ElectricTipFile,STATUS='OLD',ACTION='READ',IOSTAT=ios)
                    !------------------------------------------------------------------------------ 
                    IF(ios/=0) THEN
                        PRINT*,'Could not open file ',ElectricTipFile
                        PRINT*,'I/O Fortran Error number ',ios
                        CALL Fatal( 'NumaHeatSolve', 'Unable to load electric tips file' )
                    ELSE
                        !------------------------------------------------------------------------------ 
                        !   Read the number of points defining the tip geometry:
                        !------------------------------------------------------------------------------ 
                        READ(10,*,END=1) line
                        READ(10,*,END=1) str1, NbElectricPoints(j), str2
                        !------------------------------------------------------------------------------
                        DO i=1,NbElectricPoints(j)  
                            !------------------------------------------------------------------------------
                            !       Read the coordinates:
                            !------------------------------------------------------------------------------
                            READ(10,*,END=1) x_ElectricTip(j,i), y_ElectricTip(j,i), z_ElectricTip(j,i)
                        !------------------------------------------------------------------------------
                        END DO
                        !------------------------------------------------------------------------------
                        1 CONTINUE
                        CLOSE(10)
                    !------------------------------------------------------------------------------ 
                    END IF
                !------------------------------------------------------------------------------ 
                ELSE
                !------------------------------------------------------------------------------
                !   If the file can't be found, print an error message and stop the simulation: 
                !------------------------------------------------------------------------------
                    CALL Info('NumaHeatSolve', &
                        'Please specify electric tips file name root in input file.', Level=1 )
                    CALL Fatal( 'NumaHeatSolve', 'Unable to load electric tips file' )
                !------------------------------------------------------------------------------
                END IF ! name of the tip file found
                !------------------------------------------------------------------------------
                ! Compute the length of the tip:
                !------------------------------------------------------------------------------
                ElectricTipMesure = 0.0D0
                !------------------------------------------------------------------------------
                !   Case of point source
                !------------------------------------------------------------------------------
                IF(NbElectricPoints(j)==1) THEN
                !------------------------------------------------------------------------------
                    ElectricTipMesure = 1.0D0
                !------------------------------------------------------------------------------
                ELSE
                !------------------------------------------------------------------------------
                    DO i=1,NbElectricPoints(j)-1    
                    !------------------------------------------------------------------------------
                        ElectricTipMesure(j) = ElectricTipMesure(j) + sqrt( (x_ElectricTip(j,i+1)-x_ElectricTip(j,i))**2 +  &
                            (y_ElectricTip(j,i+1)-y_ElectricTip(j,i))**2 + &
                            (z_ElectricTip(j,i+1)-z_ElectricTip(j,i))**2 )
                    !------------------------------------------------------------------------------
                    END DO
                !------------------------------------------------------------------------------
                END IF
                !------------------------------------------------------------------------------
                ! Update the total mesure of the electric source
                !------------------------------------------------------------------------------
                TotalElectricTipsMesure = TotalElectricTipsMesure + ElectricTipMesure(j)
            !------------------------------------------------------------------------------
            END DO !j
            !------------------------------------------------------------------------------
            !   Compute the electric power distribution:
            !------------------------------------------------------------------------------
            ElectricPower = 0.0D0
            !------------------------------------------------------------------------------
            ! Go through the nodes
            !------------------------------------------------------------------------------ 
            DO i=1,LocalNodes
            !------------------------------------------------------------------------------ 
                x = Model % Nodes % x(i)
                y = Model % Nodes % y(i)
                z = Model % Nodes % z(i)
                TotalDistanceToTip = 1000000
                !------------------------------------------------------------------------------ 
                ! Go through the tips
                !------------------------------------------------------------------------------ 
                DO j=1,NbElectricTips
                    !------------------------------------------------------------------------------ 
                    ! Compute the distance to the tip
                    !------------------------------------------------------------------------------ 
                    DistanceToTip = NumaDistanceToElectricTip (x, y, z, x_ElectricTip(j,:), &
                        y_ElectricTip(j,:), z_ElectricTip(j,:), nbElectricpoints(j), &
                        xProjTip, yProjTip, zProjTip )
                    !------------------------------------------------------------------------------ 
                    ! The electric power at each node comes from the closest tip
                    !------------------------------------------------------------------------------
                    IF (TotalDistanceToTip>DistanceToTip) THEN
                    !------------------------------------------------------------------------------
                            !------------------------------------------------------------------------------
                            ! Test if we are in the neighboorhoof of the tip
                            !------------------------------------------------------------------------------
                            IF (DistanceToTip<ElectricPowerEpsilon) THEN
                                !------------------------------------------------------------------------------ 
                                ! If source = one point
                                !------------------------------------------------------------------------------ 
                                IF (nbElectricpoints(j)==1) THEN
                                    !------------------------------------------------------------------------------ 
                                    !   Constant power distibution
                                    !------------------------------------------------------------------------------ 
                                    ElectricPower(TempPerm(i)) = &
                                            (1+0.0*COS(Pi*DistanceToTip/ElectricPowerEpsilon))/&
                                            (4.0*Pi/3.0*ElectricPowerEpsilon**3)
                                    !------------------------------------------------------------------------------ 
                                    !   Smoothed power distibution
                                    !------------------------------------------------------------------------------ 
                                    ! TO BE DONE
                                !------------------------------------------------------------------------------ 
                                ELSE
                                    !------------------------------------------------------------------------------ 
                                    ! If source = closed line
                                    !------------------------------------------------------------------------------ 
                                    IF ((x_ElectricTip(j,1)==x_ElectricTip(j,nbElectricpoints(j))) .AND. &
                                    (y_ElectricTip(j,1)==y_ElectricTip(j,nbElectricpoints(j))) .AND. &
                                    (z_ElectricTip(j,1)==z_ElectricTip(j,nbElectricpoints(j)))) THEN
                                        !------------------------------------------------------------------------------ 
                                        !   Constant power distibution
                                        !------------------------------------------------------------------------------ 
                                        ElectricPower(TempPerm(i)) = &
                                            (1+0.0*COS(Pi*DistanceToTip/ElectricPowerEpsilon))/&
                                            (TotalElectricTipsMesure*Pi*ElectricPowerEpsilon**2)
                                        !------------------------------------------------------------------------------ 
                                        !   Smoothed power distibution
                                        !------------------------------------------------------------------------------ 
                                        ! TO BE DONE
                                    !------------------------------------------------------------------------------ 
                                    ELSE
                                        !------------------------------------------------------------------------------ 
                                        !   Constant power distibution
                                        !------------------------------------------------------------------------------     
                                        ElectricPower(TempPerm(i)) = &
                                            (1+0.0*COS(Pi*DistanceToTip/ElectricPowerEpsilon))/&
                                            (TotalElectricTipsMesure*Pi*ElectricPowerEpsilon**2+ &
                                            !------------------------------------------------------------------------------ 
                                            ! 0.5 coeff for particular case of electric probe 
                                            !------------------------------------------------------------------------------ 
                                            0.5*4.0*Pi/3.0*ElectricPowerEpsilon**3)
                                            !------------------------------------------------------------------------------ 
                                            !   Smoothed power distibution
                                            !------------------------------------------------------------------------------ 
                                            ! TO BE DONE
                                    !------------------------------------------------------------------------------ 
                                    END IF ! line type
                                    !------------------------------------------------------------------------------                         
                                END IF ! point or line source
                                !------------------------------------------------------------------------------ 
                            END IF !closest tip
                    !------------------------------------------------------------------------------
                    END IF !neighbourhood of the tip
                    !------------------------------------------------------------------------------ 
                END DO !j
                !------------------------------------------------------------------------------
                ! Add the electric power as a variable for vizualization only   
                !------------------------------------------------------------------------------
                ElectricPowerVisualization => ElectricPower(1:LocalNodes)
                CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
                Solver, 'Electric power', 1, &
                ElectricPowerVisualization, TempPerm )  
                !------------------------------------------------------------------------------     
            END DO !i
            !------------------------------------------------------------------------------ 
        END IF !InterpolatedElectricPower

    !------------------------------------------------------------------------------ 
    !       Read the death model in input file
    !------------------------------------------------------------------------------ 
        CellStateModel = GetInteger( Model % Simulation,'Cell State Model', Found )
        IF ( .NOT.Found ) CellStateModel = 2
        
        
        IF (CellStateModel == 2) THEN
            DeadThreshold = GetConstReal( SolverParams, 'Model 2 Dead Threshold', Found )
        END IF
        IF ( .NOT.Found ) DeadThreshold = 0.8
    !------------------------------------------------------------------------------
    !       Check if electric power has to be cut off when T>373 to approx water evaporation
    !------------------------------------------------------------------------------
        ElectricPowerCutOff = ListGetLogical( SolverParams, 'Electric Power Cut Off ', Found )
        IF ( .NOT.Found ) ElectricPowerCutOff = .TRUE.
    !------------------------------------------------------------------------------ 
    !   Read some solver options in the input file
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    !    Do some additional initialization, and go for it
    !------------------------------------------------------------------------------
    END IF ! not allocations done
!------------------------------------------------------------------------------
        dt = Timestep
        Constants => GetConstants()
        IF( IsRadiation ) THEN
           StefanBoltzmann = ListGetConstReal( Model % Constants, &
                         'Stefan Boltzmann' )
        END IF

!------------------------------------------------------------------------------
        Stabilize = GetLogical( SolverParams,'Stabilize',Found )
        IF ( .NOT.Found ) Stabilize = .TRUE.

        UseBubbles = GetLogical( SolverParams,'Bubbles',Found )
        IF ( .NOT.Found ) UseBubbles = .FALSE.

        StabilizeFlag = GetString( SolverParams, &
              'Stabilization Method',Found )

        SELECT CASE(StabilizeFlag)
        CASE('vms')
           Stabilize = .FALSE.
           UseBubbles= .FALSE.
        CASE('stabilized')
           Stabilize = .TRUE.
           UseBubbles = .FALSE.
        CASE('bubbles')
           Stabilize = .FALSE.
           UseBubbles = .TRUE.
        END SELECT

    !------------------------------------------------------------------------------ 
    !   Maximum number of nonlinear iterations
    !------------------------------------------------------------------------------ 
        NonlinearIter = GetInteger(   SolverParams, &
                         'Nonlinear System Max Iterations', Found )
        IF ( .NOT.Found ) NonlinearIter = 1

    !------------------------------------------------------------------------------ 
    !   Tolerance to be reached during nonlinear iterations
    !------------------------------------------------------------------------------ 
        NonlinearTol  = GetConstReal( SolverParams, &
                         'Nonlinear System Convergence Tolerance',    Found )
        IF ( .NOT.Found ) NonlinearTol = 3.0

        IF( IsRadiation ) THEN
           NewtonTol     = GetConstReal( SolverParams, &
                          'Nonlinear System Newton After Tolerance',  Found )
           NewtonIter    = GetInteger(   SolverParams, &
                          'Nonlinear System Newton After Iterations', Found )
        ELSE
           NewtonTol = 1.0_dp
           NewtonIter =  0
        END IF
        IF ( NewtonIter == 0) NewtonLinearization = .TRUE.

        Relax = GetCReal( SolverParams, &
                   'Nonlinear System Relaxation Factor',Found )
        IF ( .NOT.Found ) Relax = 1

        TransientAssembly = TransientSimulation
        dt0 = ListGetConstReal(SolverParams,'Steady State Transition Timestep',Found)
        IF(.NOT. Found) dt0 = ListGetConstReal(SolverParams,'Smart Heater Time Scale',Found)

        IF(Found .AND. dt > dt0) TransientAssembly = .FALSE.


    !------------------------------------------------------------------------------ 
    !   Read the dead cells criteria and characteristics in the input file
    !------------------------------------------------------------------------------ 
    !------------------------------------------------------------------------------
    !   Get body temperature from input file (for perfusion)
    !------------------------------------------------------------------------------
        BodyTemperature=310.0D0

    !------------------------------------------------------------------------------ 
    !   Check if the electric power has to be controlled in function of temperature 
    !------------------------------------------------------------------------------
        TemperatureControlledPower = GetLogical( Model % Simulation, &
           'Temperature Controlled Electric Power',Found )
        IF ( .NOT.Found ) TemperatureControlledPower = .FALSE.
    !------------------------------------------------------------------------------
        IF(TemperatureControlledPower) THEN
    !------------------------------------------------------------------------------ 
    !       Get the target temperature:
    !------------------------------------------------------------------------------     
            TargetTemperature = GetConstReal( Model % Simulation, &
               'Target Temperature For Electric Power Control', Found )
            IF ( .NOT.Found ) TargetTemperature = 373.15
    !------------------------------------------------------------------------------ 
    !       Get the PID parameters:
    !------------------------------------------------------------------------------     
            PowerControl_Kp = GetConstReal( Model % Simulation, &
                'Proportional Gain For Electric Power Control', Found )
            IF ( .NOT.Found ) PowerControl_Kp = 0.0

            PowerControl_Kd = GetConstReal( Model % Simulation, &
                'Derivative Gain For Electric Power Control', Found )
            IF ( .NOT.Found ) PowerControl_Kd = 0.0

            PowerControl_Ki = GetConstReal( Model % Simulation, &
                'Integral Gain For Electric Power Control', Found )
            IF ( .NOT.Found ) PowerControl_Ki = 0.0

            PowerControl_Integ_Length = GetInteger( Model % Simulation, &
                'Integral Length For Electric Power Control', Found )
            IF ( .NOT.Found ) PowerControl_Integ_Length = 0

            ALLOCATE(HistoryErrorToTargetTemperature(MAX(PowerControl_Integ_Length,2),LocalNodes))

            PRINT *,'PID Temperature Controlled Electric Power, with parameters:'
            PRINT *,'- PID Proportional gain = ',PowerControl_Kp            
            PRINT *,'- PID Derivative gain = ',PowerControl_Kd  
            PRINT *,'- PID Integral gain = ',PowerControl_Ki    
            PRINT *,'- PID Target Temperature (K) = ',TargetTemperature 
            PRINT *,'- PID Integration Length = ',PowerControl_Integ_Length
    !------------------------------------------------------------------------------     
    !           Initialize the corresponding variables:
    !------------------------------------------------------------------------------     
            CurrentPower = 0.0D0
            HistoryErrorToTargetTemperature = 0.0D0
            DO i=1,LocalNodes
                HistoryErrorToTargetTemperature(1,i) = TargetTemperature-Temperature(i)
            END DO
            IntegErrorToTargetTemperature(:) = HistoryErrorToTargetTemperature(1,:)*Timestep
            DerivErrorToTargetTemperature = 0.0D0
    !------------------------------------------------------------------------------     
        END IF !TemperatureControlledPower
    !------------------------------------------------------------------------------ 
    !       Compute and print allocation time:
    !------------------------------------------------------------------------------
        alloctime = CPUTime() - alloctime
        allocrealtime = RealTime() - allocrealtime
        WRITE(Message,'(a,F8.2,F8.2)') 'Allocation time (CPU,REAL): (s)', alloctime,allocrealtime
        CALL Info('NumaHeatSolve',Message,Level=4 )

    TransientHeaterControl = .FALSE.
    IF(SmartHeaterControl) THEN

        ! Mark the smart heaters 
        SmartHeaters = .FALSE.
        bf_id = 0
        DO i = 1,Model % NumberOfBodyForces
         IF( ListGetLogical( Model % BodyForces(i) % Values, &
             'Smart Heater Control', Found ) ) THEN
           SmartHeaters(i) = .TRUE.	     
           bf_id = i
         END IF
        END DO

        ! Find the BC that controls the heater 
        ! If not found assume that smart heater is related to phase change 
        MeltPoint = GetCReal( Model % BodyForces(bf_id) % Values,&
           'Smart Heater Temperature',GotMeltPoint)           
              
        SmartHeaterAverage = .FALSE.
        SmartHeaterNode = ListGetInteger( Model % BodyForces(bf_id) % Values,&
           'Smart Heater Control Node',GotIt) 
        IF(.NOT. GotIt) THEN
         RealWork => ListGetConstRealArray( Model % BodyForces(bf_id) % Values,&
             'Smart Heater Control Point',GotIt) 
         IF( GotIt ) THEN
           ControlPoint(1:3) = RealWork(1:3,1)
           
           mindist = HUGE( mindist )
           DO l=1,Model % NumberOfNodes
             IF( TempPerm(l) == 0 ) CYCLE
             
             jx = Model % Mesh % Nodes % x(l)
             jy = Model % Mesh % Nodes % y(l)
             jz = Model % Mesh % Nodes % z(l)
             
             dist = (ControlPoint(1)-jx)**2 + (ControlPoint(2)-jy)**2 + (ControlPoint(3)-jz)**2
             IF( dist < mindist ) THEN
               mindist = dist
               SmartHeaterNode = l
             END IF
           END DO
         END IF

         WRITE(Message,*) 'Found Control Point at distance:',SQRT(mindist)
         CALL Info('HeatSolve',Message)
         WRITE(Message,*) 'Control Point Index:',SmartHeaterNode
         CALL Info('HeatSolve',Message)        
        END IF
       
        IF( .NOT. GotMeltPoint .OR. SmartHeaterNode == 0) THEN
         GotIt = .FALSE.
         Found = .FALSE.
         SmartHeaterBC = 0
         
         DO i=1,Model % NumberOfBCs
           GotIt = ListGetLogical( Model % BCs(i) % Values,'Smart Heater Boundary', Found ) 
           IF(GotIt) THEN
             SmartHeaterBC = i
             EXIT
           END IF
         END DO
         IF(.NOT. GotIt) THEN
           DO i=1,Model % NumberOfBCs
             GotIt = ListGetLogical( Model % BCs(i) % Values,'Phase Change', Found ) 
             IF(GotIt) THEN
               SmartHeaterBC = i
               EXIT
             END IF
           END DO
         END IF
         IF(SmartHeaterBC == 0) THEN
           CALL Fatal('HeatSolve','Smart Heater Boundary / Phase Change is undefined')
         END IF
         
         MeltPoint = GetCReal( Model % BCs(SmartHeaterBC) % Values,&
             'Smart Heater Temperature',Found)
         IF(.NOT. Found) THEN
           DO k=1, Model % NumberOfMaterials
             MeltPoint = GetCReal( Model % Materials(k) % Values, &
                 'Melting Point', Found )
             IF(Found) EXIT
           END DO
           IF(.NOT. Found) THEN
             CALL Fatal('HeatSolver','Smart Heater Temperature / Melting Point is undefined')
           END IF
         END IF
         
         ! Find the node related to temperature control 
         SmartHeaterAverage = ListGetLogical(Solver % Values,'Smart Heater Average', Found)
         IF(.NOT. SmartHeaterAverage) THEN
           jx = -HUGE(jx)
           DO k = Model % Mesh % NumberOfBulkElements + 1, &
               Model % Mesh % NumberOfBulkElements + Model % Mesh % NumberOfBoundaryElements
             
             Element => Model % Mesh % Elements(k)
             
             IF ( Element % BoundaryInfo % Constraint == SmartHeaterBC ) THEN
               DO l=1,Element % TYPE % NumberOfNodes
                 IF ( Model % Mesh % Nodes % x(Element % NodeIndexes(l)) >= jx ) THEN
                   j = Element % NodeIndexes(l) 
                   jx = Model % Mesh % Nodes % x(Element % NodeIndexes(l))
                 END IF
               END DO
             END IF
           END DO
           SmartHeaterNode = j
         END IF
        END IF

        SmartTol  = GetConstReal( SolverParams, &
             'Smart Heater Control After Tolerance',  Found )
        IF(.NOT. Found) THEN
          SmartTolReached = .TRUE.
          SmartTol = 1.0
        END IF   
     
        PowerTimeScale = ListGetConstReal(Solver % Values, &
             'Smart Heater Time Scale',Found)

        IF(TransientSimulation .AND. dt < PowerTimeScale) THEN
           TransientHeaterControl = .TRUE.
           CALL Info( 'HeatSolve', 'Using Transient Heater Control')
        ELSE
           TransientHeaterControl = .FALSE.
           CALL Info( 'HeatSolve', 'Using Steady-state Heater Control')
        END IF
        
        IF(Solver % DoneTime /= DoneTime) THEN
           PrevPowerScaling = PowerScaling
           DoneTime = Solver % DoneTime
        END IF
    END IF
    
    IF( IntegralHeaterControl) THEN
        CALL Info( 'HeatSolve', 'Using Integral Heater Control')       
        IntegralHeaters = .FALSE.
        DO i = 1,Model % NumberOfBodyForces
           IntegralHeaters(i) = ListCheckPresent( Model % BodyForces(i) % Values, &
                'Integral Heat Source')
        END DO
    END IF

!------------------------------------------------------------------------------

    ConstantBulk = GetLogical( SolverParams, 'Constant Bulk System', Found )
    SaveBulk = ConstantBulk .OR. GetLogical( SolverParams, 'Save Bulk System', Found )
    SaveBulk = ConstantBulk .OR. GetLogical( SolverParams, 'Calculate Loads', Found )

!------------------------------------------------------------------------------

    SaveRelax = Relax
    CumulativeTime = 0.0d0

!------------------------------------------------------------------------------
    FirstTime = .TRUE.
    ALLOCATE(PrevSolution(LocalNodes))
     
    DO WHILE( CumulativeTime < Timestep-1.0d-12 .OR. .NOT. TransientSimulation )
!------------------------------------------------------------------------------
!    The first time around this has been done by the caller...
!------------------------------------------------------------------------------
    IF ( TransientSimulation .AND. .NOT.FirstTime ) &
        CALL InitializeTimestep(Solver)
    FirstTime = .FALSE.
!------------------------------------------------------------------------------
!    Save current solution
!------------------------------------------------------------------------------
    PrevSolution = Temperature(1:LocalNodes)
    IF ( TransientSimulation ) THEN
        PrevTemperature => Solver % Variable % PrevValues(:,1)
    END IF
!------------------------------------------------------------------------------

    totat = 0.0d0
    totst = 0.0d0

    !------------------------------------------------------------------------------ 
    !  Get current (physical) time
    !------------------------------------------------------------------------------
    TimeVar => VariableGet( CurrentModel % Solver % Mesh % Variables, 'Time' )
    Time = TimeVar % Values(1)  
    !------------------------------------------------------------------------------         
    !   Save Max temperature of the previous time in a file 
    !   Control that this is a new global time iteration, i.e. not only the nonlinear and 
    !  coupled system iterations have changed
    !------------------------------------------------------------------------------             
    IF ( (ControlMaxTemperature) .AND. (PrevTime < Time) ) THEN
        IF(ParEnv % PEs>1) THEN
            WRITE(char_MyPe,*) ParEnv % MyPe
            char_MyPe = ADJUSTL(char_MyPe)
            MaxTemperatureFilename = 'maxtemp'//'_'//TRIM(TestName)//'_'//TRIM(char_MyPe)//'.dat'       
        ELSE 
            MaxTemperatureFilename = 'maxtemp'//'_'//TRIM(TestName)//'.dat'         
        END IF
        OPEN(UNIT=1,FILE=MaxTemperatureFilename, FORM='FORMATTED', ACTION='WRITE', POSITION='APPEND', &
            IOSTAT=ios)
        WRITE(UNIT=1,FMT='(F13.4)', ADVANCE='no') PrevTime
        WRITE(UNIT=1,FMT='(F13.4)', ADVANCE='yes') MaxTemperature
        CLOSE(1)
    END IF
    !------------------------------------------------------------------------------         
    !   Save total power of the previous time in a file 
    !   Control that this is a new global time iteration, i.e. not only the nonlinear and 
    !  coupled system iterations have changed
    !  If parralel, only mester proc writes, after getting contribution from other procs
    !------------------------------------------------------------------------------         
    ! Get contribution from the other procs for parallel execution
    !------------------------------------------------------------------------------
    IF (ParEnv % PEs > 1) THEN
        !------------------------------------------------------------------------------
            tmpPower = TotalPower
            CALL MPI_ALLREDUCE( tmpPower, TotalPower, 1, MPI_DOUBLE_PRECISION, &
                MPI_SUM, MPI_COMM_WORLD, ierr )
        !------------------------------------------------------------------------------
    END IF
    !------------------------------------------------------------------------------
    ! Update the variable PreviousPower
    !------------------------------------------------------------------------------
    IF(TemperatureControlledPower) THEN
    !------------------------------------------------------------------------------
        PreviousPower = CurrentPower
    !------------------------------------------------------------------------------
    END IF !TemperatureControlledPower
    !------------------------------------------------------------------------------
    IF ( (ControlTotalPower) .AND. (PrevTime < Time) .AND. (ParEnv % MyPe==0) ) THEN
    !------------------------------------------------------------------------------
        TotalPowerFilename = 'totalpower'//'_'//TRIM(TestName)//'.csv'  
        !------------------------------------------------------------------------------
        ! Modify the file
        !------------------------------------------------------------------------------         
        OPEN(UNIT=1,FILE=TotalPowerFilename, FORM='FORMATTED', ACTION='WRITE', POSITION='APPEND', &
            IOSTAT=ios)
        WRITE(UNIT=1,FMT='(F13.4)', ADVANCE='no') PrevTime
        WRITE(UNIT=1,FMT='(A)', ADVANCE='no') ','
        WRITE(UNIT=1,FMT='(F16.4)', ADVANCE='yes') TotalPower
        CLOSE(1)
    !------------------------------------------------------------------------------
    END IF
    !------------------------------------------------------------------------------ 
    ! Compute the norm of the solution
    !------------------------------------------------------------------------------         

    Norm = Solver % Variable % Norm

    DO iter=1,NonlinearIter
        at  = CPUTime()
        at0 = RealTime()
        arealt = RealTime()

        CALL Info( 'HeatSolve', ' ', Level=4 )
        CALL Info( 'HeatSolve', ' ', Level=4 )
        CALL Info( 'HeatSolve', '-------------------------------------',Level=4 )
        WRITE( Message,* ) 'TEMPERATURE ITERATION', iter
        CALL Info( 'HeatSolve', Message, Level=4 )
        CALL Info( 'HeatSolve', '-------------------------------------',Level=4 )
        CALL Info( 'HeatSolve', ' ', Level=4 )
        CALL Info( 'HeatSolve', 'Starting Assembly...', Level=4 )

500    IF ( ConstantBulk .AND. ASSOCIATED(Solver % Matrix % BulkValues) ) THEN
         Solver % Matrix % Values = Solver % Matrix % BulkValues
         Solver % Matrix % RHS = Solver % Matrix % BulkRHS
         GOTO 1000
        END IF

!------------------------------------------------------------------------------
        CALL DefaultInitialize()
!------------------------------------------------------------------------------
 
        IF ( SmartHeaterControl .OR. IntegralHeaterControl ) THEN
          IF( SmartHeaterControl) ForceHeater = 0.0d0
          HeaterArea = 0.0d0
          HeaterSource = 0.0d0
          HeaterScaling = 1.0d0
          HeaterDensity = 0.0d0
          HeaterTarget = 0.0d0
          HeaterControlLocal = .FALSE.

          DO t=1,Solver % NumberOfActiveElements             
             Element => GetActiveElement(t)             
             Material => GetMaterial()

             BodyForce => GetBodyForce()
             
             IF ( .NOT. ASSOCIATED( BodyForce ) ) CYCLE
             bf_id = GetBodyForceId()
             
             IF( .NOT. (SmartHeaters(bf_id) .OR. IntegralHeaters(bf_id) ) ) CYCLE

             n = GetElementNOFNodes()

             Density(1:n) = GetReal( Material, 'Density' )
             Load(1:n) = GetReal( BodyForce, 'Heat Source' )

             s = ElementArea( Solver % Mesh, Element, n )

             IF( CurrentCoordinateSystem() == AxisSymmetric .OR. &
                  CurrentCoordinateSystem() == CylindricSymmetric ) s = 2 * PI * s

             HeaterSource(bf_id) = HeaterSource(bf_id) + s * SUM(Density(1:n) * Load(1:n)) / n
             HeaterArea(bf_id) = HeaterArea(bf_id) + s
             HeaterDensity(bf_id) = HeaterDensity(bf_id) + s * SUM( Density(1:n) ) / n
          END DO

          DO i = 1,Model % NumberOfBodyForces
             IF( IntegralHeaters(i) .OR. SmartHeaters(i) ) THEN
                HeaterDensity(i) = HeaterDensity(i) / HeaterArea(i)
             END IF
             IF(IntegralHeaters(i)) THEN
                HeaterTarget(i) = GetCReal(  Model % BodyForces(i) % Values, &
                     'Integral Heat Source', Found )
                HeaterScaling(i) = HeaterTarget(i) / HeaterSource(i)
             END IF
          END DO
        END IF

!------------------------------------------------------------------------------
        body_id = -1
        NULLIFY(Material)
!------------------------------------------------------------------------------
!      Bulk elements
!------------------------------------------------------------------------------
        CALL StartAdvanceOutput( 'HeatSolve', 'Assembly:' )
        DO t=1,Solver % NumberOfActiveElements

         CALL AdvanceOutput(t,GetNOFActive())
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where temperature 
!        should be calculated
!------------------------------------------------------------------------------
         Element => GetActiveElement(t)
         Material => GetMaterial()

         !------------------------------------------------------------------------------ 
         !       Specified heat capacity after death
         !------------------------------------------------------------------------------                 
         DeathCapacity = GetConstReal( Material, &
             'Death Heat Capacity', Found )
         IF ( .NOT.Found ) THEN
             DeathCapacity = 670.0
         END IF

!------------------------------------------------------------------------------
         IF ( Element % BodyId /= body_id ) THEN
!------------------------------------------------------------------------------
           Equation => GetEquation()
           ConvectionFlag = GetString( Equation, 'Convection', Found )

!------------------------------------------------------------------------------
           CompressibilityFlag = GetString( Material, &
                 'Compressibility Model', Found)
           IF ( .NOT.Found ) CompressibilityModel = Incompressible

           SELECT CASE( CompressibilityFlag )

             CASE( 'incompressible' )
               CompressibilityModel = Incompressible

             CASE( 'user defined' )
               CompressibilityModel = UserDefined1

             CASE( 'perfect gas', 'perfect gas equation 1' )
               CompressibilityModel = PerfectGas1

             CASE( 'thermal' )
               CompressibilityModel = Thermal

             CASE DEFAULT
               CompressibilityModel = Incompressible
           END SELECT
!------------------------------------------------------------------------------

           PhaseModel = GetString( Equation, 'Phase Change Model',Found )
           IF(.NOT. Found) PhaseModel = GetString( Material, 'Phase Change Model',Found )

           PhaseChange = Found .AND. (PhaseModel(1:4) /= 'none')
           IF ( PhaseChange ) THEN
              CheckLatentHeatRelease = GetLogical( Equation, &
                   'Check Latent Heat Release',Found )
           END IF
         END IF
!------------------------------------------------------------------------------

         n = GetElementNOFNodes()
         CALL GetElementNodes( ElementNodes )

         CALL GetScalarLocalSolution( LocalTemperature )
!------------------------------------------------------------------------------
!        Get element material parameters
!------------------------------------------------------------------------------
         !------------------------------------------------------------------------------
         ! Heat capacity: c 
         !------------------------------------------------------------------------------
         HeatCapacity(1:n) = GetReal( Material, 'Heat Capacity', Found )
         IF ( .NOT.Found ) THEN
             HeatCapacity(1:n)= 4180.0
         END IF


         CALL ListGetRealArray( Material,'Heat Conductivity',Hwrk,n, &
                      Element % NodeIndexes )

         HeatConductivity = 0.0d0
         IF ( SIZE(Hwrk,1) == 1 ) THEN
           DO i=1,3
             HeatConductivity( i,i,1:n ) = Hwrk( 1,1,1:n )
           END DO
         ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
           DO i=1,MIN(3,SIZE(Hwrk,1))
             HeatConductivity(i,i,1:n) = Hwrk(i,1,1:n)
           END DO
         ELSE
           DO i=1,MIN(3,SIZE(Hwrk,1))
             DO j=1,MIN(3,SIZE(Hwrk,2))
               HeatConductivity( i,j,1:n ) = Hwrk(i,j,1:n)
             END DO
           END DO
         END IF

         IF (ASSOCIATED(Hwrk)) DEALLOCATE( Hwrk )
         !------------------------------------------------------------------------------
         ! State of cells
         !------------------------------------------------------------------------------
         IF ( ASSOCIATED( CellStateSol ) ) THEN
             DO i=1,n
                 DO l=1,ADOFs
                     k = CellStatePerm(Element % NodeIndexes(i))
                     LocalCellState(i) = CellState((k-1)*ADOFs + l)
                 END DO
             END DO
         ELSE
             LocalCellState(n) = 0.0D0
         ENDIF
!------------------------------------------------------------------------------

         IF ( CompressibilityModel == PerfectGas1 ) THEN

           ! Read Specific Heat Ratio:
           !--------------------------
           SpecificHeatRatio = GetConstReal( Material, &
               'Specific Heat Ratio', Found )
           IF ( .NOT.Found ) SpecificHeatRatio = 5.d0/3.d0

           ! For an ideal gas, \gamma, c_p and R are really a constant
           ! GasConstant is an array only since HeatCapacity formally is:
           !-------------------------------------------------------------
           GasConstant(1:n) = ( SpecificHeatRatio - 1.d0 ) * &
               HeatCapacity(1:n) / SpecificHeatRatio

           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 1.0_dp
         ELSE IF ( CompressibilityModel == Thermal ) THEN
           ReferenceTemperature(1:n) = GetReal( Material, 'Reference Temperature' )
           HeatExpansionCoeff(1:n) = GetReal( Material, 'Heat Expansion Coefficient' )

           Density(1:n) = GetReal( Material,'Density' )
           Density(1:n) = Density(1:n) * ( 1 - HeatExpansionCoeff(1:n)  * &
                ( LocalTemperature(1:n) - ReferenceTemperature(1:n) ) )

           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) &
             PressureCoeff(1:n) = LocalTemperature(1:n) * HeatExpansionCoeff(1:n) / &
                   ( 1-HeatExpansionCoeff(1:n)*( &
                               LocalTemperature(1:n)-ReferenceTemperature(1:n)) )
         ELSE IF ( CompressibilityModel == UserDefined1 ) THEN
           IF ( ASSOCIATED( DensitySol ) ) THEN
             CALL GetScalarLocalSolution( Density, 'Density' ) 
           ELSE
             Density(1:n) = GetReal( Material,'Density' )
           END IF
           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 0.0_dp
         ELSE
           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 0.0_dp
           Density(1:n) = GetReal( Material, 'Density' )
         END IF

!------------------------------------------------------------------------------
! Take pressure deviation p_d as the dependent variable, p = p_0 + p_d
! for PerfectGas, read p_0
!------------------------------------------------------------------------------
         IF ( CompressibilityModel /= Incompressible ) THEN
           ReferencePressure = ListGetConstReal( Material, &
               'Reference Pressure', Found)
           IF ( .NOT.Found ) ReferencePressure = 0.0d0
         END IF
!------------------------------------------------------------------------------

         HeaterControlLocal = .FALSE.
         Load = 0.0D0
         Pressure = 0.0d0
         dPressuredt = 0.0d0
!------------------------------------------------------------------------------
!        Check for convection model
!------------------------------------------------------------------------------
         C1 = 1.0D0
         U = 0._dp
         V = 0._dp
         W = 0._dp

         MU = 0.0d0
         CALL GetVectorLocalSolution( MU, 'Mesh Velocity' )

         IF ( ConvectionFlag == 'constant' ) THEN

           U(1:n) = GetReal( Material, 'Convection Velocity 1', Found )
           IF ( .NOT. Found ) &
              U(1:n) = GetReal( Equation, 'Convection Velocity 1', Found )
           V(1:n) = GetReal( Material, 'Convection Velocity 2', Found )
           IF ( .NOT. Found ) &
             V(1:n) = GetReal( Equation, 'Convection Velocity 2', Found )
           W(1:n) = GetReal( Material, 'Convection Velocity 3', Found )
           IF ( .NOT. Found ) &
             W(1:n) = GetReal( Equation, 'Convection Velocity 3', Found )

         ELSE IF ( ConvectionFlag == 'computed' .AND. &
              ASSOCIATED(FlowSolution) ) THEN
           DO i=1,n
             k = FlowPerm(Element % NodeIndexes(i))
             IF ( k > 0 ) THEN
!------------------------------------------------------------------------------
               Pressure(i) = FlowSolution(NSDOFs*k) + ReferencePressure
               SELECT CASE( CompressibilityModel )
                 CASE( PerfectGas1 )
                   Density(i)  = Pressure(i) / &
                       ( GasConstant(i) * LocalTemperature(i) )
               END SELECT
               IF ( TransientSimulation ) THEN
                 dPressureDt(i) = ( FlowSolution(NSDOFs*k) - &
                     FlowSol % PrevValues(NSDOFs*k,1) ) / dt
               END IF
!------------------------------------------------------------------------------

               SELECT CASE( NSDOFs )
               CASE(3)
                 U(i) = FlowSolution( NSDOFs*k-2 )
                 V(i) = FlowSolution( NSDOFs*k-1 )
                 W(i) = 0.0D0

               CASE(4)
                 U(i) = FlowSolution( NSDOFs*k-3 )
                 V(i) = FlowSolution( NSDOFs*k-2 )
                 W(i) = FlowSolution( NSDOFs*k-1 )
               END SELECT
             ELSE
               U(i) = 0.0d0
               V(i) = 0.0d0
               W(i) = 0.0d0
             END IF
           END DO
         ELSE IF (ConvectionFlag=='computed' ) THEN
           CALL Warn( 'HeatSolver', 'Convection model specified ' //  &
                 'but no associated flow field present?' )
         ELSE
           IF ( ALL(MU==0) ) C1 = 0.0D0 
         END IF

         C0 = 0.0_dp

         BodyForce => GetBodyForce()
!------------------------------------------------------------------------------
!          Perfusion (added as suggested by Matthias Zenker)
!------------------------------------------------------------------------------
         PerfusionRate(1:n) = GetReal( BodyForce, 'Perfusion Rate', Found )

         IF ( Found ) THEN
           PerfusionRefTemperature(1:n) = GetReal( BodyForce, 'Perfusion Reference Temperature' )
           PerfusionDensity(1:n) = GetReal( BodyForce, 'Perfusion Density' )
           PerfusionHeatCapacity(1:n) = GetReal( BodyForce, 'Perfusion Heat Capacity' )
           C0(1:n) = PerfusionHeatCapacity(1:n) * PerfusionRate(1:n) * PerfusionDensity(1:n) 
         END IF

         !------------------------------------------------------------------------------
         !  No more convection in died cells...
         !------------------------------------------------------------------------------
         IF (CellsDeath) THEN
         !------------------------------------------------------------------------------
             DO i=1,n
                     !------------------------------------------------------------------------------
                     !  Adapt some parameters in function of the alive state of the cells
                     !------------------------------------------------------------------------------
                     IF (CellStateModel==1) THEN
                         IF (LocalCellState(i) == 0) THEN
                             C1(i) = 0.0D0
                             U(i) = 0.0
                             V(i) = 0.0
                             W(i) = 0.0
                         END IF
                     ELSE
                         IF (LocalCellState(i) > DeadThreshold) THEN
                             C1(i) = 0.0D0
                             U(i) = 0.0
                             V(i) = 0.0
                             W(i) = 0.0
                         END IF
                     END IF  
             END DO
             !------------------------------------------------------------------------------
             !   Death Heat Capacity in died cells
             !------------------------------------------------------------------------------
             DO i = 1,n
         !------------------------------------------------------------------------------
         ! Multiply perfusion coeff by blood density and heat capacity (HeatCapacity(1,i)=product already)
         !RMV CHECK THIS!!
         !------------------------------------------------------------------------------
                     !------------------------------------------------------------------------------
                     !  Check the alive state of the cells
                     !------------------------------------------------------------------------------
                     IF (CellStateModel==1) THEN
                         IF ((LocalCellState(i) == 0) .AND. (DeathCapacity/=0)) THEN
                             HeatCapacity( i ) = DeathCapacity
                         END IF
                     ELSE
                         IF ((LocalCellState(i) > DeadThreshold) .AND. (DeathCapacity/=0)) THEN
                             HeatCapacity( i ) = DeathCapacity
                         END IF
                     END IF
                 END DO
         !------------------------------------------------------------------------------
         END IF ! Cells Death modelling
!------------------------------------------------------------------------------
!        Check if modelling Phase Change with Eulerian approach 
!------------------------------------------------------------------------------
         PhaseSpatial = .FALSE.
         IF (  PhaseChange ) THEN
           CALL EffectiveHeatCapacity()
         ELSE
           HeatCapacity(1:n) = Density(1:n) * HeatCapacity(1:n)
         END IF

         Viscosity = 0.0d0
!------------------------------------------------------------------------------
!        Add body forces, if any
!------------------------------------------------------------------------------
         LOAD = 0.0D0
         power = 0.0D0         
         IF ( ASSOCIATED( BodyForce ) ) THEN
           bf_id = GetBodyForceId()
!------------------------------------------------------------------------------
!          Frictional viscous heating
!------------------------------------------------------------------------------
           IF ( GetLogical( BodyForce, 'Friction Heat',Found) ) THEN
              Viscosity(1:n) = GetReal( Material,'Viscosity' )
           END IF
           
           IF(TemperatureControlledPower) THEN
           !------------------------------------------------------------------------------
           !   If temperature-controlled power, get the power from computation: 
           !------------------------------------------------------------------------------
               DO i = 1,n
               !------------------------------------------------------------------------------
                   k = TempPerm(Element % NodeIndexes(i))
                   power(i) = PreviousPower(k) + &
                       PowerControl_Kp * (TargetTemperature-Temperature(k)) + & 
                       PowerControl_Ki * IntegErrorToTargetTemperature(k) + & 
                       PowerControl_Kd * DerivErrorToTargetTemperature(k)
                   !------------------------------------------------------------------------------
                   ! Control of max and min power: 
                   !------------------------------------------------------------------------------
                   power(i) = MIN(power(i),100000000.0)
                   power(i) = MAX(power(i),0.0)
                   !------------------------------------------------------------------------------
                   CurrentPower(k) = power(i)
                   !------------------------------------------------------------------------------
                   ! If T>373K, conductivity=0, modelised by power=0 in heat equation: 
                   !------------------------------------------------------------------------------
                   IF ((ElectricPowerCutOff) .AND. (Temperature(k)>373)) THEN
                       power(i) = 0.0
                   END IF
               !------------------------------------------------------------------------------ 
               END DO !i
           !------------------------------------------------------------------------------
           ELSE IF (UseElectricPower) THEN
               !------------------------------------------------------------------------------
               ! Read the electric power used as a coefficient of heat source: 
               !------------------------------------------------------------------------------
               power(1:n) = power(1:n) + & 
                   GetReal( BodyForce, 'Electric Power', Found )
               !------------------------------------------------------------------------------
               ! If T>373K, conductivity=0, modelised by power=0: 
               !------------------------------------------------------------------------------
               IF (ElectricPowerCutOff) THEN
                   DO i = 1,n
                       k = TempPerm(Element % NodeIndexes(i))
                       IF (Temperature(k)>373) THEN
                           power(i) = 0.0
                       END IF
                   END DO
               END IF
           !------------------------------------------------------------------------------
           END IF
!------------------------------------------------------------------------------
!          Given heat source
!------------------------------------------------------------------------------
           !Load(1:n) = GetReal( BodyForce, 'Heat Source', Found )
           !------------------------------------------------------------------------------
           ! Read the body force value in the input file and modify Load
           !------------------------------------------------------------------------------
           IF (InterpolatedElectricPower) THEN
               DO i = 1,n
                   k = TempPerm(Element % NodeIndexes(i))
                   LOAD(i) = LOAD(i) + power(i) * ElectricPower(k)
               END DO
           ELSE IF (UseElectricPower) THEN
               LOAD(1:n) = LOAD(1:n) + power(1:n) * & 
                   GetReal( BodyForce, 'Heat Source', Found ) 
           ELSE
               Load(1:n) = Density(1:n) *  GetReal( BodyForce, 'Heat Source', Found )
           END IF

           IF ( SmartHeaterControl .AND. NewtonLinearization .AND. SmartTolReached) THEN
              IF(  SmartHeaters(bf_id) ) THEN
               HeaterControlLocal = .TRUE.
               IF( TransientHeaterControl ) THEN
                 Load(1:n) = PrevPowerScaling * Load(1:n)
                 HeaterScaling(bf_id) = PrevPowerScaling
               END IF
             END IF
           END IF

           IF ( IntegralHeaterControl ) THEN
              IF( IntegralHeaters(bf_id) ) THEN
                 Load(1:n) = Load(1:n) * HeaterScaling(bf_id) 
              END IF
           END IF

         END IF

         CALL TotalPowerCompose( LOAD, Element, n, ElementNodes ) 

!------------------------------------------------------------------------------
! Note at this point HeatCapacity = \rho * c_p OR \rho * (c_p - R)
! and C1 = 0 (diffusion) or 1 (convection)
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!          Perfusion (added as suggested by Matthias Zenker)
!------------------------------------------------------------------------------

         IF ( Found ) THEN
                Load(1:n) = Load(1:n) + C0(1:n) * PerfusionRefTemperature(1:n)           
         END IF

!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!------------------------------------------------------------------------------
           CALL DiffuseConvectiveCompose( &
               MASS, STIFF, FORCE, LOAD, &
               HeatCapacity, C0, C1(1:n)*HeatCapacity(1:n), HeatConductivity, &
               PhaseSpatial, LocalTemperature, Enthalpy, U, V, W, &
               MU(1,1:n),MU(2,1:n),MU(3,1:n), Viscosity, Density, Pressure, &
               dPressureDt, PressureCoeff, CompressibilityModel /= Incompressible, &
               Stabilize, UseBubbles, Element, n, ElementNodes )
!------------------------------------------------------------------------------
         ELSE
!------------------------------------------------------------------------------
           CALL DiffuseConvectiveGenCompose( &
               MASS, STIFF, FORCE, LOAD, &
               HeatCapacity, C0, C1(1:n)*HeatCapacity(1:n), HeatConductivity, &
               PhaseSpatial, LocalTemperature, Enthalpy, U, V, W, &
               MU(1,1:n),MU(2,1:n),MU(3,1:n), Viscosity, Density, Pressure, &
               dPressureDt, PressureCoeff, CompressibilityModel /= Incompressible, &
               Stabilize, Element, n, ElementNodes )
!------------------------------------------------------------------------------
         END IF
!------------------------------------------------------------------------------

         IF ( HeaterControlLocal .AND. .NOT. TransientHeaterControl) THEN

           IF ( TransientAssembly .AND. .NOT. ConstantBulk ) THEN
             CALL Default1stOrderTime( MASS, STIFF, FORCE )
           END IF

           CALL UpdateGlobalEquations( Solver % Matrix, STIFF, &
               ForceHeater, FORCE, n, 1, TempPerm(Element % NodeIndexes) )
         ELSE
            Bubbles = UseBubbles .AND. .NOT.Stabilize .AND. &
            ( ConvectionFlag == 'computed' .OR. ConvectionFlag == 'constant' )
            
!------------------------------------------------------------------------------
!           If time dependent simulation add mass matrix to stiff matrix
!------------------------------------------------------------------------------
            TimeForce  = 0.0_dp
            IF ( TransientAssembly ) THEN
               IF ( ConstantBulk ) THEN
                 CALL DefaultUpdateMass( MASS )
               ELSE
                 CALL Default1stOrderTime( MASS,STIFF,FORCE )
               END IF
            ELSE IF ( Solver % NOFEigenValues>0 ) THEN
              CALL DefaultUpdateDamp(MASS)
            END IF
!------------------------------------------------------------------------------
!           Update global matrices from local matrices
!------------------------------------------------------------------------------
            IF (  Bubbles ) THEN
               CALL Condensate( N, STIFF, FORCE, TimeForce )
            END IF

            CALL DefaultUpdateEquations( STIFF, FORCE )
         END IF
!------------------------------------------------------------------------------
      END DO     !  Bulk elements
!------------------------------------------------------------------------------

      CALL DefaultFinishBulkAssembly()


      PRINT *,'TotalPower=',TotalPower
      PRINT *,'Integ_Force=',Integ_Force 
1000  CONTINUE

!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
      DO t=1, Solver % Mesh % NumberOfBoundaryElements

        !------------------------------------------------------------------------------
        ! Go through boundary elements
        !------------------------------------------------------------------------------
        Element => GetBoundaryElement(t)
        !------------------------------------------------------------------------------
        ! Check if the element is active and of suitable type
        !------------------------------------------------------------------------------
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        !TODO: clarify this
        IF ( GetElementFamily() == 1 ) CYCLE

        !------------------------------------------------------------------------------
        ! Get the number of nodes 
        !------------------------------------------------------------------------------
        n = GetElementNOFNodes()

        ! Check that the dimension of element is suitable for fluxes
        IF( .NOT. PossibleFluxElement(Element) ) CYCLE

        BC => GetBC()
        IF ( .NOT. ASSOCIATED(BC) ) CYCLE

        ! This checks whether there are any Dirichlet conditions on the 
        ! smart heater boundary. If there are the r.h.s. must be zero as 
        ! there can possibly not be any effect on temperature.
        !-----------------------------------------------------------------
        IF ( HeaterControlLocal .AND. .NOT. TransientHeaterControl) THEN
          IF( ListCheckPresent(BC, Varname) ) THEN
             nd = GetElementDOFs(Indexes)
             ForceHeater(TempPerm(Indexes(1:nd))) = 0.0_dp
          END IF
        END IF

        HeatFluxBC = GetLogical( BC, 'Heat Flux BC', Found )
        IF ( Found .AND. .NOT. HeatFluxBC ) CYCLE

        HeatGapBC = ListGetLogical( BC, 'Heat Gap', Found )
        CALL AddHeatFluxBC()

        IF ( HeatGapBC ) THEN
          CALL FindGapIndexes( Element, Indexes, n )
          SaveIndexes(1:n) = Element % NodeIndexes
          Element % NodeIndexes = Indexes(1:n)
          CALL AddHeatFluxBC()
          Element % NodeIndexes = SaveIndexes(1:n)
        END IF

      END DO   ! Neumann & Newton BCs
!------------------------------------------------------------------------------

      IF ( TransientSimulation .AND. ConstantBulk ) CALL AddGlobalTime()

      CALL DefaultFinishAssembly()
      CALL Info( 'HeatSolve', 'Assembly done', Level=4 )

      !------------------------------------------------------------------------------
      !  Dirichlet boundary conditions
      !------------------------------------------------------------------------------

      !------------------------------------------------------------------------------
      !  Dirichlet boundary conditions
      !------------------------------------------------------------------------------
      CALL NumaDefaultDirichletBCs(Solver) 


      !------------------------------------------------------------------------------
      ! Compute assembly CPU time, save current CPU time for solving CPU time below, 
      ! and save current solution norm
      !------------------------------------------------------------------------------
      at = CPUTime() - at
      arealt = RealTime() -arealt
      st = CPUTime()
      srealt = RealTime()

      PrevNorm = Norm
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------

      IF(SmartHeaterControl .AND. NewtonLinearization .AND. SmartTolReached) THEN
      
        IF(.NOT. TransientHeaterControl) THEN
 
          CALL ListAddLogical(SolverParams, &
              'Skip Compute Nonlinear Change',.TRUE.)

          Relax = GetCReal( SolverParams, &
              'Nonlinear System Relaxation Factor', Found )
          
          IF ( Found .AND. Relax /= 1.0d0 ) THEN
            CALL ListAddConstReal( Solver % Values, &
                'Nonlinear System Relaxation Factor', 1.0d0 )
          ELSE
            Relax = 1.0d0
          END IF          

          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              ForceHeater, XX, Norm, 1, Solver )
         
          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              Solver % Matrix % RHS, YY, Norm, 1, Solver )

          CALL ListAddLogical(SolverParams,'Skip Compute Nonlinear Change',.FALSE.)
        ELSE          
          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              Solver % Matrix % RHS, Temperature, Norm, 1, Solver )
          YY = Temperature
        END IF

        IF(.NOT. SmartHeaterAverage) THEN
          xave = XX(TempPerm(SmartHeaterNode))
          yave = YY(TempPerm(SmartHeaterNode))
        ELSE          
          xave = 0.0d0
          yave = 0.0d0
          j = 0
          
          DO k = Model % Mesh % NumberOfBulkElements + 1, &
              Model % Mesh % NumberOfBulkElements + Model % Mesh % NumberOfBoundaryElements            

            Element => Model % Mesh % Elements(k)            
            IF ( Element % BoundaryInfo % Constraint == SmartHeaterBC ) THEN
              l = Element % TYPE % NumberOfNodes
              j = j + l
              xave = xave + SUM( XX(TempPerm(Element % NodeIndexes)) )
              yave = yave + SUM( YY(TempPerm(Element % NodeIndexes)) )
            END IF
          END DO
          xave = xave / j
          yave = yave / j 
          CALL ListAddConstReal(Model % Simulation,'res: Smart Heater Temperature',yave)
        END IF

        IF(.NOT. TransientHeaterControl) THEN
          IF ( ASSOCIATED(Solver % Variable % NonlinValues) ) THEN
            Solver % Variable % NonlinValues = Temperature
          END IF

          PowerScaling = (MeltPoint - yave) / xave 
          Temperature = YY + PowerScaling * XX

          ! The change is computed separately for the controlled temperature field
          !-----------------------------------------------------------------------
          CALL ComputeChange(Solver,.FALSE.,LocalNodes,Temperature)
          Norm = Solver % Variable % Norm

        END IF

        IF(dt > PowerTimeScale) THEN
          IF ( Relax /= 1.0d0 ) THEN
            CALL ListAddConstReal( Solver % Values,  &
                'Nonlinear System Relaxation Factor', Relax )
          END IF
        END IF
      ELSE
!------------------------------------------------------------------------------
!     Check stepsize for nonlinear iteration
!------------------------------------------------------------------------------
        IF( DefaultLinesearch( Converged ) ) GOTO 500
        IF( Converged ) EXIT

        Norm = DefaultSolve()
      END IF


      IF( SmartHeaterControl .OR. IntegralHeaterControl) THEN
         
         CALL ListAddConstReal(Model % Simulation,'res: Heater Power Scaling',PowerScaling)
         
         CALL Info( 'HeatSolve', 'Heater Control Information', Level=4 )
         DO i=1,Model % NumberOfBodyForces
            IF( .NOT. (SmartHeaters(i) .OR. IntegralHeaters(i))) CYCLE
            IF( SmartHeaters(i) )  HeaterScaling(i) = PowerScaling

            WRITE( Message, '(A,T35,I15)' ) 'Heater for body: ', i
            CALL Info( 'HeatSolve', Message, Level=4 )
            IF(SmartHeaters(i)) WRITE( Message, '(A,T35,A)' ) 'Heater type:','Smart heater'
            IF(IntegralHeaters(i)) WRITE( Message, '(A,T35,A)' ) 'Heater type:','Integral heater'
            CALL Info( 'HeatSolve', Message, Level=4 )

            WRITE( Message,'(A,T35,ES15.4)' ) 'Heater Volume (m^3): ', HeaterArea(i)
            CALL Info( 'HeatSolve', Message, Level=4 )
            s = HeaterSource(i) * HeaterScaling(i)
            WRITE( Message,'(A,T35,ES15.4)' ) 'Heater Power (W): ', s
            CALL Info( 'HeatSolve', Message, Level=4 )

            WRITE( Message,'(A,T35,ES15.4)' ) 'Heater scaling: ', HeaterScaling(i)
            CALL Info( 'HeatSolve', Message, Level=4 )
            WRITE( Message, '(A,T35,ES15.4)' ) 'Heater Power Density (W/kg): ', s/(HeaterDensity(i) * HeaterArea(i))
            CALL Info( 'HeatSolve', Message, Level=4 )
            
            IF( SmartHeaters(i)) CALL ListAddConstReal(Model % Simulation,'res: Heater Power Density',&
                 s/(HeaterDensity(i) * HeaterArea(i)))
         END DO
      END IF

      !------------------------------------------------------------------------------         
      ! Compute Max temperature 
      !------------------------------------------------------------------------------         
      IF(ControlMaxTemperature) THEN
          MaxTemperature = 0.0d0
          DO i=1,LocalNodes
                  IF (MaxTemperature < Temperature(i)) THEN
                      MaxTemperature = Temperature(i)
                  END IF
          END DO
      END IF

      !------------------------------------------------------------------------------         
      ! Compute variables for temperature-controlled power:
      !------------------------------------------------------------------------------         
      IF (TemperatureControlledPower) THEN
      !------------------------------------------------------------------------------ 
          ! Save error between target and current temperature and update history:
          !------------------------------------------------------------------------------ 
          DO j=1,PowerControl_Integ_Length-1
              HistoryErrorToTargetTemperature(PowerControl_Integ_Length-j+1,:) = &
                  HistoryErrorToTargetTemperature(PowerControl_Integ_Length-j,:)
          END DO
          HistoryErrorToTargetTemperature(1,:) = TargetTemperature-Temperature

          HistoryErrorToTargetTemperatureVariable => HistoryErrorToTargetTemperature(1,:)

          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
          Solver, 'Error to target Temperature', 1, &
          HistoryErrorToTargetTemperatureVariable, TempPerm ) 
          !------------------------------------------------------------------------------ 
          ! Integral of error between target and current temperature:
          !------------------------------------------------------------------------------ 
          IF (PowerControl_Integ_Length /= 0) THEN
          !------------------------------------------------------------------------------
          !If specified length of integral, use history variable
          !------------------------------------------------------------------------------
              IntegErrorToTargetTemperature = 0.0D0
              DO j=1,PowerControl_Integ_Length
                  IntegErrorToTargetTemperature(:) = IntegErrorToTargetTemperature(:) + &
                      HistoryErrorToTargetTemperature(j,:)*dt
              END DO
          !------------------------------------------------------------------------------
          ELSE
          !------------------------------------------------------------------------------
          !Else, make the integral over all time-steps
          !------------------------------------------------------------------------------             
              IntegErrorToTargetTemperature(:) = IntegErrorToTargetTemperature(:) + &
                      HistoryErrorToTargetTemperature(1,:)*dt
          !------------------------------------------------------------------------------
          END IF
          !------------------------------------------------------------------------------ 
          ! Derivative of error between target and current temperature:
          !------------------------------------------------------------------------------ 
          DerivErrorToTargetTemperature(:) = (HistoryErrorToTargetTemperature(1,:)- &
              HistoryErrorToTargetTemperature(2,:))/dt
      !------------------------------------------------------------------------------ 
      END IF !TemperatureControlledPower

      !------------------------------------------------------------------------------ 
      !   Add the heat source as a variable for visualization
      !------------------------------------------------------------------------------ 
      VisualizeHeatSource = GetLogical( SolverParams,'Heat Source Visualization',Found )
      IF ( .NOT.Found ) VisualizeHeatSource = .FALSE.
      !------------------------------------------------------------------------------
      IF(VisualizeHeatSource) THEN
      !------------------------------------------------------------------------------
              HeatSource = 0.0D0
      !------------------------------------------------------------------------------ 
      !   Go through bulk elements and get heat source:
      !------------------------------------------------------------------------------
              DO t=1,Solver % NumberOfActiveElements
      !------------------------------------------------------------------------------
                  Element => GetActiveElement(t)
                  n = GetElementNOFNodes()
                  CALL GetElementNodes( ElementNodes )    

                  BodyForce => GetBodyForce()
      !------------------------------------------------------------------------------
                  IF ( ASSOCIATED( BodyForce ) ) THEN
      !------------------------------------------------------------------------------ 
                      LocalHeatSource = 0.0D0
      !------------------------------------------------------------------------------
                          !------------------------------------------------------------------------------
                          ! Read the body force value in the input file: 
                          !------------------------------------------------------------------------------
                          LocalHeatSource(1:n) = LocalHeatSource(1:n) + &
                                  GetReal( BodyForce, 'Heat Source', Found )
      !------------------------------------------------------------------------------
                          DO i=1,n
      !------------------------------------------------------------------------------
                              k = TempPerm(Element % NodeIndexes(i))
                              HeatSource(k) = LocalHeatSource(i)
      !------------------------------------------------------------------------------
                          END DO ! i
      !------------------------------------------------------------------------------ 
                  END IF
      !------------------------------------------------------------------------------ 
              END DO ! t
      !------------------------------------------------------------------------------ 
              CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
              Solver, 'Heat source', 1, &
              HeatSource, TempPerm ) 
      !------------------------------------------------------------------------------ 
      END IF

      !------------------------------------------------------------------------------ 
      !   Add the heat source as a variable for visualization
      !------------------------------------------------------------------------------ 
      VisualizePerfusion = GetLogical( SolverParams,'Perfusion Visualization',Found )
      IF ( .NOT.Found ) VisualizePerfusion = .FALSE.
      !------------------------------------------------------------------------------
      IF(VisualizePerfusion) THEN
      !------------------------------------------------------------------------------
          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
              Solver, 'Perfusion', 1, C0, TempPerm ) 
      !------------------------------------------------------------------------------ 
      END IF !VisualizePerfusion
      !------------------------------------------------------------------------------      
      ! Write some information messages
      !------------------------------------------------------------------------------      
      !------------------------------------------------------------------------------
      ! Compute solving CPU time, and total assembly and solving CPU time in the 
      ! coupled system iteration (may contain several nonlinear iterations) 
      !------------------------------------------------------------------------------


      st = CPUTIme()-st
      totat = totat + at
      totst = totst + st
        WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,', Assembly CPU Time (nonlinear it., coupled it.): (s)', at, totat
        CALL Info( 'NumaHeatSolve', Message, Level=4 )
        WRITE(Message,'(a,i4,a,F8.2)') 'iter: ',iter,', Assembly Real Time (nonlinear it): (s)', arealt
        CALL Info( 'NumaHeatSolve', Message, Level=4 )
        WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,', Solve CPU Time (nonlinear it., coupled it.):    (s)', st, totst
        CALL Info( 'NumaHeatSolve', Message, Level=4 )
        WRITE(Message,'(a,i4,a,F8.2)') 'iter: ',iter,', Solve Real Time (nonlinear it): (s)', srealt
        CALL Info( 'NumaHeatSolve', Message, Level=4 )
!------------------------------------------------------------------------------
!     If modelling phase change (and if requested by the user), check if any
!     node has jumped over the phase change interval, and if so, reduce
!     timestep and or relaxation and recompute.
!------------------------------------------------------------------------------
      IF (PhaseChange .AND. CheckLatentHeatRelease ) THEN
!------------------------------------------------------------------------------
        IF ( CheckLatentHeat() ) THEN
          Temperature(1:LocalNodes) = PrevSolution
          Norm = PrevNorm

          IF ( TransientSimulation ) THEN
            dt = dt / 2
            Solver % dt = dt
            WRITE( Message, * ) &
                  'Latent heat release check: reducing timestep to: ',dt
            CALL Info( 'HeatSolve', Message, Level=4 )
          ELSE
            Relax = Relax / 2
            CALL  ListAddConstReal( Solver % Values,  &
                 'Nonlinear System Relaxation Factor', Relax )
            WRITE( Message, * ) &
                 'Latent heat release check: reducing relaxation to: ',Relax
            CALL Info( 'HeatSolve', Message, Level=4 )
          END IF

          CYCLE
        END IF
        IF ( .NOT.TransientSimulation ) PrevSolution=Temperature(1:LocalNodes)
      END IF
!------------------------------------------------------------------------------
     
      RelativeChange = Solver % Variable % NonlinChange

      WRITE( Message, * ) 'Result Norm   : ',Norm
      CALL Info( 'HeatSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'HeatSolve', Message, Level=4 )

      IF ( RelativeChange < NewtonTol .OR. iter >= NewtonIter ) &
               NewtonLinearization = .TRUE.
      Converged =  ( Solver % Variable % NonlinConverged == 1 ) .AND. &
          ( .NOT. SmartHeaterControl .OR. SmartTolReached )
      IF( Converged ) EXIT

      IF(SmartHeaterControl) THEN
        IF ( RelativeChange < SmartTol ) THEN
          SmartTolReached = .TRUE.
          YY = Temperature
        END IF
      END IF
      
!------------------------------------------------------------------------------
    END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------
    IF(TransientHeaterControl) THEN
      PowerRelax = GetCReal(Solver % Values,'Smart Heater Relaxation Factor', GotIt)
      IF(.NOT. GotIt) PowerRelax = 1.0_dp
      PowerSensitivity = ListGetConstReal(Solver % Values,'Smart Heater Power Sensivity',GotIt)
      IF(.NOT. GotIt) PowerSensitivity = 4.0_dp
      PowerScaling = PowerScaling * (1 + PowerSensitivity * PowerRelax * (MeltPoint/yave - 1.0d0) ) 

      IF( ListGetLogical( Solver % Values,'Smart Heater Transient Speedup',GotIt ) ) THEN
        Temperature = Temperature * (1 + PowerRelax * (MeltPoint/yave - 1.0d0)   )     
      END IF
      YY = Temperature
    END IF

!------------------------------------------------------------------------------
!   Compute cumulative time done by now and time remaining
!------------------------------------------------------------------------------
    IF ( .NOT. TransientSimulation ) EXIT
    CumulativeTime = CumulativeTime + dt
    dt = Timestep - CumulativeTime

   END DO ! time interval
   Solver % dt = Timestep

!------------------------------------------------------------------------------
   CALL  ListAddConstReal( Solver % Values,  &
        'Nonlinear System Relaxation Factor', SaveRelax )
!------------------------------------------------------------------------------

   DEALLOCATE( PrevSolution )

   IF ( ListGetLogical( Solver % Values, 'Adaptive Mesh Refinement', Found ) ) &
      CALL RefineMesh( Model,Solver,Temperature,TempPerm, &
            HeatInsideResidual, HeatEdgeResidual, HeatBoundaryResidual )

CONTAINS



!------------------------------------------------------------------------------
    SUBROUTINE TotalPowerCompose( LoadVector,Element,n,Nodes )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Return element local matrices and RHS vector for diffusion-convection
!  equation: 
!
!  ARGUMENTS:
!
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT: Nodal values of RHS
!
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!******************************************************************************
        REAL(KIND=dp) :: LoadVector(:)
        INTEGER :: n
        TYPE(Nodes_t) :: Nodes
        TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
        REAL(KIND=dp) :: Basis(2*n), dBasisdx(2*n,3), ddBasisddx(n,3,3), &
            SqrtElementMetric, s, u, v, w, Force

        REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
            
        INTEGER :: i,j,k,p,q,t,N_Integ,NBasis
        LOGICAL :: stat
        TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------
!       Get numbre of basis functions       
!------------------------------------------------------------------------------ 
        NBasis = n
!------------------------------------------------------------------------------
!     Integration stuff
!------------------------------------------------------------------------------     
        IntegStuff = GaussPoints( element )

        U_Integ => IntegStuff % u
        V_Integ => IntegStuff % v
        W_Integ => IntegStuff % w
        S_Integ => IntegStuff % s
        N_Integ =  IntegStuff % n
!------------------------------------------------------------------------------
!     Loop over Gauss integration points
!------------------------------------------------------------------------------
        DO t=1,N_Integ
!------------------------------------------------------------------------------
            u = U_Integ(t)
            v = V_Integ(t)
            w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
            stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx,ddBasisddx,.FALSE. )         
            s = SqrtElementMetric * S_Integ(t)      
!------------------------------------------------------------------------------
!        Loop over basis functions
!------------------------------------------------------------------------------
            Force = SUM( LoadVector(1:n)*Basis(1:n) )
            TotalPower =  TotalPower + Force * s 
        END DO
!------------------------------------------------------------------------------
   END SUBROUTINE TotalPowerCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION NumaDistanceToElectricTip (x, y, z, xtip, ytip, ztip, nbpoints, xproj, yproj, zproj ) RESULT( Distance )
!------------------------------------------------------------------------------ 
!******************************************************************************
!
!   Compute the distance of the point (x,y,z) from a given Tip defined by using
!     linear interpolation between given points + parametric representation
!
!   ARGUMENTS:
!
!       REAL(KIND=dp) :: x,y,z
!           INPUT: Coordinates of the point at which we compute the distance
!
!       REAL(KIND=dp) :: xtip(nbpoints),ytip(nbpoints),ztip(nbpoints)
!           INPUT: Coordinates of the points used for the linear interpolation 
!
!   INTEGER :: nbpoints
!           INPUT: Number of interpolation points
!
!******************************************************************************
        REAL(KIND=dp) :: Distance, x, y, z, xtip(nbpoints), ytip(nbpoints), ztip(nbpoints), &
            xproj, yproj, zproj
        INTEGER :: nbpoints 
        !------------------------------------------------------------------------------ 
        !   Local variables
        !------------------------------------------------------------------------------ 
        INTEGER :: i,j, S_MAX 
        REAL(KIND=dp) :: s, d, x_s, y_s, z_s
        !------------------------------------------------------------------------------   
        Distance = 100000.0D0   
        !------------------------------------------------------------------------------ 
        !   Case of point source
        !------------------------------------------------------------------------------ 
        IF (nbpoints==1) THEN
        !------------------------------------------------------------------------------ 
            Distance = sqrt( (x-xtip(1))**2 + (y-ytip(1))**2 + (z-ztip(1))**2 )
            xproj =xtip(1)
            yproj =ytip(1)
            zproj =ztip(1)
        !------------------------------------------------------------------------------ 
        ELSE
            !------------------------------------------------------------------------------ 
            !   For each linear part, compute the minimal distance in function of parameter s
            !------------------------------------------------------------------------------  
            DO i=1,nbpoints-1
            !------------------------------------------------------------------------------
                s = 0.0
                S_MAX = 100
                !------------------------------------------------------------------------------
                DO j = 1,S_MAX-1
                !------------------------------------------------------------------------------
                    x_s = (xtip(i+1)-xtip(i))*s + xtip(i)
                    y_s = (ytip(i+1)-ytip(i))*s + ytip(i)
                    z_s = (ztip(i+1)-ztip(i))*s + ztip(i)

                    d = sqrt( (x-x_s)**2 + (y-y_s)**2 + (z-z_s)**2 )
                    IF (d<Distance) THEN 
                        Distance = d
                        xproj =x_s
                        yproj =y_s
                        zproj =z_s
                    END IF
                    s = j*1.0/(S_MAX-1)
                !------------------------------------------------------------------------------
                END DO
            !------------------------------------------------------------------------------
            END DO
        !------------------------------------------------------------------------------
        END IF
!------------------------------------------------------------------------------ 
    END FUNCTION NumaDistanceToElectricTip
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE AddHeatFluxBC()
!------------------------------------------------------------------------------
      CALL GetElementNodes( ElementNodes )

      HeatTransferCoeff = 0.0D0
      LOAD  = 0.0D0
!------------------------------------------------------------------------------
!     BC: -k@T/@n = \epsilon\sigma(T^4 - Text^4)
!------------------------------------------------------------------------------
      RadiationFlag = GetString( BC, 'Radiation', Found )

      IF ( Found .AND. RadiationFlag(1:4) /= 'none' ) THEN
        NodalEmissivity(1:n) = GetReal(BC, 'Emissivity', Found)
        IF(.NOT. Found) THEN
           NodalEmissivity(1:n) = GetParentMatProp( 'Emissivity' )
        END IF
        Emissivity = SUM( NodalEmissivity(1:n) ) / n

!------------------------------------------------------------------------------
        IF (  RadiationFlag(1:9) == 'idealized' ) THEN
          AText(1:n) = GetReal( BC, 'Radiation External Temperature',Found )
          IF(.NOT. Found) AText(1:n) = GetReal( BC, 'External Temperature' )
        ELSE
          CALL DiffuseGrayRadiation( Model, Solver, Element, & 
              Temperature, TempPerm, ForceVector, VisibleFraction, Text)

          IF( GetLogical( BC, 'Radiation Boundary Open', Found) ) THEN
            AText(1:n) = GetReal( BC, 'Radiation External Temperature',Found )
            IF(.NOT. Found) AText(1:n) = GetReal( BC, 'External Temperature' )
            IF( VisibleFraction >= 1.0_dp ) THEN
              Atext(1:n) = Text
            ELSE
              Atext(1:n) = ( (1 - VisibleFraction) * Atext(1:n)**4 + &
                  VisibleFraction * Text**4 ) ** 0.25_dp
            END IF
          ELSE
            AText(1:n) = Text
          END IF
        END IF
!------------------------------------------------------------------------------
!       Add our own contribution to surface temperature (and external
!       if using linear type iteration or idealized radiation)
!------------------------------------------------------------------------------
        DO j=1,n
          k = TempPerm(Element % NodeIndexes(j))
          Text = AText(j)

          IF ( .NOT. HeatGapBC .AND. NewtonLinearization ) THEN
             HeatTransferCoeff(j) = Emissivity * 4*Temperature(k)**3 * &
                               StefanBoltzmann
             LOAD(j) = Emissivity*(3*Temperature(k)**4+Text**4) * &
                               StefanBoltzmann
          ELSE
             HeatTransferCoeff(j) = Emissivity * (Temperature(k)**3 + &
             Temperature(k)**2*Text+Temperature(k)*Text**2 + Text**3) * &
                               StefanBoltzmann 
             LOAD(j) = HeatTransferCoeff(j) * Text
          END IF
        END DO
      END IF  ! of radition
!------------------------------------------------------------------------------

      Work(1:n)  = GetReal( BC, 'Heat Transfer Coefficient',Found )
      IF ( Found ) THEN
        AText(1:n) = GetReal( BC, 'External Temperature',Found )
        DO j=1,n
!------------------------------------------------------------------------------
!         BC: -k@T/@n = \alpha(T - Text)
!------------------------------------------------------------------------------
          k = TempPerm(Element % NodeIndexes(j))
          LOAD(j) = LOAD(j) + Work(j) * AText(j)
          HeatTransferCoeff(j) = HeatTransferCoeff(j) + Work(j)
        END DO
      END IF

!------------------------------------------------------------------------------
!     BC: -k@T/@n = (rho*L)*v.n 
!     Heating related to pulling is possible only in ss cases where pull velocity
!     is desrcibed.
!------------------------------------------------------------------------------

      IF( GetLogical( BC, 'Phase Change',Found ) ) THEN
         PhaseVelocity(1,1:n) = GetReal( BC,'Phase Velocity 1', Found  )
         PhaseVelocity(2,1:n) = GetReal( BC,'Phase Velocity 2', Found  )
         PhaseVelocity(3,1:n) = GetReal( BC,'Phase Velocity 3', Found  )
  
         ! Ensure that the latent heat and density come from the same side
         LatentHeat(1:n) = GetParentMatProp( 'Latent Heat', &
              UElement = Element, UParent = Parent )
         IF(.NOT. ASSOCIATED(Parent) ) THEN
           CALL Warn('HeatSolve','Parent not associated')
         ELSE
           k = GetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')
           Density(1:n) = GetReal( Model % Materials(k) % Values, 'Density' )
         END IF

         ! This could be rather put as a new type of BC into the assembly routine and 
         ! then the Normal could be taken at the proper Gaussian integration points. 
         Normal = NormalVector( Element, ElementNodes, 0.0_dp, 0.0_dp, .TRUE. )

         DO i=1,n
            LOAD(i) = LOAD(i) + &
                 LatentHeat(i) * Density(i) * SUM( Normal(1:3) * PhaseVelocity(1:3,i))
         END DO
      END IF

!------------------------------------------------------------------------------
!     BC: -k@T/@n = g
!------------------------------------------------------------------------------
      LOAD(1:n) = LOAD(1:n) +  GetReal( BC, 'Heat Flux', Found )

      InfBC = ListGetLogical( BC,'Infinity BC '//TRIM(VarName),GotIt)
      IF( InfBC ) THEN
        AText(1:n) = GetReal( BC,'Infinity BC '//TRIM(VarName)//' Offset',GotIt)
        ! currently only isotropic heat conductivity supported
        HeatConductivityIso(1:n) = GetParentMatProp('Heat Conductivity',Element,GotIt)
        IF(.NOT. GotIt) THEN
          CALL Fatal( 'HeatSolver','Could not find > Heat Conductivity < for parent!' )           
        END IF
      END IF

!------------------------------------------------------------------------------
!     Get element matrix and rhs due to boundary conditions ...
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() == Cartesian ) THEN
        CALL DiffuseConvectiveBoundary( STIFF,FORCE, &
            LOAD,HeatTransferCoeff,InfBC,HeatConductivityIso,AText(1:n),&
            Element,n,ElementNodes )
      ELSE
        IF( InfBC ) THEN
          CALL Fatal('HeatSolver','Infinity BC not implemented only for cartersian case!')
        END IF
        CALL DiffuseConvectiveGenBoundary(STIFF,FORCE,&
            LOAD,HeatTransferCoeff,Element,n,ElementNodes ) 
      END IF

!------------------------------------------------------------------------------
!     Update global matrices from local matrices
!------------------------------------------------------------------------------
      IF ( TransientAssembly .AND. .NOT. ConstantBulk ) THEN
        MASS = 0.d0
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
      END IF

      IF ( HeatGapBC ) &
        CALL AddHeatGap( Solver, Element, STIFF, TempPerm)

      CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
    END SUBROUTINE AddHeatFluxBC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE AddGlobalTime()
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,l,n
      REAL(KIND=dp) :: FORCE(1)
      REAL(KIND=dp), POINTER :: SaveValues(:) => NULL()
      SAVE STIFF, MASS, X
      REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:),MASS(:,:), X(:,:)

      IF ( .NOT.ASSOCIATED(Solver % Variable % Values, SaveValues) ) THEN
         IF ( ALLOCATED(STIFF) ) DEALLOCATE( STIFF,MASS,X )
         n = 0
         DO i=1,Solver % Matrix % NumberOfRows
           n = MAX( n,Solver % Matrix % Rows(i+1)-Solver % Matrix % Rows(i) )
         END DO
         k = SIZE(Solver % Variable % PrevValues,2)
         ALLOCATE( STIFF(1,n),MASS(1,n),X(n,k) )
 
         SaveValues => Solver % Variable % Values
      END IF

      DO i=1,Solver % Matrix % NumberOFRows
        n = 0
        DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
          n=n+1
          STIFF(1,n) = Solver % Matrix % Values(j)
          MASS(1,n)  = Solver % Matrix % MassValues(j)
          X(n,:) = Solver % Variable % PrevValues(Solver % Matrix % Cols(j),:)
        END DO
        FORCE(1) = Solver % Matrix % RHS(i)
        Solver % Matrix % Force(i,1) = FORCE(1)
        k = MIN( Solver % DoneTime, Solver % Order )
        CALL BDFLocal( n, dt, MASS, STIFF, FORCE, X, k )

        n = 0
        DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
           n=n+1
          Solver % Matrix % Values(j) = STIFF(1,n)
        END DO
        Solver % Matrix % RHS(i) = FORCE(1)
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE AddGlobalTime
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    SUBROUTINE DiffuseGrayRadiation( Model, Solver, Element,  &
      Temperature, TempPerm, ForceVector,AngleFraction, Text)
!------------------------------------------------------------------------------
      TYPE(Model_t)  :: Model
      TYPE(Solver_t) :: Solver
      TYPE(Element_t), POINTER :: Element
      INTEGER :: TempPerm(:)
      REAL(KIND=dp) :: Temperature(:), ForceVector(:)
      REAL(KIND=dp) :: AngleFraction, Text
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Area, Asum
      INTEGER :: i,j,k,l,m,ImplicitFactors
      INTEGER, POINTER :: ElementList(:)
!------------------------------------------------------------------------------
!     If linear iteration compute radiation load
!------------------------------------------------------------------------------


      Asum = 0.0d0
      IF ( .NOT. NewtonLinearization ) THEN
        Text = ComputeRadiationLoad( Model, Solver % Mesh, Element, &
                 Temperature, TempPerm, Emissivity, AngleFraction)
      ELSE   !  Full Newton-Raphson solver
!------------------------------------------------------------------------------
!       Go trough surfaces (j) this surface (i) is getting
!       radiated from.
!------------------------------------------------------------------------------

        Area  = ElementArea( Solver % Mesh, Element, n )
        ElementList => Element % BoundaryInfo % GebhardtFactors % Elements

        DO j=1,Element % BoundaryInfo % GebhardtFactors % NumberOfFactors

          RadiationElement => Solver % Mesh % Elements( ElementList(j) )

          Text = ComputeRadiationCoeff(Model,Solver % Mesh,Element,j) / ( Area )
          Asum = Asum + Text
!------------------------------------------------------------------------------
!         Gebhardt factors are given elementwise at the center
!         of the element, so take avarage of nodal temperatures
!         (or integrate over surface j)
!------------------------------------------------------------------------------

          k = RadiationElement % TYPE % NumberOfNodes
          ImplicitFactors = Element % BoundaryInfo % GebhardtFactors % NumberOfImplicitFactors
          IF(ImplicitFactors == 0) &
              ImplicitFactors = Element % BoundaryInfo % GebhardtFactors % NumberOfFactors

          IF(j <= ImplicitFactors) THEN
            
            S = (SUM( Temperature( TempPerm( RadiationElement % &
                NodeIndexes))**4 )/k )**(1.0d0/4.0d0)
!------------------------------------------------------------------------------
!         Linearization of the G_jiT^4_j term
!------------------------------------------------------------------------------
            HeatTransferCoeff(1:n) = -4 * Text * S**3 * StefanBoltzmann
            LOAD(1:n) = -3 * Text * S**4 * StefanBoltzmann
!------------------------------------------------------------------------------
!         Integrate the contribution of surface j over surface i
!         and add to global matrix
!------------------------------------------------------------------------------
            CALL IntegOverA( STIFF, FORCE, LOAD, &
                HeatTransferCoeff, Element, n, k, ElementNodes ) 
            
            IF ( TransientAssembly ) THEN
              MASS = 0.d0
              CALL Add1stOrderTime( MASS, STIFF, &
                  FORCE,dt,n,1,TempPerm(Element % NodeIndexes),Solver )
            END IF
            
            DO m=1,n
              k1 = TempPerm( Element % NodeIndexes(m) )
              DO l=1,k
                k2 = TempPerm( RadiationElement % NodeIndexes(l) )
                CALL AddToMatrixElement( StiffMatrix,k1, &
                    k2,STIFF(m,l) )
              END DO
              ForceVector(k1) = ForceVector(k1) + FORCE(m)
            END DO

          ELSE

            S = (SUM( Temperature( TempPerm( RadiationElement % &
                NodeIndexes))**4 )/k )
            
            HeatTransferCoeff(1:n) = 0.0d0
            LOAD(1:n) = Text * S * StefanBoltzmann
            
            CALL IntegOverA( STIFF, FORCE, LOAD, &
                HeatTransferCoeff, Element, n, k, ElementNodes ) 
            
            DO m=1,n
              k1 = TempPerm( Element % NodeIndexes(m) )
              ForceVector(k1) = ForceVector(k1) + FORCE(m)
            END DO
            
          END IF 

        END DO

!------------------------------------------------------------------------------
!       We have already added all external temperature contributions
!       to the matrix for the Newton type iteration
!------------------------------------------------------------------------------
        AngleFraction = Asum / Emissivity
        Text = 0.0

      END IF  !  of newton-raphson

    END SUBROUTINE DiffuseGrayRadiation
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE EffectiveHeatCapacity()
      LOGICAL :: Found, Specific
      REAL(KIND=dp), ALLOCATABLE :: dT(:)

!------------------------------------------------------------------------------
!     See if temperature gradient indside the element is large enough 
!     to use  the c_p = SQRT( (dH/dx)^2 / (dT/dx)^2 ), otherwise
!     use c_p = dH/dT, or if in time dependent simulation, use
!     c_p = (dH/dt) / (dT/dt), if requested. 
!------------------------------------------------------------------------------

      SELECT CASE(PhaseModel)
!------------------------------------------------------------------------------
        CASE( 'spatial 1' )
          PhaseChangeModel = PHASE_SPATIAL_1
!------------------------------------------------------------------------------

        CASE( 'spatial 2' )
!------------------------------------------------------------------------------
! Check if local variation of temperature is large enough to actually use the
! Spatial 2 model. Should perhaps be scaled to element size (or actually
! compute the gradient, but this will do for now...).
!------------------------------------------------------------------------------
          s = MAXVAL(LocalTemperature(1:n))-MINVAL(LocalTemperature(1:n))
          IF ( s < AEPS ) THEN
            PhaseChangeModel = PHASE_SPATIAL_1
          ELSE
            PhaseChangeModel = PHASE_SPATIAL_2
          END IF

!------------------------------------------------------------------------------
! Note that here HeatCapacity is miused for saving dT.
!------------------------------------------------------------------------------
        CASE('temporal')
          IF ( TransientSimulation )  THEN
            ALLOCATE( dT(n) )
            dT(1:n) = Temperature(TempPerm(Element % NodeIndexes)) - &
                     PrevTemperature(TempPerm(Element % NodeIndexes))

            IF ( ANY(ABS(dT(1:n)) < AEPS) ) THEN
              PhaseChangeModel = PHASE_SPATIAL_1
            ELSE
              PhaseChangeModel = PHASE_TEMPORAL
            END IF
          ELSE
             PhaseChangeModel = PHASE_SPATIAL_1
          END IF

!------------------------------------------------------------------------------
        CASE DEFAULT
          PhaseChangeModel = PHASE_SPATIAL_1

      END SELECT
!------------------------------------------------------------------------------

      PhaseSpatial = ( PhaseChangeModel == PHASE_SPATIAL_2 )
      Specific = ListCheckPresent( Material,'Specific Enthalpy')

!-----------------------------------------------------------------------------
      SELECT CASE( PhaseChangeModel )

!------------------------------------------------------------------------------
! This phase change model is available only for some type of real entries 
! that have an implemented analytical derivation rule.
!-----------------------------------------------------------------------------
      CASE( PHASE_SPATIAL_1 )
        HeatCapacity(1:n) = ListGetReal( Material, &
             'Effective Heat Capacity', n,Element % NodeIndexes, Found )
        IF ( .NOT. Found ) THEN
          IF( Specific ) THEN
            HeatCapacity(1:n) = ListGetDerivValue( Material, &
                'Specific Enthalpy', n,Element % NodeIndexes )
            HeatCapacity(1:n) = Density(1:n) * HeatCapacity(1:n)
          ELSE
            HeatCapacity(1:n) = ListGetDerivValue( Material, &
                'Enthalpy', n,Element % NodeIndexes )            
          END IF
        END IF
          
!---------------------------------------------------------------------------------------
! Note that for the 'spatial 2' model the evaluation of c_p is done in each integration
! point and thus Enthalphy and PhaseSpatial flag are used instead of HeatCapacity directly.
!-----------------------------------------------------------------------------------------
      CASE( PHASE_SPATIAL_2 )
        IF( Specific ) THEN
          Enthalpy(1:n) = ListGetReal(Material,'Specific Enthalpy',n,Element % NodeIndexes)
          Enthalpy(1:n) = Density(1:n) * Enthalpy(1:n)
        ELSE
          Enthalpy(1:n) = ListGetReal(Material,'Enthalpy',n,Element % NodeIndexes)          
        END IF
          
!------------------------------------------------------------------------------
      CASE( PHASE_TEMPORAL )
        ! When retrieving the value of enthalphy on the previous timestep 
        ! the relevant entries of the Temperature solution in the global vector
        ! are tampered in order to make the ListGetReal command work as wanted. 
        ! 1) Values at current temperature     
        !------------------------------------------------------------------------
        IF( Specific ) THEN
          Work(1:n) = ListGetReal( Material,'Specific Enthalpy',n,Element % NodeIndexes )
        ELSE
          Work(1:n) = ListGetReal( Material,'Enthalpy',n,Element % NodeIndexes )
        END IF

        ! 2) Values at previous temperature
        Temperature(TempPerm(Element % NodeIndexes)) = & 
            PrevTemperature(TempPerm(Element % NodeIndexes)) 

        IF( Specific ) THEN
          Work(1:n) = Work(1:n) - ListGetReal( Material,'Specific Enthalpy', &
              n,Element % NodeIndexes )          
          HeatCapacity(1:n) = Density(1:n) * Work(1:n) / dT(1:n)
        ELSE
          Work(1:n) = Work(1:n) - ListGetReal( Material,'Enthalpy', &
              n,Element % NodeIndexes )
          HeatCapacity(1:n) = Work(1:n) / dT(1:n)
        END IF

        ! Revert to current temperature
        Temperature(TempPerm(Element % NodeIndexes)) = & 
            PrevTemperature(TempPerm(Element % NodeIndexes)) + dT(1:n)

!------------------------------------------------------------------------------
      END SELECT
!------------------------------------------------------------------------------
    END SUBROUTINE EffectiveHeatCapacity
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    FUNCTION CheckLatentHeat() RESULT(Failure)
!------------------------------------------------------------------------------
      LOGICAL :: Failure, PhaseChange, CheckLatentHeatRelease
      INTEGER :: t, eq_id, body_id
      CHARACTER(LEN=MAX_NAME_LEN) :: PhaseModel
!------------------------------------------------------------------------------

      Failure = .FALSE.
!------------------------------------------------------------------------------
      DO t=1,Solver % Mesh % NumberOfBulkElements
!------------------------------------------------------------------------------
!       Check if this element belongs to a body where temperature 
!       has been calculated
!------------------------------------------------------------------------------
        Element => Solver % Mesh % Elements(t)

        NodeIndexes => Element % NodeIndexes
        IF ( ANY( TempPerm( NodeIndexes ) <= 0 ) ) CYCLE

        body_id = Element % Bodyid
        eq_id = ListGetInteger( Model % Bodies(body_id) % Values, &
            'Equation', minv=1, maxv=Model % NumberOfEquations )

        PhaseModel = ListGetString( Model % Equations(eq_id) % Values, &
                          'Phase Change Model',Found )

        PhaseChange = Found .AND. (PhaseModel(1:4) /= 'none')

        IF ( PhaseChange ) THEN
          CheckLatentHeatRelease = ListGetLogical(Model % Equations(eq_id) % &
                    Values, 'Check Latent Heat Release',Found )
        END IF
        IF ( .NOT. ( PhaseChange .AND. CheckLatentHeatRelease ) ) CYCLE

        n = Element % TYPE % NumberOfNodes
!------------------------------------------------------------------------------
!       Set the current element pointer in the model structure to
!       reflect the element being processed
!------------------------------------------------------------------------------
        Model % CurrentElement => Element
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!       Get element material parameters
!------------------------------------------------------------------------------
        k = ListGetInteger( Model % Bodies(body_id) % Values,'Material', &
                minv=1, maxv=Model % NumberOfMaterials )
        Material => Model % Materials(k) % Values

        PhaseChangeIntervals => ListGetConstRealArray( Material, &
                        'Phase Change Intervals' )

        DO k=1,n
          i = TempPerm( NodeIndexes(k) )
          DO j=1,SIZE(PhaseChangeIntervals,2)
            IF ( ( Temperature(i)  < PhaseChangeIntervals(1,j) .AND. &
                   PrevSolution(i) > PhaseChangeIntervals(2,j) ).OR. &
                 ( Temperature(i)  > PhaseChangeIntervals(2,j) .AND. &
                   PrevSolution(i) < PhaseChangeIntervals(1,j) )  ) THEN
              Failure = .TRUE.
              EXIT
            END IF
          END DO
          IF ( Failure ) EXIT
        END DO
        IF ( Failure ) EXIT
      END DO
!------------------------------------------------------------------------------
    END FUNCTION CheckLatentHeat
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
   SUBROUTINE IntegOverA( BoundaryMatrix, BoundaryVector, &
    LOAD, NodalAlpha, Element, n, m, Nodes )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), &
                    LOAD(:),NodalAlpha(:)

    TYPE(Nodes_t)   :: Nodes
    TYPE(Element_t) :: Element

    INTEGER :: n,  m

    REAL(KIND=dp) :: Basis(n)
    REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

    REAL(KIND=dp) :: u,v,w,s,x,y,z
    REAL(KIND=dp) :: Force,Alpha
    REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

    INTEGER :: i,t,q,p,N_Integ

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    LOGICAL :: stat
!------------------------------------------------------------------------------

    BoundaryVector = 0.0D0
    BoundaryMatrix = 0.0D0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
    DO t=1,N_Integ
        u = U_Integ(t)
        v = V_Integ(t)
        w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
        stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                  Basis,dBasisdx )

        s = SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( Nodes % x(1:n)*Basis )
         y = SUM( Nodes % y(1:n)*Basis )
         z = SUM( Nodes % z(1:n)*Basis )
         s = s * CoordinateSqrtMetric( x,y,z )
        END IF
!------------------------------------------------------------------------------
        Force = SUM( LOAD(1:n) * Basis )
        Alpha = SUM( NodalAlpha(1:n) * Basis )

        DO p=1,N
         DO q=1,M
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
                  s * Alpha * Basis(p) / m
         END DO
        END DO

        DO p=1,N
         BoundaryVector(p) = BoundaryVector(p) + s * Force * Basis(p)
        END DO
    END DO
   END SUBROUTINE IntegOverA
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE FindGapIndexes( Element, Indexes, n )
!------------------------------------------------------------------------------
      TYPE(Element_t) :: Element
      INTEGER :: n,Indexes(:)
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Parent,Left,Right
      INTEGER :: i,j,k,l
      REAL(KIND=dp) :: x0,y0,z0,x,y,z
!------------------------------------------------------------------------------
      Left  => Element % BoundaryInfo % Left
      Right => Element % BoundaryInfo % Right

      IF ( .NOT.ASSOCIATED(Left) .OR. .NOT.ASSOCIATED(Right) ) RETURN

      l = 0
      DO i=1,n
        Parent => Left
        k = Element % NodeIndexes(i)

        IF ( ANY( Parent % NodeIndexes == k ) ) &
          Parent => Right

        x0 = ElementNodes % x(i)
        y0 = ElementNodes % y(i)
        z0 = ElementNodes % z(i)
        DO j=1,Parent % TYPE % NumberOfNodes
          k = Parent % NodeIndexes(j)
          x = Solver % Mesh % Nodes % x(k) - x0
          y = Solver % Mesh % Nodes % y(k) - y0
          z = Solver % Mesh % Nodes % z(k) - z0
          IF ( x**2 + y**2 + z**2 < AEPS ) EXIT
        END DO
        Indexes(i) = k
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE FindGapIndexes
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE AddHeatGap( Solver, Element, STIFF, TempPerm )
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      REAL(KIND=dp) :: STIFF(:,:)
      INTEGER :: TempPerm(:)
      TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Parent,Left,Right
      INTEGER :: i,j,k,l, Ind(n)
      REAL(KIND=dp) :: x0,y0,z0,x,y,z
!------------------------------------------------------------------------------
      CALL FindGapIndexes( Element, Ind, n )
      DO i=1,n
        DO j=1,n
          k = TempPerm( Element % NodeIndexes(i) )
          l = TempPerm( Ind(j) )
          IF ( k > 0 .AND. l > 0 ) THEN
            CALL AddToMatrixElement( Solver % Matrix,k,l,-STIFF(i,j) )
          END IF
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE AddHeatGap
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE NumaDefaultDirichletBCs( Solver )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Set dirichlet condition on boundary and interior nodes 
!
!   ARGUMENT:
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!******************************************************************************
        TYPE(Solver_t), TARGET :: Solver
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER   :: A
        TYPE(Variable_t), POINTER :: x
        TYPE(ValueList_t), POINTER :: BC, SolverParams
        TYPE(Element_t), POINTER :: Element, Parent, Edge, Face, SaveElement
        
        REAL(KIND=dp), POINTER    :: b(:)
        REAL(KIND=dp) :: xx
        REAL(KIND=dp), ALLOCATABLE :: Work(:), STIFF(:,:)
        
        INTEGER, ALLOCATABLE :: lInd(:), gInd(:)
        INTEGER :: i,j, k, kk, l, m, n,nd, nb, mb, nn, DOF, local, numEdgeDofs,istat

        LOGICAL :: Flag,Found
        
        CHARACTER(LEN=MAX_NAME_LEN) :: name

        TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
        SAVE gInd, lInd, STIFF, Work
!------------------------------------------------------------------------------
!       Get the linear system components
!------------------------------------------------------------------------------
        A => Solver % Matrix
        x => Solver % Variable
        b => A % RHS
!------------------------------------------------------------------------------
!       Get the maximal number of DOF 
!------------------------------------------------------------------------------
        n = Solver % Mesh % MaxElementDOFs
!------------------------------------------------------------------------------
!       Allocations
!------------------------------------------------------------------------------     
        IF ( .NOT. ALLOCATED( gInd ) ) THEN
            ALLOCATE( gInd(n), lInd(n), STIFF(n,n), Work(n), stat=istat )
            IF ( istat /= 0 ) &
                CALL Fatal('DefUtils::NumaDefaultDirichletBCs','Memory allocation failed.' )
        ELSE IF ( SIZE(gInd) < n ) THEN
            DEALLOCATE( gInd, lInd, STIFF, Work )
            ALLOCATE( gInd(n), lInd(n), STIFF(n,n), Work(n), stat=istat )
            IF ( istat /= 0 ) &
                CALL Fatal('DefUtils::NumaDefaultDirichletBCs','Memory allocation failed.' )
        END IF
!------------------------------------------------------------------------------
!       Special treatment if several DOFs ?
!------------------------------------------------------------------------------
        IF ( x % DOFs > 1 ) THEN
            ! TEMP!!!
            CALL SetDirichletBoundaries( CurrentModel,A, b, x % Name,-1,x % DOFs,x % Perm )
            CALL SetInteriorDirichletConditions( CurrentModel, A, b, x % &
                Name, -1, x % DOFs, x % Perm )
        END IF
!------------------------------------------------------------------------------
!       Clear dirichlet BCs for face & edge DOFs:
!------------------------------------------------------------------------------
        DO DOF=1,x % DOFs
!------------------------------------------------------------------------------     
!           Get the name of the DOF:
!------------------------------------------------------------------------------     
            name = x % name
            IF ( x % DOFs > 1 ) name = ComponentName(name,DOF)
!------------------------------------------------------------------------------
!           Go through the boundary elements:
!------------------------------------------------------------------------------
            DO i=1,Solver % Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
                Element => GetBoundaryElement(i)
!------------------------------------------------------------------------------
!               Check if the element is active:
!------------------------------------------------------------------------------
                IF ( .NOT. ActiveBoundaryElement() ) CYCLE
!------------------------------------------------------------------------------
!               Get the BC associated to this element
!               Check if the BC exists
!------------------------------------------------------------------------------
                BC => GetBC()
                IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
                IF ( .NOT. ListCheckPresent(BC, Name) ) CYCLE 
!------------------------------------------------------------------------------
!               Get the parent element: bulk element with common face or edge.
!               Left or right
!------------------------------------------------------------------------------
                Parent => Element % BoundaryInfo % Left
                IF ( .NOT. ASSOCIATED( Parent ) ) THEN
                    Parent => Element % BoundaryInfo % Right
                END IF
                IF ( .NOT. ASSOCIATED( Parent ) ) CYCLE
!------------------------------------------------------------------------------
!               Clear dofs associated with element edges:
!------------------------------------------------------------------------------
                IF ( ASSOCIATED( Solver % Mesh % Edges ) ) THEN
!------------------------------------------------------------------------------             
                    DO j=1,Parent % TYPE % NumberOfEdges
!------------------------------------------------------------------------------                 
                        Edge => Solver % Mesh % Edges( Parent % EdgeIndexes(j) )
                        IF ( Edge % BDOFs == 0 ) CYCLE

                        n = 0
                        DO k=1,Element % TYPE % NumberOfNodes
                            DO l=1,Edge % TYPE % NumberOfNodes
                                IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) n=n+1
                            END DO
                        END DO

                        IF ( n ==  Edge % Type % NumberOfNodes ) THEN
                            DO k=1,Edge % BDOFs
                                n = Solver % Mesh % NumberofNodes + &
                                    (Parent % EdgeIndexes(j)-1)*Solver % Mesh % MaxEdgeDOFs+k
                                n = x % Perm( n )
                                IF ( n <= 0 ) CYCLE
                                n = x % DOFs*(n-1) + DOF
                                CALL CRS_ZeroRow( A, n )
                                A % RHS(n) = 0.0d0
                            END DO
                        END IF
                    END DO
                END IF
!------------------------------------------------------------------------------
!              Clear dofs associated with element faces:
!------------------------------------------------------------------------------
                IF ( ASSOCIATED( Solver % Mesh % Faces ) ) THEN
!------------------------------------------------------------------------------             
                    DO j=1,Parent % Type % NumberOfFaces
!------------------------------------------------------------------------------
                        Face => Solver % Mesh % Faces( Parent % FaceIndexes(j) )
                        IF ( Face % BDOFs == 0 ) CYCLE

                        n = 0
                        DO k=1,Element % TYPE % NumberOfNodes
                            DO l=1,Face % TYPE % NumberOfNodes
                                IF ( Face % NodeIndexes(l) == Element % NodeIndexes(k) ) n=n+1
                            END DO
                        END DO
                        IF ( n /= Face % TYPE % NumberOfNodes ) CYCLE

                        DO k=1,Face % BDOFs
                            n = Solver % Mesh % NumberofNodes + &
                            Solver % Mesh % MaxEdgeDOFs * Solver % Mesh % NumberOfEdges + &
                            (Parent % FaceIndexes(j)-1) * Solver % Mesh % MaxFaceDOFs + k
                            n = x % Perm( n )
                            IF ( n <= 0 ) CYCLE
                            n = x % DOFs*(n-1) + DOF
                            CALL CRS_ZeroRow( A, n )
                            A % RHS(n) = 0.0d0
                        END DO
                    END DO
                END IF
!------------------------------------------------------------------------------          
            END DO ! i
!------------------------------------------------------------------------------
        END DO ! DOF
!------------------------------------------------------------------------------
        CALL Info('HeatSolve: ', &
            'Setting Dirichlet boundary conditions', Level=5)
!------------------------------------------------------------------------------
!     Set Dirichlet dofs for edges and faces
!------------------------------------------------------------------------------
        DO DOF=1,x % DOFs
!------------------------------------------------------------------------------     
!           Get the name of the DOF:
!------------------------------------------------------------------------------     
            name = x % name
            IF ( x % DOFs > 1 ) name = ComponentName(name,DOF)
!------------------------------------------------------------------------------     
!           Set nodal loads:
!------------------------------------------------------------------------------     
            !CALL SetNodalLoads( CurrentModel,A, b, &
            !   Name,DOF,x % DOFs,x % Perm )
!------------------------------------------------------------------------------     
!           Set Dirichlet Conditions for Boundaries:
!------------------------------------------------------------------------------     
            CALL SetDirichletBoundaries( CurrentModel, A, b, &
                Name, DOF, x % DOFs, x % Perm )
!------------------------------------------------------------------------------
!           Set Dirichlet conditions for interior zones
!------------------------------------------------------------------------------ 
            CALL SetInteriorDirichletConditions( CurrentModel, A, b, &
                Name, DOF, x % DOFs, x % Perm )
            SaveElement => CurrentModel % CurrentElement
!------------------------------------------------------------------------------
!           Dirichlet BCs for face & edge DOFs:
!-------------------------------------------------------------------------------
            DO i=1,Solver % Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------         
                Element => GetBoundaryElement(i)
                IF ( .NOT. ActiveBoundaryElement() ) CYCLE
!------------------------------------------------------------------------------
!               Get the BC associated to this element
!               Check if the BC exists
!------------------------------------------------------------------------------
                BC => GetBC()
                IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
                IF ( .NOT. ListCheckPresent(BC, Name) ) CYCLE
!------------------------------------------------------------------------------
!               Get parent element:
!------------------------------------------------------------------------------
                Parent => Element % BoundaryInfo % Left
                IF ( .NOT. ASSOCIATED( Parent ) ) THEN
                    Parent => Element % BoundaryInfo % Right
                END IF
                IF ( .NOT. ASSOCIATED( Parent ) )   CYCLE
                IF ( .NOT. ASSOCIATED( Parent % pDefs ) ) CYCLE
!------------------------------------------------------------------------------
!               
!------------------------------------------------------------------------------
                n = Element % Type % NumberOfNodes
                DO j=1,n
                    l = Element % NodeIndexes(j)
                    Work(j)  = ListGetConstReal( BC, Name, Found, &
                        CurrentModel % Mesh % Nodes % x(l), &
                        CurrentModel % Mesh % Nodes % y(l), &
                        CurrentModel % Mesh % Nodes % z(l) )
                END DO
!------------------------------------------------------------------------------
                SELECT CASE(Parent % Type % Dimension)
!------------------------------------------------------------------------------
                    CASE(2)
!------------------------------------------------------------------------------
                    ! If no edges do not try to set boundary conditions
                    IF ( .NOT. ASSOCIATED( Solver % Mesh % Edges ) ) CYCLE

                    ! If boundary edge has no dofs move on to next edge
                    IF (Element % BDOFs <= 0) CYCLE

                    ! Number of nodes for this element
                    n = Element % TYPE % NumberOfNodes

                    ! Get indexes for boundary and values for dofs associated to them
                    CALL getBoundaryIndexes( Solver % Mesh, Element, Parent, gInd, numEdgeDofs )
                    CALL LocalBcBDOFs( BC, Element, numEdgeDofs, Name, STIFF, Work )

                    ! Contribute this boundary to global system
                    ! (i.e solve global boundary problem)
                    DO k=n+1,numEdgeDofs
                        nb = x % Perm( gInd(k) )
                        IF ( nb <= 0 ) CYCLE
                        nb = x % DOFs * (nb-1) + DOF
                        A % RHS(nb) = A % RHS(nb) + Work(k)
                        DO l=1,numEdgeDofs
                            mb = x % Perm( gInd(l) )
                            IF ( mb <= 0 ) CYCLE
                            mb = x % DOFs * (mb-1) + DOF
                            DO kk=A % Rows(nb)+DOF-1,A % Rows(nb+1)-1,x % DOFs
                                IF ( A % Cols(kk) == mb ) THEN
                                    A % Values(kk) = A % Values(kk) + STIFF(k,l)
                                    EXIT
                                END IF
                            END DO
                        END DO
                    END DO
!------------------------------------------------------------------------------                 
                    CASE(3)
!------------------------------------------------------------------------------
                    ! If no faces present do not try to set boundary conditions
                    ! @todo This should be changed to EXIT
                    IF ( .NOT. ASSOCIATED( Solver % Mesh % Faces ) ) CYCLE

                    ! Parameters of element
                    n = Element % TYPE % NumberOfNodes

                    ! Get global boundary indexes and solve dofs associated to them
                    CALL getBoundaryIndexes( Solver % Mesh, Element,  &
                        Parent, gInd, numEdgeDofs )
                    ! If boundary face has no dofs skip to next boundary element
                    IF (numEdgeDOFs == n) CYCLE

                    ! Get local solution
                    CALL LocalBcBDofs( BC, Element, numEdgeDofs, Name, STIFF, Work )

                    ! Contribute this entry to global boundary problem
                    DO k=n+1, numEdgeDOFs
                        nb = x % Perm( gInd(k) )
                        IF ( nb <= 0 ) CYCLE
                        nb = x % DOFs * (nb-1) + DOF
                        A % RHS(nb) = A % RHS(nb) + Work(k)
                        DO l=1, numEdgeDOFs
                            mb = x % Perm( gInd(l) )
                            IF ( mb <= 0 ) CYCLE
                            mb = x % DOFs * (mb-1) + DOF
                            DO kk=A % Rows(nb)+DOF-1,A % Rows(nb+1)-1,x % DOFs
                                IF ( A % Cols(kk) == mb ) THEN
                                    A % Values(kk) = A % Values(kk) + STIFF(k,l)
                                    EXIT
                                END IF
                            END DO
                        END DO
                    END DO
!------------------------------------------------------------------------------
                END SELECT
!------------------------------------------------------------------------------
        END DO ! elements
!------------------------------------------------------------------------------
        CurrentModel % CurrentElement => SaveElement
!------------------------------------------------------------------------------        
        END DO ! DOFs
!------------------------------------------------------------------------------
        CALL Info('HeatSolve: ', &
            'Dirichlet boundary conditions set', Level=5)
!------------------------------------------------------------------------------
    END SUBROUTINE NumaDefaultDirichletBCs
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
    SUBROUTINE SetInteriorDirichletConditions( Model, A, b, Name, DOF, NDOFs, Perm )
!------------------------------------------------------------------------------
!******************************************************************************
!
! Set dirichlet boundary condition for non boundary elements
!
! TYPE(Model_t) :: Model
!   INPUT: the current model structure
!
! TYPE(Matrix_t), POINTER :: A
!   INOUT: The global stiff matrix
!
! REAL(KIND=dp) :: b
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
        TYPE(Matrix_t), POINTER :: A
        REAL(KIND=dp) :: b(:)
        CHARACTER(LEN=*) :: Name 
        INTEGER :: DOF, NDOFs, Perm(:)
!------------------------------------------------------------------------------
!       Local variables
!------------------------------------------------------------------------------
        TYPE(ValueList_t), POINTER :: ValueList
        
        INTEGER, POINTER :: NodeIndexes(:)
        INTEGER, ALLOCATABLE :: IndNodes(:)

        INTEGER :: BC,i,j,n, NoNodes, NOFNodesFound, dim, Nb_Target_Spherical
        
        LOGICAL :: GotIt, NodesFound, Interior, Target_spherical, Interior1, Interior2, Stat
        REAL(KIND=dp) ::  s, min_x, min_y, max_x, max_y,min_z,max_z, dist
        REAL(KIND=dp), POINTER :: c_x_Array(:,:), c_y_Array(:,:), c_z_Array(:,:), &
            radius_Array(:,:)
!------------------------------------------------------------------------------
!       Dimension of the model
!------------------------------------------------------------------------------     
        dim=CoordinateSystemDimension()
!------------------------------------------------------------------------------
!     Go through the BCs 
!------------------------------------------------------------------------------
        DO BC=1,Model % NumberOfBCs
!------------------------------------------------------------------------------ 
!           Check if the variable is concerned by this BC
!------------------------------------------------------------------------------      
            IF( .NOT. ListCheckPresent( Model % BCs(BC) % Values,Name )) CYCLE
!------------------------------------------------------------------------------         
            NodesFound = .FALSE.
!------------------------------------------------------------------------------
!           The areas in which Dirichlet conditions have to be applied are defined either 
!           in terms of coordinates or in terms of nodes.
!           At the first calling the list of coordinates is transformed to a list of nodes 
!------------------------------------------------------------------------------
            IF(.NOT. NodesFound) THEN
!------------------------------------------------------------------------------            
                IF(.NOT. ALLOCATED(IndNodes)) ALLOCATE( IndNodes(Model % NumberOfNodes) )
                IndNodes = -1
                NoNodes=0
                GotIt = .FALSE.
!------------------------------------------------------------------------------ 
!               Check if the target zone is spherical 
!------------------------------------------------------------------------------                 
                Target_spherical = GetLogical( Model % BCs(BC) % Values,'Target Spherical',GotIt )      
                IF ( .NOT. GotIt ) Target_spherical = .FALSE.
                IF(Target_spherical) THEN
!------------------------------------------------------------------------------ 
!                   Check if several spherical zones are treated:
!------------------------------------------------------------------------------ 
                    Nb_Target_Spherical = GetInteger( Model % BCs(BC) % Values,'Nb Target Spherical',GotIt )        
                    IF ( .NOT. GotIt ) THEN 
                        Nb_Target_Spherical = 1
                    ELSE
                        IF ( Nb_Target_Spherical < 1 ) THEN 
                            PRINT*,'Nb Target Spherical = ',Nb_Target_Spherical
                            CALL Fatal( 'SetInteriorDirichletConditions', 'Check Nb Target Spherical!' )
                        END IF
                    END IF
!------------------------------------------------------------------------------ 
!                   Read the centre coordinates and radius for spherical areas in which a 
!                   Dirichlet BC has to be set
!                   These areas are then defined by
!                   (x-x_c)^2 + (y-y_c)^2 + (z-z_c)^2 < radius^2
!------------------------------------------------------------------------------                     
                    min_x=1
                    max_x=0
                    min_y=1
                    max_y=0
                    min_z=1
                    max_z=0
!------------------------------------------------------------------------------                                 
                    c_x_Array => ListGetConstRealArray( Model % BCs(BC) % Values, 'Target CentreX', GotIt )
                    IF (.NOT. GotIt) CYCLE
                    
                    c_y_Array => ListGetConstRealArray( Model % BCs(BC) % Values, 'Target CentreY', GotIt )
                    IF (.NOT. GotIt) CYCLE

                    radius_Array => ListGetConstRealArray( Model % BCs(BC) % Values, 'Target Radius', GotIt )
                    IF (.NOT. GotIt) CYCLE

                    IF (dim>2) THEN
                        c_z_Array => ListGetConstRealArray( Model % BCs(BC) % Values, 'Target CentreZ', GotIt )
                        IF (.NOT. GotIt) CYCLE
                    END IF  
 
!------------------------------------------------------------------------------                     
                ELSE
!------------------------------------------------------------------------------ 
!                   Read the upper/lower bounds of coordinates for areas in which a 
!                   Dirichlet BC has to be set
!                   These areas are then defined by
!                   min_x < x < max_x, min_y < y < max_y, min_z < z < max_z
!------------------------------------------------------------------------------                                             
                    min_x =  ListGetConstReal(  Model % BCs(BC) % Values, 'Target min x',GotIt)
                    IF (.NOT. GotIt) CYCLE
                    min_y =  ListGetConstReal(  Model % BCs(BC) % Values, 'Target min y',GotIt)
                    IF (.NOT. GotIt) CYCLE
                    max_x =  ListGetConstReal(  Model % BCs(BC) % Values, 'Target max x',GotIt)
                    IF (.NOT. GotIt) CYCLE
                    max_y =  ListGetConstReal(  Model % BCs(BC) % Values, 'Target max y',GotIt)
                    IF (.NOT. GotIt) CYCLE
                    IF (dim>2) THEN
                        min_z =  ListGetConstReal(  Model % BCs(BC) % Values, 'Target min z',GotIt)
                        IF (.NOT. GotIt) CYCLE
                        max_z =  ListGetConstReal(  Model % BCs(BC) % Values, 'Target max z',GotIt)
                        IF (.NOT. GotIt) CYCLE
                    END IF
!------------------------------------------------------------------------------
                END IF
!------------------------------------------------------------------------------ 
!               Go through the nodes of the model and check for each if it belongs 
!               to an area with Dirichlet BC
!------------------------------------------------------------------------------             
                DO i=1,Model % NumberOfNodes
!------------------------------------------------------------------------------ 
!                       Check if the element is active 
!------------------------------------------------------------------------------    
                    IF( Perm(i) == 0) CYCLE
                    Interior = .FALSE.
!------------------------------------------------------------------------------ 
!                   Comparisons between node coordinates and area bounds in spherical case
!------------------------------------------------------------------------------     
                    Interior1 = .FALSE.
!------------------------------------------------------------------------------
                    IF(Target_spherical) THEN
!------------------------------------------------------------------------------
                        DO j=1,Nb_Target_Spherical
                            dist = (Model % Mesh % Nodes % x(i)-c_x_Array(j,1))**2 + &
                                (Model % Mesh % Nodes % y(i)-c_y_Array(j,1))**2 
                            IF (dim>2) THEN 
                                dist = dist + (Model % Mesh % Nodes % z(i)-c_z_Array(j,1))**2 
                            END IF
                            Interior1= Interior1 .OR. (dist<=radius_Array(j,1)**2)
                        END DO
!------------------------------------------------------------------------------
                    END IF
!------------------------------------------------------------------------------ 
!                   Comparisons between node coordinates and area bounds in non-spherical case
!------------------------------------------------------------------------------    
                    Interior2= (Model % Mesh % Nodes % x(i)>min_x) .AND. &
                        (Model % Mesh % Nodes % x(i)<max_x)
                    Interior2= Interior2 .AND. (Model % Mesh % Nodes % y(i)>min_y) .AND. &
                        (Model % Mesh % Nodes % y(i)<max_y)
                    IF (dim>2) THEN
                        Interior2= Interior2 .AND. (Model % Mesh % Nodes % z(i)>min_z) .AND. &
                        (Model % Mesh % Nodes % z(i)<max_z)
                    END IF
!------------------------------------------------------------------------------ 
                    Interior = Interior1 .OR. Interior2                 
!------------------------------------------------------------------------------ 
!                   Update the number of nodes for the BC, and the vector of these nodes
!------------------------------------------------------------------------------
                    IF( Interior ) THEN
                        NoNodes = NoNodes+1
                        IndNodes(NoNodes) = i
                    END IF
!------------------------------------------------------------------------------
                END DO ! Model % NumberOfNodes
!--------------------------------------------------------------------------     
!               Check if all the selected nodes are active
!------------------------------------------------------------------------------         
                NOFNodesFound = 0
                DO j=1,NoNodes
                    IF ( IndNodes(j)>0 ) THEN
                        NOFNodesFound=NOFNodesFound+1
                        IndNodes(NOFNodesFound) = IndNodes(j)
                    END IF
                END DO
!------------------------------------------------------------------------------ 
!               In the first time add the found nodes to the list structure
!------------------------------------------------------------------------------ 
                IF ( NOFNodesFound > 0 ) THEN
                    CALL ListAddIntegerArray( Model % BCs(BC) % Values,'Target Nodes', &
                        NOFNodesFound, IndNodes) 
                    DEALLOCATE(IndNodes)
                    NodesFound = .TRUE.               
                END IF  
!------------------------------------------------------------------------------
            END IF ! NOT NodesFound
!------------------------------------------------------------------------------
!           If the nodes are specified in the input file, or if the coordinates 
!        have been transformed in nodes, we apply the conditions:           
!------------------------------------------------------------------------------
            IF(NodesFound) THEN          
!------------------------------------------------------------------------------     
!               Read the nodes in the input file
!------------------------------------------------------------------------------             
                NodeIndexes => ListGetIntegerArray( Model % BCs(BC) % Values,'Target Nodes')
!------------------------------------------------------------------------------     
!               Get the number of nodes concerned be Dirichlet BC
!------------------------------------------------------------------------------                     
                n = SIZE(NodeIndexes)
!------------------------------------------------------------------------------     
!               Get the fixed values of the Dirichlet BC
!------------------------------------------------------------------------------     
                ValueList => Model % BCs(BC) % Values
!------------------------------------------------------------------------------     
!               Modify the system with the Dirichlet BC
!------------------------------------------------------------------------------     
                CALL SetPointValues(A, b,Name, DOF,NDOFs,ValueList,n,NodeIndexes,Perm)
!------------------------------------------------------------------------------
            END IF ! NodesFound
!------------------------------------------------------------------------------
        END DO ! BC
!------------------------------------------------------------------------------
  END SUBROUTINE SetInteriorDirichletConditions
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------  
    SUBROUTINE SetPointValues(A,b,Name,DOF,NDOFs,ValueList,n,NodeIndexes,Perm)
!------------------------------------------------------------------------------ 
!******************************************************************************
!
! Set values related to individual points
!
! ARGUMENTS:
!
! TYPE(Matrix_t), POINTER :: A
!   INOUT: The global stiff matrix
!
! REAL(KIND=dp) :: b
!   INOUT: The global RHS vector
! 
! CHARACTER(LEN=*) :: Name
!   INPUT: name of the dof to be set
!
! INTEGER :: DOF, NDOFs
!   INPUT: The order number of the dof and the total number of DOFs for
!          this equation
!
! TYPE(ValueList_t), POINTER :: ValueList
! INPUT: Values to be set as Dirichlet BC
!
! INTEGER :: n
!   INPUT: Number of entries to be modified
!
! INTEGER :: NodeIndexes(:)
!   INPUT: List indexes of nodes modified
!
! INTEGER :: Perm(:)
!   INPUT: The node reordering info, this has been generated at the
!          beginning of the simulation for bandwidth optimization
!******************************************************************************   
        TYPE(Matrix_t), POINTER :: A
        REAL(KIND=dp) :: b(:)
        CHARACTER(LEN=*) :: Name 
        INTEGER :: n,DOF,NDOFs, Perm(:),NodeIndexes(:)
        TYPE(ValueList_t), POINTER :: ValueList
!------------------------------------------------------------------------------ 
!       Local variables
!------------------------------------------------------------------------------         
        REAL(KIND=dp) :: Work(n)        
        INTEGER :: j,k,k1,l
        REAL(KIND=dp), POINTER :: WorkA(:,:,:) => NULL()
!------------------------------------------------------------------------------
!       Get the nodes indexes in Work or WorkA in function of DOF order
!------------------------------------------------------------------------------             
        IF ( DOF > 0 ) THEN
            Work(1:n)  = ListGetReal( ValueList, Name, n, NodeIndexes, gotIt )
        ELSE
            CALL ListGetRealArray( ValueList, Name, WorkA, n, NodeIndexes, gotIt )
        END IF
!------------------------------------------------------------------------------
        IF ( gotIt ) THEN
!------------------------------------------------------------------------------
!           Go through the nodes indexes in Work or WorkA 
!------------------------------------------------------------------------------
            DO j=1,n
!------------------------------------------------------------------------------
!               Check if the index is valid
!------------------------------------------------------------------------------
                IF ( NodeIndexes(j) > SIZE(Perm) .OR. NodeIndexes(j) < 1 ) THEN
                    CALL Warn('SetDirichletBoundaries','Invalid Node Number')
                    CYCLE
                END IF
!------------------------------------------------------------------------------
!               Check if the index is valid
!------------------------------------------------------------------------------
                k = Perm(NodeIndexes(j))
                IF ( k > 0 ) THEN
!------------------------------------------------------------------------------
!                   DOF>0
!------------------------------------------------------------------------------                 
                    IF ( DOF>0 ) THEN
                        k = NDOFs * (k-1) + DOF
!------------------------------------------------------------------------------
!                       Band matrix
!------------------------------------------------------------------------------                         
                        IF ( A % FORMAT == MATRIX_SBAND ) THEN
                            CALL SBand_SetDirichlet( A,b,k,Work(j) )
!------------------------------------------------------------------------------
!                       CRS & symmetric matrix
!------------------------------------------------------------------------------
                        ELSE IF ( A % FORMAT == MATRIX_CRS .AND. A % Symmetric ) THEN 
                            CALL CRS_SetSymmDirichlet( A,b,k,Work(j) )
!------------------------------------------------------------------------------
!                       General case 
!------------------------------------------------------------------------------                         
                        ELSE
                            b(k) = Work(j)
                            CALL ZeroRow( A,k )
                            CALL SetMatrixElement( A,k,k,1.0d0 )
                        END IF
                    ELSE
!------------------------------------------------------------------------------
!                       DOF<=0
!------------------------------------------------------------------------------     
                        DO l=1,MIN( NDOFs, SIZE(Worka,1) )
                            k1 = NDOFs * (k-1) + l
!------------------------------------------------------------------------------
!                           Band matrix
!------------------------------------------------------------------------------     
                            IF ( A % FORMAT == MATRIX_SBAND ) THEN
                                CALL SBand_SetDirichlet( A,b,k1,WorkA(l,1,j) )
!------------------------------------------------------------------------------
!                           CRS & symmetric matrix
!------------------------------------------------------------------------------
                            ELSE IF ( A % FORMAT == MATRIX_CRS .AND. A % Symmetric ) THEN
                                CALL CRS_SetSymmDirichlet( A,b,k1,WorkA(l,1,j) )
!------------------------------------------------------------------------------
!                           General case 
!------------------------------------------------------------------------------
                            ELSE
                                b(k1) = WorkA(l,1,j)
                                CALL ZeroRow( A,k1 )
                                CALL SetMatrixElement( A,k1,k1,1.0d0 )
                            END IF
                        END DO ! l
!------------------------------------------------------------------------------                     
                    END IF ! DOF>0
!------------------------------------------------------------------------------
                END IF ! k>0
!------------------------------------------------------------------------------             
            END DO ! j
!------------------------------------------------------------------------------         
        END IF ! gotIt
!------------------------------------------------------------------------------
    END SUBROUTINE SetPointValues
!------------------------------------------------------------------------------










!------------------------------------------------------------------------------
  END SUBROUTINE HeatSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION HeatBoundaryResidual( Model, Edge, Mesh, Quant, Perm,Gnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
    USE DefUtils
    USE Radiation

    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: Perm(:)
    REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
    TYPE( Mesh_t ), POINTER    :: Mesh
    TYPE( Element_t ), POINTER :: Edge
!------------------------------------------------------------------------------

    TYPE(Nodes_t) :: Nodes, EdgeNodes
    TYPE(Element_t), POINTER :: Element, Bndry

    INTEGER :: i,j,k,n,l,t,DIM,Pn,En
    LOGICAL :: stat, Found

    REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

    REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

    REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), ExtTemperature(:), &
        TransferCoeff(:), EdgeBasis(:), Basis(:), x(:), y(:), z(:), &
        dBasisdx(:,:), Temperature(:), Flux(:), NodalEmissivity(:)

    REAL(KIND=dp) :: Conductivity, Emissivity, StefanBoltzmann

    REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength, gx, gy, gz

    REAL(KIND=dp) :: u, v, w, s, detJ

    REAL(KIND=dp) :: Source, Residual, ResidualNorm, Area

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    LOGICAL :: First = .TRUE., Dirichlet
    SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
    IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
    END IF

    Indicator = 0.0d0
    Gnorm     = 0.0d0

    Metric = 0.0d0
    DO i=1,3
        Metric(i,i) = 1.0d0
    END DO

    SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
    END SELECT
!    
!    ---------------------------------------------

    Element => Edge % BoundaryInfo % Left

    IF ( .NOT. ASSOCIATED( Element ) ) THEN
        Element => Edge % BoundaryInfo % Right
    ELSE IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) THEN
        Element => Edge % BoundaryInfo % Right
    END IF

    IF ( .NOT. ASSOCIATED( Element ) ) RETURN
    IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

    En = Edge % TYPE % NumberOfNodes
    Pn = Element % TYPE % NumberOfNodes

    ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )

    EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
    EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
    EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

    ALLOCATE( Nodes % x(Pn), Nodes % y(Pn), Nodes % z(Pn) )

    Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
    Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
    Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

    ALLOCATE( Temperature(Pn), Basis(Pn), ExtTemperature(En), &
        TransferCoeff(En), x(En), y(En), z(En), EdgeBasis(En), &
        dBasisdx(Pn,3), NodalConductivity(En), Flux(En), &
        NodalEmissivity(En) ) 

    DO l = 1,En
        DO k = 1,Pn
          IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) THEN
             x(l) = Element % TYPE % NodeU(k)
             y(l) = Element % TYPE % NodeV(k)
             z(l) = Element % TYPE % NodeW(k)
             EXIT
          END IF
        END DO
    END DO
!
!    Integrate square of residual over boundary element:
!    ---------------------------------------------------

    Indicator    = 0.0d0
    EdgeLength   = 0.0d0
    ResidualNorm = 0.0d0

    DO j=1,Model % NumberOfBCs
        IF ( Edge % BoundaryInfo % Constraint /= Model % BCs(j) % Tag ) CYCLE

!       IF ( .NOT. ListGetLogical( Model % BCs(j) % Values, &
!                 'Heat Flux BC', Found ) ) CYCLE

!
!       Check if dirichlet BC given:
!       ----------------------------
        s = ListGetConstReal( Model % BCs(j) % Values,'Temperature',Dirichlet )

!       Get various flux bc options:
!       ----------------------------

!       ...given flux:
!       --------------
        Flux(1:En) = ListGetReal( Model % BCs(j) % Values, &
          'Heat Flux', En, Edge % NodeIndexes, Found )

!       ...convective heat transfer:
!       ----------------------------
        TransferCoeff(1:En) =  ListGetReal( Model % BCs(j) % Values, &
          'Heat Transfer Coefficient', En, Edge % NodeIndexes, Found )

        ExtTemperature(1:En) = ListGetReal( Model % BCs(j) % Values, &
          'External Temperature', En, Edge % NodeIndexes, Found )

!       ...black body radiation:
!       ------------------------
        Emissivity      = 0.0d0
        StefanBoltzmann = 0.0d0

        SELECT CASE(ListGetString(Model % BCs(j) % Values,'Radiation',Found))
           !------------------
           CASE( 'idealized' )
           !------------------

              NodalEmissivity(1:En) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', En, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:En) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:En)) / En

              StefanBoltzMann = &
                    ListGetConstReal( Model % Constants,'Stefan Boltzmann' )

           !---------------------
           CASE( 'diffuse gray' )
           !---------------------

              NodalEmissivity(1:En) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', En, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:En) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:En)) / En

              StefanBoltzMann = &
                    ListGetConstReal( Model % Constants,'Stefan Boltzmann' )

              ExtTemperature(1:En) =  ComputeRadiationLoad( Model, &
                      Mesh, Edge, Quant, Perm, Emissivity )
        END SELECT

!       get material parameters:
!       ------------------------
        k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material', &
                    minv=1, maxv=Model % NumberOFMaterials)

        CALL ListGetRealArray( Model % Materials(k) % Values, &
               'Heat Conductivity', Hwrk, En, Edge % NodeIndexes )

        NodalConductivity( 1:En ) = Hwrk( 1,1,1:En )

!       elementwise nodal solution:
!       ---------------------------
        Temperature(1:Pn) = Quant( Perm(Element % NodeIndexes) )

!       do the integration:
!       -------------------
        EdgeLength   = 0.0d0
        ResidualNorm = 0.0d0

        IntegStuff = GaussPoints( Edge )

        DO t=1,IntegStuff % n
           u = IntegStuff % u(t)
           v = IntegStuff % v(t)
           w = IntegStuff % w(t)

           stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
               EdgeBasis, dBasisdx )

           Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              s = IntegStuff % s(t) * detJ
           ELSE
              gx = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
              gy = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
              gz = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )
      
              CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                         Symb, dSymb, gx, gy, gz )

              s = IntegStuff % s(t) * detJ * SqrtMetric
           END IF

!
!          Integration point in parent element local
!          coordinates:
!          -----------------------------------------
           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
                 Basis, dBasisdx )
!
!          Heat conductivity at the integration point:
!          --------------------------------------------
           Conductivity = SUM( NodalConductivity(1:En) * EdgeBasis(1:En) )
!
!          given flux at integration point:
!          --------------------------------
           Residual = -SUM( Flux(1:En) * EdgeBasis(1:En) )

!          convective ...:
!          ----------------
           Residual = Residual + SUM(TransferCoeff(1:En) * EdgeBasis(1:En)) * &
                     ( SUM( Temperature(1:Pn) * Basis(1:Pn) ) - &
                       SUM( ExtTemperature(1:En) * EdgeBasis(1:En) ) )

!          black body radiation...:
!          -------------------------
           Residual = Residual + &
                Emissivity * StefanBoltzmann * &
                     ( SUM( Temperature(1:Pn) * Basis(1:Pn) ) ** 4 - &
                       SUM( ExtTemperature(1:En) * EdgeBasis(1:En) ) ** 4 )

!          flux given by the computed solution, and 
!          force norm for scaling the residual:
!          -----------------------------------------
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO k=1,DIM
                 Residual = Residual + Conductivity  * &
                    SUM( dBasisdx(1:Pn,k) * Temperature(1:Pn) ) * Normal(k)

                 Gnorm = Gnorm + s * (Conductivity * &
                       SUM(dBasisdx(1:Pn,k) * Temperature(1:Pn)) * Normal(k))**2
              END DO
           ELSE
              DO k=1,DIM
                 DO l=1,DIM
                    Residual = Residual + Metric(k,l) * Conductivity  * &
                       SUM( dBasisdx(1:Pn,k) * Temperature(1:Pn) ) * Normal(l)

                    Gnorm = Gnorm + s * (Metric(k,l) * Conductivity * &
                      SUM(dBasisdx(1:Pn,k) * Temperature(1:Pn) ) * Normal(l))**2
                 END DO
              END DO
           END IF

           EdgeLength   = EdgeLength + s
           IF ( .NOT. Dirichlet ) THEN
              ResidualNorm = ResidualNorm + s * Residual ** 2
           END IF
        END DO
        EXIT
    END DO

    IF ( CoordinateSystemDimension() == 3 ) THEN
        EdgeLength = SQRT(EdgeLength)
    END IF

!    Gnorm = EdgeLength * Gnorm
    Indicator = EdgeLength * ResidualNorm

    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
    DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

    DEALLOCATE( Temperature, Basis, ExtTemperature, TransferCoeff,  &
        x, y, z, EdgeBasis, dBasisdx, NodalConductivity, Flux, &
        NodalEmissivity ) 
!------------------------------------------------------------------------------
  END FUNCTION HeatBoundaryResidual
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION HeatEdgeResidual( Model, Edge, Mesh, Quant, Perm ) RESULT( Indicator )
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t) :: Model
    INTEGER :: Perm(:)
    REAL(KIND=dp) :: Quant(:), Indicator(2)
    TYPE( Mesh_t ), POINTER    :: Mesh
    TYPE( Element_t ), POINTER :: Edge
!------------------------------------------------------------------------------

    TYPE(Nodes_t) :: Nodes, EdgeNodes
    TYPE(Element_t), POINTER :: Element, Bndry

    INTEGER :: i,j,k,l,n,t,DIM,En,Pn
    LOGICAL :: stat, Found
    REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

    REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

    REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), x(:), y(:), z(:), &
            EdgeBasis(:), Basis(:), dBasisdx(:,:), Temperature(:)

    REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength, Jump, Conductivity

    REAL(KIND=dp) :: u, v, w, s, detJ

    REAL(KIND=dp) :: Residual, ResidualNorm, Area

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    LOGICAL :: First = .TRUE.
    SAVE Hwrk, First
!------------------------------------------------------------------------------

    !    Initialize:
    !    -----------
    IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
    END IF

    SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
    END SELECT

    Metric = 0.0d0
    DO i = 1,3
        Metric(i,i) = 1.0d0
    END DO

    Grad = 0.0d0
!
!    ---------------------------------------------

    Element => Edge % BoundaryInfo % Left
    n = Element % TYPE % NumberOfNodes

    Element => Edge % BoundaryInfo % Right
    n = MAX( n, Element % TYPE % NumberOfNodes )

    ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

    En = Edge % TYPE % NumberOfNodes
    ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )

    EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
    EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
    EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

    ALLOCATE( NodalConductivity(En), EdgeBasis(En), Basis(n), &
        dBasisdx(n,3), x(En), y(En), z(En), Temperature(n) )

!    Integrate square of jump over edge:
!    -----------------------------------
    ResidualNorm = 0.0d0
    EdgeLength   = 0.0d0
    Indicator    = 0.0d0

    IntegStuff = GaussPoints( Edge )

    DO t=1,IntegStuff % n

        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
             EdgeBasis, dBasisdx )

        Normal = NormalVector( Edge, EdgeNodes, u, v, .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
           v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
           w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                      Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        ! 
        ! Compute flux over the edge as seen by elements
        ! on both sides of the edge:
        ! ----------------------------------------------
        DO i = 1,2
           SELECT CASE(i)
              CASE(1)
                 Element => Edge % BoundaryInfo % Left
              CASE(2)
                 Element => Edge % BoundaryInfo % Right
           END SELECT
!
!          Can this really happen (maybe it can...)  ?      
!          -------------------------------------------
           IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE
!
!          Next, get the integration point in parent
!          local coordinates:
!          -----------------------------------------
           Pn = Element % TYPE % NumberOfNodes

           DO j = 1,En
              DO k = 1,Pn
                 IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
                    x(j) = Element % TYPE % NodeU(k)
                    y(j) = Element % TYPE % NodeV(k)
                    z(j) = Element % TYPE % NodeW(k)
                    EXIT
                 END IF
              END DO
           END DO

           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )
!
!          Get parent element basis & derivatives at the integration point:
!          -----------------------------------------------------------------
           Nodes % x(1:Pn) = Mesh % Nodes % x(Element % NodeIndexes)
           Nodes % y(1:Pn) = Mesh % Nodes % y(Element % NodeIndexes)
           Nodes % z(1:Pn) = Mesh % Nodes % z(Element % NodeIndexes)

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
             Basis, dBasisdx )
!
!          Material parameters:
!          --------------------
           k = ListGetInteger( Model % Bodies( &
                    Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOFMaterials )

           CALL ListGetRealArray( Model % Materials(k) % Values, &
                   'Heat Conductivity', Hwrk,En, Edge % NodeIndexes )

           NodalConductivity( 1:En ) = Hwrk( 1,1,1:En )
           Conductivity = SUM( NodalConductivity(1:En) * EdgeBasis(1:En) )
!
!          Temperature at element nodal points:
!          ------------------------------------
           Temperature(1:Pn) = Quant( Perm(Element % NodeIndexes) )
!
!          Finally, the flux:
!          ------------------
           DO j=1,DIM
              Grad(j,i) = Conductivity * SUM( dBasisdx(1:Pn,j) * Temperature(1:Pn) )
           END DO
        END DO

!       Compute squre of the flux jump:
!       -------------------------------   
        EdgeLength  = EdgeLength + s
        Jump = 0.0d0
        DO k=1,DIM
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              Jump = Jump + (Grad(k,1) - Grad(k,2)) * Normal(k)
           ELSE
              DO l=1,DIM
                 Jump = Jump + &
                       Metric(k,l) * (Grad(k,1) - Grad(k,2)) * Normal(l)
              END DO
           END IF
        END DO
        ResidualNorm = ResidualNorm + s * Jump ** 2
    END DO

    IF ( CoordinateSystemDimension() == 3 ) THEN
        EdgeLength = SQRT(EdgeLength)
    END IF
    Indicator = EdgeLength * ResidualNorm

    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, x, y, z)
    DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)
    DEALLOCATE( NodalConductivity, EdgeBasis, Basis, dBasisdx, Temperature)

!------------------------------------------------------------------------------
  END FUNCTION HeatEdgeResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION HeatInsideResidual( Model, Element, Mesh, &
        Quant, Perm, Fnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
    USE DefUtils
!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: Perm(:)
    REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
    TYPE( Mesh_t ), POINTER    :: Mesh
    TYPE( Element_t ), POINTER :: Element
!------------------------------------------------------------------------------

    TYPE(Nodes_t) :: Nodes

    INTEGER :: i,j,k,l,n,t,DIM

    LOGICAL :: stat, Found, Compressible
    TYPE( Variable_t ), POINTER :: Var

    REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

    REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

    REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
    REAL(KIND=dp), ALLOCATABLE :: NodalCapacity(:)
    REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
    REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:)
    REAL(KIND=dp), ALLOCATABLE :: Velo(:,:), Pressure(:)
    REAL(KIND=dp), ALLOCATABLE :: NodalSource(:), Temperature(:), PrevTemp(:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), ddBasisddx(:,:,:)

    REAL(KIND=dp) :: u, v, w, s, detJ, Density, Capacity

    REAL(KIND=dp) :: SpecificHeatRatio, ReferencePressure, dt
    REAL(KIND=dp) :: Source, Residual, ResidualNorm, Area, Conductivity

    TYPE( ValueList_t ), POINTER :: Material

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    LOGICAL :: First = .TRUE.
    SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
    Indicator = 0.0d0
    Fnorm     = 0.0d0
!
!    Check if this eq. computed in this element:
!    -------------------------------------------
    IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

    IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
    END IF

    Metric = 0.0d0
    DO i=1,3
        Metric(i,i) = 1.0d0
    END DO

    SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
    END SELECT
!
!    Element nodal points:
!    ---------------------
    n = Element % TYPE % NumberOfNodes

    ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

    Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
    Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
    Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

    ALLOCATE( NodalDensity(n), NodalCapacity(n), NodalConductivity(n),       &
         Velo(3,n), Pressure(n), NodalSource(n), Temperature(n), PrevTemp(n), &
         Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3) )
!
!    Elementwise nodal solution:
!    ---------------------------
    Temperature(1:n) = Quant( Perm(Element % NodeIndexes) )
!
!    Check for time dep.
!    -------------------
    PrevTemp(1:n) = Temperature(1:n)
    dt = Model % Solver % dt
    IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Var => VariableGet( Model % Variables, 'Temperature', .TRUE. )
        PrevTemp(1:n) = Var % PrevValues(Var % Perm(Element % NodeIndexes),1)
    END IF
!
!    Material parameters: conductivity, heat capacity and density
!    -------------------------------------------------------------
    k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOfMaterials )

    Material => Model % Materials(k) % Values

    CALL ListGetRealArray( Material, &
                  'Heat Conductivity', Hwrk,n, Element % NodeIndexes )

    NodalConductivity( 1:n ) = Hwrk( 1,1,1:n )

    NodalDensity(1:n) = ListGetReal( Material, &
            'Density', n, Element % NodeIndexes, Found )

    NodalCapacity(1:n) = ListGetReal( Material, &
          'Heat Capacity', n, Element % NodeIndexes, Found )
!
!    Check for compressible flow equations:
!    --------------------------------------
    Compressible = .FALSE.

    IF (  ListGetString( Material, 'Compressibility Model', Found ) == &
                 'perfect gas equation 1' ) THEN

        Compressible = .TRUE.

        Pressure = 0.0d0
        Var => VariableGet( Mesh % Variables, 'Pressure', .TRUE. )
        IF ( ASSOCIATED( Var ) ) THEN
           Pressure(1:n) = &
               Var % Values( Var % Perm(Element % NodeIndexes) )
        END IF

        ReferencePressure = ListGetConstReal( Material, &
                   'Reference Pressure' )

        SpecificHeatRatio = ListGetConstReal( Material, &
                   'Specific Heat Ratio' )

        NodalDensity(1:n) =  (Pressure(1:n) + ReferencePressure) * SpecificHeatRatio / &
              ( (SpecificHeatRatio - 1) * NodalCapacity(1:n) * Temperature(1:n) )
    END IF
!
!    Get (possible) convection velocity at the nodes of the element:
!    ----------------------------------------------------------------
    k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
                minv=1, maxv=Model % NumberOFEquations )

    Velo = 0.0d0
    SELECT CASE( ListGetString( Model % Equations(k) % Values, &
                         'Convection', Found ) )

        !-----------------
        CASE( 'constant' )
        !-----------------

           Velo(1,1:n) = ListGetReal( Material, &
              'Convection Velocity 1', n, Element % NodeIndexes, Found )

           Velo(2,1:n) = ListGetReal( Material, &
              'Convection Velocity 2', n, Element % NodeIndexes, Found )

           Velo(3,1:n) = ListGetReal( Material, &
              'Convection Velocity 3', n, Element % NodeIndexes, Found )

        !-----------------
        CASE( 'computed' )
        !-----------------

           Var => VariableGet( Mesh % Variables, 'Velocity 1', .TRUE. )
           IF ( ASSOCIATED( Var ) ) THEN
              IF ( ALL( Var % Perm( Element % NodeIndexes ) > 0 ) ) THEN
                 Velo(1,1:n) = Var % Values(Var % Perm(Element % NodeIndexes))
   
                 Var => VariableGet( Mesh % Variables, 'Velocity 2', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(2,1:n) = Var % Values( &
                              Var % Perm(Element % NodeIndexes ) )
   
                 Var => VariableGet( Mesh % Variables, 'Velocity 3', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(3,1:n) = Var % Values( &
                             Var % Perm( Element % NodeIndexes ) )
              END IF
           END IF

    END SELECT

!
!    Heat source:
!    ------------
!
    k = ListGetInteger( &
         Model % Bodies(Element % BodyId) % Values,'Body Force',Found, &
                 1, Model % NumberOFBodyForces)

    NodalSource = 0.0d0
    IF ( Found .AND. k > 0  ) THEN
        NodalSource(1:n) = ListGetReal( Model % BodyForces(k) % Values, &
               'Heat Source', n, Element % NodeIndexes, Found )
    END IF

!
!    Integrate square of residual over element:
!    ------------------------------------------

    ResidualNorm = 0.0d0
    Area = 0.0d0

    IntegStuff = GaussPoints( Element )

    DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx, ddBasisddx, .TRUE., .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( Basis(1:n) * Nodes % x(1:n) )
           v = SUM( Basis(1:n) * Nodes % y(1:n) )
           w = SUM( Basis(1:n) * Nodes % z(1:n) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Capacity     = SUM( NodalCapacity(1:n) * Basis(1:n) )
        Density      = SUM( NodalDensity(1:n) * Basis(1:n) )
        Conductivity = SUM( NodalConductivity(1:n) * Basis(1:n) )
!
!       Residual of the convection-diffusion (heat) equation:
!        R = \rho * c_p * (@T/@t + u.grad(T)) - &
!            div(C grad(T)) + p div(u) - h,
!       ---------------------------------------------------
!
!       or more generally:
!
!        R = \rho * c_p * (@T/@t + u^j T_{,j}) - &
!          g^{jk} (C T_{,j}}_{,k} + p div(u) - h
!       ---------------------------------------------------
!
        Residual = -Density * SUM( NodalSource(1:n) * Basis(1:n) )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           DO j=1,DIM
!
!             - grad(C).grad(T):
!             --------------------
!
              Residual = Residual - &
                 SUM( Temperature(1:n) * dBasisdx(1:n,j) ) * &
                 SUM( NodalConductivity(1:n) * dBasisdx(1:n,j) )

!
!             - C div(grad(T)):
!             -------------------
!
              Residual = Residual - Conductivity * &
                 SUM( Temperature(1:n) * ddBasisddx(1:n,j,j) )
           END DO
        ELSE
           DO j=1,DIM
              DO k=1,DIM
!
!                - g^{jk} C_{,k}T_{j}:
!                ---------------------
!
                 Residual = Residual - Metric(j,k) * &
                    SUM( Temperature(1:n) * dBasisdx(1:n,j) ) * &
                    SUM( NodalConductivity(1:n) * dBasisdx(1:n,k) )

!
!                - g^{jk} C T_{,jk}:
!                -------------------
!
                 Residual = Residual - Metric(j,k) * Conductivity * &
                    SUM( Temperature(1:n) * ddBasisddx(1:n,j,k) )
!
!                + g^{jk} C {_jk^l} T_{,l}:
!                ---------------------------
                 DO l=1,DIM
                    Residual = Residual + Metric(j,k) * Conductivity * &
                      Symb(j,k,l) * SUM( Temperature(1:n) * dBasisdx(1:n,l) )
                 END DO
              END DO
           END DO
        END IF

!       + \rho * c_p * (@T/@t + u.grad(T)):
!       -----------------------------------
        Residual = Residual + Density * Capacity *  &
           SUM((Temperature(1:n)-PrevTemp(1:n))*Basis(1:n)) / dt

        DO j=1,DIM
           Residual = Residual + &
              Density * Capacity * SUM( Velo(j,1:n) * Basis(1:n) ) * &
                    SUM( Temperature(1:n) * dBasisdx(1:n,j) )
        END DO


        IF ( Compressible ) THEN
!
!          + p div(u) or p u^j_{,j}:
!          -------------------------
!
           DO j=1,DIM
              Residual = Residual + &
                 SUM( Pressure(1:n) * Basis(1:n) ) * &
                      SUM( Velo(j,1:n) * dBasisdx(1:n,j) )

              IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                 DO k=1,DIM
                    Residual = Residual + &
                       SUM( Pressure(1:n) * Basis(1:n) ) * &
                           Symb(j,k,j) * SUM( Velo(k,1:n) * Basis(1:n) )
                 END DO
              END IF
           END DO
        END IF

!
!       Compute also force norm for scaling the residual:
!       -------------------------------------------------
        DO i=1,DIM
           Fnorm = Fnorm + s * ( Density * &
             SUM( NodalSource(1:n) * Basis(1:n) ) ) ** 2
        END DO

        Area = Area + s
        ResidualNorm = ResidualNorm + s *  Residual ** 2
    END DO

!    Fnorm = Element % hk**2 * Fnorm
    Indicator = Element % hK**2 * ResidualNorm

    DEALLOCATE( NodalDensity, NodalCapacity, NodalConductivity,    &
         Velo, Pressure, NodalSource, Temperature, PrevTemp, Basis, &
         dBasisdx, ddBasisddx )
!------------------------------------------------------------------------------
  END FUNCTION HeatInsideResidual
!------------------------------------------------------------------------------
