PROGRAM NQMCC_ALCF_2025
! ----------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: spi=>int32, dpi=>int64, dpf=>real64
! ----------------------------------------------------------------------
  USE NQMCC_M
  USE MPI_F08
  USE OMP_LIB
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
  TYPE(MPI_COMM)  :: COMM
  TYPE(CONTROL_T) :: PARAMS
  TYPE(PHI_T)     :: PHI
  REAL(dpf)       :: STATE,FIRST,INIT,LAST,DOT
  COMPLEX(dpf)    :: CX
  INTEGER(spi)    :: I,J,RANK,SIZE,IERROR,ROOT,MAX_IJ,PARTITION,EXTRA,ISTART,IEND,IDX
  INTEGER(dpi)    :: N,FIRST_IDX,LAST_IDX,IJ
  INTEGER(dpi),DIMENSION(:), ALLOCATABLE :: MY_FIRST,MY_LAST
  COMPLEX(dpf), DIMENSION(:,:), ALLOCATABLE :: PSIJ,PSIJ2
  COMPLEX(dpf), DIMENSION(:),   ALLOCATABLE :: YLM_PROD
! ----------------------------------------------------------------------
  CALL MPI_INIT(IERROR)
  COMM=MPI_COMM_WORLD
  ROOT=0_spi
  CALL MPI_COMM_RANK(COMM, RANK, IERROR)
  CALL MPI_COMM_SIZE(COMM, SIZE, IERROR)
!-----------------------------------------------------------------------
  FIRST=0._dpf
  INIT=0._dpf
  LAST=0._dpf
!-----------------------------------------------------------------------
  IF (RANK.EQ.ROOT) THEN
    FIRST = omp_get_wtime()
    CALL PARAMS%READ_CONTROL()
  END IF
  CALL PARAMS%BROADCAST(ROOT,COMM,IERROR)
  CALL PHI%ALLOCATE_PHI(PARAMS)
! ----------------------------------------------------------------------
  ALLOCATE(PSIJ2(PARAMS%NS,PARAMS%NT),SOURCE=CMPLX(0._dpf,0._dpf,KIND=dpf))
  ALLOCATE(PSIJ(PARAMS%NS,PARAMS%NT),SOURCE=CMPLX(0._dpf,0._dpf,KIND=dpf))
  ALLOCATE(YLM_PROD(PARAMS%NYLM),SOURCE=CMPLX(1.12345_dpf,0.6789_dpf,KIND=dpf))
! ----------------------------------------------------------------------
  STATE = PARAMS%SEED
  IF (PARAMS%RW_PHI.EQ.1) THEN
    CALL PHI%READ_PHI(PARAMS%PHI_FILE)
  ELSE
    CALL PHI%INIT_PHI(PARAMS,STATE,RANK)
  END IF
! ----------------------------------------------------------------------
  MAX_IJ=PARAMS%NS*PARAMS%NT
  EXTRA = MOD(MAX_IJ,SIZE)
  PARTITION=MAX_IJ/SIZE
  IF (EXTRA.EQ.0) THEN
    ISTART=1_spi+RANK*PARTITION
    IEND=(1_spi+RANK)*PARTITION
  ELSE
    IF (RANK.LT.EXTRA) THEN
      PARTITION=PARTITION+1
      ISTART = 1_spi+RANK*PARTITION
      IEND = (1_spi+RANK)*PARTITION
    ELSE 
      ISTART = 1_spi+RANK*PARTITION+EXTRA
      IEND = (1_spi+RANK)*PARTITION+EXTRA
    END IF
  END IF
! ----------------------------------------------------------------------
  ALLOCATE(MY_FIRST(ISTART:IEND),SOURCE=0_dpi)
  ALLOCATE(MY_LAST(ISTART:IEND),SOURCE=0_dpi)
! ----------------------------------------------------------------------
  DO IDX=ISTART,IEND
    I=1_spi+(IDX-1_spi)/PARAMS%NT
    J=MOD((IDX-1_spi),PARAMS%NT)+1_spi
    MY_FIRST(IDX)=PHI%START_IDX(I,J)
    MY_LAST(IDX)=PHI%START_IDX(I,J)+PHI%NUM_ELEMENTS(I,J)-1
  END DO
! ----------------------------------------------------------------------
  IF (RANK.EQ.ROOT) THEN
    INIT = omp_get_wtime()
    PRINT *, "INIT TIME (s)",INIT-FIRST
  END IF
! ----------------------------------------------------------------------
  DO IDX=ISTART,IEND
    CX=0._dpf
    DO N=MY_FIRST(IDX),MY_LAST(IDX)
      CX=CX+YLM_PROD(PHI%YLM_IDX(N))*PHI%PHI_DAT(N)
    END DO
    I=1_spi+(IDX-1_spi)/PARAMS%NT
    J=MOD((IDX-1_spi),PARAMS%NT)+1_spi
    PSIJ2(I,J)=CX
  END DO
! ----------------------------------------------------------------------
! should be gather
  CALL MPI_REDUCE(PSIJ2,PSIJ,MAX_IJ,MPI_COMPLEX16,MPI_SUM,ROOT,COMM,IERROR)
! ----------------------------------------------------------------------
! ALL RANKS DO THIS LOOP MANY TIMES : 6*A + 1 * 3 * (A CHOOSE 2)
! TASK 1: 
!  SPLIT THIS LOOP UP FOR 1 NODE ON CPU USING MPI
!   => ONE SIDED COMM
!   => EACH CORE GETS SOME PARTITION OF I,J SPACE
! TASK 2: 
!  OFFLOAD TO GPU
! ----------------------------------------------------------------------
!  DO J=1,PARAMS%NT
!    DO I=1,PARAMS%NS
!      CX=0._dpf
!      FIRST_IDX=PHI%START_IDX(I,J)
!      LAST_IDX=PHI%START_IDX(I,J)+PHI%NUM_ELEMENTS(I,J)-1
!      DO N=FIRST_IDX,LAST_IDX
!        CX=CX+YLM_PROD(PHI%YLM_IDX(N))*PHI%PHI_DAT(N)
!      END DO
!      PSIJ(I,J)=CX
!    END DO
!  END DO
! ----------------------------------------------------------------------
  IF (RANK.EQ.ROOT) THEN
    LAST = omp_get_wtime()
    PRINT *, "PHI TIME (s)",LAST-INIT
    DOT=0._dpf
    DO J=1,PARAMS%NT
      DO I=1,PARAMS%NS
        DOT=DOT+PSIJ(I,J)%RE*PSIJ(I,J)%RE+PSIJ(I,J)%IM*PSIJ(I,J)%IM/REAL(PARAMS%NPHIM,KIND=dpf)
      END DO
    END DO
    INIT = omp_get_wtime()
    PRINT *, "DOT TIME (s)",INIT-LAST
    PRINT *, "DOT ERROR: ",ABS(DOT-26097012.832728475_dpf)
    PRINT *, "TOTAL TIME (s)",LAST-FIRST
  END IF
! ----------------------------------------------------------------------
  CALL MPI_FINALIZE(IERROR)
! ----------------------------------------------------------------------
END PROGRAM NQMCC_ALCF_2025
