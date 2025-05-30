PROGRAM NQMCC_ALCF_2025
! ----------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: spi=>int32, dpi=>int64, dpf=>real64 &
  &,COMPILER_OPTIONS,COMPILER_VERSION
! ----------------------------------------------------------------------
  USE NQMCC_M
  USE MPI_F08
  USE OMP_LIB
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
  TYPE(MPI_COMM)  :: COMM,NODE_COMM
  TYPE(CONTROL_T) :: PARAMS
  TYPE(PHI_T)     :: PHI
  TYPE(PHI_SCAT_T):: PHI_SC
  REAL(dpf)       :: STATE,FIRST,INIT,LAST,DOT
  COMPLEX(dpf)    :: CX
  INTEGER(spi)    :: I,J,RANK,SIZE,IERROR,ROOT,MAX_IJ
  INTEGER(spi)    :: NODE_RANK,NODE_SIZE,NUM_NODES,NODE_ROOT
  INTEGER(spi)    :: MY_SAMPLES,PARTITION,EXTRA,ISTART,IEND,IDX
  INTEGER(dpi)    :: S,N
  COMPLEX(dpf), DIMENSION(:,:), ALLOCATABLE :: PSIJ
  COMPLEX(dpf), DIMENSION(:),   ALLOCATABLE :: YLM_PROD,PSIJ_FLAT
  INTEGER(spi), DIMENSION(:),   ALLOCATABLE :: COUNTS,DISP
! ----------------------------------------------------------------------
  CALL MPI_INIT(IERROR)
  COMM=MPI_COMM_WORLD
  ROOT=0_spi
  CALL MPI_COMM_RANK(COMM, RANK, IERROR)
  CALL MPI_COMM_SIZE(COMM, SIZE, IERROR)
!-----------------------------------------------------------------------
  NODE_ROOT=0_spi
!-----------------------------------------------------------------------
  CALL MPI_COMM_SPLIT_TYPE(COMM, MPI_COMM_TYPE_SHARED, NODE_ROOT, MPI_INFO_NULL, NODE_COMM, IERROR)
  CALL MPI_COMM_SIZE(NODE_COMM, NODE_SIZE, IERROR)
  CALL MPI_COMM_RANK(NODE_COMM, NODE_RANK, IERROR)
  ! free the sub-communicator
  ! determine the number of allocated nodes
  NUM_NODES = SIZE / NODE_SIZE
!-----------------------------------------------------------------------
  FIRST=0._dpf
  INIT=0._dpf
  LAST=0._dpf
  DOT=0._dpf
!-----------------------------------------------------------------------
  IF (RANK.EQ.ROOT) THEN
    FIRST = omp_get_wtime()
    CALL PARAMS%READ_CONTROL()
    PRINT *, REPEAT("*",72)
    PRINT *," COMPILER: ", COMPILER_VERSION()
    PRINT *," OPTIONS: ", COMPILER_OPTIONS()
    PRINT *, REPEAT("*",72)
  END IF
  CALL PARAMS%BROADCAST(ROOT,COMM,IERROR)
  CALL PHI%ALLOCATE_PHI(PARAMS)
! ----------------------------------------------------------------------
  ALLOCATE(PSIJ(PARAMS%NS,PARAMS%NT),SOURCE=CMPLX(0._dpf,0._dpf,KIND=dpf))
  ALLOCATE(YLM_PROD(PARAMS%NYLM),SOURCE=CMPLX(0._dpf,0._dpf,KIND=dpf))
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
  IF (RANK.EQ.ROOT) THEN
    ALLOCATE(COUNTS(SIZE),SOURCE=0_spi)
    ALLOCATE(DISP(SIZE),SOURCE=0_spi)
    IF (EXTRA.EQ.0) THEN
      DO I=1,SIZE
        COUNTS(I)=PARTITION
      END DO
    ELSE
      DO I=1,SIZE
        IF (I.LE.EXTRA) THEN
          COUNTS(I)=PARTITION
        ELSE
          COUNTS(I)=PARTITION-1
        END IF
      END DO
    END IF
    DO I=2,SIZE
      DISP(I)=SUM(COUNTS(1:I-1))
    END DO
  ELSE
    ALLOCATE(COUNTS(1),SOURCE=0_spi)
    ALLOCATE(DISP(1),SOURCE=0_spi)
  END IF
! ----------------------------------------------------------------------
  MY_SAMPLES=PARAMS%NSAMPLES
! ----------------------------------------------------------------------
  ALLOCATE(PSIJ_FLAT(PARTITION),SOURCE=CMPLX(0._dpf,0._dpf,KIND=dpf))
! ----------------------------------------------------------------------
  print *, maxval(phi%ylm_idx)
  CALL PHI%SCATTER_PHI(PHI_SC,PARAMS,ISTART,IEND)
#if 1 == gpu_offload
  PRINT *, "ENTERING PHI_SC ON DEVICE"
  CALL PHI_SC%PHI_MAP_TO_DEVICE()
#endif
! ----------------------------------------------------------------------
  IF (RANK.EQ.ROOT) THEN
    INIT = omp_get_wtime()
    PRINT *, "INIT TIME (s)",INIT-FIRST
  END IF
! ----------------------------------------------------------------------
  DO S=1,MY_SAMPLES
    DO I=1,PARAMS%NYLM
      YLM_PROD(I)=CMPLX(2.3_dpf*I*10._dpf**(-12),1.23_dpf*S*10._dpf**(-12),KIND=dpf)
    END DO
! ----------------------------------------------------------------------
! ALL WALKERS DO THIS LOOP MANY TIMES FOR EACH SAMPLE : 6*A + 1 * 3 * (A CHOOSE 2)
! TASK 1: 
!  SPLIT THIS LOOP UP FOR 1 NODE ON CPU USING MPI
!   => EACH CORE GETS SOME PARTITION OF I,J SPACE
! TASK 2: 
!  OFFLOAD TO GPU
! ----------------------------------------------------------------------
#if 1 == gpu_offload
    !$omp target teams distribute parallel do &
    !$omp& map(tofrom:psij_flat) &
    !$omp& map(to:phi_sc,ylm_prod) &
    !$omp& private(n,cx) &
    !$omp& shared(psij_flat,phi,ylm_prod)
#endif
    DO IDX=1,PARTITION
      CX=0._dpf
      DO N=PHI_SC%FIRST(IDX),PHI_SC%LAST(IDX)
        CX=CX+YLM_PROD(PHI_SC%YLM_IDX_SC(N))*PHI_SC%PHI_DAT_SC(N)
      END DO
      PSIJ_FLAT(IDX)=CX
    END DO
    CALL MPI_GATHERV(PSIJ_FLAT,PARTITION,MPI_COMPLEX16,PSIJ,COUNTS,DISP,MPI_COMPLEX16,ROOT,COMM,IERROR)
! ----------------------------------------------------------------------
    IF (RANK.EQ.ROOT) THEN
      DOT=0._dpf
      DO J=1,PARAMS%NT
        DO I=1,PARAMS%NS
          DOT=DOT+PSIJ(I,J)%RE*PSIJ(I,J)%RE+PSIJ(I,J)%IM*PSIJ(I,J)%IM
        END DO
      END DO
    END IF
  END DO
  IF (RANK.EQ.ROOT) THEN
    LAST = omp_get_wtime()
    PRINT *, "SAMPLE TIME (s)",LAST-INIT
    PRINT *, "TIME PER SAMPLE (s)",(LAST-INIT)/PARAMS%NSAMPLES
    IF (PARAMS%NT.EQ.70) PRINT *, "DOT ERROR: ",ABS(DOT-3.7704504565921858E-002_dpf)
    IF (PARAMS%NT.EQ.429) PRINT *, "DOT ERROR: ",ABS(DOT-8.367575412999564E-014_dpf)
    PRINT *, "DOT: ",DOT
    PRINT *, "TOTAL TIME (s)",LAST-FIRST
  END IF
! ----------------------------------------------------------------------
  CALL MPI_COMM_FREE(NODE_COMM, IERROR)
  CALL MPI_FINALIZE(IERROR)
! ----------------------------------------------------------------------
END PROGRAM NQMCC_ALCF_2025
