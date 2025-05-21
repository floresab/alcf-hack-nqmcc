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
  INTEGER(spi)    :: I,J,RANK,SIZE,IERROR,ROOT,MAX_IJ
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
  FIRST=0._dpf
  INIT=0._dpf
  LAST=0._dpf
  DOT=0._dpf
!-----------------------------------------------------------------------
  IF (RANK.EQ.ROOT) THEN
    FIRST = omp_get_wtime()
    CALL PARAMS%READ_CONTROL()
  END IF
  CALL PARAMS%BROADCAST(ROOT,COMM,IERROR)
  CALL PHI%ALLOCATE_PHI(PARAMS)
! ----------------------------------------------------------------------
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
  IF (RANK.EQ.0) THEN
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
  CALL PHI%SCATTER_PHI(PARAMS,ISTART,IEND)
! ----------------------------------------------------------------------
  IF (RANK.EQ.ROOT) THEN
    INIT = omp_get_wtime()
    PRINT *, "INIT TIME (s)",INIT-FIRST
  END IF
! ----------------------------------------------------------------------
  DO S=1,MY_SAMPLES
! ----------------------------------------------------------------------
! ALL WALKERS DO THIS LOOP MANY TIMES FOR EACH SAMPLE : 6*A + 1 * 3 * (A CHOOSE 2)
! TASK 1: 
!  SPLIT THIS LOOP UP FOR 1 NODE ON CPU USING MPI
!   => EACH CORE GETS SOME PARTITION OF I,J SPACE
! TASK 2: 
!  OFFLOAD TO GPU
! ----------------------------------------------------------------------
    DO IDX=1,PARTITION
      CX=0._dpf
      DO N=PHI%FIRST(IDX),PHI%LAST(IDX)
        CX=CX+YLM_PROD(PHI%YLM_IDX_SC(N))*PHI%PHI_DAT_SC(N)
      END DO
      PSIJ_FLAT(IDX)=CX
    END DO
    CALL MPI_GATHERV(PSIJ_FLAT,PARTITION,MPI_COMPLEX16,PSIJ,COUNTS,DISP,MPI_COMPLEX16,ROOT,COMM,IERROR)
    !CALL MPI_GATHER(PSIJ_FLAT,PARTITION,MPI_COMPLEX16,PSIJ,PARTITION,MPI_COMPLEX16,ROOT,COMM,IERROR)
! ----------------------------------------------------------------------
    IF (RANK.EQ.0) THEN
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
    IF (PARAMS%NT.EQ.70) PRINT *, "DOT ERROR: ",ABS(DOT-4.2908172651317376_dpf)
    IF (PARAMS%NT.EQ.70) PRINT *, "DOT: ",DOT
    PRINT *, "TOTAL TIME (s)",LAST-FIRST
  END IF
! ----------------------------------------------------------------------
  CALL MPI_FINALIZE(IERROR)
! ----------------------------------------------------------------------
END PROGRAM NQMCC_ALCF_2025
