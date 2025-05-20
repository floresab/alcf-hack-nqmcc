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
  REAL(dpf)       :: STATE,FIRST,INIT,LAST
  COMPLEX(dpf)    :: CX
  INTEGER(spi)    :: I,J,RANK,SIZE,IERROR,ROOT
  INTEGER(dpi)    :: N,FIRST_IDX,LAST_IDX
  COMPLEX(dpf), DIMENSION(:,:), ALLOCATABLE :: PSIJ
  COMPLEX(dpf), DIMENSION(:),   ALLOCATABLE :: YLM_PROD
! ----------------------------------------------------------------------
  CALL MPI_INIT(IERROR)
  CALL MPI_COMM_RANK(COMM, RANK, IERROR)
  CALL MPI_COMM_SIZE(COMM, SIZE, IERROR)
  ROOT=0_spi
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
  ALLOCATE(PSIJ(PARAMS%NS,PARAMS%NT),SOURCE=CMPLX(0._dpf,0._dpf,KIND=dpf))
  ALLOCATE(YLM_PROD(PARAMS%NYLM),SOURCE=CMPLX(1.12345_dpf,0.6789_dpf,KIND=dpf))
! ----------------------------------------------------------------------
  CALL PHI%INIT_PHI(PARAMS,STATE,RANK)
  IF (RANK.EQ.ROOT) THEN
    INIT = omp_get_wtime()
    PRINT *, "INIT TIME (s)",INIT-FIRST
  END IF
! ----------------------------------------------------------------------
  !MAX_IJ=IMAX*JMAX
  !EXTRA = MOD(MAX_IJ,SIZE)
  !PARTITION=MAX_IJ/SIZE
  !IF (EXTRA.EQ.0) THEN
  !  ISTART=1_spi+RANK*PARTITION
  !  IEND=(1_spi+RANK)*PARTITION
  !ELSE
  !  IF (RANK.LT.EXTRA) THEN
  !    ISTART = 1_spi+RANK*(PARTITION+1)
  !    IEND = (1_spi+RANK)*(PARTITION+1)
  !  ELSE 
  !    ISTART = 1_spi+RANK*PARTITION+EXTRA
  !    IEND = (1_spi+RANK)*PARTITION+EXTRA
  !  END IF
  !END IF
! ----------------------------------------------------------------------
  !TOTAL_SUM=0._dpf
  !!$omp target teams distribute parallel do &
  !!$omp& map(tofrom:TOTAL_SUM) &
  !!$omp& map(to:K_DATA) &
  !!$omp& private(k) &
  !!$omp& shared(IMAX,JMAX,KMAX,K_DATA) &
  !!$omp& reduction(+:total_sum)
  !DO IDX=ISTART,IEND
  !  I=1_spi+(IDX-1_spi)/JMAX
  !  J=MOD((IDX-1_spi),JMAX)+1_spi
  !  LOCAL_SUM=0._dpf
  !  DO K=K_DATA%K_START(I,J),K_DATA%K_END(I,J)
  !    LOCAL_SUM=LOCAL_SUM+GET_VALUE(K_DATA,I,J,K)
  !  END DO
  !  TOTAL_SUM=TOTAL_SUM+LOCAL_SUM
  !END DO
! ----------------------------------------------------------------------
  DO J=1,PARAMS%NT
    DO I=1,PARAMS%NS
      CX=0._dpf
      FIRST_IDX=PHI%START_IDX(I,J)
      LAST_IDX=PHI%START_IDX(I,J)+PHI%NUM_ELEMENTS(I,J)-1
      DO N=FIRST_IDX,LAST_IDX
        CX=CX+YLM_PROD(PHI%YLM_IDX(N))*PHI%PHI_DAT(N)
      END DO
      PSIJ(I,J)=CX
    END DO
  END DO
! ----------------------------------------------------------------------
  IF (RANK.EQ.ROOT) THEN
    LAST = omp_get_wtime()
    PRINT *, "PHI TIME (s)",LAST-INIT
    PRINT *, "TOTAL TIME (s)",LAST-FIRST
  END IF
! ----------------------------------------------------------------------
  CALL MPI_FINALIZE(IERROR)
! ----------------------------------------------------------------------
END PROGRAM NQMCC_ALCF_2025
