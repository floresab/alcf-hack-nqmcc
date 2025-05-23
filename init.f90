! ----------------------------------------------------------------------
! Abraham R. Flores
! WASHINGTON UNIVERSITY IN ST. LOUIS
! ALCF HACKATHON 2025
! ----------------------------------------------------------------------
MODULE RNG_M
! ----------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: spi=>int32, dpf=>real64,stdin=>input_unit
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
  CONTAINS
! ----------------------------------------------------------------------
    SUBROUTINE RANDINT(RNG_STATE,NMAX,N)
      REAL(dpf),    INTENT(INOUT)  :: RNG_STATE
      INTEGER(spi), INTENT(IN)     :: NMAX
      INTEGER(spi), INTENT(OUT)    :: N
      REAL(dpf), PARAMETER :: M1 = 2147483647._dpf
      REAL(dpf), PARAMETER :: M = 2147483648._dpf
      REAL(dpf), PARAMETER :: A = 48271._dpf
      REAL(dpf) :: RNG
      RNG_STATE=MOD(A*RNG_STATE,M1)
      RNG=RNG_STATE/M*NMAX
      N=INT(RNG,KIND=spi)+1 ! randint FROM 1 -> Nmax
    END SUBROUTINE RANDINT
! ----------------------------------------------------------------------
END MODULE RNG_M
! ----------------------------------------------------------------------
MODULE NQMCC_M
! ----------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: spi=>int32,dpi=>int64,spf=>real32, dpf=>real64
! ----------------------------------------------------------------------
  USE MPI_F08
  USE RNG_M
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
  TYPE :: CONTROL_T
    INTEGER(dpi) :: NPHIM
    INTEGER(spi) :: NSAMPLES,NS,NT,NYLM,MAX_NZ,RW_PHI
    REAL(dpf) :: SEED
    CHARACTER(len=512) :: PHI_FILE
    CONTAINS
      PROCEDURE :: READ_CONTROL
      PROCEDURE :: BROADCAST
      PROCEDURE :: CONTROL_MAP_TO_DEVICE
  END TYPE CONTROL_T
! ----------------------------------------------------------------------
  TYPE :: PHI_T
    INTEGER(spi),DIMENSION(:),  ALLOCATABLE :: YLM_IDX
    REAL(spf),   DIMENSION(:),  ALLOCATABLE :: PHI_DAT
    INTEGER(spi),DIMENSION(:,:),ALLOCATABLE :: NUM_ELEMENTS
    INTEGER(dpi),DIMENSION(:,:),ALLOCATABLE :: START_IDX
! ----------------------------------------------------------------------
    CONTAINS
      PROCEDURE :: ALLOCATE_PHI
      PROCEDURE :: INIT_PHI
      PROCEDURE :: READ_PHI
      PROCEDURE :: SCATTER_PHI
  END TYPE PHI_T
! ----------------------------------------------------------------------
  TYPE :: PHI_SCAT_T
    INTEGER(spi),DIMENSION(:),  ALLOCATABLE :: YLM_IDX_SC
    REAL(spf),   DIMENSION(:),  ALLOCATABLE :: PHI_DAT_SC
    INTEGER(dpi),DIMENSION(:),  ALLOCATABLE :: FIRST,LAST
! ----------------------------------------------------------------------
    CONTAINS
      PROCEDURE :: ALLOCATE_PHI_SCAT
      PROCEDURE :: PHI_MAP_TO_DEVICE
  END TYPE PHI_SCAT_T
! ----------------------------------------------------------------------
  CONTAINS
! ----------------------------------------------------------------------
  SUBROUTINE  READ_CONTROL(SELF)
    CLASS(CONTROL_T),  INTENT(INOUT) :: SELF
    INTEGER(spi) :: NEWUNIT
    CHARACTER(len=13) :: FILE_NAMe
    FILE_NAME="inputs/c13.in"
    OPEN(newUNIT=NEWUNIT,file=FILE_NAME,STATUS='OLD',ACTION='READ',FORM='FORMATTED')
    !OPEN(UNIT=NEWUNIT,STATUS='OLD',ACTION='READ',FORM='FORMATTED')
    READ(NEWUNIT,*) SELF%NSAMPLES
    READ(NEWUNIT,*) SELF%SEED
    READ(NEWUNIT,*) SELF%NS
    READ(NEWUNIT,*) SELF%NT
    READ(NEWUNIT,*) SELF%NYLM
    READ(NEWUNIT,*) SELF%NPHIM
    READ(NEWUNIT,*) SELF%MAX_NZ
    READ(NEWUNIT,*) SELF%RW_PHI,SELF%PHI_FILE
    CLOSE(UNIT=NEWUNIT)
  END SUBROUTINE READ_CONTROL
! ----------------------------------------------------------------------
  SUBROUTINE  BROADCAST(SELF,ROOT,COMM,IERROR)
    CLASS(CONTROL_T),INTENT(INOUT) :: SELF
    TYPE(MPI_COMM),  INTENT(IN)    :: COMM
    INTEGER(spi),    INTENT(IN)    :: ROOT
    INTEGER(spi),    INTENT(INOUT) :: IERROR
    CALL MPI_BCAST(SELF%NSAMPLES,1,MPI_INTEGER4,ROOT,COMM,IERROR)
    CALL MPI_BCAST(SELF%SEED,1,MPI_REAL8,ROOT,COMM,IERROR)
    CALL MPI_BCAST(SELF%NS,1,MPI_INTEGER4,ROOT,COMM,IERROR)
    CALL MPI_BCAST(SELF%NT,1,MPI_INTEGER4,ROOT,COMM,IERROR)
    CALL MPI_BCAST(SELF%NYLM,1,MPI_INTEGER4,ROOT,COMM,IERROR)
    CALL MPI_BCAST(SELF%NPHIM,1,MPI_INTEGER8,ROOT,COMM,IERROR)
    CALL MPI_BCAST(SELF%MAX_NZ,1,MPI_INTEGER4,ROOT,COMM,IERROR)
    CALL MPI_BCAST(SELF%RW_PHI,1,MPI_INTEGER4,ROOT,COMM,IERROR)
    CALL MPI_BCAST(SELF%PHI_FILE,512,MPI_CHARACTER,ROOT,COMM,IERROR)
  END SUBROUTINE BROADCAST
! ----------------------------------------------------------------------
  SUBROUTINE  CONTROL_MAP_TO_DEVICE(SELF)
    CLASS(CONTROL_T),INTENT(INOUT) :: SELF
    !$omp target enter data map(to: SELF)
  END SUBROUTINE CONTROL_MAP_TO_DEVICE
! ----------------------------------------------------------------------
  SUBROUTINE  ALLOCATE_PHI(SELF,PARAMS)
    CLASS(PHI_T),   INTENT(INOUT) :: SELF
    TYPE(CONTROL_T),INTENT(IN)    :: PARAMS
! ----------------------------------------------------------------------
    ALLOCATE(SELF%YLM_IDX(PARAMS%NPHIM),SOURCE=0_spi)
    ALLOCATE(SELF%PHI_DAT(PARAMS%NPHIM),SOURCE=0._spf)
    ALLOCATE(SELF%NUM_ELEMENTS(PARAMS%NS,PARAMS%NT),SOURCE=0_spi)
    ALLOCATE(SELF%START_IDX(PARAMS%NS,PARAMS%NT),SOURCE=0_dpi)
! ----------------------------------------------------------------------
  END SUBROUTINE ALLOCATE_PHI
! ----------------------------------------------------------------------
  SUBROUTINE  INIT_PHI(SELF,PARAMS,STATE,RANK)
    CLASS(PHI_T),   INTENT(INOUT) :: SELF
    TYPE(CONTROL_T),INTENT(IN)    :: PARAMS
    REAL(dpf),      INTENT(INOUT) :: STATE
    INTEGER(spi),   INTENT(IN)    :: RANK
    INTEGER(spi) :: I,J,NX,NY
    INTEGER(dpi) :: II
! ----------------------------------------------------------------------
    DO II=1,PARAMS%NPHIM
      CALL RANDINT(STATE,PARAMS%NYLM,NY)
      SELF%YLM_IDX(II)=NY
      SELF%PHI_DAT(II)=REAL(NY,KIND=spf)/PARAMS%NYLM
    END DO
    II=1
    DO J=1,PARAMS%NT
      DO I=1,PARAMS%NS
        CALL RANDINT(STATE,PARAMS%MAX_NZ,NX)
        SELF%NUM_ELEMENTS(I,J)=NX
        SELF%START_IDX(I,J)=II
        II=II+NX
      END DO
    END DO
    IF (RANK.EQ.0) PRINT *, "HOW FULL",REAL(II,kind=dpf)/PARAMS%NPHIM
! ----------------------------------------------------------------------
  END SUBROUTINE INIT_PHI
! ----------------------------------------------------------------------
  SUBROUTINE  ALLOCATE_PHI_SCAT(SELF,MY_IJ_PART,MY_NZ)
    CLASS(PHI_SCAT_T),INTENT(INOUT) :: SELF
    INTEGER(spi),     INTENT(IN) :: MY_IJ_PART
    INTEGER(dpi),     INTENT(IN) :: MY_NZ
! ----------------------------------------------------------------------
    ALLOCATE(SELF%FIRST(MY_IJ_PART),SOURCE=0_dpi)
    ALLOCATE(SELF%LAST(MY_IJ_PART),SOURCE=0_dpi)
    ALLOCATE(SELF%PHI_DAT_SC(MY_NZ),SOURCE=0._spf)
    ALLOCATE(SELF%YLM_IDX_SC(MY_NZ),SOURCE=0_spi)
! ----------------------------------------------------------------------
  END SUBROUTINE ALLOCATE_PHI_SCAT
! ----------------------------------------------------------------------
  SUBROUTINE  SCATTER_PHI(SELF,PHI_SC,PARAMS,IJ_START,IJ_END)
    CLASS(PHI_T),     INTENT(INOUT) :: SELF
    CLASS(PHI_SCAT_T),INTENT(INOUT) :: PHI_SC
    TYPE(CONTROL_T),  INTENT(IN)    :: PARAMS
    INTEGER(spi),     INTENT(IN)    :: IJ_START,IJ_END
! ----------------------------------------------------------------------
    INTEGER(dpi) :: MY_NZ,FIRST_IDX,LAST_IDX,IPHI
    INTEGER(spi) :: I,J,IDX,MY_IJ_PART,IJ
! ----------------------------------------------------------------------
    MY_NZ=0_dpi
    DO IDX=IJ_START,IJ_END
      J=1_spi+(IDX-1_spi)/PARAMS%NS
      I=MOD((IDX-1_spi),PARAMS%NS)+1_spi
      MY_NZ=MY_NZ+SELF%NUM_ELEMENTS(I,J)
    END DO
! ----------------------------------------------------------------------
    MY_IJ_PART=IJ_END-IJ_START+1_spi
    print *, MY_IJ_PART,MY_NZ,REAL(MY_NZ,kind=dpf)/PARAMS%NPHIM*100,"%"
    CALL PHI_SC%ALLOCATE_PHI_SCAT(MY_IJ_PART,MY_NZ)
! ----------------------------------------------------------------------
    IDX=0_spi
    IJ=0_spi
    IPHI=0_dpi
    DO IDX=IJ_START,IJ_END
      IJ=IJ+1_spi
      J=1_spi+(IDX-1_spi)/PARAMS%NS
      I=MOD((IDX-1_spi),PARAMS%NS)+1_spi
      PHI_SC%FIRST(IJ)=IPHI+1
      IPHI=IPHI+SELF%NUM_ELEMENTS(I,J)
      PHI_SC%LAST(IJ)=IPHI
      FIRST_IDX=SELF%START_IDX(I,J)
      LAST_IDX=SELF%START_IDX(I,J)+SELF%NUM_ELEMENTS(I,J)-1
      PHI_SC%PHI_DAT_SC(PHI_SC%FIRST(IJ):PHI_SC%LAST(IJ))=SELF%PHI_DAT(FIRST_IDX:LAST_IDX)
      PHI_SC%YLM_IDX_SC(PHI_SC%FIRST(IJ):PHI_SC%LAST(IJ))=SELF%YLM_IDX(FIRST_IDX:LAST_IDX)
    END DO
! ----------------------------------------------------------------------
    IF (ALLOCATED(SELF%YLM_IDX)) DEALLOCATE(SELF%YLM_IDX)
    IF (ALLOCATED(SELF%PHI_DAT)) DEALLOCATE(SELF%PHI_DAT)
    IF (ALLOCATED(SELF%NUM_ELEMENTS)) DEALLOCATE(SELF%NUM_ELEMENTS)
    IF (ALLOCATED(SELF%START_IDX)) DEALLOCATE(SELF%START_IDX)
! ----------------------------------------------------------------------
  END SUBROUTINE SCATTER_PHI
! ----------------------------------------------------------------------
  SUBROUTINE READ_PHI(SELF,FILE_NAME,IS_COMMANDER)
    CLASS (PHI_T),   INTENT(INOUT) :: SELF
    CHARACTER(len=*),INTENT(IN) :: FILE_NAME
    LOGICAL,OPTIONAL,INTENT(IN) :: IS_COMMANDER
    INTEGER(spi) :: IN_UNIT
    LOGICAL :: FILE_EXISTS,IPRINT
    IF (PRESENT(IS_COMMANDER)) THEN
      IPRINT=IS_COMMANDER
    ELSE
      IPRINT=.TRUE.
    END IF
    INQUIRE(FILE=FILE_NAME,EXIST=FILE_EXISTS)
    IF (.NOT.FILE_EXISTS) STOP ("BAD FILE IN READ PHI")
    IF (IPRINT) WRITE(*,*) "READING PHI FROM: ",TRIM(FILE_NAME)
    OPEN(NEWUNIT=IN_UNIT,FILE=FILE_NAME,FORM='UNFORMATTED',ACTION="READ")
    READ(IN_UNIT) SELF%YLM_IDX
    READ(IN_UNIT) SELF%NUM_ELEMENTS
    READ(IN_UNIT) SELF%PHI_DAT
    READ(IN_UNIT) SELF%START_IDX
    CLOSE(IN_UNIT)
    IF (IPRINT) WRITE(*,*) "READING PHI DONE"
  END SUBROUTINE READ_PHI
! ----------------------------------------------------------------------
  SUBROUTINE  PHI_MAP_TO_DEVICE(SELF)
    CLASS(PHI_SCAT_T),INTENT(INOUT) :: SELF
    !$omp target enter data map(to: self &
    !$omp& ,self%ylm_idx_sc(:) &
    !$omp& ,self%phi_dat_sc(:) &
    !$omp& ,self%first(:) &
    !$omp& ,self%last(:))
  END SUBROUTINE PHI_MAP_TO_DEVICE
! ----------------------------------------------------------------------
END MODULE NQMCC_M
! ----------------------------------------------------------------------
