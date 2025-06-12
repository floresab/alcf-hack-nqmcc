!This is the mpi logic for generating the walker sub-comm and the wavefuction comm
PROGRAM nQMCC_COMM_SPLIT
! ----------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: spi=>int32, dpf=>real64,stdin=>input_unit
! ----------------------------------------------------------------------
  USE MPI_F08
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
  BLOCK
! ----------------------------------------------------------------------
    TYPE(MPI_COMM) :: MPI_GLOBAL,MPI_WALKER,MPI_WAVEFUNCTION,MPI_NODE
! ----------------------------------------------------------------------
    INTEGER(spi) :: GLOBAL_SIZE,GLOBAL_RANK,GLOBAL_ROOT
    INTEGER(spi) :: WALKER_SIZE,WALKER_RANK,WALKER_ROOT
    INTEGER(spi) :: WF_SIZE,WF_RANK,WF_ROOT
    INTEGER(spi) :: NODE_SIZE,NODE_RANK,NODE_ROOT
! ----------------------------------------------------------------------
    LOGICAL      :: IS_WALKER
    INTEGER(spi) :: WF_ID,WF_KEY,IERROR
    INTEGER(spi) :: WALKER_ID,WALKER_KEY
    INTEGER(spi) :: NODE_ID,NUM_NODES,NUM_WALKERS_PER_NODE
    REAL(dpf)    :: NUM_WALKERS_PER_NODE_INV
    INTEGER(spi) :: CSUM,NUM_HELPERS_PER_WALKER,WF_SUM
! ----------------------------------------------------------------------
    GLOBAL_ROOT=0_spi
    WALKER_ROOT=0_spi
    WF_ROOT=0_spi
    NODE_ROOT=0_spi
! ----------------------------------------------------------------------
    IERROR=0_spi
    CALL MPI_INIT(IERROR)
    MPI_GLOBAL=MPI_COMM_WORLD
    CALL MPI_COMM_SIZE(MPI_GLOBAL,GLOBAL_SIZE,IERROR)
    CALL MPI_COMM_RANK(MPI_GLOBAL,GLOBAL_RANK,IERROR)
! ----------------------------------------------------------------------
    CALL MPI_COMM_SPLIT_TYPE(MPI_GLOBAL, MPI_COMM_TYPE_SHARED, NODE_ROOT, MPI_INFO_NULL, MPI_NODE, IERROR)
    CALL MPI_COMM_SIZE(MPI_NODE,NODE_SIZE,IERROR)
    CALL MPI_COMM_RANK(MPI_NODE,NODE_RANK,IERROR)
    CALL MPI_COMM_FREE(MPI_NODE,IERROR)
    NODE_ID   = GLOBAL_RANK / NODE_SIZE
    NUM_NODES = GLOBAL_SIZE / NODE_SIZE
! ----------------------------------------------------------------------
    IF (GLOBAL_RANK.EQ.GLOBAL_ROOT) THEN
! ----------------------------------------------------------------------
      OPEN(UNIT=STDIN,STATUS='OLD',ACTION='READ',FORM='FORMATTED')
      READ(STDIN,*) NUM_WALKERS_PER_NODE
      CLOSE(UNIT=STDIN)
      PRINT *, "NUMBER OF NODES: ",NUM_NODES
! ----------------------------------------------------------------------
      IF (NUM_WALKERS_PER_NODE.EQ.0) THEN            !trivial case each rank is a walker
         WF_SIZE=1_spi
         WALKER_SIZE=NUM_NODES*NODE_SIZE
         PRINT *, "NUM_WALKERS_PER_NODE: ",NODE_SIZE
      ELSE IF (NUM_WALKERS_PER_NODE.GT.0) THEN !non tivial multiple walkers per node
         WF_SIZE=NODE_SIZE/NUM_WALKERS_PER_NODE
         WALKER_SIZE=NUM_NODES*NUM_WALKERS_PER_NODE
         PRINT *, "NUM_WALKERS_PER_NODE: ",NUM_WALKERS_PER_NODE
      ELSE                                     !non trivial multiple nodes per walker
         NUM_WALKERS_PER_NODE_INV=REAL(1._dpf/ABS(NUM_WALKERS_PER_NODE),KIND=dpf)
         WF_SIZE=INT(NODE_SIZE/NUM_WALKERS_PER_NODE_INV,KIND=spi)
         WALKER_SIZE=INT(NUM_NODES*NUM_WALKERS_PER_NODE_INV,KIND=spi)
         PRINT *, "NUM_WALKERS_PER_NODE: ",NUM_WALKERS_PER_NODE_INV
      END IF
! ----------------------------------------------------------------------
      PRINT *, "NUM_HELPERS_PER_WALKER: ",WF_SIZE
      PRINT *, "TOTAL WALKERS: ",WALKER_SIZE
! ----------------------------------------------------------------------
    END IF
! ----------------------------------------------------------------------
    CALL MPI_BCAST(NUM_WALKERS_PER_NODE,1,MPI_INTEGER4,GLOBAL_ROOT,MPI_GLOBAL,IERROR)
! ----------------------------------------------------------------------
! trival case: num_walkers_per_node == node_size => walker_comm == mpi_comm_world
! ----------------------------------------------------------------------
    IF (NUM_WALKERS_PER_NODE.EQ.0) THEN
      IS_WALKER=.TRUE.
      MPI_WALKER=MPI_GLOBAL
      WALKER_SIZE=GLOBAL_SIZE
      WALKER_RANK=GLOBAL_RANK
      WF_ID=GLOBAL_RANK
      WF_RANK=0
! ----------------------------------------------------------------------
! non-trival case: num_walkers_per_node != node_size
! ----------------------------------------------------------------------
    ELSE IF (NUM_WALKERS_PER_NODE.GT.0) THEN
      NUM_HELPERS_PER_WALKER=NODE_SIZE/NUM_WALKERS_PER_NODE
      WALKER_SIZE=NUM_NODES*NUM_WALKERS_PER_NODE
      WF_ID=GLOBAL_RANK/NUM_HELPERS_PER_WALKER
      WF_KEY=GLOBAL_RANK
      IF (MOD(NODE_RANK,NUM_HELPERS_PER_WALKER).EQ.0) THEN
        IS_WALKER=.TRUE.
        WALKER_ID=0
        WALKER_KEY=GLOBAL_RANK
      ELSE
        IS_WALKER=.FALSE.
        WALKER_ID=MPI_UNDEFINED
        WALKER_KEY=0
      END IF
! ----------------------------------------------------------------------
      CALL MPI_COMM_SPLIT(MPI_GLOBAL,WALKER_ID, WALKER_KEY, MPI_WALKER      , IERROR)
      CALL MPI_COMM_SPLIT(MPI_GLOBAL,WF_ID    , WF_KEY    , MPI_WAVEFUNCTION, IERROR)
      CALL MPI_COMM_RANK(MPI_WAVEFUNCTION,WF_RANK,IERROR)
! ----------------------------------------------------------------------
! non-trival case: num_walkers_per_node != node_size and each walker needs more than 1 node
! ----------------------------------------------------------------------
    ELSE
      NUM_WALKERS_PER_NODE_INV=REAL(1._dpf/ABS(NUM_WALKERS_PER_NODE),KIND=dpf) !USE NEGATIVE TO INDICATE INVERSE
      NUM_HELPERS_PER_WALKER=INT(NODE_SIZE/NUM_WALKERS_PER_NODE_INV,KIND=spi)
      WALKER_SIZE=INT(NUM_NODES*NUM_WALKERS_PER_NODE_INV,KIND=spi)
      WF_ID=GLOBAL_RANK/NUM_HELPERS_PER_WALKER
      WF_KEY=GLOBAL_RANK
      IF (MOD(GLOBAL_RANK,NUM_HELPERS_PER_WALKER).EQ.0) THEN
        IS_WALKER=.TRUE.
        WALKER_ID=0
        WALKER_KEY=GLOBAL_RANK
      ELSE
        IS_WALKER=.FALSE.
        WALKER_ID=MPI_UNDEFINED
        WALKER_KEY=0
      END IF
! ----------------------------------------------------------------------
      CALL MPI_COMM_SPLIT(MPI_GLOBAL,WALKER_ID, WALKER_KEY, MPI_WALKER      , IERROR)
      CALL MPI_COMM_SPLIT(MPI_GLOBAL,WF_ID    , WF_KEY    , MPI_WAVEFUNCTION, IERROR)
      CALL MPI_COMM_RANK(MPI_WAVEFUNCTION,WF_RANK,IERROR)
! ----------------------------------------------------------------------
    END IF
! ----------------------------------------------------------------------
! example reduction over walkers => SHOULD BE ZERO
! ----------------------------------------------------------------------
    IF (IS_WALKER) CALL MPI_REDUCE(WALKER_ID, CSUM, 1, MPI_INTEGER4, MPI_SUM, WALKER_ROOT, MPI_WALKER)
    IF (GLOBAL_RANK.EQ.GLOBAL_ROOT) PRINT *, "COLOR SUM: ",CSUM
! ----------------------------------------------------------------------
! example reduction over wavefunction => sum(global rank inside wavefunction helper for a given walker)
! ----------------------------------------------------------------------
    WF_SUM=0_spi
    IF (NUM_WALKERS_PER_NODE.NE.0) THEN
      CALL MPI_REDUCE(GLOBAL_RANK, WF_SUM, 1, MPI_INTEGER4, MPI_SUM, WF_ROOT, MPI_WAVEFUNCTION)
    ELSE !TRIVAL CASE
      WF_SUM=GLOBAL_RANK
    END IF
! ----------------------------------------------------------------------
    101 FORMAT(A,I4,A,I4,A,L,A,I4,A,I4,A,I8)
    WRITE(*,101) '[MPI PROCESS ', GLOBAL_RANK,' ON NODE ',NODE_ID,'] &
           & AM I A WALKER? '      , IS_WALKER, &
           & " | I AM PART OF WF " , WF_ID, &
           & " | I AM WF HELPER "  , WF_RANK, &
           & " | MY WF SUM IS "    , WF_SUM
! ----------------------------------------------------------------------
    !IF (MPI_WALKER.NE.MPI_GLOBAL) CALL MPI_COMM_FREE(MPI_WALKER,IERROR)
    CALL MPI_FINALIZE()
! ----------------------------------------------------------------------
  END BLOCK
! ----------------------------------------------------------------------
END PROGRAM nQMCC_COMM_SPLIT
