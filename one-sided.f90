program start
!$ use omp_lib
#ifdef MPI
  use mpi_f08
#endif
  implicit none
  !
  integer :: iam, nRanks, nNodes, nodeSize, subCommSize
#ifdef MPI
  type(MPI_Comm) :: MPI_COMM_NODE
  type(MPI_Comm) :: MPI_COMM_SUBCOMM
  integer :: mpiError, color
  double precision  :: timerStart, timerEnd
  type(MPI_Win)                     :: window
  integer (kind=MPI_ADDRESS_KIND)   :: buffer_size, localDisplacement, remoteDisplacement
  real(kind=8), target, allocatable :: array(:)
  integer(kind=4) :: i, allocStat, iRank, jRank, data_size, nElements
  integer(kind=4) :: nItems = 1000000
#endif
  real(kind=8) :: dev_value, ref_value

#ifdef MPI
  call MPI_Init(mpiError)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, mpiError)

  ! create a node-based sub-communicator to determe the number of ranks per node, nodeSize
  call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_NODE, mpiError)
  call MPI_Comm_size(MPI_COMM_NODE, nodeSize, mpiError)
  ! free the sub-communicator
  call MPI_Comm_Free(MPI_COMM_NODE, mpiError)

  ! determine the number of allocated nodes
  nNodes = nRanks / nodeSize
#else
  iam = 0
  nRanks = 1
  nNodes = 1
  nodeSize = 1
  subCommSize = 1
  iRank = 0
#endif

#ifdef MPI
  ! create 3 sub-communicators per node
  color = iam / (nodeSize / 3)
  call MPI_Comm_split(MPI_COMM_WORLD, color, 0, MPI_COMM_SUBCOMM, mpiError)
  call MPI_Comm_rank(MPI_COMM_SUBCOMM, iRank, mpiError)
  call MPI_Comm_size(MPI_COMM_SUBCOMM, subCommSize, mpiError)
#endif

if(iam == 0) then
  write(*,'(/a,i0)') "Number of nodes:             ", nNodes
  write(*,'(a,i0)')  "Number of ranks per node:    ", nodeSize
  write(*,'(a,i0)')  "Total number of ranks:       ", nRanks
  write(*,'(a,i0)')  "Number of sub-comms:         ", nRanks / subCommSize
  write(*,'(a,i0)')  "Number of ranks in sub-comm: ", subCommSize
endif

  ! allocate array on the device
  !$omp allocators allocate (allocator(omp_target_device_mem_alloc):array)
  allocate(array(0 : nItems * subCommSize), stat=allocStat)

  ! In this example, each sub-communicator performs identical computation
  ! In real life, each sub-communicator may work on a dedicated portion of the total workload
  ! The total workload may be split among the individual nodes
  ! Within the node, the assigned workload may be further split among sub-communicators
  ! Within sub-communicator, each rank computes its own portion of the work

  ! Following is a mock problem for the ranks of the sub-communicator:
  ! Sum up elements of the array
  ! 1111...2222...3333...
  ! the array has subCommSize groups of nItems elements
  ! where first nItems elements of the array have the computed value of 1, 
  !       second nItems elements of the array have the computer value of 2, 
  !       and so on

  ! Solution:
  ! initialize array on the device with the data computed on the rank
  ! assume that iRank+1 is the result of computation on the rank iRank in the sub-communicator
  ! assign the value iRank+1 to nItems on the rank iRank
  ! for example: on iRank==1, the array has the values 0000...2222...0000... 
  ! where zeros are default and 2 is the computed value
  !$omp target teams distribute parallel do has_device_addr(array) map(to: nItems, iRank)
  do i = nItems * iRank, nItems * (iRank+1) - 1
    array(i) = dble(iRank + 1)  ! initialize nItems elements of the array
  enddo

  ! turn the timer on to measure the communication time
#ifdef MPI
  ! make sure all ranks reach this point; timer needs a barrier
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timerStart = MPI_WTime()
#endif

  ! at this stage, each rank in the sub-communicator has only the data generated on this rank
  ! use one-sided communication to exchange data between ranks in the sub-communicator
  ! when it is done, the array on each rank in the sub-communicator will obtain the form
  ! 1111...2222...3333...
#ifdef MPI
  ! iRank reads from or writes to all other ranks jRank in MPI_COMM_SUBCOMM
  ! uncomment either USE_GET or USE_PUT
!#define USE_GET
#define USE_PUT
  ! one-sided GPU-to-GPU communication
  call MPI_Type_size(MPI_DOUBLE_PRECISION, data_size)
  buffer_size = nItems * subCommSize * data_size            ! size in bytes
  call MPI_Win_Create(array, buffer_size, data_size, MPI_INFO_NULL, MPI_COMM_SUBCOMM, window, mpiError)
  call MPI_Win_fence(0, window)

  do jRank = 0, subCommSize-1
    if(iRank /= jRank) then
      do i = 0, nItems-1
        ! transfer nElements in a single Get/Put
        ! we can potentially transfer all nItems in a single Get/Put without using the do-loop
        nElements = 1  ! transfer 1 array element in the mock exercise
#if defined(USE_GET)
        ! jRank holds the data for iRank in the sub-communicator
        ! iRank reads the data from each jRank in the sub-communicator
        localDisplacement = jRank * nItems + i ! where to save the data locally
        remoteDisplacement = jRank * nItems + i   ! where to take the data from
        call MPI_Get(array(localDisplacement), nElements, MPI_DOUBLE_PRECISION, &
                jRank, remoteDisplacement, nElements, MPI_DOUBLE_PRECISION, window, mpiError)
#else defined(USE_PUT)
        ! iRank has the data for other ranks in the sub-communicator
        ! iRank writes the data to each jRank in the sub-communicator
        localDisplacement = iRank * nItems + i ! where to get the local data
        remoteDisplacement = iRank * nItems + i   ! where to save the data on a remote rank
        call MPI_Put(array(localDisplacement), nElements, MPI_DOUBLE_PRECISION, &
                jRank, remoteDisplacement, nElements, MPI_DOUBLE_PRECISION, window, mpiError)
#endif
      enddo
    endif
  enddo

  call MPI_Win_fence(0, window)
  call MPI_Win_free(window, mpiError)
#endif

  ! turn the timer off
#ifdef MPI
  ! make sure all ranks reach this point
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timerEnd = MPI_WTime()

  if(iam == 0) write(*,'(/a,f7.1,a)') "Communication time: ", timerEnd - timerStart, " s"
#endif

  ! sum up array elements and return dev_value to the host
  !$omp target teams distribute parallel do has_device_addr(array) map(to: nItems, subCommSize) &
  !$omp reduction(+:dev_value) map(from: dev_value)
  do i = 0, nItems * subCommSize - 1
    dev_value = dev_value + array(i)
  enddo

  ! compute a reference value on the host
  ! shortcut: ref_value = factorial(subCommSize) * nItems
  !$omp parallel do reduction(+:ref_value)
  do i = 0, nItems * subCommSize - 1
    ref_value = ref_value + dble(i / nItems + 1)
  enddo

  ! at this point, all ranks in MPI_COMM_WORLD in this example should have an identical result, 
  ! so we check the correctness only on rank 0
  if(iam == 0) write(*,'(/a,2f20.2/)') "Computed and reference values: ", dev_value, ref_value

  ! deallocate the device-resident array
  deallocate(array)

#ifdef MPI
  call MPI_Comm_Free(MPI_COMM_SUBCOMM, mpiError)
  call MPI_Finalize(mpiError)
#endif
end program start
