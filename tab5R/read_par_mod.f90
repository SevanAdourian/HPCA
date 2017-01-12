subroutine read_parameters(parameters_filename)

  !use parameters
  !use wavenumber
  !use constants
  implicit none

  character (len=256)	:: model_filename
  real			:: htol
  real			:: facmin, facmax, HH
  real, allocatable	:: kv(:)
  real			:: kmin, dk, kmax
  integer		:: nk
  real, parameter		:: Id2(2,2) = reshape(source = (/ 1., 0., 0., 1. /), &
       shape = (/ 2, 2 /))
  real, parameter		:: Ze2(2,2) = reshape(source = (/ 0., 0., 0., 0. /), &
       shape = (/ 2, 2 /))
  real, parameter		:: DET_MIN  = 1.e-20
  real, parameter		:: M_PI = 3.14159265
  
  character (len = 256)	:: parameters_filename

  ! local variables
  integer, parameter	:: UNIT_PARAMETER = 8

  open (status = "old", unit = UNIT_PARAMETER, file = parameters_filename)

  ! Model
  read(UNIT_PARAMETER, err = 20, end = 20, fmt = *) model_filename

  ! Integration parameters
  read(UNIT_PARAMETER, err = 20, end = 20, fmt = *) facmin
  read(UNIT_PARAMETER, err = 20, end = 20, fmt = *) facmax
  read(UNIT_PARAMETER, err = 20, end = 20, fmt = *) HH

  close(UNIT_PARAMETER)
  return

20 print *, "Error: reading the parameters file"

end subroutine read_parameters

subroutine read_model(model_filename)

  implicit none

  integer, parameter    :: MXMD=10
  integer                       :: N, Np, M
  real,         allocatable     :: vp0(:), vs0(:), rh0(:), th0(:)
  real,         allocatable     ::  vp(:),  vs(:),  rh(:),  th(:)
  real,         allocatable     :: lambda(:), mu(:)
  real,         allocatable     :: delta(:), gamma(:)
  
  external store_struct
  character (len=256)	:: model_filename

  ! local variables

  integer, parameter	:: UNIT_MODEL = 9
  integer			:: j
  real			:: factor

  open (status = "old", unit = UNIT_MODEL, file = model_filename)
  read (UNIT_MODEL, end = 10, FMT = *)  N, factor

  allocate(vp0(N), vs0(N), rh0(N), th0(N))
  do j=1,N
     read (UNIT_MODEL, end = 10, fmt = *)  rh0(j), vp0(j), vs0(j), th0(j)
     call store_struct(vp0(j), vs0(j), rh0(j), th0(j))
  end do

  close(UNIT_MODEL)

  vp0 = vp0*factor
  vs0 = vs0*factor
  rh0 = rh0*factor
  th0 = th0*factor

  print *, vs0
  return

10 print *, "Error reading the model file: ", model_filename
  stop 

end subroutine read_model
