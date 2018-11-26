
!---------------------------------------------------------------------
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! USE in any subroutine or function or other module.
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!---------------------------------------------------------------------
module grav_types

! integer,parameter    :: dp = SELECTED_REAL_KIND(15,307)

type Tparameter_grid   ! regular grid defining the density as a function of position
                       ! Note: usually each density grid is defined over the entire propagation
                       ! grid, but only some nodes actually influence the corresponding region

    double precision :: h,lat,long      ! geospatial coordinates of point
    double precision :: density         ! density at point
    character        :: layer*2
    logical          :: moho

end type Tparameter_grid
!-------------------------------------------------------------------------------------------------------

type Telement

    integer          :: id                ! id # of element based on crustal type
    character        :: type*2
    logical          :: crust    ! element is above 30km
    logical          :: active   ! set to true for the nodes that actually
    logical          :: sheet    ! sheet element
    double precision, dimension (8) :: x,y,z
    double precision, pointer :: rho
    double precision :: thickness

end type Telement
!--------------------------------------------------------------!

type Telemento

    integer          :: id                ! id # of element based on crustal type
    character        :: type*2
    logical          :: crust    ! element is above 30km
    logical          :: active   ! set to true for the nodes that actually
    logical          :: sheet    ! sheet element
    double precision, dimension (8) :: x,y,z
    double precision :: rho,thickness

end type Telemento
!--------------------------------------------------------------!
type Tobservation

    integer                      :: id                    ! identifies the receiver
    double precision             :: lat,long            ! position
    double precision             :: x,y,z                 ! position
    double precision             :: bouguer_anom
    double precision             :: dgr                   ! calculated bouguer anomaly
    double precision             :: ux,uy,uz              ! unit vectors of point

end type Tobservation

!--------------------------------------------------------------!
type Tcrust
    double precision     :: rhoc
    integer, dimension(2):: elmnts
    logical              :: double

end type Tcrust


end module grav_types

!--------------------------------------------------------------!

module global_grav

    use grav_types

    integer                                         :: n_paragrids,n_observations,n_elements,c_elmnts
    integer                                         :: nlat,nlong,ncol
    double precision                                :: refrho,rmserr,refmoho
    double precision,parameter                      :: pi=3.1415926535
    double precision, parameter                     :: earth_radius = 6371.0
    double precision, parameter                     :: ugc = 6.67408/10**(3) ! after converting to mgal
    double precision, dimension(:,:,:),allocatable  :: x,y,z,density
    type(Tcrust), allocatable, dimension(:),target  :: cr
    character(len=2), dimension(:,:,:),allocatable  :: layer_type
    logical, dimension(:,:,:),allocatable           :: MohoP
    double precision, dimension(:,:),allocatable    :: specG         ! result of SpecGrav3d for each element-observation pair
    type(Tparameter_grid),dimension(:),allocatable  :: paragrid
    type(Telement),dimension(:),allocatable         :: element
    type(Telemento),dimension(:),allocatable         :: elemento
    type(Tobservation),dimension(:),allocatable     :: observation


end module global_grav


