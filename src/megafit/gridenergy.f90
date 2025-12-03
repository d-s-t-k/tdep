#include "precompilerdefinitions"
module gridenergy
!! Evaluate forceconstants and associated quantities at arbitrary points
use konstanter, only: flyt,r8,i8,lo_huge,lo_hugeint,lo_pi,lo_twopi,lo_imag,lo_status,lo_exitcode_param,lo_exitcode_symmetry,&
                      lo_tol,lo_sqtol,lo_pressure_HartreeBohr_to_GPa,lo_pressure_GPa_to_HartreeBohr,lo_Hartree_to_eV,&
                      lo_volume_bohr_to_A,lo_volume_A_to_bohr,lo_freqtol,lo_kb_Hartree
use gottochblandat, only: walltime,tochar,lo_progressbar_init,lo_progressbar,lo_looptimer,lo_sqnorm,lo_mean,&
                          lo_planck,open_file,lo_does_file_exist,lo_flattentensor,lo_linspace,qsort,lo_return_unique,&
                          lo_linear_least_squares
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_IN_PLACE
use lo_memtracker, only: lo_mem_helper
use geometryfunctions, only: lo_linesegment,lo_convex_hull_2d
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_firstorder,  only: lo_forceconstant_firstorder
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder,  only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_jij_secondorder, only: lo_jij_secondorder
use type_forcemap, only: lo_forcemap
use type_qpointmesh, only: lo_qpoint_mesh,lo_fft_mesh,lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use type_blas_lapack_wrappers, only: lo_dgels
use hdf5_wrappers, only: lo_hdf5_helper,lo_h5_store_data,lo_h5_store_attribute,HID_T,H5F_ACC_TRUNC_F,&
                         h5close_f,h5open_f,h5fclose_f,h5fopen_f,h5fcreate_f,h5gclose_f,h5gopen_f,h5gcreate_f

use type_gridsim, only: lo_gridsim
use type_equation_of_state, only: lo_eos,lo_eos_1d,lo_eos_2d,lo_eos_birch_murnaghan,lo_eos_vinet,lo_eos_2d_birch_murnaghan
use type_polynomial_interpolation, only: lo_grid_interpolation
use helperobjects, only: megafit_secondorder_constraints

implicit none

private
public :: lo_gridenergy

!> quasiharmonic volume grid
type lo_gridenergy_vol
    !> how many volumes
    integer :: nv=-lo_hugeint
    integer :: nt=-lo_hugeint
    !> volume axis
    real(flyt), dimension(:), allocatable :: volume
    !> temperature axis
    real(flyt), dimension(:), allocatable :: temperature
    !> energies
    real(flyt), dimension(:,:), allocatable :: U,U0,fph,ah3,ah4
end type

!> volume-temperature grid
type lo_gridenergy_vol_temp
    !> how many volumes
    integer :: nv=-lo_hugeint
    !> how many temperatures
    integer :: nt=-lo_hugeint
    !> volume axis
    real(flyt), dimension(:,:), allocatable :: volume
    !> temperature axis
    real(flyt), dimension(:), allocatable :: temperature
    !> energies
    real(flyt), dimension(:,:), allocatable :: U,U0,fph,ah3,ah4
    real(flyt), dimension(:,:), allocatable :: qh_fph,qh_ah3,qh_ah4
end type

!> volume-temperature-eta grid
type lo_gridenergy_vol_temp_eta
    !> how many volumes
    integer :: nv=-lo_hugeint
    !> how many temperatures
    integer :: nt=-lo_hugeint
    !> how many eta
    integer :: neta=-lo_hugeint
    !> volume axis
    real(flyt), dimension(:), allocatable :: volume
    !> temperature axis
    real(flyt), dimension(:), allocatable :: temperature
    !> eta axis axis
    real(flyt), dimension(:), allocatable :: eta
    !> energies
    real(flyt), dimension(:,:,:), allocatable :: U,U0,fph,ah3,ah4
    real(flyt), dimension(:,:,:), allocatable :: qh_fph,qh_ah3,qh_ah4
end type

!> a-c lattice parameter grid (for hexagonal/tetragonal QHA)
type lo_gridenergy_ac
    !> number of a values
    integer :: na=-lo_hugeint
    !> number of c values
    integer :: nc=-lo_hugeint
    !> number of temperatures for evaluation
    integer :: nt=-lo_hugeint
    !> a axis (Angstrom)
    real(flyt), dimension(:), allocatable :: a_values
    !> c axis (Angstrom)
    real(flyt), dimension(:), allocatable :: c_values
    !> temperature axis
    real(flyt), dimension(:), allocatable :: temperature
    !> energies F(a,c,T)
    real(flyt), dimension(:,:,:), allocatable :: U,U0,fph,ah3,ah4,Ftot
    !> minimum finding results
    real(flyt), dimension(:), allocatable :: a_min, c_min, F_min  ! minimum a,c,F at each T
    !> polynomial coefficients for the 2D fit
    real(flyt), dimension(:,:), allocatable :: poly_coeffs  ! coefficients at each T
end type

!> a-c-T lattice parameter grid with T-dependent force constants
type lo_gridenergy_act
    !> number of a values
    integer :: na=-lo_hugeint
    !> number of c values  
    integer :: nc=-lo_hugeint
    !> number of temperatures (for FC interpolation)
    integer :: nt=-lo_hugeint
    !> a axis (Angstrom)
    real(flyt), dimension(:), allocatable :: a_values
    !> c axis (Angstrom)
    real(flyt), dimension(:), allocatable :: c_values
    !> temperature axis (K)
    real(flyt), dimension(:), allocatable :: temperature
    !> energies F(a,c,T) - all at same T as FC evaluation
    real(flyt), dimension(:,:,:), allocatable :: U,U0,fph,ah3,ah4,Ftot
    !> minimum finding results
    real(flyt), dimension(:), allocatable :: a_min, c_min, F_min  ! minimum a,c at each T
    !> polynomial coefficients for the 2D fit at each T
    real(flyt), dimension(:,:), allocatable :: poly_coeffs
end type


!> evaluate the free energy across the mesh
type lo_gridenergy
    !> how many dimensions
    integer :: ndim
    !> gridpoints per dimension
    integer, dimension(:), allocatable :: pts_per_dim
    !> min per dimension
    real(flyt), dimension(:), allocatable :: min_per_dim
    !> min max per dimension
    real(flyt), dimension(:), allocatable :: max_per_dim
    !> evaluate quasiharmonic energyies
    logical :: quasiharmonic=.false.
    !> what kind of grid do we have
    integer :: gridtype
    !> Do it slightly differently depending on the number of dimenisons
    !type(lo_gridenergy_2d) :: grid_2d
    type(lo_gridenergy_vol) :: grid_V
    type(lo_gridenergy_vol_temp) :: grid_VT
    type(lo_gridenergy_vol_temp_eta) :: grid_VTeta
    type(lo_gridenergy_ac) :: grid_AC
    type(lo_gridenergy_act) :: grid_ACT
    contains
        !> create everything
        procedure :: generate
end type

! how oftern to report the slow loops
real(flyt), parameter :: timereport=30.0_flyt
integer, parameter :: pm_vgrid=1
integer, parameter :: pm_vtgrid=2
integer, parameter :: pm_vtetagrid=3
integer, parameter :: pm_acgrid=4
integer, parameter :: pm_actgrid=5

contains

#include "gridenergy_aux.f90"
#include "gridenergy_vt_grid.f90"
#include "gridenergy_vteta_grid.f90"
#include "gridenergy_impossible.f90"

!> evaluate the irreducible representation at a certain point
subroutine generate(ge,gs,map,qgrid_harm,qgrid_anharm,quasiharmonic,dumpgrid,mw,mem)
    !> interpolation grid
    class(lo_gridenergy), intent(out) :: ge
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> q-point mesh dimensions
    integer, dimension(3) :: qgrid_harm
    !> support q-grid mesh dimensions
    integer, dimension(3) :: qgrid_anharm
    !> should I calculate the quasiharmonic as well
    logical, intent(in) :: quasiharmonic
    !> should I dump the input files for the entire grid
    logical, intent(in) :: dumpgrid
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    character(len=3), dimension(gs%ndim) :: pointspacing
    real(flyt), dimension(:,:), allocatable :: gridcoord,gridenergy
    real(flyt), dimension(:), allocatable :: adapt_minvol,adapt_maxvol
    real(flyt) :: pressurestep,timer
    integer, dimension(:,:), allocatable :: gridind
    integer :: verbosity,evalmode

    ! read the input file and figure out some stuff
    init: block
        integer :: i,u
        if ( mw%talk ) then
            write(*,*) ''
            write(*,*) 'INTERPOLATING FREE ENERGY'
            timer=walltime()
            verbosity=1
        else
            timer=walltime()
            verbosity=-1
        endif
        ! grab dimensions and grid stuff from file
        ge%ndim=gs%ndim
        u=open_file('in','infile.evalpoints')
            read(u,*) evalmode
            select case(evalmode)
            case(1)
                ! generic way of definig the mesh
                lo_allocate(ge%min_per_dim(ge%ndim))
                lo_allocate(ge%max_per_dim(ge%ndim))
                lo_allocate(ge%pts_per_dim(ge%ndim))
                ge%min_per_dim=0.0_flyt
                ge%max_per_dim=0.0_flyt
                ge%pts_per_dim=0
                do i=1,ge%ndim
                    read(u,*) ge%pts_per_dim(i),pointspacing(i)
                    read(u,*) ge%min_per_dim(i),ge%max_per_dim(i)
                enddo
                read(u,*) pressurestep
                if ( ge%ndim .eq. 1 .and. gs%info%dim_volume .gt. 0 ) then
                    write(*,*) 'onlyqh'
                endif
                ! Convert to atomic units
                i=gs%info%dim_volume
                if ( i .gt. 0 ) then
                    ge%min_per_dim(i)=ge%min_per_dim(i)*lo_volume_A_to_bohr
                    ge%max_per_dim(i)=ge%max_per_dim(i)*lo_volume_A_to_bohr
                endif
                pressurestep=pressurestep*lo_pressure_GPa_to_HartreeBohr
            case(3)
                ! specialized that only works in 2D now, and with V-T as the grid
                lo_allocate(ge%min_per_dim(ge%ndim))
                lo_allocate(ge%max_per_dim(ge%ndim))
                lo_allocate(ge%pts_per_dim(ge%ndim))
                ge%min_per_dim=0.0_flyt
                ge%max_per_dim=0.0_flyt
                ge%pts_per_dim=0

                i=gs%info%dim_temperature
                read(u,*) ge%pts_per_dim(i),pointspacing(i)
                read(u,*) ge%min_per_dim(i),ge%max_per_dim(i)
                i=gs%info%dim_volume
                read(u,*) ge%pts_per_dim(i),pointspacing(i)
                read(u,*) pressurestep
                ! Instead of reading the volume limits from file, get them per temperature
                ! from the convex hull.
                allocate(adapt_minvol( ge%pts_per_dim(gs%info%dim_temperature) ))
                allocate(adapt_maxvol( ge%pts_per_dim(gs%info%dim_temperature) ))
                adapt_minvol=0.0_flyt
                adapt_minvol=0.0_flyt
                i=gs%info%dim_temperature
                call volume_limits_from_convex_hull(gs%grid_coordinates,gs%info%dim_volume,gs%info%dim_temperature,&
                    ge%min_per_dim(i),ge%max_per_dim(i),ge%pts_per_dim(i),pointspacing(i),adapt_minvol,adapt_maxvol,verbosity)
                ! Convert to atomic units
                pressurestep=pressurestep*lo_pressure_GPa_to_HartreeBohr
            case(4)
                ! a-c lattice parameter grid mode (for hexagonal/tetragonal QHA)
                ! Expects 2D: a and c as the grid dimensions
                if ( gs%ndim .ne. 2 ) then
                    call lo_stop_gracefully(['evalmode=4 requires 2D grid with a and c dimensions'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                endif
                if ( gs%info%dim_a .le. 0 .or. gs%info%dim_c .le. 0 ) then
                    call lo_stop_gracefully(['evalmode=4 requires a and c dimension names in infile.simulations'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                endif
                lo_allocate(ge%min_per_dim(ge%ndim))
                lo_allocate(ge%max_per_dim(ge%ndim))
                lo_allocate(ge%pts_per_dim(ge%ndim))
                ge%min_per_dim=0.0_flyt
                ge%max_per_dim=0.0_flyt
                ge%pts_per_dim=0
                do i=1,ge%ndim
                    read(u,*) ge%pts_per_dim(i),pointspacing(i)
                    read(u,*) ge%min_per_dim(i),ge%max_per_dim(i)
                enddo
                ! Force the grid type to a-c grid
                ge%gridtype=pm_acgrid
            case(5)
                ! a-c-T lattice parameter grid mode with T-dependent force constants
                ! Expects 3D: a, c, and T as the grid dimensions
                if ( gs%ndim .ne. 3 ) then
                    call lo_stop_gracefully(['evalmode=5 requires 3D grid with a, c, and T dimensions'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                endif
                if ( gs%info%dim_a .le. 0 .or. gs%info%dim_c .le. 0 .or. gs%info%dim_temperature .le. 0 ) then
                    call lo_stop_gracefully(['evalmode=5 requires a, c, and T dimension names in infile.simulations'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                endif
                lo_allocate(ge%min_per_dim(ge%ndim))
                lo_allocate(ge%max_per_dim(ge%ndim))
                lo_allocate(ge%pts_per_dim(ge%ndim))
                ge%min_per_dim=0.0_flyt
                ge%max_per_dim=0.0_flyt
                ge%pts_per_dim=0
                do i=1,ge%ndim
                    read(u,*) ge%pts_per_dim(i),pointspacing(i)
                    read(u,*) ge%min_per_dim(i),ge%max_per_dim(i)
                enddo
                ! Force the grid type to a-c-T grid
                ge%gridtype=pm_actgrid
            case default
                call lo_stop_gracefully(['NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            end select
        close(u)
        if ( mw%talk ) then
            write(*,*) '              ndim:',ge%ndim
        endif
        ! Figure out what kind of grid we have.
        select case(gs%ndim)
        case(pm_vgrid)
            if ( gs%info%dim_volume .gt. 0 ) then
                ! V grid
                ge%gridtype=pm_vgrid !  1
            else
                call lo_stop_gracefully(['1-D NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif
        case(pm_vtgrid)
            if ( gs%info%dim_temperature .gt. 0 .and. gs%info%dim_volume .gt. 0 ) then
                ! V-T grid
                ge%gridtype=pm_vtgrid !  1
                if ( mw%talk ) then
                    write(*,*) '          gridtype:',ge%gridtype,'(V-T grid)'
                    write(*,*) '   temperature-dim:',gs%info%dim_temperature
                    write(*,*) '        volume-dim:',gs%info%dim_volume
                endif
            else
                call lo_stop_gracefully(['NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                stop
            endif
        case(pm_vtetagrid)
            if ( gs%info%dim_temperature .gt. 0 .and. gs%info%dim_volume .gt. 0 .and. gs%info%dim_eta .gt. 0 ) then
                ! V-T-eta grid
                ge%gridtype=pm_vtetagrid !2
                if ( mw%talk ) then
                    write(*,*) '          gridtype:',ge%gridtype,'(V-T-eta grid)'
                    write(*,*) '   temperature-dim:',gs%info%dim_temperature
                    write(*,*) '        volume-dim:',gs%info%dim_volume
                    write(*,*) '           eta-dim:',gs%info%dim_eta
                endif
            else
                call lo_stop_gracefully(['NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif
        end select
        
        ! Override gridtype if evalmode=4 was used (a-c grid)
        if ( evalmode .eq. 4 ) then
            ge%gridtype=pm_acgrid
            if ( mw%talk ) then
                write(*,*) '          gridtype:',ge%gridtype,'(a-c lattice parameter grid)'
                write(*,*) '             a-dim:',gs%info%dim_a
                write(*,*) '             c-dim:',gs%info%dim_c
            endif
        endif
        ! Override gridtype if evalmode=5 was used (a-c-T grid with T-dependent FCs)
        if ( evalmode .eq. 5 ) then
            ge%gridtype=pm_actgrid
            if ( mw%talk ) then
                write(*,*) '          gridtype:',ge%gridtype,'(a-c-T lattice parameter grid with T-dependent FCs)'
                write(*,*) '             a-dim:',gs%info%dim_a
                write(*,*) '             c-dim:',gs%info%dim_c
                write(*,*) '             T-dim:',gs%info%dim_temperature
            endif
        endif
    end block init

    ! Set the flat-ish grid
    setgrid: block
        real(flyt), dimension(:), allocatable :: dumvol
        real(flyt) :: f0,f1
        integer :: npts,i,j,k,l,ii,jj,kk

        npts=product(ge%pts_per_dim)
        lo_allocate(gridcoord(ge%ndim,npts))
        lo_allocate(gridenergy(8,npts))
        lo_allocate(gridind(ge%ndim,npts))
        gridcoord=0.0_flyt
        gridenergy=0.0_flyt
        gridind=0

        ! build the grid
        select case(ge%gridtype)
        case(pm_vgrid) ! V-grid (QHA mode: interpolate FC in volume, evaluate at multiple T)
            ge%grid_V%nv=ge%pts_per_dim( gs%info%dim_volume )
            ! For pure QHA, we need a temperature range for evaluation
            ! Read it from the next lines in infile.evalpoints
            ge%grid_V%nt=100  ! Default number of temperature points
            
            allocate( ge%grid_V%volume( ge%grid_V%nv ) )
            allocate( ge%grid_V%temperature( ge%grid_V%nt ) )
            allocate( ge%grid_V%U     ( ge%grid_V%nv, ge%grid_V%nt ) )
            allocate( ge%grid_V%U0    ( ge%grid_V%nv, ge%grid_V%nt ) )
            allocate( ge%grid_V%fph   ( ge%grid_V%nv, ge%grid_V%nt ) )
            allocate( ge%grid_V%ah3   ( ge%grid_V%nv, ge%grid_V%nt ) )
            allocate( ge%grid_V%ah4   ( ge%grid_V%nv, ge%grid_V%nt ) )
            ge%grid_V%volume=0.0_flyt
            ge%grid_V%temperature=0.0_flyt
            ge%grid_V%U     =0.0_flyt
            ge%grid_V%U0    =0.0_flyt
            ge%grid_V%fph   =0.0_flyt
            ge%grid_V%ah3   =0.0_flyt
            ge%grid_V%ah4   =0.0_flyt

            ! fix volume axis
            f0=ge%min_per_dim( gs%info%dim_volume )
            f1=ge%max_per_dim( gs%info%dim_volume )
            select case( pointspacing(gs%info%dim_volume) )
            case('lin') ! linearly spaced in volume
                call lo_linspace(f0,f1,ge%grid_V%volume)
            case('den') ! linearly spaced in density
                call lo_linspace(1.0_flyt/f1,1.0_flyt/f0,ge%grid_V%volume)
                ge%grid_V%volume=1.0_flyt/ge%grid_V%volume
            end select

            ! Temperature grid for QHA evaluation (0 to 2000K default)
            call lo_linspace(1.0_flyt,2000.0_flyt,ge%grid_V%temperature)

            ! Set up the grid coordinates (only volume dimension)
            npts=ge%grid_V%nv
            deallocate(gridcoord)
            deallocate(gridenergy)
            deallocate(gridind)
            lo_allocate(gridcoord(ge%ndim,npts))
            lo_allocate(gridenergy(8,npts))
            lo_allocate(gridind(ge%ndim,npts))
            gridcoord=0.0_flyt
            gridenergy=0.0_flyt
            gridind=0
            do i=1,ge%grid_V%nv
                gridind(1,i)=i
                gridcoord(gs%info%dim_volume,i)=ge%grid_V%volume(i)
            enddo

            if ( mw%talk ) then
                write(*,*) '           volumes: ',tochar(minval(ge%grid_V%volume*lo_volume_bohr_to_A)),' -> ',tochar(maxval(ge%grid_V%volume*lo_volume_bohr_to_A)),' with ',tochar(ge%grid_V%nv),' points'
                write(*,*) '      temperatures: ',tochar(minval(ge%grid_V%temperature)),' -> ',tochar(maxval(ge%grid_V%temperature)),' with ',tochar(ge%grid_V%nt),' points (for QHA evaluation)'
            endif
        case(pm_vtgrid) ! V-T grid
            ge%grid_VT%nv=ge%pts_per_dim( gs%info%dim_volume )
            ge%grid_VT%nt=ge%pts_per_dim( gs%info%dim_temperature )
            allocate( ge%grid_VT%U     ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%U0    ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%fph   ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%ah3   ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%ah4   ( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%qh_fph( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%qh_ah3( ge%grid_VT%nv, ge%grid_VT%nt ) )
            allocate( ge%grid_VT%qh_ah4( ge%grid_VT%nv, ge%grid_VT%nt ) )
            ge%grid_VT%U     =0.0_flyt
            ge%grid_VT%U0    =0.0_flyt
            ge%grid_VT%fph   =0.0_flyt
            ge%grid_VT%ah3   =0.0_flyt
            ge%grid_VT%ah4   =0.0_flyt
            ge%grid_VT%qh_fph=0.0_flyt
            ge%grid_VT%qh_ah3=0.0_flyt
            ge%grid_VT%qh_ah4=0.0_flyt

            lo_allocate(ge%grid_VT%volume( ge%grid_VT%nv,ge%grid_VT%nt ))
            lo_allocate(ge%grid_VT%temperature( ge%grid_VT%nt ))
            ge%grid_VT%volume=0.0_flyt
            ge%grid_VT%temperature=0.0_flyt
            lo_allocate(dumvol( ge%grid_VT%nv ))

            ! fix volume axis
            select case(evalmode)
            case(1)
                ! same volume for all temperatures
                f0=ge%min_per_dim( gs%info%dim_volume )
                f1=ge%max_per_dim( gs%info%dim_volume )
                select case( pointspacing(gs%info%dim_volume) )
                case('lin') ! linearly spaced in volume
                    call lo_linspace(f0,f1,dumvol)
                case('den') ! linearly spaced in density
                    call lo_linspace(1.0_flyt/f1,1.0_flyt/f0,dumvol)
                    dumvol=1.0_flyt/dumvol
                end select
                do i=1,ge%grid_VT%nt
                    ge%grid_VT%volume(:,i)=dumvol
                enddo
            case(3)
                ! different volumes for different temperatures
                do i=1,ge%grid_VT%nt
                    f0=adapt_minvol(i)
                    f1=adapt_maxvol(i)
                    select case( pointspacing(gs%info%dim_volume) )
                    case('lin') ! linearly spaced in volume
                        call lo_linspace(f0,f1,dumvol)
                    case('den') ! linearly spaced in density
                        call lo_linspace(1.0_flyt/f1,1.0_flyt/f0,dumvol)
                        dumvol=1.0_flyt/dumvol
                    end select
                    ge%grid_VT%volume(:,i)=dumvol
                enddo
            case default
            end select

            ! fix temperature axis
            f0=ge%min_per_dim( gs%info%dim_temperature )
            f1=ge%max_per_dim( gs%info%dim_temperature )
            select case( pointspacing(gs%info%dim_temperature) )
            case('lin') ! linearly spaced in temperature
                call lo_linspace(f0,f1,ge%grid_VT%temperature)
            case('log') ! log=spaced
                call logspace(f0,f1,ge%grid_VT%temperature)
            end select
            if ( mw%talk ) then
                write(*,*) '      temperatures: ',tochar(minval(ge%grid_VT%temperature)),' -> ',tochar(maxval(ge%grid_VT%temperature)),' with ',tochar(ge%grid_VT%nt),' points'
                write(*,*) '           volumes: ',tochar(minval(ge%grid_VT%volume*lo_volume_bohr_to_A)),' -> ',tochar(maxval(ge%grid_VT%volume*lo_volume_bohr_to_A)),' with ',tochar(ge%grid_VT%nv),' points'
            endif
            l=0
            ii=gs%info%dim_volume
            jj=gs%info%dim_temperature
            do i=1,ge%grid_VT%nv
            do j=1,ge%grid_VT%nt
                l=l+1
                gridind(:,l)=[i,j]
                gridcoord(ii,l)=ge%grid_VT%volume( i,j )
                gridcoord(jj,l)=ge%grid_VT%temperature( j )
            enddo
            enddo
        case(pm_vtetagrid) ! V-T-eta grid
            ge%grid_VTeta%nv=ge%pts_per_dim( gs%info%dim_volume )
            ge%grid_VTeta%nt=ge%pts_per_dim( gs%info%dim_temperature )
            ge%grid_VTeta%neta=ge%pts_per_dim( gs%info%dim_eta )
            allocate( ge%grid_VTeta%U     ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%U0    ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%fph   ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%ah3   ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%ah4   ( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%qh_fph( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%qh_ah3( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            allocate( ge%grid_VTeta%qh_ah4( ge%grid_VTeta%nv, ge%grid_VTeta%nt, ge%grid_VTeta%neta ) )
            ge%grid_VTeta%U     =0.0_flyt
            ge%grid_VTeta%U0    =0.0_flyt
            ge%grid_VTeta%fph   =0.0_flyt
            ge%grid_VTeta%ah3   =0.0_flyt
            ge%grid_VTeta%ah4   =0.0_flyt
            ge%grid_VTeta%qh_fph=0.0_flyt
            ge%grid_VTeta%qh_ah3=0.0_flyt
            ge%grid_VTeta%qh_ah4=0.0_flyt

            allocate(ge%grid_VTeta%volume( ge%grid_VTeta%nv ))
            allocate(ge%grid_VTeta%temperature( ge%grid_VTeta%nt ))
            allocate(ge%grid_VTeta%eta( ge%grid_VTeta%neta ))
            ! fix volume axis
            f0=ge%min_per_dim( gs%info%dim_volume )
            f1=ge%max_per_dim( gs%info%dim_volume )
            select case( pointspacing(gs%info%dim_volume) )
            case('lin') ! linearly spaced in volume
                call lo_linspace(f0,f1,ge%grid_VTeta%volume)
            case('den') ! linearly spaced in density
                call lo_linspace(1.0_flyt/f1,1.0_flyt/f0,ge%grid_VTeta%volume)
                ge%grid_VTeta%volume=1.0_flyt/ge%grid_VTeta%volume
            end select
            ! fix temperature axis
            f0=ge%min_per_dim( gs%info%dim_temperature )
            f1=ge%max_per_dim( gs%info%dim_temperature )
            select case( pointspacing(gs%info%dim_temperature) )
            case('lin') ! linearly spaced in temperature
                call lo_linspace(f0,f1,ge%grid_VTeta%temperature)
            case('log') ! log=spaced
                call logspace(f0,f1,ge%grid_VTeta%temperature)
            end select
            ! fix eta axis
            f0=ge%min_per_dim( gs%info%dim_eta )
            f1=ge%max_per_dim( gs%info%dim_eta )
            select case( pointspacing(gs%info%dim_eta) )
            case('lin') ! linearly spaced in eta
                call lo_linspace(f0,f1,ge%grid_VTeta%eta)
            case('log') ! log=spaced
                call logspace(f0,f1,ge%grid_VTeta%eta)
            end select
            ! and get the grid coordinates
            l=0
            ii=gs%info%dim_volume
            jj=gs%info%dim_temperature
            kk=gs%info%dim_eta
            do i=1,ge%grid_VTeta%nv
            do j=1,ge%grid_VTeta%nt
            do k=1,ge%grid_VTeta%neta
                l=l+1
                gridind(:,l)=[i,j,k]
                gridcoord(ii,l)=ge%grid_VTeta%volume( i )
                gridcoord(jj,l)=ge%grid_VTeta%temperature( j )
                gridcoord(kk,l)=ge%grid_VTeta%eta( k )
            enddo
            enddo
            enddo

            if ( mw%talk ) then
                write(*,*) '      temperatures: ',tochar(minval(ge%grid_VTeta%temperature)),' -> ',tochar(maxval(ge%grid_VTeta%temperature)),' with ',tochar(ge%grid_VTeta%nt),' points'
                write(*,*) '           volumes: ',tochar(minval(ge%grid_VTeta%volume)),' -> ',tochar(maxval(ge%grid_VTeta%volume)),' with ',tochar(ge%grid_VTeta%nv),' points'
                write(*,*) '               eta: ',tochar(minval(ge%grid_VTeta%eta)),' -> ',tochar(maxval(ge%grid_VTeta%eta)),' with ',tochar(ge%grid_VTeta%neta),' points'
            endif
        case(pm_acgrid) ! a-c lattice parameter grid
            ge%grid_AC%na=ge%pts_per_dim( gs%info%dim_a )
            ge%grid_AC%nc=ge%pts_per_dim( gs%info%dim_c )
            ge%grid_AC%nt=100  ! Default number of temperature points for evaluation
            
            ! Allocate storage
            allocate( ge%grid_AC%a_values( ge%grid_AC%na ) )
            allocate( ge%grid_AC%c_values( ge%grid_AC%nc ) )
            allocate( ge%grid_AC%temperature( ge%grid_AC%nt ) )
            allocate( ge%grid_AC%U     ( ge%grid_AC%na, ge%grid_AC%nc, ge%grid_AC%nt ) )
            allocate( ge%grid_AC%U0    ( ge%grid_AC%na, ge%grid_AC%nc, ge%grid_AC%nt ) )
            allocate( ge%grid_AC%fph   ( ge%grid_AC%na, ge%grid_AC%nc, ge%grid_AC%nt ) )
            allocate( ge%grid_AC%ah3   ( ge%grid_AC%na, ge%grid_AC%nc, ge%grid_AC%nt ) )
            allocate( ge%grid_AC%ah4   ( ge%grid_AC%na, ge%grid_AC%nc, ge%grid_AC%nt ) )
            allocate( ge%grid_AC%Ftot  ( ge%grid_AC%na, ge%grid_AC%nc, ge%grid_AC%nt ) )
            allocate( ge%grid_AC%a_min ( ge%grid_AC%nt ) )
            allocate( ge%grid_AC%c_min ( ge%grid_AC%nt ) )
            allocate( ge%grid_AC%F_min ( ge%grid_AC%nt ) )
            ! For 4th order 2D polynomial: (a^0 to a^4) x (c^0 to c^4) = 25 terms
            allocate( ge%grid_AC%poly_coeffs( 25, ge%grid_AC%nt ) )
            
            ge%grid_AC%a_values=0.0_flyt
            ge%grid_AC%c_values=0.0_flyt
            ge%grid_AC%temperature=0.0_flyt
            ge%grid_AC%U     =0.0_flyt
            ge%grid_AC%U0    =0.0_flyt
            ge%grid_AC%fph   =0.0_flyt
            ge%grid_AC%ah3   =0.0_flyt
            ge%grid_AC%ah4   =0.0_flyt
            ge%grid_AC%Ftot  =0.0_flyt
            ge%grid_AC%a_min =0.0_flyt
            ge%grid_AC%c_min =0.0_flyt
            ge%grid_AC%F_min =0.0_flyt
            ge%grid_AC%poly_coeffs=0.0_flyt

            ! fix a axis (convert from Angstrom to Bohr)
            f0=ge%min_per_dim( gs%info%dim_a )
            f1=ge%max_per_dim( gs%info%dim_a )
            select case( pointspacing(gs%info%dim_a) )
            case('lin')
                call lo_linspace(f0,f1,ge%grid_AC%a_values)
            case default
                call lo_linspace(f0,f1,ge%grid_AC%a_values)
            end select
            
            ! fix c axis (convert from Angstrom to Bohr)
            f0=ge%min_per_dim( gs%info%dim_c )
            f1=ge%max_per_dim( gs%info%dim_c )
            select case( pointspacing(gs%info%dim_c) )
            case('lin')
                call lo_linspace(f0,f1,ge%grid_AC%c_values)
            case default
                call lo_linspace(f0,f1,ge%grid_AC%c_values)
            end select

            ! Temperature grid for QHA evaluation (0 to 2000K default)
            call lo_linspace(1.0_flyt,2000.0_flyt,ge%grid_AC%temperature)

            ! Set up the grid coordinates (only a and c dimensions)
            npts=ge%grid_AC%na * ge%grid_AC%nc
            deallocate(gridcoord)
            deallocate(gridenergy)
            deallocate(gridind)
            lo_allocate(gridcoord(ge%ndim,npts))
            lo_allocate(gridenergy(8,npts))
            lo_allocate(gridind(ge%ndim,npts))
            gridcoord=0.0_flyt
            gridenergy=0.0_flyt
            gridind=0
            l=0
            ii=gs%info%dim_a
            jj=gs%info%dim_c
            do i=1,ge%grid_AC%na
            do j=1,ge%grid_AC%nc
                l=l+1
                gridind(:,l)=[i,j]
                gridcoord(ii,l)=ge%grid_AC%a_values(i)
                gridcoord(jj,l)=ge%grid_AC%c_values(j)
            enddo
            enddo

            if ( mw%talk ) then
                write(*,*) '         a values: ',tochar(minval(ge%grid_AC%a_values)),' -> ',tochar(maxval(ge%grid_AC%a_values)),' with ',tochar(ge%grid_AC%na),' points (Angstrom)'
                write(*,*) '         c values: ',tochar(minval(ge%grid_AC%c_values)),' -> ',tochar(maxval(ge%grid_AC%c_values)),' with ',tochar(ge%grid_AC%nc),' points (Angstrom)'
                write(*,*) '      temperatures: ',tochar(minval(ge%grid_AC%temperature)),' -> ',tochar(maxval(ge%grid_AC%temperature)),' with ',tochar(ge%grid_AC%nt),' points (for QHA evaluation)'
            endif
        case(pm_actgrid) ! a-c-T lattice parameter grid with T-dependent FCs
            ge%grid_ACT%na=ge%pts_per_dim( gs%info%dim_a )
            ge%grid_ACT%nc=ge%pts_per_dim( gs%info%dim_c )
            ge%grid_ACT%nt=ge%pts_per_dim( gs%info%dim_temperature )
            
            ! Allocate storage
            allocate( ge%grid_ACT%a_values( ge%grid_ACT%na ) )
            allocate( ge%grid_ACT%c_values( ge%grid_ACT%nc ) )
            allocate( ge%grid_ACT%temperature( ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%U     ( ge%grid_ACT%na, ge%grid_ACT%nc, ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%U0    ( ge%grid_ACT%na, ge%grid_ACT%nc, ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%fph   ( ge%grid_ACT%na, ge%grid_ACT%nc, ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%ah3   ( ge%grid_ACT%na, ge%grid_ACT%nc, ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%ah4   ( ge%grid_ACT%na, ge%grid_ACT%nc, ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%Ftot  ( ge%grid_ACT%na, ge%grid_ACT%nc, ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%a_min ( ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%c_min ( ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%F_min ( ge%grid_ACT%nt ) )
            allocate( ge%grid_ACT%poly_coeffs( 25, ge%grid_ACT%nt ) )
            
            ge%grid_ACT%a_values=0.0_flyt
            ge%grid_ACT%c_values=0.0_flyt
            ge%grid_ACT%temperature=0.0_flyt
            ge%grid_ACT%U     =0.0_flyt
            ge%grid_ACT%U0    =0.0_flyt
            ge%grid_ACT%fph   =0.0_flyt
            ge%grid_ACT%ah3   =0.0_flyt
            ge%grid_ACT%ah4   =0.0_flyt
            ge%grid_ACT%Ftot  =0.0_flyt
            ge%grid_ACT%a_min =0.0_flyt
            ge%grid_ACT%c_min =0.0_flyt
            ge%grid_ACT%F_min =0.0_flyt
            ge%grid_ACT%poly_coeffs=0.0_flyt

            ! fix a axis
            f0=ge%min_per_dim( gs%info%dim_a )
            f1=ge%max_per_dim( gs%info%dim_a )
            call lo_linspace(f0,f1,ge%grid_ACT%a_values)
            
            ! fix c axis
            f0=ge%min_per_dim( gs%info%dim_c )
            f1=ge%max_per_dim( gs%info%dim_c )
            call lo_linspace(f0,f1,ge%grid_ACT%c_values)

            ! fix temperature axis
            f0=ge%min_per_dim( gs%info%dim_temperature )
            f1=ge%max_per_dim( gs%info%dim_temperature )
            call lo_linspace(f0,f1,ge%grid_ACT%temperature)

            ! Set up the grid coordinates
            npts=ge%grid_ACT%na * ge%grid_ACT%nc * ge%grid_ACT%nt
            deallocate(gridcoord)
            deallocate(gridenergy)
            deallocate(gridind)
            lo_allocate(gridcoord(ge%ndim,npts))
            lo_allocate(gridenergy(8,npts))
            lo_allocate(gridind(ge%ndim,npts))
            gridcoord=0.0_flyt
            gridenergy=0.0_flyt
            gridind=0
            l=0
            ii=gs%info%dim_a
            jj=gs%info%dim_c
            kk=gs%info%dim_temperature
            do i=1,ge%grid_ACT%na
            do j=1,ge%grid_ACT%nc
            do k=1,ge%grid_ACT%nt
                l=l+1
                gridind(:,l)=[i,j,k]
                gridcoord(ii,l)=ge%grid_ACT%a_values(i)
                gridcoord(jj,l)=ge%grid_ACT%c_values(j)
                gridcoord(kk,l)=ge%grid_ACT%temperature(k)
            enddo
            enddo
            enddo

            if ( mw%talk ) then
                write(*,*) '         a values: ',tochar(minval(ge%grid_ACT%a_values)),' -> ',tochar(maxval(ge%grid_ACT%a_values)),' with ',tochar(ge%grid_ACT%na),' points (Angstrom)'
                write(*,*) '         c values: ',tochar(minval(ge%grid_ACT%c_values)),' -> ',tochar(maxval(ge%grid_ACT%c_values)),' with ',tochar(ge%grid_ACT%nc),' points (Angstrom)'
                write(*,*) '      temperatures: ',tochar(minval(ge%grid_ACT%temperature)),' -> ',tochar(maxval(ge%grid_ACT%temperature)),' with ',tochar(ge%grid_ACT%nt),' points (T-dependent FCs)'
            endif
        case default
            call lo_stop_gracefully(['NOT DONE'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
        end select
    end block setgrid

    ! if ( dumpgrid .and. mw%talk ) then
    ! dgrid: block
    !     type(lo_crystalstructure) :: p
    !     type(lo_forceconstant_secondorder) :: fc
    !     type(lo_forceconstant_thirdorder) :: fct
    !     type(lo_forceconstant_fourthorder) :: fcf
    !     type(lo_jij_secondorder) :: jij
    !     real(flyt), dimension(:,:), allocatable :: pairconstraints
    !     integer :: i,j,l,nconstr
    !
    !     write(*,*) '... dumping input files across the entire grid'
    !
    !     select case(ge%gridtype)
    !     case(pm_vtgrid)
    !         l=0
    !         do i=1,ge%grid_VT%nv
    !         do j=1,ge%grid_VT%nt
    !             l=l+1
    !             ! Get the structure
    !             call gs%structure%interpolate(gridcoord(:,l),p,gs%info%dim_volume)
    !             call lo_secondorder_rot_herm_huang( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
    !             if ( nconstr .gt. 0 ) then
    !                 call gs%eval(map,gridcoord(:,l),pairconstraints)
    !             else
    !                 call gs%eval(map,gridcoord(:,l))
    !             endif
    !             p%info%title='interp struct V= '//tochar(ge%grid_VT%volume(i,j)*lo_volume_bohr_to_A,ndecimals=11)//' T= '//tochar(ge%grid_VT%temperature(j),ndecimals=10)
    !             call p%writetofile('uc_'//tochar(i)//'_'//tochar(j),1)
    !
    !             ! Now dump all the component guys
    !             if ( gs%info%secondorder ) then
    !                 call map%get_secondorder_forceconstant(p,fc,mem,-1)
    !                 call fc%writetofile(p,'fc2_'//tochar(i)//'_'//tochar(j))
    !             endif
    !             if ( gs%info%thirdorder ) then
    !                 call map%get_thirdorder_forceconstant(p,fct)
    !                 call fct%writetofile(p,'fc3_'//tochar(i)//'_'//tochar(j))
    !             endif
    !             if ( gs%info%fourthorder ) then
    !                 call map%get_fourthorder_forceconstant(p,fcf)
    !                 call fcf%writetofile(p,'fc4_'//tochar(i)//'_'//tochar(j))
    !             endif
    !             ! if ( gs%info%magnetic_pair_interactions ) then
    !             !     call map%get_secondorder_jij(p,jij)
    !             !     call jij%writetofile(p,'jij_'//tochar(i)//'_'//tochar(j))
    !             ! endif
    !         enddo
    !         enddo
    !     case default
    !         write(*,*) 'FIXME DUMPGRID'
    !         stop
    !     end select
    !
    ! end block dgrid
    ! endif

    ! Evaluate stuff - handle V-grid (QHA) case specially
    select case(ge%gridtype)
    case(pm_vgrid)
        call evaluate_vgrid_qha(ge,gs,map,qgrid_harm,qgrid_anharm,mw,mem,verbosity)
    case(pm_acgrid)
        call evaluate_acgrid_qha(ge,gs,map,qgrid_harm,qgrid_anharm,mw,mem,verbosity)
    case(pm_actgrid)
        call evaluate_actgrid(ge,gs,map,qgrid_harm,qgrid_anharm,mw,mem,verbosity)
    case default
        call evaluate_vt_grid(ge,gs,map,gridcoord,gridind,gridenergy,qgrid_harm,qgrid_anharm,quasiharmonic,mw,mem,verbosity)
    end select

    ! then finalize it
    select case(ge%gridtype)
    case(pm_vgrid)
        call Vfinalize( ge%grid_V,gs,map,'outfile.interpolated_free_energy.hdf5',qgrid_harm,mw,mem )
    case(pm_acgrid)
        call ACfinalize( ge%grid_AC,gs,map,'outfile.interpolated_free_energy.hdf5',qgrid_harm,mw,mem )
    case(pm_actgrid)
        call ACTfinalize( ge%grid_ACT,gs,map,'outfile.interpolated_free_energy.hdf5',qgrid_harm,mw,mem )
    case(pm_vtgrid)
        call VTfinalize( ge%grid_VT,gs,map,'outfile.interpolated_free_energy.hdf5',qgrid_harm,quasiharmonic,pressurestep,dumpgrid,mw,mem )
    case(pm_vtetagrid)
        call VTetafinalize( ge%grid_VTeta,gs,map,'outfile.interpolated_free_energy.hdf5',qgrid_harm,quasiharmonic,pressurestep,mw,mem )
    end select
end subroutine

!> Evaluate energies for V-grid (pure QHA mode)
subroutine evaluate_vgrid_qha(ge,gs,map,qgrid_harm,qgrid_anharm,mw,mem,verbosity)
    !> grid energy
    type(lo_gridenergy), intent(inout) :: ge
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> q-grid for harmonic
    integer, dimension(3), intent(in) :: qgrid_harm
    !> q-grid for anharmonic
    integer, dimension(3), intent(in) :: qgrid_anharm
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> verbosity
    integer, intent(in) :: verbosity

    class(lo_qpoint_mesh), allocatable :: qp
    type(lo_phonon_dispersions) :: dr
    type(lo_crystalstructure) :: p
    type(lo_forceconstant_secondorder) :: fc
    real(flyt), dimension(:,:), allocatable :: pairconstraints
    real(flyt), dimension(1) :: depvar
    real(flyt) :: t0,volume,fph
    integer :: iv,it,nconstr

    if ( mw%talk ) then
        t0=walltime()
        write(*,*) ''
        write(*,*) 'Evaluating QHA free energy (volume-only mode)'
        call lo_progressbar_init()
    endif

    ! Loop over volumes
    do iv=1,ge%grid_V%nv
        if ( mod(iv,mw%n) .ne. mw%r ) cycle

        volume=ge%grid_V%volume(iv)
        depvar(1)=volume

        ! Get structure at this volume
        call gs%structure%interpolate(depvar,p,gs%info%dim_volume)
        call p%classify('wedge',timereversal=.true.)

        ! Get force constants at this volume
        call megafit_secondorder_constraints( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
        if ( nconstr .gt. 0 ) then
            call gs%eval(map,depvar,pairconstraints)
        else
            call gs%eval(map,depvar)
        endif
        call map%get_secondorder_forceconstant(p,fc,mem,-1)

        ! Get static energy from EOS
        select type(eos=>gs%eos)
        class is(lo_eos_1d)
            ge%grid_V%U(iv,:)=eos%energy_from_volume( volume )
        class default
            ge%grid_V%U(iv,:)=0.0_flyt
        end select

        ! Get delta U0 from interpolation
        call gs%energy%interpolate( depvar, ge%grid_V%U0(iv,1) )
        ge%grid_V%U0(iv,:)=ge%grid_V%U0(iv,1)

        ! Generate q-mesh and dispersions once per volume
        call lo_generate_qmesh(qp,p,qgrid_harm,'fft',timereversal=.true.,headrankonly=.false.,mw=mw,mem=mem,verbosity=-1)
        call dr%generate(qp,fc,p,mw=mw,mem=mem,verbosity=-1)

        ! Check for unstable modes
        if ( dr%omega_min .lt. lo_freqtol ) then
            ge%grid_V%fph(iv,:) = 123456789.0_flyt
            ge%grid_V%ah3(iv,:) = 0.0_flyt
            ge%grid_V%ah4(iv,:) = 0.0_flyt
        else
            ! Now evaluate phonon free energy at each temperature (this is fast - just Bose-Einstein)
            do it=1,ge%grid_V%nt
                ge%grid_V%fph(iv,it) = dr%phonon_free_energy(ge%grid_V%temperature(it))
            enddo
            ge%grid_V%ah3(iv,:) = 0.0_flyt
            ge%grid_V%ah4(iv,:) = 0.0_flyt
        endif

        if ( mw%talk .and. iv .lt. ge%grid_V%nv ) then
            call lo_progressbar(' ... QHA free energy',iv,ge%grid_V%nv,walltime()-t0)
        endif
    enddo

    if ( mw%talk ) call lo_progressbar(' ... QHA free energy',ge%grid_V%nv,ge%grid_V%nv,walltime()-t0)

    ! Collect results across MPI ranks
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_V%U,  ge%grid_V%nv*ge%grid_V%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_V%U0, ge%grid_V%nv*ge%grid_V%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_V%fph,ge%grid_V%nv*ge%grid_V%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_V%ah3,ge%grid_V%nv*ge%grid_V%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_V%ah4,ge%grid_V%nv*ge%grid_V%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
end subroutine

!> Evaluate energies for V-T and V-T-eta grids
subroutine evaluate_vt_grid(ge,gs,map,gridcoord,gridind,gridenergy,qgrid_harm,qgrid_anharm,quasiharmonic,mw,mem,verbosity)
    !> grid energy
    type(lo_gridenergy), intent(inout) :: ge
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> grid coordinates
    real(flyt), dimension(:,:), intent(in) :: gridcoord
    !> grid indices
    integer, dimension(:,:), intent(in) :: gridind
    !> energies
    real(flyt), dimension(:,:), intent(inout) :: gridenergy
    !> q-grid for harmonic
    integer, dimension(3), intent(in) :: qgrid_harm
    !> q-grid for anharmonic
    integer, dimension(3), intent(in) :: qgrid_anharm
    !> evaluate quasiharmonic
    logical, intent(in) :: quasiharmonic
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> verbosity
    integer, intent(in) :: verbosity

    ! Original evalen block code
    evalen: block
        type(lo_mpi_helper) :: ml
        real(flyt), dimension(gs%ndim) :: dcrd
        real(flyt) :: t0,timer_anharmonic,temperature,volume,eta
        integer :: i,l,ii,jj,kk,npts

        ! Split communicator to serial things
        call mw%split(ml,mw%r,__FILE__,__LINE__)

        npts=size(gridcoord,2)
        ! First intepolate U0 and static energy
        if ( mw%talk ) then
            t0=walltime()
            call lo_progressbar_init()
        endif
        do i=1,npts
            ! to make it MPI parallel
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            ! interpolate the internal energy
            call gs%energy%interpolate( gridcoord(:,i),gridenergy(2,i) )
            ! add the energy from the provided equation of state
            select type(eos=>gs%eos)
            class is(lo_eos_1d)
                volume=gridcoord(gs%info%dim_volume,i)
                gridenergy(1,i)=eos%energy_from_volume( volume )
            class is(lo_eos_2d)
                volume=gridcoord(gs%info%dim_volume,i)
                eta=gridcoord(gs%info%dim_eta,i)
                gridenergy(1,i)=eos%energy_from_volume_eta( volume,eta )
            class default
                ! with no equation of state, do nothing.
                gridenergy(1,i)=0.0_flyt
            end select
            if ( mw%talk .and. i .lt. npts ) call lo_progressbar(' ... interpolating U0',i,npts,walltime()-t0)
        enddo
        if ( mw%talk ) call lo_progressbar(' ... interpolating U0',npts,npts,walltime()-t0)

        ! Now the phonon free energy
        if ( mw%talk ) then
            t0=walltime()
            call lo_progressbar_init()
        endif
        do i=1,npts
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            temperature=gridcoord( gs%info%dim_temperature,i )
            call phonon_free_energy_for_single_point( gs,map,gridcoord(:,i),qgrid_harm,temperature,ml,mem,gridenergy(3,i) )
            if ( quasiharmonic ) then
                dcrd=gridcoord(:,i)
                dcrd( gs%info%dim_temperature )=0.0_flyt
                call phonon_free_energy_for_single_point( gs,map,dcrd,qgrid_harm,temperature,ml,mem,gridenergy(6,i) )
            endif
            if ( mw%talk .and. i .lt. npts ) call lo_progressbar(' ... phonon free energy',i,npts,walltime()-t0)
        enddo
        if ( mw%talk ) call lo_progressbar(' ... phonon free energy',npts,npts,walltime()-t0)

        ! And destroy the temporary communicators
        call ml%free(__FILE__,__LINE__)

        ! Add it up over ranks, what we have so far
        call mpi_allreduce(MPI_IN_PLACE,gridenergy,8*npts,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error )

        ! Now the anharmonic free energy. This is parallel internally instead.
        if ( map%have_fc_triplet ) then
        if ( map%have_fc_quartet ) then
        if ( qgrid_anharm(1) .gt. 0 ) then
            if ( mw%talk ) then
                timer_anharmonic=walltime()
                t0=walltime()
            endif
            do i=1,npts
                temperature=gridcoord( gs%info%dim_temperature,i )
                call anharmonic_free_energy_for_single_point( gs,map,gridcoord(:,i),qgrid_anharm,temperature,&
                                                              gridenergy(4,i),gridenergy(5,i),mw,mem )
                if ( quasiharmonic ) then
                   dcrd=gridcoord(:,i)
                   dcrd( gs%info%dim_temperature )=0.0_flyt
                   call anharmonic_free_energy_for_single_point( gs,map,dcrd,qgrid_anharm,temperature,&
                                                                 gridenergy(7,i),gridenergy(8,i),mw,mem )
                endif
                if ( mw%talk ) then
                if ( walltime()-t0 .gt. timereport ) then
                    call lo_looptimer('... anharmonic free energy',timer_anharmonic,walltime(),i,npts)
                    t0=walltime()
                endif
                endif
            enddo
        endif
        endif
        endif

        ! Store it organized
        select case(ge%gridtype)
        case(pm_vtgrid) ! V-T grid
            do l=1,npts
                ii=gridind(1,l)
                jj=gridind(2,l)
                ge%grid_VT%U     ( ii,jj )=gridenergy(1,l)
                ge%grid_VT%U0    ( ii,jj )=gridenergy(2,l)
                ge%grid_VT%fph   ( ii,jj )=gridenergy(3,l)
                ge%grid_VT%ah3   ( ii,jj )=gridenergy(4,l)
                ge%grid_VT%ah4   ( ii,jj )=gridenergy(5,l)
                ge%grid_VT%qh_fph( ii,jj )=gridenergy(6,l)
                ge%grid_VT%qh_ah3( ii,jj )=gridenergy(7,l)
                ge%grid_VT%qh_ah4( ii,jj )=gridenergy(8,l)
            enddo
        case(pm_vtetagrid) ! V-T-eta grid
            do l=1,npts
                ii=gridind(1,l)
                jj=gridind(2,l)
                kk=gridind(3,l)
                ge%grid_VTeta%U     ( ii,jj,kk )=gridenergy(1,l)
                ge%grid_VTeta%U0    ( ii,jj,kk )=gridenergy(2,l)
                ge%grid_VTeta%fph   ( ii,jj,kk )=gridenergy(3,l)
                ge%grid_VTeta%ah3   ( ii,jj,kk )=gridenergy(4,l)
                ge%grid_VTeta%ah4   ( ii,jj,kk )=gridenergy(5,l)
                ge%grid_VTeta%qh_fph( ii,jj,kk )=gridenergy(6,l)
                ge%grid_VTeta%qh_ah3( ii,jj,kk )=gridenergy(7,l)
                ge%grid_VTeta%qh_ah4( ii,jj,kk )=gridenergy(8,l)
            enddo
        end select
    end block evalen
end subroutine

!> Finalize V-grid (QHA) output
subroutine Vfinalize(gr,gs,map,filename,qgrid,mw,mem)
    !> grid
    type(lo_gridenergy_vol), intent(in) :: gr
    !> gridsim
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> filename
    character(len=*), intent(in) :: filename
    !> q-grid
    integer, dimension(3), intent(in) :: qgrid
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:,:), allocatable :: Ftot
    integer :: iv,it

    if ( mw%talk ) then
        write(*,*) ''
        write(*,*) 'Writing QHA results to ',trim(filename)

        call h5%init(__FILE__,__LINE__)
        call h5%open_file('write',trim(filename))
        call h5%open_group('write','grid_QHA')

        ! Store axes
        call h5%store_data(gr%volume*lo_volume_Bohr_to_A, h5%group_id,'volumes',enhet='A^3/atom')
        call h5%store_data(gr%temperature,                h5%group_id,'temperatures',enhet='K')

        ! Store energies
        call h5%store_data(gr%U*lo_Hartree_to_eV,   h5%group_id,'static_internal_energy',enhet='eV/atom')
        call h5%store_data(gr%U0*lo_Hartree_to_eV,  h5%group_id,'delta_U0',enhet='eV/atom')
        call h5%store_data(gr%fph*lo_Hartree_to_eV, h5%group_id,'phonon_free_energy',enhet='eV/atom')

        ! Total Helmholtz free energy
        allocate(Ftot(gr%nv,gr%nt))
        Ftot = gr%U + gr%U0 + gr%fph
        call h5%store_data(Ftot*lo_Hartree_to_eV, h5%group_id,'Helmholtz_free_energy',enhet='eV/atom')
        deallocate(Ftot)

        call h5%close_group()
        call h5%close_file()

        write(*,*) '... done writing QHA output'
    endif
end subroutine

!> Evaluate energies for a-c grid (hexagonal/tetragonal QHA)
subroutine evaluate_acgrid_qha(ge,gs,map,qgrid_harm,qgrid_anharm,mw,mem,verbosity)
    !> grid energy
    type(lo_gridenergy), intent(inout) :: ge
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> q-grid for harmonic
    integer, dimension(3), intent(in) :: qgrid_harm
    !> q-grid for anharmonic
    integer, dimension(3), intent(in) :: qgrid_anharm
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> verbosity
    integer, intent(in) :: verbosity

    class(lo_qpoint_mesh), allocatable :: qp
    type(lo_phonon_dispersions) :: dr
    type(lo_crystalstructure) :: p
    type(lo_forceconstant_secondorder) :: fc
    type(lo_forceconstant_thirdorder) :: fct
    type(lo_forceconstant_fourthorder) :: fcf
    real(flyt), dimension(:,:), allocatable :: pairconstraints
    real(flyt), dimension(2) :: depvar
    real(flyt) :: t0,fph,f3,f4
    integer :: ia,ic,it,nconstr,ipt,npts
    logical :: have_anharmonic

    if ( mw%talk ) then
        t0=walltime()
        write(*,*) ''
        write(*,*) 'Evaluating QHA free energy (a-c lattice parameter grid)'
        call lo_progressbar_init()
    endif

    npts = ge%grid_AC%na * ge%grid_AC%nc
    
    ! Check if we have anharmonic force constants
    have_anharmonic = map%have_fc_triplet .and. map%have_fc_quartet .and. (qgrid_anharm(1) .gt. 0)
    if ( mw%talk .and. have_anharmonic ) then
        write(*,*) '... including anharmonic (3rd+4th order) contributions'
    endif

    ! Loop over a,c pairs
    ipt=0
    do ia=1,ge%grid_AC%na
    do ic=1,ge%grid_AC%nc
        ipt=ipt+1
        if ( mod(ipt-1,mw%n) .ne. mw%r ) cycle

        depvar(gs%info%dim_a)=ge%grid_AC%a_values(ia)
        depvar(gs%info%dim_c)=ge%grid_AC%c_values(ic)

        ! Get structure at this a,c
        call gs%structure%interpolate(depvar,p,-1)  ! -1 means no special volume handling
        call p%classify('wedge',timereversal=.true.)

        ! Get force constants at this a,c
        call megafit_secondorder_constraints( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
        if ( nconstr .gt. 0 ) then
            call gs%eval(map,depvar,pairconstraints)
        else
            call gs%eval(map,depvar)
        endif
        call map%get_secondorder_forceconstant(p,fc,mem,-1)

        ! Get static energy (interpolated)
        call gs%energy%interpolate( depvar, ge%grid_AC%U(ia,ic,1) )
        ge%grid_AC%U(ia,ic,:) = ge%grid_AC%U(ia,ic,1)

        ! Get delta U0 from interpolation (if available)
        ge%grid_AC%U0(ia,ic,:) = 0.0_flyt

        ! Generate q-mesh and dispersions once per (a,c) point
        call lo_generate_qmesh(qp,p,qgrid_harm,'fft',timereversal=.true.,headrankonly=.false.,mw=mw,mem=mem,verbosity=-1)
        call dr%generate(qp,fc,p,mw=mw,mem=mem,verbosity=-1)

        ! Check for unstable modes
        if ( dr%omega_min .lt. lo_freqtol ) then
            ge%grid_AC%fph(ia,ic,:) = 123456789.0_flyt
            ge%grid_AC%ah3(ia,ic,:) = 0.0_flyt
            ge%grid_AC%ah4(ia,ic,:) = 0.0_flyt
        else
            ! Evaluate phonon free energy at each temperature
            do it=1,ge%grid_AC%nt
                ge%grid_AC%fph(ia,ic,it) = dr%phonon_free_energy(ge%grid_AC%temperature(it))
            enddo
            
            ! Evaluate anharmonic contributions if available
            if ( have_anharmonic ) then
                call map%get_thirdorder_forceconstant(p,fct)
                call map%get_fourthorder_forceconstant(p,fcf)
                do it=1,ge%grid_AC%nt
                    select type(qp); type is(lo_fft_mesh)
                        call anharmonic_free_energy(p,fct,fcf,qp,dr,ge%grid_AC%temperature(it),f3,f4,mw,mem,verbosity=-1)
                    end select
                    ge%grid_AC%ah3(ia,ic,it) = f3
                    ge%grid_AC%ah4(ia,ic,it) = f4
                enddo
            else
                ge%grid_AC%ah3(ia,ic,:) = 0.0_flyt
                ge%grid_AC%ah4(ia,ic,:) = 0.0_flyt
            endif
        endif

        if ( mw%talk .and. ipt .lt. npts ) then
            call lo_progressbar(' ... a-c QHA free energy',ipt,npts,walltime()-t0)
        endif
    enddo
    enddo

    if ( mw%talk ) call lo_progressbar(' ... a-c QHA free energy',npts,npts,walltime()-t0)

    ! Collect results across MPI ranks
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_AC%U,  ge%grid_AC%na*ge%grid_AC%nc*ge%grid_AC%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_AC%U0, ge%grid_AC%na*ge%grid_AC%nc*ge%grid_AC%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_AC%fph,ge%grid_AC%na*ge%grid_AC%nc*ge%grid_AC%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_AC%ah3,ge%grid_AC%na*ge%grid_AC%nc*ge%grid_AC%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_AC%ah4,ge%grid_AC%na*ge%grid_AC%nc*ge%grid_AC%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)

    ! Compute total free energy (including anharmonic)
    ge%grid_AC%Ftot = ge%grid_AC%U + ge%grid_AC%U0 + ge%grid_AC%fph + ge%grid_AC%ah3 + ge%grid_AC%ah4
end subroutine

!> Finalize a-c grid output with 4th order polynomial fitting
subroutine ACfinalize(gr,gs,map,filename,qgrid,mw,mem)
    !> grid
    type(lo_gridenergy_ac), intent(inout) :: gr
    !> gridsim
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> filename
    character(len=*), intent(in) :: filename
    !> q-grid
    integer, dimension(3), intent(in) :: qgrid
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:,:,:), allocatable :: Ftot_eV
    real(r8), dimension(:,:), allocatable :: design_matrix
    real(r8), dimension(:), allocatable :: F_vec, coeffs
    real(r8) :: a_norm, c_norm, a_mean, c_mean, a_scale, c_scale
    real(r8) :: a_opt, c_opt, F_opt
    integer :: ia, ic, it, k, l, npts, ncoeffs, u
    integer :: ia_min, ic_min, pa, pc

    ncoeffs = 25  ! 5x5 for 4th order polynomial in 2D
    npts = gr%na * gr%nc

    if ( mw%talk ) then
        write(*,*) ''
        write(*,*) 'Fitting 4th order polynomial to free energy surface'

        ! Normalize a and c for numerical stability
        a_mean = 0.5_r8*(minval(gr%a_values) + maxval(gr%a_values))
        c_mean = 0.5_r8*(minval(gr%c_values) + maxval(gr%c_values))
        a_scale = 0.5_r8*(maxval(gr%a_values) - minval(gr%a_values))
        c_scale = 0.5_r8*(maxval(gr%c_values) - minval(gr%c_values))
        if (a_scale .lt. 1.0e-10_r8) a_scale = 1.0_r8
        if (c_scale .lt. 1.0e-10_r8) c_scale = 1.0_r8

        allocate(design_matrix(npts, ncoeffs))
        allocate(F_vec(npts))
        allocate(coeffs(ncoeffs))

        ! Build design matrix once
        l = 0
        do ia = 1, gr%na
        do ic = 1, gr%nc
            l = l + 1
            a_norm = (gr%a_values(ia) - a_mean) / a_scale
            c_norm = (gr%c_values(ic) - c_mean) / c_scale
            k = 0
            do pa = 0, 4
            do pc = 0, 4 - pa  ! Keep total order <= 4
                k = k + 1
                design_matrix(l, k) = (a_norm**pa) * (c_norm**pc)
            enddo
            enddo
            ! Fill remaining columns with zeros (if any)
            do while (k < ncoeffs)
                k = k + 1
                design_matrix(l, k) = 0.0_r8
            enddo
        enddo
        enddo

        ! For each temperature, fit and find minimum
        do it = 1, gr%nt
            ! Build F vector for this temperature (skip unstable points)
            l = 0
            do ia = 1, gr%na
            do ic = 1, gr%nc
                l = l + 1
                F_vec(l) = gr%Ftot(ia, ic, it)
            enddo
            enddo

            ! Check for unstable points
            if (any(F_vec .gt. 1.0e8_r8)) then
                gr%a_min(it) = -1.0_r8
                gr%c_min(it) = -1.0_r8
                gr%F_min(it) = 1.0e10_r8
                gr%poly_coeffs(:, it) = 0.0_r8
                cycle
            endif

            ! Least squares fit
            call lo_linear_least_squares(design_matrix, F_vec, coeffs)
            gr%poly_coeffs(:, it) = coeffs

            ! Find minimum by grid search (simple approach)
            F_opt = 1.0e10_r8
            a_opt = gr%a_values(1)
            c_opt = gr%c_values(1)
            do ia = 1, gr%na
            do ic = 1, gr%nc
                if (gr%Ftot(ia, ic, it) .lt. F_opt) then
                    F_opt = gr%Ftot(ia, ic, it)
                    a_opt = gr%a_values(ia)
                    c_opt = gr%c_values(ic)
                    ia_min = ia
                    ic_min = ic
                endif
            enddo
            enddo
            gr%a_min(it) = a_opt
            gr%c_min(it) = c_opt
            gr%F_min(it) = F_opt
        enddo

        deallocate(design_matrix, F_vec, coeffs)

        ! Write output to HDF5
        write(*,*) 'Writing a-c QHA results to ',trim(filename)

        call h5%init(__FILE__,__LINE__)
        call h5%open_file('write',trim(filename))
        call h5%open_group('write','grid_ac_QHA')

        ! Store axes
        call h5%store_data(gr%a_values, h5%group_id,'a_values',enhet='Angstrom')
        call h5%store_data(gr%c_values, h5%group_id,'c_values',enhet='Angstrom')
        call h5%store_data(gr%temperature, h5%group_id,'temperatures',enhet='K')

        ! Store energies
        call h5%store_data(gr%U*lo_Hartree_to_eV,   h5%group_id,'static_internal_energy',enhet='eV/atom')
        call h5%store_data(gr%U0*lo_Hartree_to_eV,  h5%group_id,'delta_U0',enhet='eV/atom')
        call h5%store_data(gr%fph*lo_Hartree_to_eV, h5%group_id,'phonon_free_energy',enhet='eV/atom')
        call h5%store_data(gr%ah3*lo_Hartree_to_eV, h5%group_id,'anharmonic_free_energy_3rd',enhet='eV/atom')
        call h5%store_data(gr%ah4*lo_Hartree_to_eV, h5%group_id,'anharmonic_free_energy_4th',enhet='eV/atom')
        call h5%store_data(gr%Ftot*lo_Hartree_to_eV, h5%group_id,'Helmholtz_free_energy',enhet='eV/atom')

        ! Store minimum finding results
        call h5%store_data(gr%a_min, h5%group_id,'a_equilibrium',enhet='Angstrom')
        call h5%store_data(gr%c_min, h5%group_id,'c_equilibrium',enhet='Angstrom')
        call h5%store_data(gr%F_min*lo_Hartree_to_eV, h5%group_id,'F_equilibrium',enhet='eV/atom')

        ! Store polynomial coefficients
        call h5%store_data(gr%poly_coeffs, h5%group_id,'polynomial_coefficients',enhet='')

        call h5%close_group()
        call h5%close_file()

        ! Also write a simple text file with equilibrium a(T), c(T)
        u = open_file('out','outfile.ac_equilibrium.dat')
        write(u,'(A)') '# Temperature (K)    a_eq (Angstrom)    c_eq (Angstrom)    F_min (eV/atom)'
        do it = 1, gr%nt
            if (gr%a_min(it) .gt. 0.0_r8) then
                write(u,'(4E20.10)') gr%temperature(it), gr%a_min(it), gr%c_min(it), gr%F_min(it)*lo_Hartree_to_eV
            endif
        enddo
        close(u)

        write(*,*) '... done writing a-c QHA output'
    endif
end subroutine

!> Evaluate energies for a-c-T grid with T-dependent force constants
subroutine evaluate_actgrid(ge,gs,map,qgrid_harm,qgrid_anharm,mw,mem,verbosity)
    !> grid energy
    type(lo_gridenergy), intent(inout) :: ge
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> q-grid for harmonic
    integer, dimension(3), intent(in) :: qgrid_harm
    !> q-grid for anharmonic
    integer, dimension(3), intent(in) :: qgrid_anharm
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> verbosity
    integer, intent(in) :: verbosity

    class(lo_qpoint_mesh), allocatable :: qp
    type(lo_phonon_dispersions) :: dr
    type(lo_crystalstructure) :: p
    type(lo_forceconstant_secondorder) :: fc
    type(lo_forceconstant_thirdorder) :: fct
    type(lo_forceconstant_fourthorder) :: fcf
    real(flyt), dimension(:,:), allocatable :: pairconstraints
    real(flyt), dimension(3) :: depvar
    real(flyt) :: t0,fph,f3,f4,temperature
    integer :: ia,ic,it,nconstr,ipt,npts
    logical :: have_anharmonic

    if ( mw%talk ) then
        t0=walltime()
        write(*,*) ''
        write(*,*) 'Evaluating free energy (a-c-T grid with T-dependent FCs)'
        call lo_progressbar_init()
    endif

    npts = ge%grid_ACT%na * ge%grid_ACT%nc * ge%grid_ACT%nt
    
    ! Check if we have anharmonic force constants
    have_anharmonic = map%have_fc_triplet .and. map%have_fc_quartet .and. (qgrid_anharm(1) .gt. 0)
    if ( mw%talk .and. have_anharmonic ) then
        write(*,*) '... including anharmonic (3rd+4th order) contributions'
    endif

    ! Loop over all (a,c,T) points
    ipt=0
    do ia=1,ge%grid_ACT%na
    do ic=1,ge%grid_ACT%nc
    do it=1,ge%grid_ACT%nt
        ipt=ipt+1
        if ( mod(ipt-1,mw%n) .ne. mw%r ) cycle

        temperature = ge%grid_ACT%temperature(it)
        depvar(gs%info%dim_a)=ge%grid_ACT%a_values(ia)
        depvar(gs%info%dim_c)=ge%grid_ACT%c_values(ic)
        depvar(gs%info%dim_temperature)=temperature

        ! Get structure at this (a,c,T)
        call gs%structure%interpolate(depvar,p,-1)
        call p%classify('wedge',timereversal=.true.)

        ! Get force constants at this (a,c,T) - FCs are now T-dependent!
        call megafit_secondorder_constraints( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
        if ( nconstr .gt. 0 ) then
            call gs%eval(map,depvar,pairconstraints)
        else
            call gs%eval(map,depvar)
        endif
        call map%get_secondorder_forceconstant(p,fc,mem,-1)

        ! Get static energy (interpolated)
        call gs%energy%interpolate( depvar, ge%grid_ACT%U(ia,ic,it) )

        ! Get delta U0 from interpolation (if available)
        ge%grid_ACT%U0(ia,ic,it) = 0.0_flyt

        ! Generate q-mesh and dispersions
        call lo_generate_qmesh(qp,p,qgrid_harm,'fft',timereversal=.true.,headrankonly=.false.,mw=mw,mem=mem,verbosity=-1)
        call dr%generate(qp,fc,p,mw=mw,mem=mem,verbosity=-1)

        ! Check for unstable modes
        if ( dr%omega_min .lt. lo_freqtol ) then
            ge%grid_ACT%fph(ia,ic,it) = 123456789.0_flyt
            ge%grid_ACT%ah3(ia,ic,it) = 0.0_flyt
            ge%grid_ACT%ah4(ia,ic,it) = 0.0_flyt
        else
            ! Evaluate phonon free energy at this temperature
            ge%grid_ACT%fph(ia,ic,it) = dr%phonon_free_energy(temperature)
            
            ! Evaluate anharmonic contributions if available
            if ( have_anharmonic ) then
                call map%get_thirdorder_forceconstant(p,fct)
                call map%get_fourthorder_forceconstant(p,fcf)
                select type(qp); type is(lo_fft_mesh)
                    call anharmonic_free_energy(p,fct,fcf,qp,dr,temperature,f3,f4,mw,mem,verbosity=-1)
                end select
                ge%grid_ACT%ah3(ia,ic,it) = f3
                ge%grid_ACT%ah4(ia,ic,it) = f4
            else
                ge%grid_ACT%ah3(ia,ic,it) = 0.0_flyt
                ge%grid_ACT%ah4(ia,ic,it) = 0.0_flyt
            endif
        endif

        if ( mw%talk .and. ipt .lt. npts ) then
            call lo_progressbar(' ... a-c-T free energy',ipt,npts,walltime()-t0)
        endif
    enddo
    enddo
    enddo

    if ( mw%talk ) call lo_progressbar(' ... a-c-T free energy',npts,npts,walltime()-t0)

    ! Collect results across MPI ranks
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_ACT%U,  ge%grid_ACT%na*ge%grid_ACT%nc*ge%grid_ACT%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_ACT%U0, ge%grid_ACT%na*ge%grid_ACT%nc*ge%grid_ACT%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_ACT%fph,ge%grid_ACT%na*ge%grid_ACT%nc*ge%grid_ACT%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_ACT%ah3,ge%grid_ACT%na*ge%grid_ACT%nc*ge%grid_ACT%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(MPI_IN_PLACE,ge%grid_ACT%ah4,ge%grid_ACT%na*ge%grid_ACT%nc*ge%grid_ACT%nt,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)

    ! Compute total free energy (including anharmonic)
    ge%grid_ACT%Ftot = ge%grid_ACT%U + ge%grid_ACT%U0 + ge%grid_ACT%fph + ge%grid_ACT%ah3 + ge%grid_ACT%ah4
end subroutine

!> Finalize a-c-T grid output with polynomial fitting
subroutine ACTfinalize(gr,gs,map,filename,qgrid,mw,mem)
    !> grid
    type(lo_gridenergy_act), intent(inout) :: gr
    !> gridsim
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> filename
    character(len=*), intent(in) :: filename
    !> q-grid
    integer, dimension(3), intent(in) :: qgrid
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:,:), allocatable :: design_matrix
    real(r8), dimension(:), allocatable :: F_vec, coeffs
    real(r8) :: a_norm, c_norm, a_mean, c_mean, a_scale, c_scale
    real(r8) :: a_opt, c_opt, F_opt
    integer :: ia, ic, it, k, l, npts, ncoeffs, u
    integer :: ia_min, ic_min, pa, pc

    ncoeffs = 25  ! 5x5 for 4th order polynomial in 2D
    npts = gr%na * gr%nc

    if ( mw%talk ) then
        write(*,*) ''
        write(*,*) 'Fitting 4th order polynomial to free energy surface (a-c-T)'

        ! Normalize a and c for numerical stability
        a_mean = 0.5_r8*(minval(gr%a_values) + maxval(gr%a_values))
        c_mean = 0.5_r8*(minval(gr%c_values) + maxval(gr%c_values))
        a_scale = 0.5_r8*(maxval(gr%a_values) - minval(gr%a_values))
        c_scale = 0.5_r8*(maxval(gr%c_values) - minval(gr%c_values))
        if (a_scale .lt. 1.0e-10_r8) a_scale = 1.0_r8
        if (c_scale .lt. 1.0e-10_r8) c_scale = 1.0_r8

        allocate(design_matrix(npts, ncoeffs))
        allocate(F_vec(npts))
        allocate(coeffs(ncoeffs))

        ! Build design matrix once
        l = 0
        do ia = 1, gr%na
        do ic = 1, gr%nc
            l = l + 1
            a_norm = (gr%a_values(ia) - a_mean) / a_scale
            c_norm = (gr%c_values(ic) - c_mean) / c_scale
            k = 0
            do pa = 0, 4
            do pc = 0, 4 - pa
                k = k + 1
                design_matrix(l, k) = (a_norm**pa) * (c_norm**pc)
            enddo
            enddo
            do while (k < ncoeffs)
                k = k + 1
                design_matrix(l, k) = 0.0_r8
            enddo
        enddo
        enddo

        ! For each temperature, fit and find minimum
        do it = 1, gr%nt
            l = 0
            do ia = 1, gr%na
            do ic = 1, gr%nc
                l = l + 1
                F_vec(l) = gr%Ftot(ia, ic, it)
            enddo
            enddo

            ! Check for unstable points
            if (any(F_vec .gt. 1.0e8_r8)) then
                gr%a_min(it) = -1.0_r8
                gr%c_min(it) = -1.0_r8
                gr%F_min(it) = 1.0e10_r8
                gr%poly_coeffs(:, it) = 0.0_r8
                cycle
            endif

            ! Least squares fit
            call lo_linear_least_squares(design_matrix, F_vec, coeffs)
            gr%poly_coeffs(:, it) = coeffs

            ! Find minimum by grid search
            F_opt = 1.0e10_r8
            a_opt = gr%a_values(1)
            c_opt = gr%c_values(1)
            do ia = 1, gr%na
            do ic = 1, gr%nc
                if (gr%Ftot(ia, ic, it) .lt. F_opt) then
                    F_opt = gr%Ftot(ia, ic, it)
                    a_opt = gr%a_values(ia)
                    c_opt = gr%c_values(ic)
                endif
            enddo
            enddo
            gr%a_min(it) = a_opt
            gr%c_min(it) = c_opt
            gr%F_min(it) = F_opt
        enddo

        deallocate(design_matrix, F_vec, coeffs)

        ! Write output to HDF5
        write(*,*) 'Writing a-c-T results to ',trim(filename)

        call h5%init(__FILE__,__LINE__)
        call h5%open_file('write',trim(filename))
        call h5%open_group('write','grid_act')

        ! Store axes
        call h5%store_data(gr%a_values, h5%group_id,'a_values',enhet='Angstrom')
        call h5%store_data(gr%c_values, h5%group_id,'c_values',enhet='Angstrom')
        call h5%store_data(gr%temperature, h5%group_id,'temperatures',enhet='K')

        ! Store energies
        call h5%store_data(gr%U*lo_Hartree_to_eV,   h5%group_id,'static_internal_energy',enhet='eV/atom')
        call h5%store_data(gr%U0*lo_Hartree_to_eV,  h5%group_id,'delta_U0',enhet='eV/atom')
        call h5%store_data(gr%fph*lo_Hartree_to_eV, h5%group_id,'phonon_free_energy',enhet='eV/atom')
        call h5%store_data(gr%ah3*lo_Hartree_to_eV, h5%group_id,'anharmonic_free_energy_3rd',enhet='eV/atom')
        call h5%store_data(gr%ah4*lo_Hartree_to_eV, h5%group_id,'anharmonic_free_energy_4th',enhet='eV/atom')
        call h5%store_data(gr%Ftot*lo_Hartree_to_eV, h5%group_id,'Helmholtz_free_energy',enhet='eV/atom')

        ! Store minimum finding results
        call h5%store_data(gr%a_min, h5%group_id,'a_equilibrium',enhet='Angstrom')
        call h5%store_data(gr%c_min, h5%group_id,'c_equilibrium',enhet='Angstrom')
        call h5%store_data(gr%F_min*lo_Hartree_to_eV, h5%group_id,'F_equilibrium',enhet='eV/atom')

        ! Store polynomial coefficients
        call h5%store_data(gr%poly_coeffs, h5%group_id,'polynomial_coefficients',enhet='')

        call h5%close_group()
        call h5%close_file()

        ! Write text file with equilibrium a(T), c(T)
        u = open_file('out','outfile.act_equilibrium.dat')
        write(u,'(A)') '# Temperature (K)    a_eq (Angstrom)    c_eq (Angstrom)    F_min (eV/atom)'
        do it = 1, gr%nt
            if (gr%a_min(it) .gt. 0.0_r8) then
                write(u,'(4E20.10)') gr%temperature(it), gr%a_min(it), gr%c_min(it), gr%F_min(it)*lo_Hartree_to_eV
            endif
        enddo
        close(u)

        write(*,*) '... done writing a-c-T output'
    endif
end subroutine

!> anharmonic free energy at a single point
subroutine anharmonic_free_energy_for_single_point( gs,map,depvar,qgrid,temperature,free_energy_thirdorder,free_energy_fourthorder,mw,mem )
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> where to evaluate
    real(flyt), dimension(:), intent(in) :: depvar
    !> q-mesh density
    integer, dimension(3), intent(in) :: qgrid
    !> temperature
    real(flyt), intent(in) :: temperature
    !> free energy
    real(flyt), intent(out) :: free_energy_thirdorder,free_energy_fourthorder
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    class(lo_qpoint_mesh), allocatable :: qp
    type(lo_phonon_dispersions) :: dr
    type(lo_crystalstructure) :: p
    type(lo_forceconstant_secondorder) :: fc
    type(lo_forceconstant_thirdorder) :: fct
    type(lo_forceconstant_fourthorder) :: fcf
    real(flyt), dimension(:,:), allocatable :: pairconstraints
    integer :: nconstr

    free_energy_thirdorder  = 0.0_flyt
    free_energy_fourthorder = 0.0_flyt
    ! setup stuffs, first get a structure
    call gs%structure%interpolate(depvar,p,gs%info%dim_volume)
    call p%classify('wedge',timereversal=.true.)
    ! forceconstants, first the constraints
    call megafit_secondorder_constraints( map,p,pairconstraints,nconstr,.true.,.true.,.true. )
    if ( nconstr .gt. 0 ) then
        call gs%eval(map,depvar,pairconstraints)
    else
        call gs%eval(map,depvar)
    endif
    call map%get_secondorder_forceconstant(p,fc,mem,-1)
    call map%get_thirdorder_forceconstant(p,fct)
    call map%get_fourthorder_forceconstant(p,fcf)
    ! q-points
    call lo_generate_qmesh(qp,p,qgrid,'fft',timereversal=.true.,headrankonly=.false.,mw=mw,mem=mem,verbosity=-1)
    ! dispersions
    call dr%generate(qp,fc,p,mw=mw,mem=mem,verbosity=-1)
    ! check for unstable modes right away
    if ( dr%omega_min .lt. lo_freqtol ) then
        ! kinda dumb, but I have no better idea right now.
        free_energy_thirdorder  = 123456789.0_flyt
        free_energy_fourthorder = 123456789.0_flyt
        return
    endif
    select type(qp); type is(lo_fft_mesh)
        call anharmonic_free_energy(p,fct,fcf,qp,dr,temperature,free_energy_thirdorder,free_energy_fourthorder,mw,mem,verbosity=-1)
    end select
end subroutine

!> phonon free energy at a single point
subroutine phonon_free_energy_for_single_point(gs,map,depvar,qgrid,temperature,mw,mem,fph)
    !> simulation grid
    type(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> where to evaluate
    real(flyt), dimension(:), intent(in) :: depvar
    !> q-mesh density
    integer, dimension(3), intent(in) :: qgrid
    !> temperature
    real(flyt), intent(in) :: temperature
    !> mpi communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> free energy
    real(flyt), intent(out) :: fph


    class(lo_qpoint_mesh), allocatable :: qp
    type(lo_phonon_dispersions) :: dr
    type(lo_crystalstructure) :: p
    type(lo_forceconstant_secondorder) :: fc
    real(flyt), dimension(:,:), allocatable :: pairconstraints
    integer :: nconstr

    ! setup stuffs, first get a structure
    call gs%structure%interpolate(depvar,p,gs%info%dim_volume)
    call p%classify('wedge',timereversal=.true.)
    ! constraints
    call megafit_secondorder_constraints(map,p,pairconstraints,nconstr,.true.,.true.,.true.)
    ! evaluate forceconstants
    if ( nconstr .gt. 0 ) then
        call gs%eval(map,depvar,pairconstraints)
    else
        call gs%eval(map,depvar)
    endif
    ! forceconstant at this point
    call map%get_secondorder_forceconstant(p,fc,mem,-1)
    ! q-points
    call lo_generate_qmesh(qp,p,qgrid,'monkhorst',timereversal=.true.,headrankonly=.false.,mw=mw,mem=mem,verbosity=-1)
    ! phonon dispersions
    call dr%generate(qp,fc,p,mw=mw,mem=mem,verbosity=-1)
    ! and the free energy
    fph=dr%phonon_free_energy(temperature)
end subroutine

end module
