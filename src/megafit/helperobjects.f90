module helperobjects
!! small helper objects that made the other files too messy.
use konstanter, only: r8,lo_hugeint,lo_sqtol
use gottochblandat, only: lo_chop
use type_blas_lapack_wrappers, only: lo_dgesvd
use type_crystalstructure, only: lo_crystalstructure
use type_forcemap, only: lo_forcemap, lo_secondorder_rot_herm_huang
implicit none

private
public :: lo_sparsematrix
public :: reduce_equations
public :: megafit_secondorder_constraints

!> special case of sparse matrix
type lo_sparsematrix
    !> total number of elements
    integer :: n=-lo_hugeint
    !> dimensions
    integer :: nrow=-lo_hugeint
    integer :: ncol=-lo_hugeint
    !> row index
    integer, dimension(:), allocatable :: rowind
    !> column index
    integer, dimension(:), allocatable :: colind
    !> value
    real(r8), dimension(:), allocatable :: val
    contains
        procedure :: init=>init_sparsething
        procedure :: destroy=>destroy_sparsething
end type

contains

!> initialize my sparse guy. Highly specialized.
subroutine init_sparsething(A,nrow,ncol)
    !> specialized sparse guy
    class(lo_sparsematrix), intent(inout) :: A
    !> rows, columns
    integer, intent(in) :: ncol,nrow

    A%nrow=nrow
    A%ncol=ncol
    A%n=ncol
    allocate(A%rowind(A%n))
    allocate(A%colind(A%n))
    allocate(A%val(A%n))
    A%rowind=0.0_r8
    A%colind=0.0_r8
    A%val=0.0_r8
end subroutine

!> destroy the sparse guy
subroutine destroy_sparsething(A)
    !> specialized sparse guy
    class(lo_sparsematrix), intent(inout) :: A

    if ( allocated(A%rowind) ) deallocate(A%rowind)
    if ( allocated(A%colind) ) deallocate(A%colind)
    if ( allocated(A%val)    ) deallocate(A%val)
end subroutine

!> take a bunch of equations and SVD them, to get the irreducible amount
subroutine reduce_equations(allequations,redeq,nredeq)
    !> all equations
    real(r8), dimension(:,:), intent(in) :: allequations
    !> the reduced set
    real(r8), dimension(:,:), allocatable, intent(out) :: redeq
    !> how many in the reduced set
    integer, intent(out) :: nredeq

    real(r8), dimension(:,:), allocatable :: m,u,v
    real(r8), dimension(:), allocatable :: s
    integer :: i,l,nu,ne,ns
    !
    nu=size(allequations,2)
    ne=size(allequations,1)
    ns=min(nu,ne)
    allocate(m(nu,ne))
    m=transpose(allequations)
    ! remove tiny tiny numbers
    m=lo_chop(m,lo_sqtol)
    ! SVD this
    allocate(s(ns))
    allocate(u(nu,nu))
    allocate(v(ne,ne))
    call lo_dgesvd(m,s,u,v)
    l=0
    do i=1,ns
        if ( s(i) .gt. lo_sqtol ) then
            l=l+1
        endif
    enddo
    nredeq=l
    if ( l .gt. 0 ) then
        ! return the equations with non-zero singular vectors
        allocate(redeq(nu,l))
        redeq=u(:,1:l)
        ! remove tiny numbers
        redeq=lo_chop(redeq,lo_sqtol)
    endif
    ! cleanup
    deallocate(m,s,u,v)
end subroutine

!> Wrapper for lo_secondorder_rot_herm_huang with homogeneous RHS (zeros)
!> This provides the old API for megafit code compatibility
subroutine megafit_secondorder_constraints(map, uc, eq2, neq2, rotational, huang, hermitian)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> constraints matrix (output)
    real(r8), dimension(:,:), allocatable, intent(out) :: eq2
    !> number of equations (output)
    integer, intent(out) :: neq2
    !> which constraints to apply
    logical, intent(in) :: rotational, huang, hermitian

    real(r8), dimension(:), allocatable :: vD
    real(r8), dimension(:), allocatable :: hermitian_rhs, huang_rhs, rotational_rhs

    ! Allocate zero RHS arrays with proper sizes
    allocate(hermitian_rhs(uc%na * 9))
    allocate(huang_rhs(81))
    allocate(rotational_rhs(uc%na * 27))
    hermitian_rhs = 0.0_r8
    huang_rhs = 0.0_r8
    rotational_rhs = 0.0_r8

    ! Call the actual subroutine with zero RHS
    call lo_secondorder_rot_herm_huang(map, uc, eq2, vD, neq2, &
                                       rotational, huang, hermitian, &
                                       hermitian_rhs, huang_rhs, rotational_rhs)

    ! Cleanup
    deallocate(hermitian_rhs, huang_rhs, rotational_rhs)
    if (allocated(vD)) deallocate(vD)
end subroutine

end module
