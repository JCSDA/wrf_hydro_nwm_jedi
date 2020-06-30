! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wrf_hydro_nwm_jedi_covariance_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration

!use Shallow_Water_kind,     only : r8kind
use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_state_mod
use wrf_hydro_nwm_jedi_state_utils_mod, only: wrf_hydro_nwm_jedi_state

implicit none

!> Fortran derived type to hold configuration data for the background/model covariance
type :: wrf_hydro_nwm_jedi_covar
  ! Climatological variance of U, V, H
  real :: normfactor
  integer      :: rad      ! Radius, in grid points, of covariance cutoff

  ! Gaussian decay mask for covariance
  real, allocatable :: cov_mask(:,:)
end type wrf_hydro_nwm_jedi_covar

! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------

!> Setup for the model's 3d error covariance matrices (B and Q_i)

!> This routine queries the configuration for the parameters that define the
!! covariance matrix, and stores the relevant values in the
!! error covariance structure.

subroutine wrf_hydro_nwm_jedi_covar_setup(self, geom, c_conf)

  implicit none
  type(wrf_hydro_nwm_jedi_covar),  intent(inout) :: self    !< Covariance structure
  type(c_ptr),  intent(   in) :: c_conf  !< Configuration
  type(wrf_hydro_nwm_jedi_geometry), intent(   in) :: geom    !< Geometry

  type(fckit_configuration)  :: f_conf
  integer                    :: i, j
  real              :: distnorm

  f_conf = fckit_configuration(c_conf)

  ! Get field normalization factors from configuration
  call f_conf%get_or_die("normfactor", self%normfactor)

  ! Get the grid point radius for covariance cutoff
  call f_conf%get_or_die("rad", self%rad)

  ! Precompute covariance localization mask based on Gaspari-Cohn 99 function
  allocate(self%cov_mask(0:self%rad, 0:self%rad))
  do j = 0, self%rad
    do i = 0, self%rad
      distnorm = sqrt(real(i**2 + j**2)) / self%rad
      if (distnorm < 0.5) then
        self%cov_mask(i,j) = -8.0 * distnorm**5.0 + 8.0 * distnorm**4.0        &
                           +  5.0 * distnorm**3.0 - 20.0 / 3.0 * distnorm**2.0 &
                           +  1.0
      else if (distnorm < 1.0) then
        self%cov_mask(i,j) = 8.0 / 3.0 * distnorm**5.0 - 8.0 * distnorm**4.0   &
                           + 5.0 * distnorm**3.0 + 20.0 / 3.0 * distnorm**2.0  &
                           - 10.0 * distnorm + 4.0 - 1.0 / (3.0 * distnorm)
      else
        self%cov_mask(i,j) = TINY(0.0)
      end if
    end do
  end do

end subroutine wrf_hydro_nwm_jedi_covar_setup

! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_covar_delete(self)

  implicit none
  type(wrf_hydro_nwm_jedi_covar), intent(inout) :: self  !< Covariance structure

end subroutine wrf_hydro_nwm_jedi_covar_delete

! ------------------------------------------------------------------------------

!> Multiply increment by C, where C is a covariance matrix

subroutine wrf_hydro_nwm_jedi_covar_mult(self, xin, xout)

  use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

  implicit none
  type(wrf_hydro_nwm_jedi_covar), intent(   in) :: self
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: xin
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: xout
  ! type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  
  real, pointer             :: in_u(:,:), in_v(:,:), in_h(:,:)
  real, pointer             :: out_u(:,:), out_v(:,:), out_h(:,:)
  integer                           :: xps, xpe, yps, ype
  real                      :: local_out(3), global_out(3)
  integer                           :: nneighbors, nown
  integer                           :: nx, ny
  integer                           :: i, j, ki, kj
  type(fckit_mpi_comm)              :: f_comm

  ! Get geometry bounds
 ! geom = xout%get_geometry()
  ! xps = geom%get_xps()
  ! xpe = geom%get_xpe()
  ! yps = geom%get_yps()
  ! ype = geom%get_ype()
  ! nx = geom%get_nx()
  ! ny = geom%get_ny()
  nx = 1
  ny = 1

  ! Get pointers to xout fields
  call xout%get_u_ptr(out_u)
  ! call xout%get_v_ptr(out_v)
  ! call xout%get_h_ptr(out_h)

  ! Get pointers to xin fields
  call xin%get_u_ptr(in_u)

  ! Compute each element of xout
  do j = 1, ny
    do i = 1, nx

      ! Initialize neighbor counters - to be used later to minimize communication
      nneighbors = (min(nx, i + self%rad) - max(1, i - self%rad) + 1) * (min(ny, j + self%rad) - max(1, j - self%rad) + 1)
      nown = 0

      ! for this element of xout, compute our local portion of the dot product
      local_out = 0.0
      do kj = max(1, j - self%rad), min(ny, j + self%rad)
        if (kj >= yps .AND. kj <= ype) then
          do ki = max(1, i - self%rad), min(nx, i + self%rad)
            if (ki >= xps .AND. ki <= xpe) then
              nown = nown + 1
              local_out(1) = local_out(1) + self%unormfactor * self%cov_mask(abs(i - ki), abs(j - kj)) * in_u(ki, kj)
              local_out(2) = local_out(2) + self%vnormfactor * self%cov_mask(abs(i - ki), abs(j - kj)) * in_v(ki, kj)
              local_out(3) = local_out(3) + self%hnormfactor * self%cov_mask(abs(i - ki), abs(j - kj)) * in_h(ki, kj)
            end if
          end do
        end if
      end do

      ! ! Add up all the partial sums
      ! call f_comm%allreduce(local_out, global_out, fckit_mpi_sum())

      ! ! If we own this element, update it
      ! if (i >= xps .AND. i <= xpe .AND. j >= yps .AND. j <= ype) then
      !   out_u(i,j) = global_out(1)
      !   out_v(i,j) = global_out(2)
      !   out_h(i,j) = global_out(3)
      ! end if
    end do
  end do

end subroutine wrf_hydro_nwm_jedi_covar_mult

! ------------------------------------------------------------------------------

end module wrf_hydro_nwm_jedi_covariance_mod
