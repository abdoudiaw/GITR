
module share

    implicit none
    public
  
    integer, parameter :: Max_Num=2000000
! define pi
real, parameter :: pi=3.14159265359
real, parameter :: m_p = 1.67262158e-27
real, parameter :: e = 1.60217646e-19
real, parameter :: eps0 = 8.854187817e-12
real, parameter :: k = 1.3806503e-23
    real(8), dimension(Max_Num) :: x, y, z, xu, yu, zu, vx, vy, vz, ax, ay, az, ex, ey, ez, bx, by, bz,  vx0, vy0, vz0
    real(8), dimension(Max_Num) :: local_x, local_y, local_z, local_vx, local_vy, local_vz
    real(8) :: rsincm, dlinrs, edgel, atoe0, tfact, vfact, zbar3, fermit, gamma
    real(8) :: densi, zbar, pmass, tempi, tempe, dense, dt, tempeV, tempea
    real(8) :: tener, vener, Ewald_Factor
    real(8) :: t(3), s(3), v_minus(3), v_prime(3), v_plus(3), cross_v_minus_t(3), xdir(3), ydir(3), zdir(3)

    real(8) :: v_minus_x, v_minus_y, v_minus_z

    real(8) :: sigma, vxsum, vysum, vzsum, r1, r2, xx, yy, zz, distsq
    real(8) :: sheath_factor, q_over_m, Bmagnitude, alpha, larmor_over_debye_radius
    real(8) :: dens_back, zbar_back, pmass_back, tempi_back, tempe_back
    real(8) :: pmass_imp, tempi_imp, zbar_imp, Debye, larmor_radius, vti, omega_ci
    integer :: npart, steps, snaps, nequi, which, initvx, save_run, second_ion
    integer :: loop, i, j
real (8) :: xmax, xmin, ymax, ymin, zmax, zmin
integer :: npartActive
logical, allocatable :: isActive(:)
  
  end module share

  
  
