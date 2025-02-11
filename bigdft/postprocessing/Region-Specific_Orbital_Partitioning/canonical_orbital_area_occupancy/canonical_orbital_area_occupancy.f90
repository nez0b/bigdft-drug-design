!> @file
!! @author
!!    Copyright (C) 2007-2013 BigDFT group. This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt.
!!    For the list of contributors, see ~/AUTHORS

!===============================!
!> Orbitals localization with given range.
!===============================!

program toy_model
  ! use fyaml
  use box
  use at_domain
  use Poisson_Solver!, except_dp => dp, except_gp => gp
  use BigDFT_API
  ! -----newly added:
  use io
  use mpif_module
  use module_bigdft_mpi
  !use mpi
  !use wrapper_MPI
  !--------------------
  use bigdft_run
  use dynamic_memory
  use compression
  use yaml_output
  use module_input_dicts
  use module_input_keys
  !-------------------
  !use module_razero
  use f_zero_module
  !-----------------
  use module_atoms
  use module_dpbox,       only: denspot_distribution,dpbox_free,dpbox_set
  use rhopotential,       only: full_local_potential
  use locregs_init,       only: lr_set
  use locreg_operations 
  use communications_base,  only: deallocate_comms
  use communications_init,  only: orbitals_communicators
  use communications,       only: transpose_v, untranspose_v
  !-----------------
  !use module_base
  use PSbase
  !---------------
  use module_types
  use psp_projectors_base
  use psp_projectors
  use pseudopotentials
  use orbitalbasis
  use ao_inguess, only: lmax_ao
  use dictionaries
  use at_domain, domain_renamed => domain !, only: domain, domain_new
  use numerics, only: onehalf,pi
  use f_refcnts
  
  implicit none
  type(input_variables)        :: inputs
  type(atoms_data)             :: atoms
  type(local_zone_descriptors) :: Lzd
  type(orbitals_data)          :: orbs, orbsv
  type(comms_cubic)            :: comms
  type(workarr_sumrho)         :: wisf
  type(xc_info)                :: xc
  type(rho_descriptors)        :: rhodsc
  type(denspot_distribution)   :: dpcom
  type(GPU_pointers)           :: GPU
  type(workarr_locham)         :: wrk_lh
  type(coulomb_operator)       :: pkernel
  type(paw_objects)            :: paw
  type(DFT_PSP_projectors)     :: nlpsp
  type(ket)                    :: psi_it, psi_itv
  type(orbital_basis)          :: psi_ob, psi_obv
  type(DFT_PSP_projector_iter) :: psp_it
  type(atomic_proj_matrix)     :: prj
  type(dictionary), pointer :: user_inputs, options, dict
  type(domain_renamed) :: dom
  type(f_reference_counter) :: pot_ref_count

  logical  :: dosome, alive, init_projectors_completely
! logical  :: paw = .false.
! integer, parameter :: dp=8
! integer, parameter :: gp=8
  integer, parameter :: wp=8
  real(dp) :: nrm, epot_sum
  real(dp) :: psoffset, ekin, ekin_sum
  real(gp), dimension(3)            :: shift
  real(wp) :: place_holder
  real(wp), dimension(:),   pointer :: psi, psiv, psir, psir_i,psir_j, f
  real(dp), dimension(:),   pointer :: rhor,rho_ion
  real(wp), dimension(:),   allocatable :: pot_ion
  real(wp), dimension(:),   pointer :: potential
  real(wp), dimension(:)  , pointer :: hpsi_ptr
  real(gp), dimension(:,:), pointer :: rxyz,wann_func
  real(gp), dimension(:,:), pointer :: rxyz_old
  real(wp), dimension(:,:), pointer :: ovrlp, grid_xyz
  real(dp), dimension(:,:), pointer :: rho_p => null() !needs to be nullified
  real(wp), dimension(:,:), allocatable :: pot_tmp

  ! Interpolating Scaling Function (daub_to_ISF)
  real(wp), dimension(:,:), allocatable :: psirr
  ! Spin Polarized Format of psirr
  real(wp), dimension(:,:), allocatable :: phirr, phirr_v
  ! phirr_overlap = phirr(:,i)*phirr(:,j)
  real(wp), dimension(:), allocatable :: phirr_overlap, grid_r, u_matrix
  ! Spin Polarized Format of Projector
  real(wp), dimension(:,:), allocatable :: proj, hproj !Projector Array
  integer :: nnproj, iproj, iat          !Projector Index, Atomic Index
  integer :: num_proj     !Number of Different Species Atomic Projector
  integer :: nmproj       !Number of Projector Considered for each Atom

  real(wp), dimension(:)  , allocatable :: tpsi, tpsi_o, tpsi_v, hpsi
  real(wp), dimension(:)  , allocatable :: tho
  real(dp), dimension(:,:), allocatable :: E_local, E_nonlocal, E_kin
  real(dp), dimension(:,:), allocatable :: output1, output2, output3
  real(dp)                              :: eproj, ektmp1, ektmp2, Enl

  integer(kind=8), dimension(:), allocatable :: indwvl_i,indwvl_f
  integer :: ierr, iproc, nproc, nwarnings, ii,jj,nn, n_virt, num_args, n_occ
  integer :: iorb, nvirtu, nvirtd, ispsi, ilr, ilr_orb, npot, mspin, cspin
  integer :: ityp, nc, ishift, nwarning, dry_run, orbtot, nv, istate, n_u_mat, n_plot
  integer :: orbdimocc, orbdimvir, norb_start,norb_end,  orbocc,orbvir
  integer(kind=4)  :: OMP_get_max_threads, OMP_get_thread_num, OMP_get_num_threads
  integer :: max_ind, x, y, z,target_index

  character(len=100) :: cmd, name, arg
  character*32 con, tmp
  integer :: ip ,iq ,ir ,is , ihpq, ihpqrs, istat, i,j,k, re_index1, re_index2,loop_b
  integer :: ipt,iqt,irt,ist, nhpq, nhpqrs
  real(kind=8) :: hpq,hpqrs,hpq_overlap,hpq_abs_overlap
  logical, dimension(3) :: peri
  integer, dimension(3) :: n_grids_for_cube_file
  real(wp), dimension(:,:), allocatable :: target_grid_positions

  ! type(yaml) :: config
  real(wp), allocatable :: aij(:)
  integer :: n, lda, lwork, info, lwork_v
  real(wp), allocatable :: a(:,:), w_(:), work(:), a_virt(:,:), w_v(:), work_v(:)
  
  character(len=10) :: unit
  integer, allocatable :: indices(:)
  real(8), allocatable :: radii(:)
  integer :: n_items
  character(len=256) :: line
  logical :: found_unit
  integer :: idx
  real(8) :: rad
  integer :: n_lines

  ! PSolver
  ! hpqrs 
  real(dp), dimension(:), allocatable :: rhopq, rhors

!-------------------------------------------------------------------------------------------------------------------------------------
  
  open(124,file="toy_model.log") 
  open(  6,file="toy_model.out")

  !-----------------------------------------------
  ! initializes the mpi_environment and bigdft
  !-----------------------------------------------
  call f_lib_initialize() ; nullify(options)
  call bigdft_init(options) ; call dict_free(options)

  iproc = bigdft_mpi%iproc
  nproc = bigdft_mpi%nproc
  call dict_init(user_inputs)
  call user_dict_from_files(user_inputs, 'input', 'posinp', bigdft_mpi)
  call inputs_from_dict(inputs, atoms, user_inputs)

  if (iproc == 0) call print_general_parameters(inputs,atoms,'input')
  call dict_free(user_inputs)
  GPU%OCLconv = .false.

  call system_properties(iproc,nproc,inputs,atoms,orbs)

  Lzd = default_Lzd() ; Lzd%hgrids=(/ inputs%hx, inputs%hy, inputs%hz /)
  call lr_set(Lzd%Glr, iproc, GPU%OCLconv, .true., inputs%crmult, inputs%frmult, &
              Lzd%hgrids,atoms%astruct%rxyz,atoms,.true.,.false.)
  call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,comms)
  call xc_init(xc, inputs%ixc, XC_ABINIT, inputs%nspin)
  call dpbox_set(dpcom,Lzd%Glr%mesh,xc,iproc,nproc,MPI_COMM_WORLD, inputs%SIC%approach, inputs%nspin)

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

  if(orbs%nspin .eq. 1) cspin = 2 ; if(orbs%nspin .eq. 1) mspin = 1
  if(orbs%nspin .eq. 2) cspin = 1 ; if(orbs%nspin .eq. 2) mspin = 2
  orbocc = orbs%norb
  orbvir = mspin*inputs%norbv
  orbtot = orbs%norb + mspin*inputs%norbv

  allocate(orbs%eval(orbs%norb*orbs%nkpts)) ; call f_zero(orbs%eval)
  allocate(psi(  max(orbs%npsidim_orbs, orbs%npsidim_comp))) ; psi=0._wp
  allocate(psiv( 2*inputs%norbv * (Lzd%Glr%wfd%nvctr_c + 7*Lzd%Glr%wfd%nvctr_f))) ; psiv=0._wp
  allocate(rxyz_old(3, atoms%astruct%nat)) ; rxyz_old=0._gp
  

  !------------------i--------------------------------------------------!
  ! Read occupied state wavefunctions from disk and store them in psi. !
  ! orbs%norb - number of occupied orbitals                            !
  !--------------------------------------------------------------------!
  ! orbs%norbp = 20

  call check_linear_and_create_Lzd(iproc,nproc,inputs%linear,Lzd,atoms,orbs,inputs%nspin,atoms%astruct%rxyz)
  
  call readmywaves(iproc,trim(inputs%dir_output) // "wavefunction", &
  WF_FORMAT_PLAIN,orbs,lzd%glr,atoms,rxyz_old,atoms%astruct%rxyz,psi)

  if(nproc>1) call fmpi_allreduce(orbs%eval(1), orbs%norb*orbs%nkpts, op=FMPI_SUM)

  !--------------------------------------------------------------------!
  ! Read virtual  state wavefunctions from disk and store them in psi. !
  ! inputs%norbv - number of virtual orbitals                          !
  ! orbsv%norb - number of total virtual orbitals (2*inputs%norbv)     !
  !--------------------------------------------------------------------!

  nvirtu = abs(inputs%norbv) ; nvirtd = nvirtu
  call orbitals_descriptors(iproc, nproc, nvirtu+nvirtd, nvirtu, nvirtd, &
                            orbs%nspin, orbs%nspinor, orbs%nkpts, orbs%kpts, orbs%kwgts, orbsv, LINEAR_PARTITION_NONE)


  nullify(orbsv%eval)
  orbsv%eval = f_malloc_ptr(orbsv%norb*orbsv%nkpts,id='orbsv%eval')
  call check_linear_and_create_Lzd(iproc,nproc,inputs%linear,Lzd,atoms,orbsv,inputs%nspin,atoms%astruct%rxyz)
  
  call readmywaves(iproc, "data/virtuals",WF_FORMAT_PLAIN, orbsv, Lzd%Glr, atoms, rxyz_old, atoms%astruct%rxyz,psiv)
  if(nproc>1) call fmpi_allreduce(orbsv%eval(1), orbsv%norb*orbsv%nkpts, op=FMPI_SUM)

  orbdimocc = (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%norb
  orbdimvir = (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbsv%norb

  allocate(psir_i( Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3)));psir_i=0._wp
  allocate(f( Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3)));f=0._wp
  allocate(phirr(Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3) ,0:orbs%norb*cspin-1)) ; phirr=0._wp  
  allocate(phirr_v(Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3) ,0:orbsv%norb-1)) ; phirr_v=0._wp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
 ! call local_potential_dimensions(Lzd, orbs, dpcom%ngatherarr(0,1))
  dict => dict_new()

  dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(domain_geocode(atoms%astruct%dom)),&
            alpha_bc=onehalf*pi,beta_ac=onehalf*pi,gamma_ab=onehalf*pi,&
            acell=(/Lzd%Glr%nboxi(2,1), Lzd%Glr%nboxi(2,2), Lzd%Glr%nboxi(2,3)/)*&
                  (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/))

  pkernel=pkernel_init(iproc,nproc,dict,&
      dom,(/Lzd%Glr%nboxi(2,1), Lzd%Glr%nboxi(2,2), Lzd%Glr%nboxi(2,3)/),&
       (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/))

  call dict_free(dict)
  call pkernel_set(pkernel,verbose=.false.)
  
  allocate(pot_ion(Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*dpcom%n3p)) ; pot_ion=0._gp
  nullify(rho_ion) ; allocate(rho_ion(Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*dpcom%n3p)) ; rho_ion=0._gp
  write(124,*) " call createIonicPotential"
  
  call createIonicPotential(iproc, (iproc==0),atoms, atoms%astruct%rxyz,&
   inputs%elecfield, dpcom, pkernel, pot_ion, rho_ion, psoffset)
  
  write(124,*) " call full_local_potential"
 ! allocate the potential in the full box
  call full_local_potential(iproc,nproc,orbs,Lzd,dpcom,xc,pot_ion,potential,pot_ref_count)

  write(124,*) " call initialize_work_arrays_sumrho"

  call initialize_work_arrays_sumrho(Lzd%Glr , .true., wisf)

  ! Open the input file
  open(unit=10, file="localization.input", status="old", action="read")

  ! Initialize
  n_items = 0
  found_unit = .false.

  ! First Pass: Count index_and_radii entries
  do
    read(10, '(A)', iostat=i) line
    if (i /= 0) exit  ! End of file

    if (index(to_lowercase(trim(line)), "unit:") > 0 .and. .not. found_unit) then
      ! Parse the unit line
      unit = trim(adjustl(line(index(line, ":") + 1:)))  ! Extract the unit value
      found_unit = .true.
    elseif (index(to_lowercase(trim(line)), "index and radii:") > 0) then
      ! Count entries for index_and_radii
      do
        read(10, '(A)', iostat=i) line
        if (i /= 0 .or. trim(line) == "") exit  ! Stop if EOF or empty line
        n_items = n_items + 1
      end do
    end if
  end do

  rewind(10)  ! Reset file pointer

  n_virt = orbsv%norb/2
  n_occ = orbs%norb

  ! Allocate arrays
  allocate(indices(n_items))
  allocate(radii(n_items)) ; radii=0._wp
  allocate(aij(n_occ*n_occ)) ; aij=0._wp

  ! Second Pass: Read indices and radii
  n_lines = 0
  rewind(10)  ! Reset the file pointer to the beginning
  do
    read(10, '(A)', iostat=i) line
    if (i /= 0) exit  ! End of file
  
    if (index(to_lowercase(trim(line)), "index and radii:") > 0) then
      do
        read(10, '(A)', iostat=i) line
        if (i /= 0 .or. trim(line) == "") exit  ! Stop if EOF or empty line
        ! if (index(to_lowercase(trim(line)), "end") > 0) exit

  
        ! Parse index and radius
        read(line, *) idx, rad
        n_lines = n_lines + 1
        indices(n_lines) = idx
        radii(n_lines) = rad

        ! Convert radii to Bohr if needed
        if (trim(to_lowercase(unit)) == "bohr") then
          radii(n_lines) = radii(n_lines) / Lzd%Glr%mesh%hgrids(1)  ! Assuming hgrids(1)==hgrids(2)==hgrids(3)
        else
          radii(n_lines) = radii(n_lines) * 1.8897259886 / Lzd%Glr%mesh%hgrids(1)  ! Convert Angstrom to Bohr
        end if
      end do
    end if
  end do
  

  close(10)

  allocate(target_grid_positions(3, n_items)) ; target_grid_positions=0._wp




  peri = bc_periodic_dims(geocode_to_bc(domain_geocode(Lzd%Glr%mesh%dom)))

  do i=1,3
    if (peri(i)) then
      n_grids_for_cube_file(i) = Lzd%Glr%mesh%ndims(i)
    else
      n_grids_for_cube_file(i) = Lzd%Glr%mesh%ndims(i) - 31
    end if
  end do

  do target_index=1,n_items
    do i=1,3
      target_grid_positions(i, target_index)=14.0+atoms%astruct%rxyz(i,indices(target_index))/Lzd%Glr%mesh%hgrids(i)
    end do
  end do

  call create_function_3D(Lzd%Glr%mesh, target_grid_positions, radii, f)

  do i=1,orbs%norb
    if(orbs%nspin .eq. 1) then
      call daub_to_isf(Lzd%Glr,wisf,psi((i-1)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1) ,psir_i)
      phirr(:,2*(i-1)+0)=psir_i
      phirr(:,2*(i-1)+1)=psir_i
    else if(orbs%nspin .eq. 2) then
      call daub_to_isf(Lzd%Glr,wisf,psi((i-1)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1) ,psir_i)
      ! Block Spin Format
      !phirr(:,i-1)=psir_i
      ! Interleaved Spin Format
      ii=re_index1(i-1,orbs%norb*cspin+orbsv%norb,orbs%norb*cspin)
      phirr(:,ii)=psir_i
    end if
  end do
  
  ! Virtual Orbital (Spin Orbital Notation in BigDFT)
  do i=1,orbsv%norb
    call daub_to_isf(Lzd%Glr,wisf,psiv((i-1)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1) ,psir_i)
    ! Block Spin Format
    !phirr(:,orbs%norb*cspin+i-1)=psir_i
    ! Interleaved Spin Format
    ii=re_index1(i-1,orbsv%norb,0)
    phirr_v(:,ii)=psir_i
  end do

!-------------------------------------------------------------------------------------------------------------------------------------
  open(14121, file="area_occupancy_occ.out")
  open(141210, file="area_occupancy_virt.out")
  do j=1,n_occ
    write(14121, *) sum(phirr(:,2*(j-1))*f(:)*phirr(:,2*(j-1)))
  end do

  do j=1,n_virt
    write(141210, *) sum(phirr_v(:,2*(j-1))*f(:)*phirr_v(:,2*(j-1)))
  end do
  close(14121)
  close(141210)

  write(*,*) "area_occupancy_occ.out and area_occupancy_virt.out written successfully."
!-------------------------------------------------------------------------------------------------------------------------------------





 !call free_full_potential(dpcom%mpi_env%nproc,0,xc,potential)
   if (f_associated(pot_ref_count)) then
     call f_free_ptr(potential)
     call f_ref_free(pot_ref_count)
   end if
   nullify(potential)

  if (nproc>1) call fmpi_allreduce(epot_sum,1,op=FMPI_SUM)
  deallocate(pot_ion, rho_ion)
  call deallocate_work_arrays_sumrho(wisf)
    close(124)
!-------------------------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------------------------------------


    !----------------------
    ! Free allocated space.
    !----------------------
    call deallocate_comms(comms)
    call deallocate_locreg_descriptors(Lzd%Glr)
    call deallocate_Lzd_except_Glr(Lzd)
    call deallocate_orbs(orbs)
    call deallocate_atoms_data(atoms)
    call xc_end(xc)
    call dpbox_free(dpcom)
    call free_input_variables(inputs)
    call free_DFT_PSP_projectors(nlpsp)
    call pkernel_free(pkernel)
    
    !--------------------------------------
    !wait all processes before finalisation
    !--------------------------------------
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    call f_lib_finalize()
    close(6)  

stop

contains

  subroutine create_function_3D(mesh, target_centers, radii, f)
    ! Inputs
    ! integer, dimension(3), intent(in) :: grid_size  ! Grid dimensions (n_x, n_y, n_z)
    type(cell), intent(in) :: mesh
    real(wp), dimension(:), intent(in) :: radii  ! Radii for each sphere, unit = 1 grid length
    real(wp), dimension(3, size(radii)), intent(in) :: target_centers  ! Center coordinates (x, y, z) for each sphere
    ! real, dimension(:), intent(out) :: f  ! 1D array representing the 3D function
    real(wp), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(out) :: f
    ! Local variables
    integer :: n_x, n_y, n_z  ! Extracted grid sizes
    integer :: x, y, z, t, index, num_targets
    real :: dist_squared, radius_squared

    ! Extract grid sizes
    n_x = mesh%ndims(1)
    n_y = mesh%ndims(2)
    n_z = mesh%ndims(3)

    ! Number of target spheres
    num_targets = size(radii)

    ! Fill the 1D array `f` with values based on the 3D coordinates (x, y, z)
    do z = 1, n_z  ! Loop over z
      do y = 1, n_y  ! Loop over y
          do x = 1, n_x  ! Loop over x
            ! Calculate the 1D index corresponding to (x, y, z) in Fortran column-major order
            ! index = (z - 1) * n_x * n_y + (y - 1) * n_x + x

            ! Assume point is outside all spheres (default -1)
            f(x, y, z) = 1.0

            ! Check if point (x, y, z) is within any sphere
            do t = 1, num_targets
                radius_squared = radii(t)**2
                dist_squared = (real(x) - target_centers(1, t))**2 + &
                              (real(y) - target_centers(2, t))**2 + &
                              (real(z) - target_centers(3, t))**2

                if (dist_squared <= radius_squared) then
                  f(x, y, z) = -1.0  ! Inside at least one sphere
                  exit  ! No need to check further targets for this point
                end if
            end do
          end do
      end do
    end do

  end subroutine create_function_3D

  function int_to_str(i) result(s)
    integer, intent(in) :: i
    character(len=20) :: s
    write(s, '(I0)') i
  end function int_to_str

  function to_lowercase(str) result(lower_str)
    implicit none
    character(len=*), intent(in) :: str
    character(len=len(str)) :: lower_str
    integer :: i
  
    do i = 1, len(str)
      select case (str(i:i))
      case ('A':'Z')
        lower_str(i:i) = char(iachar(str(i:i)) + 32)
      case default
        lower_str(i:i) = str(i:i)
      end select
    end do
  end function to_lowercase
  
end program toy_model

integer function re_index1(i,itot,iocc)
implicit none
  integer :: i
  integer :: itot,iocc
  if(i .lt. iocc) then
    if(i .lt. iocc/2) then
      re_index1 =  i*2
    else
      re_index1 = (i-(iocc/2))*2 + 1
    end if
  else
    if(i .lt. (iocc+(itot-iocc)/2)) then
      re_index1 = iocc + (i-iocc)*2
    else
      re_index1 = 2*i - itot + 1
    end if
  end if
return
end
!------------------------------------------------------------
! -- Interleaved Spin Format    0 2 4 1 3 5    6 8 7 9  -----
!                                      transfer to           
! --       Block Spin Format    0 1 2 3 4 5    6 7 8 9  -----
!------------------------------------------------------------
! Number of Spin Up Electrons = Number of Spin Down Electrons
!------------------------------------------------------------
integer function re_index2(i,itot,iocc)
implicit none
  integer :: i
  integer :: itot,iocc
  if(i .lt. iocc) then
     if(mod(i,2) .eq. 0) then
        re_index2 = i/2
     else
        re_index2 = (i + iocc - 1)/2
     end if
  else
    if(mod(i,2) .eq. 0) then
      re_index2 = (i + iocc)/2 
    else
      re_index2 = (i + itot - 1)/2
    end if
  end if
return
end
