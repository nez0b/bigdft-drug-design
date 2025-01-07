!> @file
!! @author
!!    Copyright (C) 2007-2013 BigDFT group. This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt.
!!    For the list of contributors, see ~/AUTHORS

!===============================!
!> Generate .cube files using *_u.mat file.
!===============================!

program toy_model
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
  real(wp), dimension(:),   pointer :: w
  real(wp), dimension(:),   pointer :: psi, psiv, psir, psir_i,psir_j
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
  real(wp), dimension(:,:), allocatable :: phirr
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
  integer :: ierr, iproc, nproc, nwarnings, ii,jj,nn, n_virt, num_args
  integer :: iorb, nvirtu, nvirtd, ispsi, ilr, ilr_orb, npot, mspin, cspin
  integer :: ityp, nc, ishift, nwarning, dry_run, orbtot, nv, istate, n_u_mat, n_plot
  integer :: orbdimocc, orbdimvir, norb_start,norb_end,  orbocc,orbvir
  integer(kind=4)  :: OMP_get_max_threads, OMP_get_thread_num, OMP_get_num_threads
  integer :: max_ind, x, y, z

  character(len=100) :: cmd, name, arg
  character*32 con, tmp
  integer :: ip ,iq ,ir ,is , ihpq, ihpqrs, istat, i,j,k, re_index1, re_index2,loop_b
  integer :: ipt,iqt,irt,ist, nhpq, nhpqrs
  real(kind=8) :: hpq,hpqrs,hpq_overlap,hpq_abs_overlap

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


  ! Default values
  n_plot = 1
  name = "default_sys"

  call get_command_argument(1, name)
  call get_command_argument(2, arg)
  read(arg, *) n_plot

  n_virt = orbsv%norb/2
  n_u_mat = orbs%norb+n_virt
  allocate(u_matrix(n_u_mat*n_u_mat)) ; u_matrix=0._wp

  open(1111,file=trim(name)//'_u.mat')
  do i=1,4
    read(1111,*)
  end do
  do i=1, n_u_mat*n_u_mat
    read(1111,*) u_matrix(i), place_holder
  end do
  close(1111)


  allocate(wann_func((Lzd%Glr%wfd%nvctr_c + 7*Lzd%Glr%wfd%nvctr_f), n_plot)) ; wann_func=0._wp
  do j=1,n_plot
    do loop_b=1, orbs%norb
      wann_func(:,j) = &
      wann_func(:,j) + &
      u_matrix(loop_b+(j-1)*(orbs%norb+n_virt))*psi((loop_b-1)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1:(loop_b)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f))
    end do

    do loop_b=1, n_virt/2
      wann_func(:,j) = &
      wann_func(:,j) + &
      u_matrix(loop_b+orbs%norb+(j-1)*(orbs%norb+n_virt))*psiv((loop_b-1)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1:(loop_b)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f))
    end do
    call plot_wf(.false.,'wann_func' // int_to_str(j),1,atoms,1.0_wp,Lzd%Glr, &
    & Lzd%hgrids,atoms%astruct%rxyz, &
    & wann_func(:,j))
  end do

  write(*,*) n_plot, '.cube files generated.'

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
  function int_to_str(i) result(s)
    integer, intent(in) :: i
    character(len=20) :: s
    write(s, '(I0)') i
  end function int_to_str

end program toy_model
