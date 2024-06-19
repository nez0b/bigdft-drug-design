!> @file
!! @author
!!    Copyright (C) 2007-2013 BigDFT group. This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt.
!!    For the list of contributors, see ~/AUTHORS
!
!==================================================================================
! JP NOTES:
!   * input files: wavefunction files in data/ folder
!                  input.yaml
!                  input-hpqrs (112)   ! SKChou modifiy; See below
!    output files: toy_model.log (124)
!                  toy_model.out (6)
!
!   * The number of virtual orbitals must be less than occupied orbitals due to 
!     the dimension of virtual orbital is identical as occupied orbitals

!  gp = wp = dp = 8
!==================================================================================
! SKChou NOTES:                                                            May 2021
!
! phirr(:,:) - Spin Orbital Index ( V.S. psirr(:,:) )
!            - Interleaved Spin Format 
! * calculate orthogonalization
! * calculate hpqrs
!
! hpq   input  file: hpq.inp
! hpqrs input  file: hpqrs.inp
!
! hpq   output file: hpq.out
! hpqrs output file: hpqrs.out
!==================================================================================
! SKChou NOTES:                                                           June 2021
!
! Pseudopotential nonlocal part
!
! Consider the number of projectors (nproj)
!
!  proj: Spin Polarized Format Array (scalar products with projectors and orbitals)
! hproj: Spin Polarized Format Array (after application of the hamiltonian on proj)
!
!  proj: ( 4 , (# of orbitals)*nproj )
! hproj: ( 4 , (# of orbitals)*nproj ),
! where "4" is the nlpsp%(h)cproj size but only 1st element has a value.
!
! Array Index:
! [1st projector with occupied and then virtual orbitals; 2nd projector … ; 3rd …]
!
!
! Virtual orbital # in the nonlocal part are iterated in the block spin format:
! [1st projector #1 #2 #3 ... and then #1 #2 #3 ...; 2nd projector ...]
!
! nspin=1:
! Due to the spin polarized format proj and hproj, it needs to do re_index1
! from block spin format to interleaved spin format. 
! Since in E_nonlocal and Output3 array for nspin=1 are in interleaved spin format.
!
! nspin=2:
! Calculations are store E_nonlocal array which is in block spin formaat.
! Transform to Output3 array which is in interleaved spin format. [by re_index1]
!
!
! For atom system, the projector array is different than molecule system.
!  proj( nproj , # of orbitals )
! hproj( nproj , # of orbitals )
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SKChou NOTES:                                                      August 2, 2021 
!
! The Gift for myslef. 
! I found the unified projector structure! 
! 
! num_proj : For each atom with projector(s) counts 1
! nmproj   : Number of Projector Considered for each Atom 
!
! (h)proj(4*MAXVAL(nlpsp%projs%mproj), (# or orbitals)*num_proj), 
! where 
!   "4" is default size for an atom with nproj <= 4
!   MAXVAL(nlpsp%projs%mproj): the maximum single atomic projector which alters
!                              the size of psp_it%parent%(h)cproj 
!   (# or orbitals)*num_proj : the number of iterations
!
! The virtual part array structure is the same as before.
!==================================================================================
! SKChou NOTES:                                                           June 2021
!
! Some segmentation fault caused by unknown reason
!
! rm toy_model.x toy_model.log toy_model.out
!
! move close(6) and close(124) to suitable line 
!
! Consistent in the integer type !
!
!==================================================================================
! SKChou NOTES:                                                           June 2021
! 
! Number of Virtual Orbitals Limis Bug:
!
! BIGDFT_RUNTIME_ERROR:
!   Action:
!     Check the exact unrolling of runtime operations, likely something has been 
!     initialized/finalized twice
!
! The subroutine call twice !
! ==> call orbital_basis_release(psi_ob)
!
! The subroutine call twice !
! ==> if(associated(nlpsp%iagamma)) call f_zero(nlpsp%gamma_mmp)
! 
! For virtual orbitals, the subroutine needs not to call twice!
! call createProjectorsArrays(Lzd%Glr, atoms%astruct%rxyz,atoms,orbsv,&
!      inputs%crmult,inputs%frmult,inputs%projection,dry_run,nlpsp,init_projectors_completely)
!
! !!!!!!!!!!
! write(6,*) 
!==================================================================================
! SKChou NOTES:                                                       December 2021
!
! After Linux was reinstalled, BigDFT also updated the newest version (1.9.1;1.9.2).
! The usage of pkernel is changed:
!
! # change 
! !use at_domain, only: domain_geocode
! use at_domain, domain_renamed => domain !, only: domain, domain_new
!
! # add
! use numerics, only: onehalf,pi
!
! # add
! dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(domain_geocode(atoms%astruct%dom)),&
!            alpha_bc=onehalf*pi,beta_ac=onehalf*pi,gamma_ab=onehalf*pi,&
!            acell=(/Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i/)*&
!                  (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/))
!
! # change
! !pkernel = pkernel_init(iproc,nproc,dict, atoms%astruct%geocode, &
! !          (/Lzd%Glr%d%n1i, Lzd%Glr%d%n2i, Lzd%Glr%d%n3i/), (/inputs%hx/2._gp,inputs%hy/2._gp,inputs%hz/2._gp/)
! pkernel=pkernel_init(iproc,nproc,dict,&
!      !domain_geocode(atoms%astruct%dom),(/Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i/),&
!      dom,(/Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i/),&
!      (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/))
!==================================================================================

!===============================!
!> Toy program to use BigDFT API!
!===============================!

program toy_model
  ! use omp_lib
  use dynamic_memory_base, only: f_memcpy
  use public_enums
  use module_fragments
  use f_enums
  use module_f_malloc
  use module_interfaces, only: system_initialization, readmywaves_linear_new
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
  type(local_zone_descriptors) :: Lzdc, lzdgb
  type(orbitals_data)          :: orbsc
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
  type(DFT_wavefunction) :: tmb
  type(f_enumerator) ::  inputpsi
  type(DFT_local_fields) :: denspot
  type(system_fragment), dimension(:), pointer :: ref_frags

  logical  :: dosome, alive, init_projectors_completely
! logical  :: paw = .false.
! integer, parameter :: dp=8
! integer, parameter :: gp=8
  integer, parameter :: wp=8
  real(dp) :: nrm, epot_sum
  real(dp) :: psoffset, ekin, ekin_sum
  real(gp), dimension(3)            :: shift
  real(wp), dimension(:),   pointer :: w
  real(wp), dimension(:),   pointer :: psir, psir_i
  real(dp), dimension(:),   pointer :: rhor,rho_ion
  real(wp), dimension(:),   allocatable :: pot_ion
  real(wp), dimension(:),   pointer :: potential
  real(wp), dimension(:)  , pointer :: hpsi_ptr
  real(gp), dimension(:,:), pointer :: rxyz
  real(gp), dimension(:,:), pointer :: rxyz_old
  real(wp), dimension(:,:), pointer :: ovrlp, grid_xyz
  real(dp), dimension(:,:), pointer :: rho_p => null() !needs to be nullified
  real(wp), dimension(:,:), allocatable :: pot_tmp
  real(kind=8),dimension(:,:),pointer :: locregcenters
  integer,dimension(:,:,:),pointer :: frag_env_mapping
  real(kind=8),dimension(:),pointer :: psi_g

  ! Interpolating Scaling Function (daub_to_ISF)
  real(wp), dimension(:,:), allocatable :: psirr
  ! Spin Polarized Format of psirr
  real(wp), dimension(:,:), allocatable :: phirr
  ! phirr_overlap = phirr(:,i)*phirr(:,j)
  real(wp), dimension(:), allocatable :: phirr_overlap, grid_r
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
  integer :: ierr, iproc, nproc, nwarnings, ii,jj,nn
  integer :: iorb, nvirtu, nvirtd, ispsi, ilr, ilr_orb, npot, mspin, cspin
  integer :: ityp, nc, ishift, nwarning, dry_run, orbtot, nv, istate
  integer :: orbdimocc, orbdimvir, norb_start,norb_end,  orbocc,orbvir
  integer(kind=4)  :: OMP_get_max_threads, OMP_get_thread_num, OMP_get_num_threads
  integer :: x, y, z
  integer :: input_wf_format, max_nbasis_env, input_mat_format
  integer,dimension(:),allocatable :: in_which_locreg_
  integer :: sdim, ldim, istp

  character*32 con, tmp
  integer :: ip ,iq ,ir ,is ,ip_ ,iq_ , is_, ihpq, ihpqrs, istat, i,j,k, re_index1, re_index2
  integer :: ipt,iqt,irt,ist, nhpq, nhpqrs, increment
  real(kind=8) :: hpq,hpqrs,hpq_overlap,hpq_abs_overlap,start_time,end_time

  ! PSolver
  ! hpqrs 
  real(dp), dimension(:), allocatable :: rhopq, rhors, rhopqxphirrir

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

  lzdgb = default_Lzd() ; lzdgb%hgrids=(/ inputs%hx, inputs%hy, inputs%hz /)
  call lr_set(lzdgb%Glr, iproc, GPU%OCLconv, .true., inputs%crmult, inputs%frmult, &
              lzdgb%hgrids,atoms%astruct%rxyz,atoms,.true.,.false.)


!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------
  inputpsi = inputs%inputPsiId
  locregcenters=f_malloc_ptr((/3,atoms%astruct%nat/),id=' locregcenters')
  locregcenters = atoms%astruct%rxyz

  call system_initialization(iproc,nproc,.true.,inputpsi,input_wf_format,.false.,inputs,atoms,atoms%astruct%rxyz,GPU%OCLconv,&
  orbsc,tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp,tmb%orbs,Lzdc,tmb%Lzd,nlpsp,&
  comms,ref_frags,denspot=denspot,locregcenters=locregcenters)

  call lr_set(Lzdc%Glr, iproc, GPU%OCLconv, .true., inputs%crmult, inputs%frmult, &
              Lzdc%hgrids,atoms%astruct%rxyz,atoms,.true.,.false.)

  call check_linear_and_create_Lzd(iproc,nproc,inputs%linear,Lzdc,atoms,tmb%orbs,inputs%nspin,atoms%astruct%rxyz)

  tmb%npsidim_orbs = tmb%orbs%npsidim_orbs
  tmb%npsidim_comp = tmb%orbs%npsidim_comp

  if(tmb%orbs%nspin .eq. 1) cspin = 2 ; if(tmb%orbs%nspin .eq. 1) mspin = 1
  if(tmb%orbs%nspin .eq. 2) cspin = 1 ; if(tmb%orbs%nspin .eq. 2) mspin = 2
  orbocc = tmb%orbs%norb
  orbvir = 0
  orbtot = tmb%orbs%norb
  allocate( E_kin(      orbtot,orbtot) ) ; E_kin=0.d0
  allocate( E_local(    orbtot,orbtot) ) ; E_local=0.d0
  allocate( E_nonlocal( orbtot,orbtot) ) ; E_nonlocal=0.d0

  allocate(tmb%orbs%eval(tmb%orbs%norb*tmb%orbs%nkpts)) ; call f_zero(tmb%orbs%eval)
  allocate(tmb%psi(  max(tmb%orbs%npsidim_orbs, tmb%orbs%npsidim_comp)+1 )) ; tmb%psi=0._gp
  allocate(rxyz_old(3, atoms%astruct%nat)) ; rxyz_old=0._gp

  !------------------i--------------------------------------------------!
  ! Read occupied state wavefunctions from disk and store them in psi. !
  ! orbs%norb - number of occupied orbitals                            !
  !--------------------------------------------------------------------!
  allocate(ref_frags(1))
  call init_fragments(inputs, tmb%orbs, atoms%astruct, ref_frags)

  max_nbasis_env = max(0,ref_frags(1)%nbasis_env)
  frag_env_mapping = f_malloc0_ptr((/inputs%frag%nfrag,max_nbasis_env,3/),id='frag_env_mapping')
  input_mat_format = inputs%lin%output_mat_format 

  call readmywaves_linear_new(iproc,nproc,'data/','minBasis',&
        1,input_mat_format,atoms,tmb,atoms%astruct%rxyz,ref_frags,&
        inputs%frag,inputs%lin%fragment_calculation,&
        .false.,.false.,&
        max_nbasis_env,frag_env_mapping)
  !--------------------------------------------------------------------!
  if(nproc>1) call fmpi_allreduce(tmb%orbs%eval(1), tmb%orbs%norb*tmb%orbs%nkpts, op=FMPI_SUM)

  orbdimocc = (tmb%Lzd%Glr%wfd%nvctr_c+7*tmb%Lzd%Glr%wfd%nvctr_f)*tmb%orbs%norb
  orbdimvir = 0

  !--------------------------------------
  ! local tmb%psi to global psi_g
  !--------------------------------------

  in_which_locreg_ = f_malloc(tmb%orbs%norb,id='in_which_locreg_')
  call f_memcpy(src=tmb%orbs%inwhichlocreg, dest=in_which_locreg_)
  psi_g = f_malloc_ptr(array_dim(lzdc%glr)*tmb%orbs%norb, id='psi_g')
  call f_zero(psi_g)
  istp = 1
  do i=1,tmb%orbs%norb
    if(tmb%orbs%nspin .eq. 1) then
      ilr = in_which_locreg_(i)
      sdim=array_dim(tmb%lzd%llr(ilr))
      ldim=array_dim(tmb%lzd%glr)
      call lpsi_to_global2(sdim, ldim, tmb%orbs%norb, 1, tmb%lzd%glr, &
                          tmb%lzd%llr(ilr), tmb%psi(istp:), psi_g(1+(i-1)*ldim:i*ldim))
      istp = istp + sdim
    end if
  end do

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

    call system('echo "kinetic calculation"')
    !--------------------------------------
    ! kinetic energy calculation
    !--------------------------------------
    write(124,*) "================"
    write(124,*) " kinetic energy "
    write(124,*) "================"
    allocate( tpsi_o(orbdimocc), indwvl_i(orbtot) ) ; tpsi_o = 0._dp
    allocate(indwvl_f(orbtot))
  
    !------------------------------------------------------
    ! loop on the localisation regions (occupied orbitals)
    !------------------------------------------------------
    allocate(tpsi(orbdimocc)) ; tpsi = 0._dp
    write(124,20) tmb%orbs%norb, tmb%orbs%nspinor, mspin, orbdimocc, tmb%Lzd%Glr%wfd%nvctr_c+7*tmb%Lzd%Glr%wfd%nvctr_f

    call initialize_work_arrays_locham(Lzdc%nlr,Lzdc%Llr,tmb%orbs%nspinor,wrk_lh)
    ekin = 0.d0 ; ekin_sum = 0.0_gp 
    dosome=.false.   ! check if this localisation region is used by one of the orbitals
    do iorb=1,tmb%orbs%norb ; dosome = (tmb%orbs%inwhichlocreg(iorb+tmb%orbs%isorb) == 1) ; if(dosome) exit ; end do
    if (.not. dosome) error stop "dosome is false. No kinetic energy calculation."

    call initialize_work_arrays_locham(Lzdc%Llr(1),tmb%orbs%nspinor,wrk_lh)
    ispsi = 1
    loop_orbs: do iorb=1,tmb%orbs%norb
      indwvl_i(iorb) = ispsi
!+++++++++++++error++++++++++++here+++++++++++++++++++
      !call psi_to_tpsi(Lzd%hgrids, orbs%kpts(1,orbs%iokpt(iorb)),orbs%nspinor,Lzd%Llr(ilr), psi(ispsi), wrk_lh, tpsi(ispsi), ekin)
      
      call psi_to_tpsi(tmb%orbs%kpts(1,tmb%orbs%iokpt(iorb)),tmb%orbs%nspinor,Lzdc%Llr(1), psi_g(ispsi), wrk_lh,tpsi(ispsi), ekin)
      ekin_sum = ekin_sum + tmb%orbs%kwgts(tmb%orbs%iokpt(iorb))*tmb%orbs%occup(iorb+tmb%orbs%isorb)*ekin
      ispsi  = ispsi + array_dim(lzdc%glr)*tmb%orbs%nspinor
      indwvl_f(iorb) = ispsi
      write(124,22) iorb, ekin, indwvl_i(iorb),indwvl_f(iorb)
    end do loop_orbs
    write(124,*) " occupied orbitals ... DONE"
    write(124,'("total kinetic energy = ",f15.8)') ekin_sum ; write(124,*)
    tpsi_o=tpsi
    call deallocate_work_arrays_locham(wrk_lh)
    deallocate(tpsi)

  20 format("   # occupied orbitals:",i3,", nspinor:",i3,", mspin: ",i4,", psi_dim:",3i9)
  21 format("   # virtual  orbitals:",i3,", nspinor:",i3,", mspin: ",i4,", psi_dim:",3i9)
  22 format(" orbital ",i3,"  ekin = ",f15.8,2x,"index of wavelet:",2i10)

    do i=1,   tmb%orbs%norb       ; do j=1,  tmb%orbs%norb
      E_kin(i,j) = sum( psi_g( indwvl_i(i):indwvl_f(i)-1) * tpsi_o(indwvl_i(j):indwvl_f(j)-1) ) 
    end do ; end do
 
    deallocate(tpsi_o, indwvl_i, indwvl_f)
    call system('echo "kinetic calculation ... DONE"')

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

  call system('echo "potential calculation (local)"')
  !---------------------------------------------------------------------
  !          energy of the local potential, unit in Hartree         -
  !---------------------------------------------------------------------
  write(124,*) 
  write(124,*) "================================================"
  write(124,*) " energy of the LOCAL potential, unit in Hartree "
  write(124,*) "       -- psir(i)*potential*psir(j) --       "
  write(124,*) "================================================"
  
  allocate(psir(  tmb%Lzd%Glr%nboxi(2,1)*tmb%Lzd%Glr%nboxi(2,2)*tmb%Lzd%Glr%nboxi(2,3)))
  allocate(psir_i( tmb%Lzd%Glr%nboxi(2,1)*tmb%Lzd%Glr%nboxi(2,2)*tmb%Lzd%Glr%nboxi(2,3)))
  allocate(psirr(  tmb%Lzd%Glr%nboxi(2,1)*tmb%Lzd%Glr%nboxi(2,2)*tmb%Lzd%Glr%nboxi(2,3),0:tmb%orbs%norb-1));psirr=0._wp
 
  allocate(phirr(tmb%Lzd%Glr%nboxi(2,1)*tmb%Lzd%Glr%nboxi(2,2)*tmb%Lzd%Glr%nboxi(2,3) ,0:tmb%orbs%norb*cspin-1)) ; phirr=0._wp
  allocate(phirr_overlap(tmb%Lzd%Glr%nboxi(2,1)*tmb%Lzd%Glr%nboxi(2,2)*tmb%Lzd%Glr%nboxi(2,3))) ; phirr_overlap=0._wp
  allocate(grid_xyz(tmb%Lzd%Glr%nboxi(2,1)*tmb%Lzd%Glr%nboxi(2,2)*tmb%Lzd%Glr%nboxi(2,3), 3)) ; grid_xyz=0._wp
  allocate(grid_r(tmb%Lzd%Glr%nboxi(2,1)*tmb%Lzd%Glr%nboxi(2,2)*tmb%Lzd%Glr%nboxi(2,3))) ; grid_r=0._wp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
 ! call local_potential_dimensions(Lzd, orbs, dpcom%ngatherarr(0,1))
  call xc_init(xc, inputs%ixc, XC_ABINIT, inputs%nspin)
  call dpbox_set(dpcom,Lzdc%Glr%mesh,xc,iproc,nproc,MPI_COMM_WORLD, inputs%SIC%approach, inputs%nspin)
  dict => dict_new()
  
  dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(domain_geocode(atoms%astruct%dom)),&
            alpha_bc=onehalf*pi,beta_ac=onehalf*pi,gamma_ab=onehalf*pi,&
            acell=(/tmb%Lzd%Glr%nboxi(2,1), tmb%Lzd%Glr%nboxi(2,2), tmb%Lzd%Glr%nboxi(2,3)/)*&
                  (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/))

  pkernel=pkernel_init(iproc,nproc,dict,&
      dom,(/tmb%Lzd%Glr%nboxi(2,1), tmb%Lzd%Glr%nboxi(2,2), tmb%Lzd%Glr%nboxi(2,3)/),&
       (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/))

  call dict_free(dict)
  call pkernel_set(pkernel,verbose=.false.)
  !nullify(pot_ion) ;
  allocate(pot_ion(tmb%Lzd%Glr%nboxi(2,1)*tmb%Lzd%Glr%nboxi(2,2)*dpcom%n3p)) ; pot_ion=0._gp
  nullify(rho_ion) ; allocate(rho_ion(tmb%Lzd%Glr%nboxi(2,1)*tmb%Lzd%Glr%nboxi(2,2)*dpcom%n3p)) ; rho_ion=0._gp
  write(124,*) " call createIonicPotential"
  
  call createIonicPotential(iproc, (iproc==0),atoms, atoms%astruct%rxyz,&
   inputs%elecfield, dpcom, pkernel, pot_ion, rho_ion, psoffset)
  
 ! allocate the potential in the full box
  ! Lzdc%linear = .true.
  call full_local_potential(iproc,nproc,orbsc,Lzdc,dpcom,xc,pot_ion,potential,pot_ref_count)


  call initialize_work_arrays_sumrho(Lzdc%Glr , .true., wisf)

!------------------------------------------------------------------------------------------

!!!

!!!
  do i=1,tmb%orbs%norb
    if(tmb%orbs%nspin .eq. 1) then
      ldim=array_dim(tmb%lzd%glr)
      call daub_to_isf(lzdgb%Glr,wisf,psi_g(1+(i-1)*ldim:i*ldim),psir_i)
      phirr(:,2*(i-1)+0)=psir_i
      phirr(:,2*(i-1)+1)=psir_i
    else if(tmb%orbs%nspin .eq. 2) then
      call daub_to_isf(tmb%Lzd%Glr,wisf,tmb%psi((i-1)*(tmb%Lzd%Glr%wfd%nvctr_c+7*tmb%Lzd%Glr%wfd%nvctr_f)+1) ,psir_i)
      ! Block Spin Format
      !phirr(:,i-1)=psir_i
      ! Interleaved Spin Format
      ii=re_index1(i-1,tmb%orbs%norb*cspin,tmb%orbs%norb*cspin)
      phirr(:,ii)=psir_i
    end if
  end do




!------------------------------------------------------------------------------------------
  if(tmb%orbs%nspin .eq. 1) then
    !$omp parallel default(private) shared(orbtot, E_local, phirr, potential)
    !$omp do
    do i=1,orbtot
      do j=1,orbtot
        E_local(i,j) = 0!sum(phirr(:,2*i-2)*potential(:)*phirr(:,2*j-2))
      end do
    end do 
    !$omp end do
    !$omp end parallel
  else
    error stop "line 564: nspin/=1 is not supported."
  end if


 !call free_full_potential(dpcom%mpi_env%nproc,0,xc,potential)
   if (f_associated(pot_ref_count)) then
     call f_free_ptr(potential)
     call f_ref_free(pot_ref_count)
   end if
   nullify(potential)

  if (nproc>1) call fmpi_allreduce(epot_sum,1,op=FMPI_SUM)
  deallocate(pot_ion, rho_ion, psir)
  call deallocate_work_arrays_sumrho(wisf)
  call system('echo "potential calculation (local) ... DONE"')


!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

    call system('echo "potential calculation (nonlocal)"')
  !---------------------------------------------------------------------
  !          energy of the nonlocal potential, unit in Hartree         -
  !---------------------------------------------------------------------

    write(124,*) 
    write(124,*) "==============="
    write(124,*) " nonlocal term "
    write(124,*) "==============="

    !-------------------
    !  occupied orbitals
    !-------------------
    allocate(hpsi(orbdimocc)) 
    hpsi_ptr => ob_ket_map(hpsi,psi_it) 
    if(orbdimocc > 0) call f_zero(hpsi(1),orbdimocc)

    init_projectors_completely = .true. ; dry_run = 0 ; nwarnings = 0

    ! cpmult; fpmult
    call createProjectorsArrays(Lzdc%Glr, atoms%astruct%rxyz,atoms,orbsc,&
                                inputs%crmult,inputs%frmult,inputs%projection,dry_run,nlpsp,init_projectors_completely)

    ! Number of Atomic Species Projector    
    num_proj=0
    do iat=1,atoms%astruct%nat
      ! Atom with projector(s) counts 1 
      if(nlpsp%projs(iat)%mproj .ne. 0) then
        num_proj=num_proj+1
      end if
    end do

    ! Allocate proj and hrpoj
    allocate( proj(4*MAXVAL(nlpsp%projs%mproj), 0:((tmb%orbs%norb*cspin)*num_proj)-1)) ;  proj=0._dp
    allocate(hproj(4*MAXVAL(nlpsp%projs%mproj), 0:((tmb%orbs%norb*cspin)*num_proj)-1)) ; hproj=0._dp


    ! nn=1:orbs%norb
    nn=1
    ! Number of Projector Index
    iproj=1
    ! Projector Index
    if(tmb%orbs%nspin .eq. 1) then
      nnproj=tmb%orbs%norb
    else if(tmb%orbs%nspin .eq. 2) then
      nnproj=tmb%orbs%norb*cspin
    end if

     ! call orbital_basis_associate(psi_ob, orbs=tmb%orbs, phis_wvl=tmb%psi, Lzd=tmb%Lzd, id='nonlocalham')
    call orbital_basis_associate(psi_ob, orbs=orbsc, phis_wvl=psi_g, Lzd=Lzdc, id='nonlocalham')
    if(associated(nlpsp%iagamma)) call f_zero(nlpsp%gamma_mmp)
    psi_it = orbital_basis_iterator(psi_ob)
    loop_kpt: do while(ket_next_kpt(psi_it))
      loop_lr: do while(ket_next_locreg(psi_it, ikpt=psi_it%ikpt))
        call DFT_PSP_projectors_iter_new(psp_it, nlpsp)
        loop_proj: do while (DFT_PSP_projectors_iter_next(psp_it, ilr=psi_it%ilr, lr=psi_it%lr, glr=Lzdc%glr))
          call DFT_PSP_projectors_iter_ensure(psp_it, psi_it%kpoint, 0, nwarnings, Lzdc%Glr)
          loop_psi_kpt: do while(ket_next(psi_it, ikpt=psi_it%ikpt, ilr=psi_it%ilr))          
            call DFT_PSP_projectors_iter_apply(psp_it, psi_it, atoms, eproj, hpsi=hpsi)          

            ! Number of Projector Considered for each Atomic Species
            nmproj=nlpsp%projs(iproj)%mproj
            if(nmproj .eq. 0) then
              iproj=iproj+1
              nmproj=nlpsp%projs(iproj)%mproj
            end if

            if(tmb%orbs%nspin .eq. 1) then
   
              ! Interleaved Spin Format 
              proj(1:nmproj,2*(nn-1)+0)=psp_it%parent%cproj(1:nmproj)
              proj(1:nmproj,2*(nn-1)+1)=psp_it%parent%cproj(1:nmproj)
              hproj(1:nmproj,2*(nn-1)+0)=psp_it%parent%hcproj(1:nmproj)
              hproj(1:nmproj,2*(nn-1)+1)=psp_it%parent%hcproj(1:nmproj)

              ! Next Projector 
              if(nn==nnproj) then
                nn=nn
                nnproj=nnproj+(tmb%orbs%norb)
                iproj=iproj+1
              end if

            else if(tmb%orbs%nspin .eq. 2) then

              ! Block Spin Format
              proj(1:nmproj,nn-1)=psp_it%parent%cproj(1:nmproj)
              hproj(1:nmproj,nn-1)=psp_it%parent%hcproj(1:nmproj)

              ! Next Projector 
              if(nn-1==nnproj-1) then
                nn=nn
                nnproj=nnproj+(tmb%orbs%norb*cspin)
                iat=iat+1
              end if
            end if

            ! Next Orbital
            nn=nn+1

          end do loop_psi_kpt
        end do loop_proj
      end do loop_lr
    end do loop_kpt

    deallocate(hpsi)

    call orbital_basis_release(psi_ob)

    call system('echo "potential calculation (nonlocal) ... DONE"')


    !! Number of Projector
    !do k=1,nlpsp%nproj
    !  ! Next Projector Index
    !  nn=(orbs%norb*cspin+orbsv%norb)*(k-1)
    !  if(orbs%nspin .eq. 1) then
    !    do i=1,orbs%norb+inputs%norbv
    !      write(*,*) i, proj(1,(i-1)*2+nn),proj(1,(i-1)*2+1+nn)
    !    end do
    !  else if(orbs%nspin .eq. 2) then
    !    do i=1,orbs%norb*cspin+orbsv%norb
    !      write(*,*) i,nn, proj(1,(i-1)+nn)
    !    end do
    !  end if
    !end do


    ! Close File
    close(124)
!-------------------------------------------------------------------------------------------------------------------------------------

  ! Pseudopotential nonlocal part  

  E_nonlocal=0.d0
  if(tmb%orbs%nspin .eq. 1) then
    ! Number of Projector
    do k=1,num_proj
      ! Next Projector Index
      nn=(tmb%orbs%norb*cspin)*(k-1)
      do i=1,tmb%orbs%norb
        do j=1,tmb%orbs%norb
          E_nonlocal(i,j) = E_nonlocal(i,j) + sum(proj(:,(i-1)*2+nn)*hproj(:,(j-1)*2+nn))
        end do
      end do
    end do
  else if(tmb%orbs%nspin .eq. 2) then
    ! Number of Projector
    do k=1,num_proj
      ! Next Projector Index
      nn=(tmb%orbs%norb*cspin)*(k-1)
      do i=1,tmb%orbs%norb*cspin
        do j=1,tmb%orbs%norb*cspin
          E_nonlocal(i,j) = E_nonlocal(i,j) + sum(proj(:,(i-1)+nn)*hproj(:,(j-1)+nn))
        end do
      end do
    end do
  end if

  ! Check with .yaml - Enl:
  !Enl=0.d0
  !if(orbs%nspin .eq. 1) then
  !  do i=1,orbs%norb
  !    Enl=Enl+E_nonlocal(i,i)
  !  end do
  !  Enl=Enl*2.d0
  !else if(orbs%nspin .eq. 2) then
  !  do i=1,orbs%norb*cspin
  !    Enl=Enl+E_nonlocal(i,i)  
  !  end do
  !end if
  !write(*,*) "Enl:", Enl

  ! Note:
  ! E_nonlocal uses Block Spin Format
  ! output3 uses Interleaved Spin Format  

  deallocate(rxyz_old, tmb%psi, psir_i)
  deallocate(proj,hproj)

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------
  !----------------
  !   data output 
  !----------------
  write(6,*) 
  write(6,'("number of occupied orbitals - orbs%norb  : ",5i4)') tmb%orbs%norb
  write(6,'("    input parameter norbv - inputs%norbv : ",5i4)') inputs%norbv
  write(6,'("number of    total orbitals -     orbtot : ",5i4)') orbtot
  write(6,'("                            orbs%nspinor : ",5i4)') tmb%orbs%nspinor
  write(6,'("                                   mspin : ",5i4)') mspin
  write(6,'("                              OMP number : ",5i4)') OMP_get_thread_num(), OMP_get_num_threads(), OMP_get_max_threads()
  write(6,'("                          orbs%npsidim_orbs : ", i8)') tmb%orbs%npsidim_orbs
  write(6,'(" dimension of occupied orbitals (orbdimocc) : ", i8)') orbdimocc
  write(6,'(" dimension of  virtual orbitals (orbdimvir) : ", i8)') orbdimvir
  
  
  allocate( output1(0:orbtot-1,0:orbtot-1) ) ; output1=0._dp
  allocate( output2(0:orbtot-1,0:orbtot-1) ) ; output2=0._dp
  allocate( output3(0:orbtot-1,0:orbtot-1) ) ; output3=0._dp
  
  if(tmb%orbs%nspin .eq. 2) then
    do i=1,orbtot ; do j=1,orbtot

      !! Open System
      !! Block Spin Format
      !!ii=i-1
      !!jj=j-1
      ! Closed System: # of Spin Up = # of Spin Down
      ! Interleaved Spin Format
      ii = re_index1(i-1,orbtot,orbocc)
      jj = re_index1(j-1,orbtot,orbocc)

    output1(ii,jj) =      E_kin(i,j)
    output2(ii,jj) =    E_local(i,j)
    output3(ii,jj) = E_nonlocal(i,j)
!   write(4567,'(2i3,2x,2i3,2x,f15.8)') i-1,j-1, ii,jj, E_kin(ii,jj)+E_local(ii,jj)+E_nonlocal(ii,jj)
    end do ; end do
   
    write(6,*)  
    write(6,*) ; write(6,*) " =====  kinetic energy: psi(i)*tpsi(j) =====",orbtot
    do j=0,orbtot-1; do i=0,orbtot-1
      write(6,1,advance="NO") output1(i,j) ; if(i .eq. orbocc/2) write(6,"(A)",advance="NO") "  " ; end do
      write(6,"(A)") " " ; if(j .eq. orbocc/2) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  local potential: psir(i)*potential(ij)*psir(j)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1
      write(6,1,advance="NO") output2(i,j) ; if(i .eq. orbocc/2) write(6,"(A)",advance="NO") "  " ; end do
      write(6,"(A)") " " ; if(j .eq. orbocc/2) write(6,"(A)") "  " ; end do
    write(6,*) ; write(6,*) " =====  nonlocal potential: psir(i)*hpsir(j)  ===== "

    do j=0,orbtot-1; do i=0,orbtot-1 
      write(6,1,advance="NO") output3(i,j) ; if(i .eq. orbocc/2) write(6,"(A)",advance="NO") "  " ; end do
      write(6,"(A)") " " ; if(j .eq. orbocc/2) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  hpq:  sum(kinetic, V_local, V_nonlocal)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1 
      write(6,1,advance="NO") output1(i,j) + output2(i,j) + output3(i,j)  
      if(i .eq. orbocc/2) write(6,"(A)",advance="NO") "  "
      end do ; write(6,"(A)") " " ; if(j .eq. orbocc/2) write(6,"(A)") "  " ; end do

  else

    do i=0,orbtot-1 ; do j=0,orbtot-1
     output1(i,j) =      E_kin(i+1,j+1)
     output2(i,j) =    E_local(i+1,j+1)
     output3(i,j) = E_nonlocal(i+1,j+1)
    end do ; end do
   
    write(6,*)  
    write(6,*) ; write(6,*) " =====  kinetic energy: psi(i)*tpsi(j) ====="
    do j=0,orbtot-1; do i=0,orbtot-1; write(6,1,advance="NO") output1(i,j); if(i .eq. tmb%orbs%norb) write(6,"(A)",advance="NO") "  " &
      ; end do
     write(6,"(A)") " " ; if(j .eq. tmb%orbs%norb) write(6,"(A)") "  " ; end do 

    write(6,*) ; write(6,*) " =====  local potential: psir(i)*potential(ij)*psir(j)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1; write(6,1,advance="NO") output2(i,j); if(i .eq. tmb%orbs%norb) write(6,"(A)",advance="NO") "  " &
      ; end do
      write(6,"(A)") " " ; if(j .eq. tmb%orbs%norb) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  nonlocal potential: psir(i)*hpsir(j)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1; write(6,1,advance="NO") output3(i,j); if(i .eq. tmb%orbs%norb) write(6,"(A)",advance="NO") "  " &
      ; end do
      write(6,"(A)") " " ; if(j .eq. tmb%orbs%norb) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  hpq:  sum(kinetic, V_local, V_nonlocal)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1; write(6,1,advance="NO") output1(i,j) + output2(i,j) + output3(i,j)  
      if(i .eq. tmb%orbs%norb) write(6,"(A)",advance="NO") "  "
      end do ; write(6,"(A)") " " ; if(j .eq. tmb%orbs%norb) write(6,"(A)") "  " ; end do
    
  end if
  close(6)
! ------------------------------------------------------------------------------------------------------------------------------------


  !write(*,*)
  !write(*,*) "PHIRR, Orthogonalization"
  !do i=0,orbs%norb*cspin+orbsv%norb-1
  !  do j=0,orbs%norb*cspin+orbsv%norb-1
  !    !write(*,'(f12.5)',advance="NO"), sum(phirr(:,i)*phirr(:,j))
  !  end do
  !  !write(*,*)
  !end do

  ! Orbital Number Sign Flip
  !write(*,*)
  !do i=0,orbs%norb*cspin+orbsv%norb-1,2
  !  if(sum(phirr(:,i)*phirr(:,i+1)) .lt. 0._gp) then
  !    write(*,'(2x,"Orbital Number Sign Flip: ", i2)', advance="NO"), i
  !    phirr(:,i)=-phirr(:,i)
  !  end if
  !end do
 

! ------------------------------------------------------------------------------------------------------------------------------------
  
    ! Output hpq Index & hpq
    open(05132021,file="hpq.out")

    ! Read hpq Index input file
    open(0513,file="hpq.inp")
    nhpq=0
    do 
      read(0513,*,iostat=istat) tmp
      if(istat .ne. 0) exit
      nhpq=nhpq+1
    end do
    close(0513)

    open(0513,file="hpq.inp")
    do ihpq=1,nhpq
      read(0513,*) ip,iq
      ipt=ip
      iqt=iq
      if(tmb%orbs%nspin .eq. 1) then
        ip=int(ip/2)
        iq=int(iq/2)
      else if(tmb%orbs%nspin .eq. 2) then
        ip=ip
        iq=iq
      end if

      hpq=output1(ip,iq)+output2(ip,iq)+output3(ip,iq)

      write(05132021,'(2i4,f19.12)') ipt,iqt,hpq
    end do
    close(0513)
    close(05132021)


 !1 format(f16.12)
 1 format(f10.6)
 2 format(5x,2i3,3x,f20.16,6x,2i3)
  deallocate(E_local, E_nonlocal, E_kin)
  deallocate(output1, output2, output3)

!-------------------------------------------------------------------------------------------------------------------------------------

! SKChou 
! PSolver - Calculate hpqrs

!-------------------------------------------------------------------------------------------------------------------------------------
    open(6,file="toy_model.out")
    ! Electrostat_Solver should be called under open( ,file);
    ! otherwise, there is a additional information aboule PSolver printed out on the screen

    open(05132021,file='hpqrs.out')
 
    allocate(rhopqxphirrir(Lzdc%Glr%nboxi(2,1)*Lzdc%Glr%nboxi(2,2)*Lzdc%Glr%nboxi(2,3)))
    allocate(rhopq(Lzdc%Glr%nboxi(2,1)*Lzdc%Glr%nboxi(2,2)*Lzdc%Glr%nboxi(2,3)))
    allocate(rhors(Lzdc%Glr%nboxi(2,1)*Lzdc%Glr%nboxi(2,2)*Lzdc%Glr%nboxi(2,3)))

    if(tmb%orbs%nspin .eq. 1) increment = 2
    if(tmb%orbs%nspin .eq. 2) increment = 1 

    do ip=0,orbtot*2-1,increment
      do iq=ip,orbtot*2-1,2
        rhopq(:)=phirr(:,ip)*phirr(:,iq)
        rhopq(:)=rhopq(:)/pkernel%mesh%volume_element
        call Electrostatic_Solver(pkernel,rhopq)

        ir=ip
        rhopqxphirrir(:) = rhopq(:)*phirr(:,ir)
        !$omp parallel default(private) shared(orbtot, ip, iq, ir, phirr, rhopqxphirrir)
        !$omp do ordered
        do is=iq,orbtot*2-1,2
          hpqrs=sum(rhopqxphirrir(:)*phirr(:,is))

          !$omp ordered
          write( * ,'(4i4,f19.12)') ip,iq,ir,is, hpqrs
          write(05132021,'(4(i4),f19.12)') ip,iq,ir,is, hpqrs
          !$omp end ordered
        end do
        !$omp end do
        !$omp end parallel

        do ir=ip+increment,orbtot*2-1,increment
          rhopqxphirrir(:) = rhopq(:)*phirr(:,ir)
          !$omp parallel default(private) shared(orbtot, ip, iq, ir, phirr, rhopqxphirrir)
          !$omp do ordered
          do is=ir,orbtot*2-1,2
            ! rhors(:)=phirr(:,ir)*phirr(:,is)
            hpqrs=sum(rhopqxphirrir(:)*phirr(:,is))
            
            !$omp ordered
            write( * ,'(4i4,f19.12)') ip,iq,ir,is, hpqrs
            write(05132021,'(4(i4),f19.12)') ip,iq,ir,is, hpqrs
            !$omp end ordered
          end do
          !$omp end do
          !$omp end parallel
        end do
      end do
    end do

    close(05132021)
    deallocate(psirr,rhopq,rhors)
    deallocate(phirr)
    
!-------------------------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------------------------------------


    !----------------------
    ! Free allocated space.
    !----------------------
    call deallocate_comms(comms)
    call deallocate_locreg_descriptors(tmb%Lzd%Glr)
    call deallocate_Lzd_except_Glr(tmb%Lzd)
    call deallocate_orbs(tmb%orbs)
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
  subroutine unravel_index(idx, d1, d2, d3, i, j, k)
    integer, intent(in) :: idx, d1, d2, d3
    integer, intent(out) :: i, j, k
    integer :: temp

    temp = idx - 1
    i = mod(temp / (d2*d3), d1) + 1
    j = mod(temp / d3, d2) + 1
    k = mod(temp, d3) + 1
  end subroutine unravel_index

end program toy_model

!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------
! --       Block Spin Format    0 1 2 3 4 5    6 7 8 9  -----
!                                      transfer to           
! -- Interleaved Spin Format    0 2 4 1 3 5    6 8 7 9  -----
!------------------------------------------------------------
! Number of Spin Up Electrons = Number of Spin Down Electrons
!------------------------------------------------------------
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
