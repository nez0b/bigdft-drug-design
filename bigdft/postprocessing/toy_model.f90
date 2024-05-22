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
  type(ket)                    :: psi_it
  type(orbital_basis)          :: psi_ob
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
  real(wp), dimension(:),   pointer :: w
  real(wp), dimension(:),   pointer :: psi, psiv, psir, psir_i,psir_j
  real(dp), dimension(:),   pointer :: rhor,rho_ion
  real(wp), dimension(:),   allocatable :: pot_ion
  real(wp), dimension(:),   pointer :: potential
  real(wp), dimension(:)  , pointer :: hpsi_ptr
  real(gp), dimension(:,:), pointer :: rxyz
  real(gp), dimension(:,:), pointer :: rxyz_old
  real(wp), dimension(:,:), pointer :: ovrlp
  real(dp), dimension(:,:), pointer :: rho_p => null() !needs to be nullified
  real(wp), dimension(:,:), allocatable :: pot_tmp

  ! Interpolating Scaling Function (daub_to_ISF)
  real(wp), dimension(:,:), allocatable :: psirr
  ! Spin Polarized Format of psirr
  real(wp), dimension(:,:), allocatable :: phirr
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

  character*32 con, tmp
  integer :: ip ,iq ,ir ,is , ihpq, ihpqrs, istat, i,j,k, re_index1, re_index2
  integer :: ipt,iqt,irt,ist, nhpq, nhpqrs
  real(kind=8) :: hpq,hpqrs

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
  write(*,*) "========orbs==============="
  write(*,*) orbs%npsidim_orbs

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
  allocate( E_kin(      orbtot,orbtot) ) ; E_kin=0.d0
  allocate( E_local(    orbtot,orbtot) ) ; E_local=0.d0
  allocate( E_nonlocal( orbtot,orbtot) ) ; E_nonlocal=0.d0

  allocate(orbs%eval(orbs%norb*orbs%nkpts)) ; call f_zero(orbs%eval)
  allocate(psi(  max(orbs%npsidim_orbs, orbs%npsidim_comp)+1 )) ; psi=0._gp
  allocate(psiv( 2*inputs%norbv * (Lzd%Glr%wfd%nvctr_c + 7*Lzd%Glr%wfd%nvctr_f)+1 )) ; psiv=0._gp
  allocate(rxyz_old(3, atoms%astruct%nat)) ; rxyz_old=0._gp

  !------------------i--------------------------------------------------!
  ! Read occupied state wavefunctions from disk and store them in psi. !
  ! orbs%norb - number of occupied orbitals                            !
  !--------------------------------------------------------------------!
  call check_linear_and_create_Lzd(iproc,nproc,inputs%linear,Lzd,atoms,orbs,inputs%nspin,atoms%astruct%rxyz)
  call system('echo "========orbs====-3==========="')
  call readmywaves(iproc,trim(inputs%dir_output) // "wavefunction", &
  WF_FORMAT_PLAIN,orbs,lzd%glr,atoms,rxyz_old,atoms%astruct%rxyz,psi)
  call system('echo "========orbs====-2==========="')
  !--------------------------------------------------------------------!
  !--------------------------OLD　VERSION----------------------------!
  !call readmywaves(iproc, "data/wavefunction", WF_FORMAT_PLAIN, orbs, Lzd%Glr%d%n1, Lzd%Glr%d%n2, Lzd%Glr%d%n3, &
  !                 inputs%hx, inputs%hy, inputs%hz, atoms, rxyz_old, atoms%astruct%rxyz, Lzd%Glr%wfd, psi)

  !--------------------------NEW　VERSION   MOD by PJW-----------------!

 ! if (paw%usepaw) then
 !    call readmywaves(trim(in%dir_output) // "wavefunction", &
 !           orbs,lzd%glr,atoms,rxyz_old,rxyz,psi,pawrhoij=paw%pawrhoij)
     !call readmywaves( "data/wavefunction", orbs,lzd%glr,atoms,rxyz_old,rxyz,psi,pawrhoij=paw%pawrhoij)
!     call input_wf_disk_paw(iproc, nproc, atoms, GPU, Lzd, orbs, psi, denspot, nlpsp, paw)
!  else
     !call readmywaves(trim(in%dir_output) // "wavefunction",orbs,lzd%glr,atoms,rxyz_old,rxyz,psi)
     !write(*,*)'======================', iproc!WF_FORMAT_PLAIN == 1
     !call readmywaves(trim(in%dir_output) // "wavefunction", orbs,lzd%glr,atoms,rxyz_old,rxyz,psi)
     !call readmywaves("wavefunction",orbs,lzd%glr,atoms,rxyz_old,rxyz,psi)
     
    
 !    call readmywaves(iproc,"data/wavefunction",WF_FORMAT_PLAIN,orbs,lzd%glr,atoms,rxyz_old,atoms%astruct%rxyz,psi)
 ! end if
! call readmywaves(iproc, "data/wavefunction", WF_FORMAT_PLAIN, orbs, Lzd%Glr, &
!                  inputs%hx, inputs%hy, inputs%hz, atoms, rxyz_old, atoms%astruct%rxyz, Lzd%Glr%wfd, psi)
  !--------------------------------------------------------------------!
  if(nproc>1) call fmpi_allreduce(orbs%eval(1), orbs%norb*orbs%nkpts, op=FMPI_SUM)
! do i=1,orbs%norb ; write(1000+i,'("# orb: ",i4)') i ; write(1000+i,'(f20.12)') psi ; end do

  !--------------------------------------------------------------------!
  ! Read virtual  state wavefunctions from disk and store them in psi. !
  ! inputs%norbv - number of virtual orbitals                          !
  ! orbsv%norb - number of total virtual orbitals (2*inputs%norbv)     !
  !--------------------------------------------------------------------!
  nullify(orbsv%eval)
  orbsv%eval = f_malloc_ptr(orbsv%norb*orbsv%nkpts,id='orbsv%eval')
  nvirtu = abs(inputs%norbv) ; nvirtd = nvirtu
  call orbitals_descriptors(iproc, nproc, nvirtu+nvirtd, nvirtu, nvirtd, &
                            orbs%nspin, orbs%nspinor, orbs%nkpts, orbs%kpts, orbs%kwgts, orbsv, LINEAR_PARTITION_NONE)
  call check_linear_and_create_Lzd(iproc,nproc,inputs%linear,Lzd,atoms,orbsv,inputs%nspin,atoms%astruct%rxyz)
  !--------------------------OLD　VERSION----------------------------!
  ! call readmywaves(iproc, "data/virtuals",     WF_FORMAT_PLAIN, orbsv, Lzd%Glr%d%n1, Lzd%Glr%d%n2, Lzd%Glr%d%n3, &
  !                 inputs%hx, inputs%hy, inputs%hz, atoms, rxyz_old, atoms%astruct%rxyz, Lzd%Glr%wfd, psiv)
  !--------------------------NEW　VERSION----------------------------!
  call system('echo "========orbs====-1==========="')
  call readmywaves(iproc, "data/virtuals",WF_FORMAT_PLAIN, orbsv, Lzd%Glr, atoms, rxyz_old, atoms%astruct%rxyz,psiv)
  call system('echo "========orbs====0==========="')
  if(nproc>1) call fmpi_allreduce(orbsv%eval(1), orbsv%norb*orbsv%nkpts, op=FMPI_SUM)
! do i=1,mspin*inputs%norbv ; write(2000+i,'("# orb: ",i4)') i ;  write(2000+i,'(f20.12)') psiv ; end do

  orbdimocc = (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%norb
  orbdimvir = (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbsv%norb

  call system('echo "========orbs====1==========="')
  write(*,*) orbs%npsidim_orbs
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
    norb_start = orbs%norb+1
    norb_end   = orbs%norb+mspin*inputs%norbv
    allocate( tpsi_o(orbdimocc+1), indwvl_i(orbtot) ) ; tpsi_o = 0._dp
    allocate( tpsi_v(orbdimvir+1), indwvl_f(orbtot) ) ; tpsi_v = 0._dp
  
    !------------------------------------------------------
    ! loop on the localisation regions (occupied orbitals)
    !------------------------------------------------------
    allocate(tpsi(orbdimocc+1)) ; tpsi = 0._dp
    write(124,20) orbs%norb, orbs%nspinor, mspin, orbdimocc, Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
    !--------------------------OLD　VERSION----------------------------!
    !call initialize_work_arrays_locham(Lzd%nlr,Lzd%Llr,orbs%nspinor,.true.,wrk_lh)
    !--------------------------NEW　VERSION----------------------------!

    call initialize_work_arrays_locham(Lzd%nlr,Lzd%Llr,orbs%nspinor,wrk_lh)
    ekin = 0.d0 ; ekin_sum = 0.0_gp 
    loop_lr_kin: do ilr=1,Lzd%nlr
      dosome=.false.   ! check if this localisation region is used by one of the orbitals
      do iorb=1,orbs%norb ; dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr) ; if(dosome) exit ; end do
      if (.not. dosome) cycle loop_lr_kin
      !--------------------------OLD　VERSION----------------------------!
      !call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,.false.,wrk_lh)
      !--------------------------NEW　VERSION----------------------------!
      call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,wrk_lh)
      ispsi = 1
      loop_orbs: do iorb=1,orbs%norb
        indwvl_i(iorb) = ispsi
        ilr_orb = orbs%inwhichlocreg(iorb+orbs%isorb)
        ilr_orb = 1
        if (ilr_orb /= ilr) then
          ispsi = ispsi + (Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor
          cycle loop_orbs
        end if
!+++++++++++++error++++++++++++here+++++++++++++++++++
        !call psi_to_tpsi(Lzd%hgrids, orbs%kpts(1,orbs%iokpt(iorb)),orbs%nspinor,Lzd%Llr(ilr), psi(ispsi), wrk_lh, tpsi(ispsi), ekin)
        
  !++++++++error randomly happen here+++++++++++++++++++
        write(*,*)"debug----------bg",iorb,"/",orbs%norb
        call psi_to_tpsi(orbs%kpts(1,orbs%iokpt(iorb)),orbs%nspinor,Lzd%Llr(ilr), psi(ispsi), wrk_lh,tpsi(ispsi), ekin)
        write(*,*) "debug---------ed"

        ekin_sum = ekin_sum + orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
        ektmp1 = sum(psi(1:ispsi)*tpsi(1:ispsi))
        ispsi  = ispsi + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
        ektmp2 = sum(psi(1:ispsi)*tpsi(1:ispsi))  
        indwvl_f(iorb) = ispsi
        write(124,22) iorb, ekin, indwvl_i(iorb),indwvl_f(iorb)
      end do loop_orbs
    end do loop_lr_kin
    write(124,*) " occupied orbitals ... DONE"
    write(124,'("total kinetic energy = ",f15.8)') ekin_sum ; write(124,*)
    tpsi_o=tpsi
    call deallocate_work_arrays_locham(wrk_lh)
    deallocate(tpsi)
 
    !------------------------------------------------------
    ! loop on the localisation regions ( virtual orbitals)
    !------------------------------------------------------
    allocate(tpsi(orbdimvir)) ; tpsi = 0._dp
    write(124,21) inputs%norbv*mspin, orbs%nspinor, mspin, orbdimvir, Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
    !--------------------------OLD　VERSION----------------------------!
    !call initialize_work_arrays_locham(Lzd%nlr,Lzd%Llr,orbs%nspinor,.true.,wrk_lh)
    !--------------------------NEW　VERSION----------------------------! 
    call initialize_work_arrays_locham(Lzd%nlr,Lzd%Llr,orbs%nspinor,wrk_lh)
    ekin = 0.d0 ; ekin_sum = 0.0_gp ; tpsi=0._dp
    loop_lr_kinv: do ilr=1,Lzd%nlr
      dosome=.false.   ! check if this localisation region is used by one of the orbitals
      do iorb=1,inputs%norbv ; dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr) ; if(dosome) exit ; end do
      if (.not. dosome) cycle loop_lr_kinv
      !--------------------------OLD　VERSION----------------------------!
      !call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,.false.,wrk_lh)
      !--------------------------NEW　VERSION----------------------------!
      call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,wrk_lh)
      ispsi = 1
      loop_orbsv: do iorb=1,mspin*inputs%norbv
        indwvl_i(iorb+orbs%norb) = ispsi
        ! ilr_orb = orbs%inwhichlocreg(iorb+orbs%isorb)
        ilr_orb = 1
        if (ilr_orb /= ilr) then
          ispsi = ispsi + (Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor
          cycle loop_orbsv
        end if
        !--------------------------OLD　VERSION----------------------------!
        ! call psi_to_tpsi(Lzd%hgrids, orbs%kpts(1,orbs%iokpt(iorb)),orbs%nspinor,Lzd%Llr(ilr), psiv(ispsi), wrk_lh, tpsi(ispsi), ekin)
        !--------------------------NEW　VERSION----------------------------!
        call psi_to_tpsi(orbsv%kpts(1,orbsv%iokpt(iorb)),orbsv%nspinor, Lzd%Llr(ilr), psiv(ispsi), wrk_lh, tpsi(ispsi), ekin)
        ekin_sum = ekin_sum + ekin
        ektmp1 = sum(psiv(1:ispsi)*tpsi(1:ispsi))
        ispsi  = ispsi + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
        ektmp2 = sum(psiv(1:ispsi)*tpsi(1:ispsi))
        indwvl_f(iorb+orbs%norb) = ispsi
        write(124,22) iorb+orbs%norb, ekin, indwvl_i(iorb+orbs%norb),indwvl_f(iorb+orbs%norb)
      end do loop_orbsv
    end do loop_lr_kinv
    write(124,*) " virtual  orbitals ... DONE"
    write(124,'("total kinetic energy = ",f15.8)') ekin_sum ; write(124,*)
    write(*,*) "========orbs======2========="
    write(*,*) orbs%npsidim_orbs
    tpsi_v=tpsi
    call deallocate_work_arrays_locham(wrk_lh)
    deallocate(tpsi)
  20 format("   # occupied orbitals:",i3,", nspinor:",i3,", mspin: ",i4,", psi_dim:",3i9)
  21 format("   # virtual  orbitals:",i3,", nspinor:",i3,", mspin: ",i4,", psi_dim:",3i9)
  22 format(" orbital ",i3,"  ekin = ",f15.8,2x,"index of wavelet:",2i10)

    do i=1,   orbs%norb       ; do j=1,  orbs%norb
      E_kin(i,j) = sum( psi( indwvl_i(i):indwvl_f(i)-1) * tpsi_o(indwvl_i(j):indwvl_f(j)-1) ) 
    end do ; end do
    do i=1,   orbs%norb       ; do j=norb_start,norb_end
      E_kin(i,j) = sum( psi( indwvl_i(i):indwvl_f(i)-1) * tpsi_v(indwvl_i(j):indwvl_f(j)-1) )
    end do ; end do
    do i=norb_start,norb_end  ; do j=1,  orbs%norb
      E_kin(i,j) = sum( psiv(indwvl_i(i):indwvl_f(i)-1) * tpsi_o(indwvl_i(j):indwvl_f(j)-1) )
    end do ; end do
    do i=norb_start,norb_end  ; do j=norb_start,norb_end 
      E_kin(i,j) = sum( psiv(indwvl_i(i):indwvl_f(i)-1) * tpsi_v(indwvl_i(j):indwvl_f(j)-1) )
    end do ; end do
 
    deallocate(tpsi_o, tpsi_v, indwvl_i, indwvl_f)
    call system('echo "kinetic calculation ... DONE"')

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

  call system('echo "potential calculation (local)"')
  !---------------------------------------------------------------------
  !          energy of the    local potential, unit in Hartree         -
  !---------------------------------------------------------------------
  write(124,*) 
  write(124,*) "================================================"
  write(124,*) " energy of the LOCAL potential, unit in Hartree "
  write(124,*) "       -- psir(i)*potential*psir(j) --       "
  write(124,*) "================================================"
  
  allocate(psir(  Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3)))
  allocate(psir_i( Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3)))
  allocate(psir_j( Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3)))
  allocate(psirr(  Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3),0:orbs%norb+orbsv%norb-1));psirr=0._wp
 
 
  allocate(phirr(Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3) ,0:orbs%norb*cspin+orbsv%norb-1)) ; phirr=0._wp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
 ! call local_potential_dimensions(Lzd, orbs, dpcom%ngatherarr(0,1))
  dict => dict_new()

  call system('echo "03"')
  dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(domain_geocode(atoms%astruct%dom)),&
            alpha_bc=onehalf*pi,beta_ac=onehalf*pi,gamma_ab=onehalf*pi,&
            acell=(/Lzd%Glr%nboxi(2,1), Lzd%Glr%nboxi(2,2), Lzd%Glr%nboxi(2,3)/)*&
                  (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/))

  pkernel=pkernel_init(iproc,nproc,dict,&
      dom,(/Lzd%Glr%nboxi(2,1), Lzd%Glr%nboxi(2,2), Lzd%Glr%nboxi(2,3)/),&
       (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/))

  call system('echo "04"')
  call dict_free(dict)
  call pkernel_set(pkernel,verbose=.false.)
  !nullify(pot_ion) ;
  allocate(pot_ion(Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*dpcom%n3p)) ; pot_ion=0._gp
  nullify(rho_ion) ; allocate(rho_ion(Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*dpcom%n3p)) ; rho_ion=0._gp
  write(124,*) " call createIonicPotential"
  
  call system('echo "05"')
  call createIonicPotential(iproc, (iproc==0),atoms, atoms%astruct%rxyz,&
   inputs%elecfield, dpcom, pkernel, pot_ion, rho_ion, psoffset)
  call system('echo "06"')
  
  write(124,*) " call full_local_potential"
 ! allocate the potential in the full box
  call full_local_potential(iproc,nproc,orbs,Lzd,dpcom,xc,pot_ion,potential,pot_ref_count)

  write(124,*) " call initialize_work_arrays_sumrho"

  write(*,*) '????????????begin?????????????????'
  call initialize_work_arrays_sumrho(Lzd%Glr , .true., wisf)

  epot_sum = 0._dp
  do i=1,orbtot
    do j=1,orbtot
      if(i .le. orbs%norb) then
        call daub_to_isf(Lzd%Glr,wisf,  psi( (i-1)          *(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1) ,psir_i)
      else
        call daub_to_isf(Lzd%Glr,wisf, psiv( (i-1-orbs%norb)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1) ,psir_i)
      end if

      if(j .le. orbs%norb) then
        call daub_to_isf(Lzd%Glr,wisf,  psi( (j-1)          *(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1) ,psir_j)
        if(i .eq. 1) psirr(:,j-1) = psir_j(:)
      else
        call daub_to_isf(Lzd%Glr,wisf, psiv( (j-1-orbs%norb)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1) ,psir_j)
        if(i .eq. 1) psirr(:,j-1) = psir_j(:)
      end if

      do k=1, Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3)
         epot_sum = epot_sum + psir_i(k)*potential(k)*psir_j(k)
      end do
    ! write(124,"(f19.9)",advance = "NO") epot_sum  
      E_local(i,j) = epot_sum 
      epot_sum = 0.0_gp
    end do
  ! write(124,"(A)") " "
  end do ! write(124,*)
  write(*,*) '???????????end?????????????????????????????????'

  ! /////////////////////////////////////////////////////////////////////////////////////////////////
  ! /////////////////////////////////////////////////////////////////////////////////////////////////
  ! /////////////////////////////////////////////////////////////////////////////////////////////////
  ! Occupied Orbital
  ! Spin Orbital Notation

  write(*,*) "========orbs====3==========="
  write(*,*) orbs%npsidim_orbs
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
    ii=re_index1(orbs%norb*cspin+i-1,orbs%norb*cspin+orbsv%norb,orbs%norb*cspin)
    phirr(:,ii)=psir_i
  end do
  ! /////////////////////////////////////////////////////////////////////////////////////////////////

  !do i=0,orbs%norb+orbsv%norb-1  
  !  write(3000+i,*)   Lzd%Glr%d%n1i,   Lzd%Glr%d%n2i,   Lzd%Glr%d%n3i
  !  write(3000+i,*) inputs%hx/2._gp, inputs%hy/2._gp, inputs%hz/2._gp
  !  !write(3000+i,'(f15.8)') psirr(:,i)  
  !  write(3000+i,*) psirr(:,i)  
  !end do

  !do i=0,orbs%norb*cspin+orbsv%norb-1  
  !  write(4000+i,*)   Lzd%Glr%d%n1i,   Lzd%Glr%d%n2i,   Lzd%Glr%d%n3i
  !  write(4000+i,*) inputs%hx/2._gp, inputs%hy/2._gp, inputs%hz/2._gp
  !  write(4000+i,*) phirr(:,i)  
  !end do

  ! re_index1 (Block to Interleaved Spin Format)
  !do i=0,orbs%norb*cspin+orbsv%norb-1
  !  ii=re_index1(i,orbs%norb*cspin+orbsv%norb,orbs%norb*cspin)
  !  write(*,*) i, ii
  !end do

  ! Local Potential Energy
  !do i=0,orbs%norb*cspin+orbsv%norb-1 
  !  do j=0,orbs%norb*cspin+orbsv%norb-1 
  !    epot_sum=0.0_gp
  !    do k=1,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i
  !       epot_sum = epot_sum + phirr(k,i)*potential(k)*phirr(k,j)
  !    end do
  !    write(*,"(f19.9)",advance="NO"), epot_sum  
  !  end do
  !  write(*,*)
  !end do 

  ! /////////////////////////////////////////////////////////////////////////////////////////////////
  ! /////////////////////////////////////////////////////////////////////////////////////////////////
  ! /////////////////////////////////////////////////////////////////////////////////////////////////


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
    call createProjectorsArrays(Lzd%Glr, atoms%astruct%rxyz,atoms,orbs,&
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
    allocate( proj(4*MAXVAL(nlpsp%projs%mproj), 0:((orbs%norb*cspin+orbsv%norb)*num_proj)-1)) ;  proj=0._dp
    allocate(hproj(4*MAXVAL(nlpsp%projs%mproj), 0:((orbs%norb*cspin+orbsv%norb)*num_proj)-1)) ; hproj=0._dp


    ! nn=1:orbs%norb
    nn=1
    ! Number of Projector Index
    iproj=1
    ! Projector Index
    if(orbs%nspin .eq. 1) then
      nnproj=orbs%norb
    else if(orbs%nspin .eq. 2) then
      nnproj=orbs%norb*cspin
    end if

    call orbital_basis_associate(psi_ob, orbs=orbs, phis_wvl=psi, Lzd=Lzd, id='nonlocalham')
    if(associated(nlpsp%iagamma)) call f_zero(nlpsp%gamma_mmp)
    psi_it = orbital_basis_iterator(psi_ob)
    loop_kpt: do while(ket_next_kpt(psi_it))
      loop_lr: do while(ket_next_locreg(psi_it, ikpt=psi_it%ikpt))
        call DFT_PSP_projectors_iter_new(psp_it, nlpsp)
        loop_proj: do while (DFT_PSP_projectors_iter_next(psp_it, ilr=psi_it%ilr, lr=psi_it%lr, glr=Lzd%glr))
          call DFT_PSP_projectors_iter_ensure(psp_it, psi_it%kpoint, 0, nwarnings, Lzd%Glr)
          loop_psi_kpt: do while(ket_next(psi_it, ikpt=psi_it%ikpt, ilr=psi_it%ilr))          
            call system('echo "07"')
            call DFT_PSP_projectors_iter_apply(psp_it, psi_it, atoms, eproj, hpsi=hpsi)          
            call system('echo "08"')

            ! Number of Projector Considered for each Atomic Species
            nmproj=nlpsp%projs(iproj)%mproj
            if(nmproj .eq. 0) then
              iproj=iproj+1
              nmproj=nlpsp%projs(iproj)%mproj
            end if

            if(orbs%nspin .eq. 1) then
   
              ! Interleaved Spin Format 
              proj(1:nmproj,2*(nn-1)+0)=psp_it%parent%cproj(1:nmproj)
              proj(1:nmproj,2*(nn-1)+1)=psp_it%parent%cproj(1:nmproj)
              hproj(1:nmproj,2*(nn-1)+0)=psp_it%parent%hcproj(1:nmproj)
              hproj(1:nmproj,2*(nn-1)+1)=psp_it%parent%hcproj(1:nmproj)

              ! Next Projector 
              if(nn==nnproj) then
                nn=nn+inputs%norbv
                nnproj=nnproj+(orbs%norb+inputs%norbv)
                iproj=iproj+1
              end if

            else if(orbs%nspin .eq. 2) then

              ! Block Spin Format
              proj(1:nmproj,nn-1)=psp_it%parent%cproj(1:nmproj)
              hproj(1:nmproj,nn-1)=psp_it%parent%hcproj(1:nmproj)

              ! Next Projector 
              if(nn-1==nnproj-1) then
                nn=nn+orbsv%norb
                nnproj=nnproj+(orbs%norb*cspin+orbsv%norb)
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


    !-------------------
    !  virtual orbitals
    !-------------------
    !---- PJW: fixed hpsi/psi_it mapping mismatch-----comment out below line 816-817ish
    call orbital_basis_associate(psi_ob, orbs=orbsv, phis_wvl=psiv, Lzd=Lzd, id='nonlocalham')
    psi_it = orbital_basis_iterator(psi_ob)
    !-----------------------------------------------------------------
    allocate(hpsi(orbdimvir))
    hpsi_ptr => ob_ket_map(hpsi,psi_it) ;
    if(orbdimvir > 0) call f_zero(hpsi(1),orbdimvir)


    ! nn=1:orbsv%norb
    nn=1
    ! Number of Projector Index
    iproj=1
    ! Projector Index
    nnproj=orbs%norb*cspin+orbsv%norb

    write(*,*) "========orbs======4========="
    write(*,*) orbs%npsidim_orbs
    write(6,*)

    !call orbital_basis_associate(psi_ob, orbs=orbsv, phis_wvl=psiv, Lzd=Lzd, id='nonlocalham')
    !psi_it = orbital_basis_iterator(psi_ob)
    loop_kpt_v: do while(ket_next_kpt(psi_it))
      loop_lr_v: do while(ket_next_locreg(psi_it, ikpt=psi_it%ikpt))
        call DFT_PSP_projectors_iter_new(psp_it, nlpsp)
        loop_proj_v: do while (DFT_PSP_projectors_iter_next(psp_it, ilr=psi_it%ilr, lr=psi_it%lr, glr=Lzd%glr))
          call DFT_PSP_projectors_iter_ensure(psp_it, psi_it%kpoint, 0, nwarnings, Lzd%Glr)
          loop_psi_kpt_v: do while(ket_next(psi_it, ikpt=psi_it%ikpt, ilr=psi_it%ilr))          
            call system('echo "09"')
            call DFT_PSP_projectors_iter_apply(psp_it, psi_it, atoms, eproj, hpsi=hpsi)          
            call system('echo "10"')



            ! Number of Projector Considered for each Atomic Species
            nmproj=nlpsp%projs(iproj)%mproj
            if(nmproj .eq. 0) then
              iproj=iproj+1
              nmproj=nlpsp%projs(iproj)%mproj
            end if

            if(orbs%nspin .eq. 1) then
              ! Interleaved Spin Format
              ii=(orbs%norb*cspin+orbsv%norb)*(iproj-1) + &
                 re_index1(mod(orbs%norb*cspin+nn-1,orbs%norb*cspin+orbsv%norb),orbs%norb*cspin+orbsv%norb,orbs%norb*cspin)
            else if(orbs%nspin .eq. 2) then
              ! Block Spin Format
              ii=orbs%norb*cspin+nn-1
            end if

            proj(1:nmproj,ii)=psp_it%parent%cproj(1:nmproj)
            hproj(1:nmproj,ii)=psp_it%parent%hcproj(1:nmproj)

            ! Next Projector
            if(orbs%norb*cspin+nn-1==nnproj-1) then
              nn=nn+orbs%norb*cspin
              nnproj=nnproj+(orbs%norb*cspin+orbsv%norb)
              iproj=iproj+1
            end if

            ! Next Orbital
            nn=nn+1

          end do loop_psi_kpt_v
        end do loop_proj_v
      end do loop_lr_v
    end do loop_kpt_v
    call orbital_basis_release(psi_ob) ; deallocate(hpsi)

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
  if(orbs%nspin .eq. 1) then
    ! Number of Projector
    do k=1,num_proj
      ! Next Projector Index
      nn=(orbs%norb*cspin+orbsv%norb)*(k-1)
      do i=1,orbs%norb+inputs%norbv
        do j=1,orbs%norb+inputs%norbv
          E_nonlocal(i,j) = E_nonlocal(i,j) + sum(proj(:,(i-1)*2+nn)*hproj(:,(j-1)*2+nn))
        end do
      end do
    end do
  else if(orbs%nspin .eq. 2) then
    ! Number of Projector
    do k=1,num_proj
      ! Next Projector Index
      nn=(orbs%norb*cspin+orbsv%norb)*(k-1)
      do i=1,orbs%norb*cspin+orbsv%norb
        do j=1,orbs%norb*cspin+orbsv%norb
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

  deallocate(rxyz_old, psi, psir_i, psir_j)
  deallocate(proj,hproj)

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------
  write(*,*) "========orbs====5==========="
  write(*,*) orbs%npsidim_orbs
  !----------------
  !   data output 
  !----------------
  write(6,*) 
  write(6,'("number of occupied orbitals - orbs%norb  : ",5i4)') orbs%norb
  write(6,'("number of  virtual orbitals - orbsv%norb : ",5i4)') orbsv%norb
  write(6,'("    input parameter norbv - inputs%norbv : ",5i4)') inputs%norbv
  write(6,'("number of    total orbitals -     orbtot : ",5i4)') orbtot
  write(6,'("                            orbs%nspinor : ",5i4)') orbs%nspinor
  write(6,'("                                   mspin : ",5i4)') mspin
  write(6,'("                              OMP number : ",5i4)') OMP_get_thread_num(), OMP_get_num_threads(), OMP_get_max_threads()
  write(6,'("                          orbs%npsidim_orbs : ", i8)') orbs%npsidim_orbs
  write(6,'(" dimension of occupied orbitals (orbdimocc) : ", i8)') orbdimocc
  write(6,'(" dimension of  virtual orbitals (orbdimvir) : ", i8)') orbdimvir
  
  
  allocate( output1(0:orbtot-1,0:orbtot-1) ) ; output1=0._dp
  allocate( output2(0:orbtot-1,0:orbtot-1) ) ; output2=0._dp
  allocate( output3(0:orbtot-1,0:orbtot-1) ) ; output3=0._dp
  
  if(orbs%nspin .eq. 2) then
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
    do j=0,orbtot-1; do i=0,orbtot-1; write(6,1,advance="NO") output1(i,j); if(i .eq. orbs%norb) write(6,"(A)",advance="NO") "  " &
      ; end do
     write(6,"(A)") " " ; if(j .eq. orbs%norb) write(6,"(A)") "  " ; end do 

    write(6,*) ; write(6,*) " =====  local potential: psir(i)*potential(ij)*psir(j)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1; write(6,1,advance="NO") output2(i,j); if(i .eq. orbs%norb) write(6,"(A)",advance="NO") "  " &
      ; end do
      write(6,"(A)") " " ; if(j .eq. orbs%norb) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  nonlocal potential: psir(i)*hpsir(j)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1; write(6,1,advance="NO") output3(i,j); if(i .eq. orbs%norb) write(6,"(A)",advance="NO") "  " &
      ; end do
      write(6,"(A)") " " ; if(j .eq. orbs%norb) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  hpq:  sum(kinetic, V_local, V_nonlocal)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1; write(6,1,advance="NO") output1(i,j) + output2(i,j) + output3(i,j)  
      if(i .eq. orbs%norb) write(6,"(A)",advance="NO") "  "
      end do ; write(6,"(A)") " " ; if(j .eq. orbs%norb) write(6,"(A)") "  " ; end do
    
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

  write(*,*)

! ------------------------------------------------------------------------------------------------------------------------------------
  write(*,*)

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
      if(orbs%nspin .eq. 1) then
        ip=int(ip/2)
        iq=int(iq/2)
      else if(orbs%nspin .eq. 2) then
        ip=ip
        iq=iq
      end if

      hpq=output1(ip,iq)+output2(ip,iq)+output3(ip,iq)

      write(*,'(2i4,f19.12)') ipt,iqt,hpq
      write(05132021,'(2i4,f19.12)') ipt,iqt,hpq
    end do
    close(0513)
    close(05132021)

  write(*,*)

 !1 format(f16.12)
 1 format(f10.6)
 2 format(5x,2i3,3x,f20.16,6x,2i3)
  deallocate(E_local, E_nonlocal, E_kin)
  deallocate(output1, output2, output3)

!-------------------------------------------------------------------------------------------------------------------------------------

! SKChou 
! PSolver - Calculate hpqrs

!-------------------------------------------------------------------------------------------------------------------------------------
    !open(6,file="toy_model_1.out")
    ! Electrostat_Solver should be called under open( ,file);
    ! otherwise, there is a additional information aboule PSolver printed out on the screen

    open(05132021,file='hpqrs.out')
 
    open(0513,file="hpqrs.inp")
    nhpqrs=0
    do 
      read(0513,*,iostat=istat) tmp
      if(istat .ne. 0) exit
      nhpqrs=nhpqrs+1
    end do
    close(0513)
    
    allocate(rhopq(Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3)))
    allocate(rhors(Lzd%Glr%nboxi(2,1)*Lzd%Glr%nboxi(2,2)*Lzd%Glr%nboxi(2,3)))
    open(0513,file="hpqrs.inp")
    do ihpqrs=1,nhpqrs
      read(0513,*) ip,iq,ir,is

      ! Spin Orbital Index
      ! Interleaved Spin Format
      rhopq(:)=phirr(:,ip)*phirr(:,iq)
      rhors(:)=phirr(:,ir)*phirr(:,is)
      rhors(:)=rhors(:)/pkernel%mesh%volume_element

      write(*,*) "=================hererere===begin=========================="
      call Electrostatic_Solver(pkernel,rhors)
      write(*,*) "=================hererere===end=========================="
      hpqrs=sum(rhopq(:)*rhors(:))

      write( * ,'(4i4,f19.12)') ip,iq,ir,is, hpqrs
      write(05132021,'(4(i4),f19.12)') ip,iq,ir,is, hpqrs
    end do
    close(0513)
    close(05132021)
    deallocate(psirr,rhopq,rhors)
    deallocate(phirr)

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
