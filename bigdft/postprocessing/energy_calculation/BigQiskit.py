#!/usr/bin/env python
# coding: utf-8

"""
							SKChou, March 2024

  OneEInt: Construct One Body Integral hpq

  TwoEInt: Construct Two Body Integral hpqrs

  Qiskit Electronic Integral Convention

  Exact Diagonalization

  VQE

  Active Space

  MP2 Active Space

"""


# # BigDFT - Qiskit Module
import numpy as np


# Read hpq and Construct hpq Matrix
def OneEhpqM(fOBIhpq, norbs, threshold):
  """
  Function: For first time read and truncate hpq and save it in .npy file
  Input: 
    hpq Index & Value
    Number of Spin Orbitals (occupied and Virtual)
  Output: 
    hpq Matrix (permutation symmetry restoration)
  """
  OBI=fOBIhpq.readlines()
  OneBodyIntegrals=np.zeros((norbs,norbs))
  for i in range(len(OBI)):
    p=int(OBI[i].split()[0])
    q=int(OBI[i].split()[1])
    hpq=float(OBI[i].split()[2])
    if abs(hpq) <= threshold:
              continue
    OneBodyIntegrals[p][q]=hpq
    # Symmetry
    OneBodyIntegrals[q][p]=hpq
  np.save('OneBodyIntegral/hpq',OneBodyIntegrals)
  return OneBodyIntegrals


def TwoEhpqrsM(fTBIhpqrs, norbs, threshold):
  """
  Function: For first time read and truncate hpqrs and save it in .npy file
  Input:
    hpqrs Index & Value
  Output:
    hpqrs Tensor (permutation symmetry restoration)
  """
  TBI=fTBIhpqrs.readlines()
  TwoBodyIntegrals=np.zeros((norbs,norbs,norbs,norbs))
  for i in range(len(TBI)):
    p=int(TBI[i].split()[0])
    q=int(TBI[i].split()[1])
    r=int(TBI[i].split()[2])
    s=int(TBI[i].split()[3])
    hpqrs=float(TBI[i].split()[4])
    if abs(hpqrs) <= threshold:
              continue

    TwoBodyIntegrals[p][q][r][s]=TwoBodyIntegrals[p][q][r+1][s+1]=TwoBodyIntegrals[p+1][q+1][r][s]=TwoBodyIntegrals[p+1][q+1][r+1][s+1]=hpqrs
    # Permutation Symmetry and turn RHF to UHF
    TwoBodyIntegrals[p][q][s][r]=TwoBodyIntegrals[p][q][s+1][r+1]=TwoBodyIntegrals[p+1][q+1][s][r]=TwoBodyIntegrals[p+1][q+1][s+1][r+1]=hpqrs
    TwoBodyIntegrals[q][p][r][s]=TwoBodyIntegrals[q][p][r+1][s+1]=TwoBodyIntegrals[q+1][p+1][r][s]=TwoBodyIntegrals[q+1][p+1][r+1][s+1]=hpqrs
    TwoBodyIntegrals[q][p][s][r]=TwoBodyIntegrals[q][p][s+1][r+1]=TwoBodyIntegrals[q+1][p+1][s][r]=TwoBodyIntegrals[q+1][p+1][s+1][r+1]=hpqrs
    TwoBodyIntegrals[r][s][p][q]=TwoBodyIntegrals[r][s][p+1][q+1]=TwoBodyIntegrals[r+1][s+1][p][q]=TwoBodyIntegrals[r+1][s+1][p+1][q+1]=hpqrs
    TwoBodyIntegrals[r][s][q][p]=TwoBodyIntegrals[r][s][q+1][p+1]=TwoBodyIntegrals[r+1][s+1][q][p]=TwoBodyIntegrals[r+1][s+1][q+1][p+1]=hpqrs
    TwoBodyIntegrals[s][r][p][q]=TwoBodyIntegrals[s][r][p+1][q+1]=TwoBodyIntegrals[s+1][r+1][p][q]=TwoBodyIntegrals[s+1][r+1][p+1][q+1]=hpqrs
    TwoBodyIntegrals[s][r][q][p]=TwoBodyIntegrals[s][r][q+1][p+1]=TwoBodyIntegrals[s+1][r+1][q][p]=TwoBodyIntegrals[s+1][r+1][q+1][p+1]=hpqrs
  np.save('TwoBodyIntegral/hpqrs',TwoBodyIntegrals)
  return TwoBodyIntegrals


# Nucleus-Nuclues Repulsion Energys
def Eion(fEion):
    EionList = fEion.readline().split()
    EionValue = []
    for item in EionList:
        EionValue.append(float(item))
    return float(EionValue[-1])

# BigDFT GSE
def DFTGSE(fGSE):
    GSEList = fGSE.readline().split()
    GSEValue = []
    for item in GSEList:
        GSEValue.append(float(item))
    return float(GSEValue[-1])

# BigDFT Energy
def BigDFTEnergy(fEnergy):
    EnergyList = fEnergy.readlines()
    EnergyValue = []
    for item in EnergyList:
        EnergyValue.append(float(item))
    return EnergyValue


from BigDFT.Logfiles import Logfile
def OrbEnergy(fLogYaml,norb,norbv,nspin):
  """
  Read BigDFT Orbital Energy

  Args:
    fLogYaml: filename of log.yaml
    norb:  Number of Occupied Orbital
    norbv: Number of Virtual  Orbital
    nspin: 1: No Spin 2: Spin Polarization
  Return:
    OE: Orbital Energy
  """

  log=Logfile(fLogYaml)
  rawOE=list(log.evals[0][0])

  # Orbital Energy
  OE=[]

  if nspin==1:
    for i in range(norb+norbv):
      OE.append(rawOE[i])

  if nspin==2:
    for i in range(norb+norbv):
      OE.append(rawOE[2*i])
      OE.append(rawOE[2*i+1])

  return OE

# Hartree-Fock Energy
def HartreeFockE(noc, OneBodyIntegral, TwoBodyIntegral):
    hpq=0.0
    for a in range(noc):
        hpq=hpq+OneBodyIntegral[a][a]
        
    hpqrs=0.0
    for a in range(noc):
        for b in range(noc):
            # Chemistry Notation
            # Coulomb Integral
            Jab=TwoBodyIntegral[a][a][b][b]
            # Exchange Integral
            Kab=TwoBodyIntegral[a][b][b][a]
            hpqrs=hpqrs+Jab-Kab
    hpqrs=hpqrs/2.0

    E0=hpq+hpqrs            
    return E0


from qiskit_nature.second_q.operators import PolynomialTensor, ElectronicIntegrals
def QiskitElectronicIntegrals(one_body_integrals, two_body_integrals):
    """
    Args:
      one_body_integrals: one-electron integrals
      two_body_integrals: two-electron integrals (hpqrs a_p^\dag a_r^\dag a_s a_q) (S8 Integral)

    Return:
      integral: Instance of Qiskit ElectronicIntegrals Class 
    """

    hpq=one_body_integrals
    hpqrs=two_body_integrals
    h1_a=hpq[np.ix_(*[range(0,i,2) for i in hpq.shape])] #the alpha-spin one-body coefficients
    h2_aa=hpqrs[np.ix_(*[range(0,i,2) for i in hpqrs.shape])] #the alpha-alpha-spin two-body coefficients
    h1_b=hpq[np.ix_(*[range(1,i,2) for i in hpq.shape])] #the beta-spin one-body coefficients
    h2_bb=hpqrs[np.ix_(*[range(1,i,2) for i in hpqrs.shape])] #the beta-beta-spin two-body coefficients
    h2_ba=hpqrs[np.ix_(*([range(1,i,2) for i in hpq.shape]+[range(0,i,2) for i in hpq.shape]))] #the beta-alpha-spin two-body coefficients
    h2_aa=np.einsum('ijkl->iklj', h2_aa)
    h2_bb=np.einsum('ijkl->iklj', h2_bb)
    h2_ba=np.einsum('ijkl->iklj', h2_ba)
    alpha = PolynomialTensor({"+-": h1_a, "++--": h2_aa})
    beta = PolynomialTensor({"+-": h1_b, "++--": h2_bb})
    beta_alpha = PolynomialTensor({"++--": h2_ba})

    #integrals = ElectronicIntegrals(alpha) #RHF
    integrals = ElectronicIntegrals(alpha,beta,beta_alpha)

    return integrals


from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit_nature.second_q.mappers import JordanWignerMapper,ParityMapper,BravyiKitaevMapper

def QubitH(num_particles, one_body_integrals, two_body_integrals, map_type='parity', threshold=1e-08, freeze_list=[], remove_list=[]):

    # Closed-Shell
    num_ele=(num_particles//2,num_particles//2)

    # Qiskit Electronic Integrals
    integrals = QiskitElectronicIntegrals(one_body_integrals, two_body_integrals)

    # Second-Quantized Electronic Hamiltonian
    EleHam = ElectronicEnergy(integrals,constants=None)

    # Active Space Transformation
    # spin orbital list to spatial orbital list
    freeze_list=[i//2 for i in freeze_list if i%2==0]
    remove_list=[i//2 for i in remove_list if i%2==0]
    nso=one_body_integrals.shape[0]
    nmo=nso//2
    total_num_electrons=num_ele[0]+num_ele[1]
    total_num_orbitals=nmo
    shift=0.0
    if len(freeze_list) != 0 or len(remove_list) !=0:
      total_hamiltonian=EleHam
      num_ele=(num_particles//2-len(freeze_list),num_particles//2-len(freeze_list))
      ASn=num_ele[0]+num_ele[1]
      ASm=nmo-len(freeze_list)-len(remove_list)
      active_orbitals_list=[i for i in list(range(nmo)) if i not in freeze_list+remove_list]
      transformer=ActiveSpaceTransformer(ASn, ASm, active_orbitals=active_orbitals_list)
      transformer.prepare_active_space(total_num_electrons,total_num_orbitals)
      reduced_hamiltonian=transformer.transform_hamiltonian(total_hamiltonian)
      shift=reduced_hamiltonian.constants['ActiveSpaceTransformer']
      EleHam=reduced_hamiltonian


    # Qubit Mapper
    if map_type.upper() == 'JORDAN_WIGNER':
        mapper = JordanWignerMapper()
    elif map_type.upper() == 'BRAVYI_KITAEV':
        mapper = BravyiKitaevMapper()
    elif map_type.upper() == 'PARITY':
        #mapper = ParityMapper() # Parity Mapping
        mapper = ParityMapper(num_ele) # Two Qubit Reduction

    # Qubit Hamiltonian
    qubit_op = mapper.map(EleHam.second_q_op())

    # Chops the real and imaginary parts of the operator coefficients
    qubit_op=qubit_op.chop(threshold)

    return qubit_op, num_ele, shift


from qiskit_algorithms import NumPyMinimumEigensolver
def ExactGSE(qubit_op):
    solver_NumPy = NumPyMinimumEigensolver()
    result_NumPy = solver_NumPy.compute_minimum_eigenvalue(qubit_op)
    eigenvalue=result_NumPy.eigenvalue
    #eigenstate=result_NumPy.eigenstate
    return eigenvalue, result_NumPy



from qiskit_algorithms import VQE
from qiskit.primitives import Estimator
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit.circuit.library import *
from qiskit_algorithms.optimizers import *

def HeuristicVQE(qubitOp,num_particles,map_type='parity',optimizer="SLSQP",maxeval=100000, \
                 ansatz="RealAmplitudes",entanglement='linear',reps=4,initial_point=None,disp=False,plot=False, \
                 rotation_blocks=None,entanglement_blocks=None,su2_gates=None,initial_state='HartreeFock'):

    """
    Qiskit Variational Quantum Eigensolver (VQE)

    Args:
        qubitOp: Qubit Operator
        map_type: Qubit Mappying Type
        optimizer: Optimizer
        maxeval: Maximum Functional Evaluation
        ansatz: Ansatz (Variational Form)
        entanglement: Entanglement Type
        reps: Repetition
        initial_point: Initial Parameters
        disp: Print Circuit Information
        plot: Plot VQE Convergence
        # TwoLocal Circuit Parameters:
        rotation_blocks: The gates used in the rotation layer.
        entanglement_blocks: The gates used in the entanglement layer.
        # EfficientSU2 Circuit Parameters:
        su2_gates: The SU(2) single qubit gates to apply in single qubit gate layers.

    Returns:
        E_vqe: Ground State Energy of VQE
        vqe_backend: The Result of VQE
    """


    # Optimizer
    if optimizer.upper()=="SLSQP":
        optimizer=SLSQP(maxiter=maxeval,disp=disp)
    elif optimizer.upper()=="L_BFGS_B":
        optimizer=L_BFGS_B(maxfun=maxeval,maxiter=15000)
    elif optimizer.upper()=="COBYLA":
        optimizer=COBYLA(maxiter=maxeval)
    elif optimizer.upper()=="NFT":
        optimizer=NFT(maxfev=maxeval,disp=disp)
    # Global Optimizer
    elif optimizer.upper()=="CRS":
        optimizer=CRS(max_evals=maxeval)
    elif optimizer.upper()=="DIRECT_L":
        optimizer=DIRECT_L(max_evals=maxeval)
    elif optimizer.upper()=="DIRECT_L_RAND":
        optimizer=DIRECT_L_RAND(max_evals=maxeval)
    elif optimizer.upper()=="ESCH":
        optimizer=ESCH(max_evals=maxeval)
    elif optimizer.upper()=="ISRES":
        optimizer=ISRES(max_evals=maxeval)
    else:
        raise RuntimeError('Optimizer Unavailable')

    # Closed-Shell
    num_particles=(num_particles//2,num_particles//2)

    num_qubits=qubitOp.num_qubits
    num_spatial_orbitals=num_qubits//2

    # Qubit Mapper
    if map_type.upper() == 'JORDAN_WIGNER':
        mapper = JordanWignerMapper()
    elif map_type.upper() == 'BRAVYI_KITAEV':
        mapper = BravyiKitaevMapper()
    elif map_type.upper() == 'PARITY':
        #mapper = ParityMapper() # Parity Mapping
        mapper = ParityMapper(num_particles) # Two Qubit Reduction

    if map_type.upper()=='PARITY':
        # default
        # restore to the original system
        qubit_reduction=True
        #num_spin_orbitals=num_qubits+2
        num_spatial_orbitals=num_spatial_orbitals+1
    else:
        qubit_reduction=False
        #num_spin_orbitals=num_qubits
        num_spatial_orbitals=num_spatial_orbitals

    if initial_state.upper() == 'HARTREEFOCK' or 'HARTREE_FOCK' or 'HF':
        initial_state=HartreeFock(num_spatial_orbitals,num_particles,mapper)
    else:
        initial_state=None

    # Ansatz (Variational Form)
    if ansatz.upper()=="UCCSD":
        var_form=UCCSD(num_spatial_orbitals=num_spatial_orbitals,num_particles=num_particles,qubit_mapper=mapper,initial_state=initial_state)
    elif ansatz.upper()=="TWOLOCAL":
        var_form=TwoLocal(num_qubits,rotation_blocks=rotation_blocks,entanglement_blocks=entanglement_blocks,
                          entanglement=entanglement,reps=reps,initial_state=initial_state)
    elif ansatz.upper()=="EFFICIENTSU2":
        var_form=EfficientSU2(num_qubits,su2_gates=su2_gates,entanglement=entanglement,reps=reps,initial_state=initial_state)
    elif ansatz.upper()=="EXCITATIONPRESERVING":
        var_form=ExcitationPreserving(num_qubits,mode='iswap',entanglement=entanglement,reps=reps,initial_state=initial_state)
    elif ansatz.upper()=="REALAMPLITUDES":
        var_form=RealAmplitudes(num_qubits,entanglement=entanglement,reps=reps,initial_state=initial_state)
    else:
        raise RuntimeError('Ansatz (Variational Form) Unavailable')


    # Initial Point
    if isinstance(initial_point,(list,np.ndarray)):
        pass
    elif isinstance(initial_point,str):
        if initial_point.upper()=="NONE":
            initial_point=None
        elif initial_point.upper()=="ZERO":
            initial_point=np.zeros(var_form.num_parameters)
        elif initial_point.upper()=="RANDOM": #[0,1)
            initial_point=np.random.random(var_form.num_parameters)
    else: # isinstance(initial_point, type(None))
        initial_point=None # Qiskit VQE Default: [-2pi,2pi)

    # callback
    Count_List = []
    Value_List = []
    def store_intermediate_result(count, parameters, mean, std):
        Count_List.append(count)
        Value_List.append(mean)

    vqe_solver = VQE(Estimator(), var_form, optimizer)
    vqe_solver.initial_point = initial_point
    vqe_backend = vqe_solver.compute_minimum_eigenvalue(qubitOp)
    E_vqe=vqe_backend.eigenvalue

   
    # Circuit Information
    Ansatz=type(var_form).__name__
    opt_cir=vqe_backend.optimal_circuit
    cir_depth=opt_cir.decompose().depth()
    if Ansatz=='UCCSD':
        cir_depth=opt_cir.decompose().depth() # NOT the "depth"!
    num_para=var_form.num_parameters # len(vqe.optimal_params)
    fun_eval=vqe_backend.cost_function_evals
    opt_time=vqe_backend.optimizer_time
    if disp==True:
        print("Number of Qubit: %d"%num_qubits)
        print("Optimizer: %s"%type(optimizer).__name__)
        print("Ansatz: %s"%Ansatz)
        if Ansatz=="TwoLocal":
            print("Rotation Block: ",rotation_blocks)
            print("Entanglement Block: ",entanglement_blocks)
        if Ansatz!='UCCSD':
            print("Entanglement: %s"%var_form.entanglement)
            print("Repetition: %d"%var_form.reps)
        print("Circuit Depth: %d"%cir_depth)
        print("Gate Count:",opt_cir.decompose().count_ops())
        print("Number of Parameter: %d"%num_para)
        print("Function Evaluation: %d"%fun_eval)
        print(f"Optimizer Time: {opt_time/60:.2f} min")
    if plot==True:
        import pylab
        fname=f"{type(optimizer).__name__}_{Ansatz}_q{num_qubits}"
        if Ansatz!='UCCSD':
            fname=f"{fname}_{var_form.entanglement}_r{var_form.reps}"
        pylab.plot(Count_List,Value_List,label=fname)
        pylab.xlabel('Evaluation Count')
        pylab.ylabel('Energy (Hartree)')
        pylab.title('VQE Energy Convergence')
        pylab.legend()
        pylab.tight_layout()
        pylab.savefig(f"{fname}_d{cir_depth}_p{num_para}_e{fun_eval}")

    return E_vqe, vqe_backend


#######################################################################


def ActiveSpaceList(norb,n,m):
  """
    Args:
      norb: Number of Occupied Orbital
      n: Number of Active Electron
      m: Number of Active Orbital

    Return:
      Active Space List (0-based Interleaved Spin Orbital Index)
  """

  # Closed System
  if n%2==0 :
    # Active Space Number of Occupied Orbital
    ASnorb=n//2
    # Active Space Number of Virtual  Orbital
    ASnorbv=m-ASnorb
  # Open System
  ##   Restricted Open System (nspin=1)
  ## Unrestricted Open System (nspin=2)
  else:
    sys.exit("Unavailable Open System!\n")

  # Number of Active Space Spin Orbital Range Index
  iAS=(norb-ASnorb)*2
  fAS=(norb+ASnorbv)*2

  AS_list=[]
  for i in range(iAS,fAS):
    AS_list.append(i)

  return AS_list


def AS_to_FZ_RM_list(norbc, norbs, AS_list):
  """
    Args:
      norbc: Number of Occupied Spin Orbital
      norbs: Number of Occupied and Virtual Spin Orbital
      ASList: Active Space List (0-based Interleaved Spin Orbital Index)

    Returns:
      freeze_list (0-based Interleaved Spin Orbital Index)
      remove_list (0-based Interleaved Spin Orbital Index)
  """

  # remove_list [ virtual orbital ]
  remove_list=[]
  for i in range(norbc+1,norbs):
    if i not in AS_list:
      remove_list.append(i)

  # freeze_list [ occupied orbital ]
  freeze_list=[]
  for i in range(norbc):
    if i not in AS_list:
      freeze_list.append(i)

  return freeze_list, remove_list


def MP2Delta(norb,norbv,TBI,orbital_energies,remove_list=[],orbital_list=[]):
    """
    Args:
        norb: Number of Occupied Orbital
        norbv: Number of Virtual Orbital
        TBI: Two Body Integral (Chemist Notation; Interleaved Spin Format)
        orbital_energies: Orbital Energies[norb]
        remove_list: Remove Occupied and Virtual Spin Orbital Index List (Interleaved Spin Format)
        orbital_list: Occupied and Virtual Spin Orbital Index List (Interleaved Spin Format)
    Return:
        MP2 Delta Energy
    """

    # Number of Occupied Spin Orbital
    nsorb=norb*2
    # Number of Virtual  Spin Orbital
    nsorbv=norbv*2
    # Number of Spin Orbital
    num_so=nsorb+nsorbv

    # Orbital Energy in Spin Orbital Basis
    o_e=np.zeros(num_so)
    for i in range(len(orbital_energies)):
        o_e[2*i+0]=orbital_energies[i]
        o_e[2*i+1]=orbital_energies[i]

    # Chemist to Physicist Notation
    TBI=np.einsum('ijkm->ikmj', TBI)

    occorb_list=[]
    virorb_list=[]

    if len(orbital_list)>0 :
        #if orbital_list is not empty
        for iOrb in orbital_list:
            if iOrb < nsorb:
                occorb_list.append(iOrb)
            else:
                virorb_list.append(iOrb)

    mp2_delta = 0.0
    
    # Occupied Spin Orbital Index
    for i in range(nsorb):
        for j in range(nsorb):
            if j<=i:
                continue

            # Remove Occupied Orbital
            if len(remove_list)>0 :
                #if remove_list is not empty
                if i in remove_list or j in remove_list:
                    continue

            # Consider Occupied Orbital
            if len(occorb_list)>0 :
                #if occorb_list is not empty
                if i in occorb_list and j in occorb_list:
                    pass
                elif i in occorb_list and j not in occorb_list:
                    pass
                elif i not in occorb_list and j in occorb_list:
                    pass
                elif i not in occorb_list and j not in occorb_list:
                    continue

            # Virtual Spin Orbital Index
            for a in range(nsorb,nsorb+nsorbv):
                for b in range(nsorb,nsorb+nsorbv):
                    if b<=a:
                        continue

                    # Remove Virtual Orbital
                    if len(remove_list)>0 :
                        #if remove_list is not empty
                        if a in remove_list or b in remove_list:
                            continue

                    # Consider Virtual Orbital
                    if len(virorb_list) :
                        #if virorb_list is not empty
                        if a in virorb_list and b in virorb_list:
                            pass
                        elif a in virorb_list and b not in virorb_list:
                            pass
                        elif a not in virorb_list and b in virorb_list:
                            pass
                        elif a not in virorb_list and b not in virorb_list:
                            continue

                    # MP2 Formula (Physicist Notation)
                    num = (TBI[i][j][a][b]-TBI[i][j][b][a])**2
                    denom = o_e[i] + o_e[j] - o_e[a] - o_e[b]
                    mp2_delta += num/denom

                    #if abs(num/denom) > 1e-6:
                    #    print("%2d_%2d_%2d_%2d: %.6e" % (i,j,a,b, num/denom))

    return mp2_delta


def MP2ASMO(norb,norbv,TBI,OE,thresholdPCT=40.0):
    """
    MP2 Energy Choose Active Space MO

    Args:
      norb:         Number of Occupied Orbital
      norbv:        Number of Virtual  Orbital
      TBI:          Two Body Integral (Intverleaved Spin Format)
      OE:           Orbital Energy (Intverleaved Spin Format)
      thresholdPCT: Active Space MO Threshold in Percentage

    Returns:
      OccMP2:       Occupied Orbital MP2 Correlation Energy
      OccMP2PCT:    Occupied Orbital MP2 Correlation Percentage
      VirMP2:       Virtual  Orbital MP2 Correlation Energy
      VirMP2PCT:    Virtual  Orbital MP2 Correlation Percentage
      freeze_list:  MP2 Frozen Occupied Spin Orbital Index List (Interleaved Spin Format)
      remove_list:  MP2 Remove Virtual  Spin Orbital Index List (Interleaved Spin Format)
    """

    MP2Correlation=MP2Delta(norb,norbv,TBI,OE)

    OccMP2=[]
    OccMP2PCT=[]
    for i in range(norb):
      OccIdx=i*2
      orbital_list=[OccIdx,OccIdx+1]
      OccMP2Delta=MP2Delta(norb,norbv,TBI,OE,orbital_list=orbital_list)
      OccMP2.append(OccMP2Delta)
      OccMP2DeltaPCT=(OccMP2Delta)*100/MP2Correlation
      OccMP2PCT.append(OccMP2DeltaPCT)

    VirMP2=[]
    VirMP2PCT=[]
    for i in range(norbv):
      VirIdx=(norb+i)*2
      orbital_list=[VirIdx,VirIdx+1]
      VirMP2Delta=MP2Delta(norb,norbv,TBI,OE,orbital_list=orbital_list)
      VirMP2.append(VirMP2Delta)
      VirMP2DeltaPCT=(VirMP2Delta)*100/MP2Correlation
      VirMP2PCT.append(VirMP2DeltaPCT)

    freeze_list=[]
    for i in range(norb):
      if OccMP2PCT[i] < thresholdPCT:
        freeze_list.append(2*i+0)
        freeze_list.append(2*i+1)

    remove_list=[]
    for i in range(norbv):
      if VirMP2PCT[i] < thresholdPCT:
        remove_list.append(2*(norb+i)+0)
        remove_list.append(2*(norb+i)+1)

    return OccMP2, OccMP2PCT, VirMP2, VirMP2PCT, freeze_list, remove_list


