'''
                            SKChou,   March 2024

  Qiskit:

  Read hpq and hpqrs

  Active Space Selection

  Construct qubit Hamiltonian

  Calculate ground state energy
  
  240726: add RHF symmetric revision
  250120: add imaginary part

'''

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

import numpy as np
import BigQiskit

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# The Route for BigDFT Info and One- and Two-Electron Integrals
Data="."
# Molecule
Molecule="109_2ASP"
# Number of Occupied Orbital
norb=124
# Number of Virtual Orbital
norbv=10
# XC Energy Functional
ixc="PBE"
# Spin Polarization
nspin=1
# Initial Point
inip=0
# Number of Point
nopt=1
# Fisrt time read hpq&hpqrs or not
#FirstRead_HpqHpqrs=True
FirstRead_HpqHpqrs=False

# Active Space
'''count from 0'''
freeze_list_conv = [(i * 2, i * 2 + 1) for i in range(norb) if i not in {119,116,7}]#range(norb) if i not in {88, 97, 103, 108, 113}]
freeze_list = [x for tup in freeze_list_conv for x in tup] #flattened the list
print('freeze_list=')
print(freeze_list)

remove_list_conv = [((norb + i) * 2, (norb + i) * 2 + 1) for i in range(norbv) if i not in {3,9}]#range(norbv) if i not in {0, 1, 3, 4, 5}]
remove_list = [x for tup in remove_list_conv for x in tup]  #flattened the list
print('remove_list=')
print(remove_list)

''' count from 1
virtual: 9, 4, 8, 5, 10
wavefunction: 120, 124, 122, 117, 14
'''

Specify_Freeze_Remove_List=False
#Specify_Freeze_Remove_List=True

AS=False 
#AS=True ; (n,m)=(0,0)

# MP2 Active Space MO
MP2AS=False
#MP2AS=True; AS=False
#MP2 Threshold Pecentage
MP2thresholdPCT=14

# Exact Diag.
Exact=True
#Exact=False

# VQE
#VQE=False
VQE=True
#optimizer="L_BFGS_B"
optimizer="SLSQP"
maxeval=5e6
ansatz="UCCSD"
#ansatz="RealAmplitudes"
entanglement='linear'
reps=1
initial_point='zero'
disp=False
plot=False
info=False
# TwoLocal:
rotation_blocks="ry"
entanglement_blocks="cz"

# Save Hamiltonian Coeff. in .npy
SaveHam=True
# Save Qubit Hamiltonian
SaveQubitHam=True

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

print("BigDFT Qiskit")
print(f"Data Directory: {Data}")
print("Molecule: %s" % Molecule)
print("Number of Occupied  Orbitals: %d" % norb)
print("Number of Virtual   Orbitals: %d" % norbv)
# Number of Particle
# Closed System
norbp=norb*2
# Number of Total Spin Orbital
norbs=(norb+norbv)*2
print("Number of Molecule Spin Orbitals: %d" % norbs)
# Number of Occupied Spin Orbital
norbc=norb*2

print("DFT - XC Energy Functional: %s" % ixc)

if nspin==1:
  print("No Spin Polarization (Restricted)")
if nspin==2:
  print("Spin Polarization (Unrestricted)")

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Qubit Mapping Type 
map_type='parity'
#map_type='jordan_wigner'
print("Qubit Mapping Type:",map_type.upper())
# Chop Hamiltonian Coefficient
ChopHamCoeff=1e-06
print("Hamiltonian Coeff. Threshold: %.1e" % ChopHamCoeff)
# Chop Pauli Terms 
threshold=1e-05
print("Pauli Terms Threshold: %.1e" % threshold)

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
if Specify_Freeze_Remove_List == True:
  freeze_list=np.ndarray.tolist(np.load(Data+'/MP2List/FreezeList.npy'))
  remove_list=np.ndarray.tolist(np.load(Data+'/MP2List/RemoveList.npy'))


try:
  add_freeze_list
except NameError:
  add_freeze_list=[]
try:
  add_remove_list
except NameError:
  add_remove_list=[]

# Spin Orbital Index
# Interleaved Spin Format
if AS==True:
  AS_list=BigQiskit.ActiveSpaceList(norb,n,m)  
  freeze_list,remove_list=BigQiskit.AS_to_FZ_RM_list(norbc, norbs, AS_list)
  print("Active Space: CASCI[%d,%d]" % (n,m))
  print("remove_list:",remove_list)
  print("freeze_list:",freeze_list)

  if add_freeze_list != []:
    print("Additional freeze_list:",add_freeze_list)
    freeze_list+=add_freeze_list
  if add_remove_list != []:
    print("Additional remove_list:",add_remove_list)
    remove_list+=add_remove_list

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
fnAnsatz=f"{optimizer}_{ansatz}"
#if VQE==True and disp==False:
if VQE==True:
  print("VQE Optimizer: %s"%optimizer)
  print("VQE Ansatz: %s"%ansatz)
  if ansatz=="TwoLocal":
    print("Rotation Block:",rotation_blocks)
    print("Entanglement Block",entanglement_blocksa)
    fnAnsatz+=f"_{rotation_blocks}_{entanglement_blocks}"
  if ansatz!='UCCSD':
    print("Entanglement: %s"%entanglement)
    print("Repetition: %d"%reps)
    fnAnsatz+=f"_{entanglement}_r{reps}"

print(flush=True)
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# DFT
DFTE=[]
# Hartree Fock 
HFE=[]
# MP2
MP2E=[]
# Exact Diagonalization
ExactE=[]
# VQE
VQEE=[]
# VQE Parameter
VQEpara=[]

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Ground State
GSE=[]

# MP2 Info
MP2Info_MP2=[]
MP2Info_MP2PCT=[]
MP2Info_fzlist=[]
MP2Info_rmlist=[]

  # ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
if FirstRead_HpqHpqrs == True :
  # One Body Integral First Read
  import time
  start = time.time()
  fOBI = open(Data+'/OneBodyIntegral/'+Molecule+'.hpq','r')
  OneBodyIntegrals=BigQiskit.OneEhpqM(fOBI,norbs,ChopHamCoeff)
  fOBI.close()
  end = time.time()
  print("hpq read in seconds: ",(end-start),flush=True)

  # Two Body Integral First Read
  start = time.time()
  fTBI = open(Data+'/TwoBodyIntegral/'+Molecule+'.hpqrs','r')
  TwoBodyIntegrals=BigQiskit.TwoEhpqrsM(fTBI,norbs,ChopHamCoeff)
  fTBI.close()
  end = time.time()
  print("hpqrs read in seconds: ",(end-start),flush=True)
else :
  # One Body Integral Call .npy
  OneBodyIntegrals=np.load(Data+'/OneBodyIntegral/hpq.npy')

  # Two Body Integral Call .npy
  TwoBodyIntegrals=np.load(Data+'/TwoBodyIntegral/hpqrs.npy')

# Nuclear Repulsion Energy
fEion = open(Data+'/BigDFT_Energy/Eion/Eion','r')
shift = BigQiskit.Eion(fEion)
fEion.close()

# BigDFT DFT Energy
fDFTGSE = open(Data+'/BigDFT_Energy/GSE/GSE','r')
BigDFTGSE = BigQiskit.DFTGSE(fDFTGSE)
fDFTGSE.close()

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# The Qiskit convention is converted in BigQiskit.QubitH module
QiskitOneBodyIntegrals = OneBodyIntegrals
QiskitTwoBodyIntegrals = TwoBodyIntegrals   

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# MP2
# BigDFT log.yaml 
LogYaml=Data+"/"+"log-"+Molecule+".yaml"
# Orbital Energy
OE=BigQiskit.OrbEnergy(LogYaml,norb,norbv,nspin)

if MP2AS==True:
  # MP2: Active Space freeze_list & remove_list
  OccMP2,OccMP2PCT,VirMP2,VirMP2PCT,freeze_list,remove_list \
  =BigQiskit.MP2ASMO(norb,norbv,TwoBodyIntegrals,OE,MP2thresholdPCT)

  MP2=OccMP2+VirMP2
  MP2PCT=OccMP2PCT+VirMP2PCT
  MP2Info_MP2.append(MP2)
  MP2Info_MP2PCT.append(MP2PCT)
  MP2Info_fzlist.append(freeze_list)
  MP2Info_rmlist.append(remove_list)

MP2_delta=BigQiskit.MP2Delta(norb,norbv,TwoBodyIntegrals,OE,remove_list=freeze_list+remove_list)

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Qubit Hamiltonian
qubitOp,num_particles,freezeE=BigQiskit.QubitH(norbp,QiskitOneBodyIntegrals, QiskitTwoBodyIntegrals, map_type, threshold, freeze_list, remove_list)
core_energy=shift+freezeE

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Hartree-Fock Energy
HFEnergy = BigQiskit.HartreeFockE(norbp, OneBodyIntegrals, TwoBodyIntegrals)    
HFEnergy = HFEnergy + shift    

# MP2
MP2Energy = HFEnergy + MP2_delta

# Exact Diagonalization
ExactEnergy=0.0
if Exact==True:
  ExactEnergy, NumPyBackend = BigQiskit.ExactGSE(qubitOp)
  ExactEnergy = ExactEnergy+core_energy

# VQE
if VQE==True:
  # Read Initial Parameter
  fnOP=f"Param_{Molecule}_q{qubitOp.num_qubits}_{fnAnsatz}"
  #fnOP=fnOP+'_E'+str(i).zfill(2) now we only have one spot
  if type(initial_point)==str and initial_point.upper()=="READ":
    initial_point=np.loadtxt(f"{fnOP}.txt")      
  # VQE Convergence INFO
  if info==True:
    import logging
    fnINFO=fnOP.replace('Param','INFO')
    logging.basicConfig(filename=f'{fnINFO}.log',level=logging.INFO)
    logging.getLogger('qiskit.aqua.algorithms.minimum_eigen_solvers.vqe').setLevel(logging.INFO)
  import time
  start=time.time()
  VQEEnergy, VQEBackend = BigQiskit.HeuristicVQE(qubitOp,num_particles[0]+num_particles[1],map_type,optimizer,maxeval,ansatz,entanglement, \
                          reps,initial_point,disp,plot,rotation_blocks,entanglement_blocks)
  VQEEnergy+=core_energy
  end=time.time()
  print("VQE time in second: ",(end-start),flush=True)
  if info==True:
    with open(f'{fnINFO}.log','a') as fINFO:
        fINFO.write('\nCore Energy: %.12f'%core_energy)

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Electron Correlation Energy
eCorrelationE = ExactEnergy - HFEnergy    

DFTE.append(BigDFTGSE) 
HFE.append(HFEnergy)
MP2E.append(MP2Energy)
ExactE.append(ExactEnergy)
if VQE==True:
  VQEE.append(VQEEnergy)    

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if VQE==True:
  if Exact==True:
    print("DFT-%s=% .8f+%.8fj"%(ixc,BigDFTGSE.real,BigDFTGSE.imag),";", "HF=% .8f+%.8fj"%(HFEnergy.real,HFEnergy.imag), \
    "MP2=% .8f+%.8fj"%(MP2Energy.real,MP2Energy.imag), "Exact=% .8f+%.20fj"%(ExactEnergy.real,ExactEnergy.imag), \
    "CorrE=% .8f+%.8fj"%(eCorrelationE.real,eCorrelationE.imag), "VQE=% .8f+%.8fj"%(VQEEnergy.real,VQEEnergy.imag), \
    "VQE-Exact=%.8f+%.20fj"%(VQEEnergy.real-ExactEnergy.real,VQEEnergy.imag-ExactEnergy.imag))
  else:
    eCorrelationE = VQEEnergy - HFEnergy
    print("DFT-%s=% .8f+%.8fj"%(ixc,BigDFTGSE.real,BigDFTGSE.imag),";", "HF=% .8f+%.8fj"%(HFEnergy.real,HFEnergy.imag), \
    "MP2=% .8f+%.8fj"%(MP2Energy.real,MP2Energy.imag), "VQE=% .8f+%.8fj"%(VQEEnergy.real,VQEEnergy.imag), \
    "CorrE=% .8f+%.8fj"%(eCorrelationE.real,eCorrelationE.imag)) 
  np.savetxt(f"{fnOP}.txt",VQEBackend.optimal_point)
else:
  if Exact==True:
    print("DFT-%s=% .8f+%.8fj"%(ixc,BigDFTGSE.real,BigDFTGSE.imag),";", "HF=% .8f+%.8fj"%(HFEnergy.real,HFEnergy.imag), \
    "MP2=% .8f+%.8fj"%(MP2Energy.real,MP2Energy.imag), "Exact=% .8f+%.20fj"%(ExactEnergy.real,ExactEnergy.imag), \
    "CorrE=% .8f+%.8fj"%(eCorrelationE.real,eCorrelationE.imag))
  else:
    print("DFT-%s=% .8f+%.8fj"%(ixc,BigDFTGSE.real,BigDFTGSE.imag),";", "HF=% .8f+%.8fj"%(HFEnergy.real,HFEnergy.imag), \
    "MP2=% .8f+%.8fj"%(MP2Energy.real,MP2Energy.imag))

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Save Hamiltonian in .npy format
if SaveHam==True:
  """
  hpq=OneBodyIntegrals
  hpqrs=TwoBodyIntegrals
  fnhpq=f"hpq_finish" #hpq_E00 not sure what is it, name it finish temporarily
  fnhpqrs=f"hpqrs_finish"
  # .npy
  np.save(fnhpq,hpq)
  np.save(fnhpqrs,hpqrs)
  """
  # Save Energy
  with open('E_Eion.txt','a') as fEion:
      fEion.write("%.10f\r\n"%shift)
  with open('E_Frozen.txt','a') as fFrozen:
      fFrozen.write("%.10f\r\n"%freezeE)
  with open('E_Exact.txt','a') as fExact:
      fExact.write("%.10f+%.20fj\r\n"%(ExactEnergy.real,ExactEnergy.imag))
  # Save Bond Length
  """we now don't have different bond length
  fPosList=open(Data+"/"+Molecule+"_PosList",'r')
  PosList=fPosList.readlines()
  Distance=[]
  for item in PosList:
    Distance.append(float(item))
  fPosList.close()
  with open('PosList.txt','a') as fPos:
      fPos.write("%.10f\r\n"%Distance[i])"""


# Save Qubit Hamiltonian
if SaveQubitHam==True:
  print("\nQubit Hamiltonian\n")
  print("Number of Puali Terms: ",len(qubitOp))
  print(qubitOp)
  print()
  PaulisLabel=[]
  PaulisCoeff=[]
  for i in range(len(qubitOp)):
    PaulisCoeff.append(qubitOp.coeffs[i].real)
    PaulisLabel.append(qubitOp.paulis[i])
  QubitHam=dict(zip(PaulisLabel,PaulisCoeff))
  np.save(f'QubitHam_{Molecule}_q{qubitOp.num_qubits}_{map_type.upper()}_E%02d'%i, QubitHam)
  ##QubitHam=np.load(f"QubitHam_{Molecule}_q{qubitOp.num_qubits}_E%02d.npy"%i,allow_pickle='True')
  ##QubitHam=QubitHam.item() # numpy.ndarray to dictionary
  ##PauliLabel=list(QubitHam.keys())
  ##PauliCoeff=list(QubitHam.values())


# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Print MP2 Info
if MP2AS==True:
  import os
  MPfolder=os.path.exists(Data+'/MP2List')
  if not MPfolder:
    os.makedirs(Data+'/MP2List')
  fMP2Eng=open(Data+'/MP2List/MP2Eng','w')

  print("MP2 Active Space MO Energy Percentage")
  print("MP2 Threshold Percentage: %.2f%%" % MP2thresholdPCT,file=fMP2Eng)
  
  for i in range(nopt):#useless cause only one point of molecule distance
    for iOrb in range(norb):
      print("OccMO[%d] MP2: %.6f"%(iOrb,MP2Info_MP2[i][iOrb]),"(%5.2f %%)"%MP2Info_MP2PCT[i][iOrb],file=fMP2Eng)
    for iOrb in range(norbv):
      print("VirMO[%d] MP2: %.6f"%(iOrb,MP2Info_MP2[i][norb+iOrb]),"(%5.2f %%)"%MP2Info_MP2PCT[i][norb+iOrb],file=fMP2Eng)
    fMP2Eng.close()
    
    #save mp2 list
    print("MP2 freeze_list:",MP2Info_fzlist[i])
    np.save(Data+'/MP2List/FreezeList',MP2Info_fzlist)
    print("MP2 remove_list:",MP2Info_rmlist[i])
    np.save(Data+'/MP2List/RemoveList',MP2Info_rmlist)
    Norb=norb-len(MP2Info_fzlist[i])//2
    Norbv=norbv-len(MP2Info_rmlist[i])//2
    print("MP2 CASCI[%d,%d]-%s" % (Norb*2,Norb+Norbv,ixc))
    print()

print("complete the whole script ╮/(＞▽<)人(>▽＜)╭")
