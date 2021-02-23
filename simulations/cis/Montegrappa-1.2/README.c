

----- FORMAT --------------------

FILES

filename1.c
filename2.c
....


* filename1.c

  Brief Description
 
  func1
  func2
  func3
  ...


* filename2.c 

...
...


--------------------------------


MONTEGRAPPA ==================================================


FILES


adjust_st.c
do_mhistogram.c
energy.c
geometry.c
grappino.c
io.c
local_move.c 
mc.c
memory1.c
memory.c
misc.c
montegrappa.c
MPIfunc.c
optimizepot.c
pdb.c
peptides.c
potential.c
rotamers.c
stempering.c



* adjust_st.c

ManageRestart
AdjustSTempering
AddNewTemperature
OptimalWeights
EstimatedJumpProb
MaximizeTempProb
ReduceTemperatures
AddNewTemperature2
EstimatedJumpSingleProbability
CheckHistogramsAdjust
NormalizeHistograms
CheckHistograms
AddHistogramPile
CheckHistogramsAdjust2
ExtrapolateLinear
CheckHistogramsAdjust3
FilterHistograms
CutHistograms
CalculateG

* do_mhistogram.c

MHistogram
DensityOfStates
SolveEquations
PrintThermodynamics
CalculateThermodynamics
GMHequation
print_state
FatalError
FreeEnergy
SigmaEnergy
AverageEnergy
GMHequation_df
GMHequation_fdf

* energy.c 

GO_Paris
GO_Angles
Ram_Dihedrals
ContactMap
SetGoTypes
ReadTypes
DisulphideBonds
PrintOpGoFile

* geometry.c

Flip
MoveHead
MoveTail
Cartesian2Spherical
SPherical2Cartesian
RotateVector
PivotForward
PivotBackward
AddSidechain
Dist2
Angle
Dihedral
FlipFragment
CopyPolymer
CopyAllPolymers
DisplaceCom


*io.c

Error
ReadPolymer
SetLookbackTables
FindKeyword
Printpolymer
PrintPDBStream
ReadMcParms
ResetStuff

ReadParLLU
ReadParL
ReadParF
ReadParS
ReadParN

ReadPotential
PrintVector
PrintAngles
PrintStructure


* mc.c

Do_MC
MoveBackboneFlip
CopyResiduePositions
MoveBackbonePIvot
MoveMultiplePivot
Metropolis
UpdateMonomer
UpdateMonomerRange
MoveSideChain
SoftExit


*memory1.c

AlloDoubleMat
FreeDoubleMat
AlloInt
AlloDouble
AlloRestart
FreeRestart



*memory.c
AlloPolymer
InitTables

FreePolymer
FreeTables
AlloPotential
FreePotential
AlloIntMatrix
AlloDoubleMatrix
AlloInt
AlloDouble


* misc.c

irand
frand
FastSqrt
FastSin
FastCos
FastAcos
FastExp
Norm2
MatrixVectorProduct
Abs
DAbs
gasdev
ran1
gauss
CopyDoubleMatrix
CopyVector


* MPIFunc.c

Crete_vector_datatype
Create_rot_datatype
Create_side_datatype
Create_back_datatype
Crate_pot_datatype
Create_parms_datat
send_parms
send_pol
send_double_matrix
send_int_matrix
ExchangePol
Allo_op_input
send_op_input
send_op


* optimizepot.c

InitializeOptimizePot
AlloOptimizepot
ReadOPRestrains
FreeOpInput
OP_GetREstrain
OP_AddEnergy
OP_SamplePotential
OP_function
Chi2
OP_funcionDiff


* pdb.c

ReadPDB
PDB2CACB
PDB2NCAC
CreateBackboneFromPDB
CreateSidecainFromPDB
FindSidechain
Dist
SetRotarmersSimilarToPDB
DumbRMSD2
CopyPDB
SimplfyPDB
ReadPropensity

* peptides.c

MakePeptide

* potential.c
TotalEnergy
EnergyPair
EnergyMonomer
EnergyMonomerRange
AddShell
ResetContactsMonomer
GetEnergyMonomer
GetEnergyMonomerRange
PrintContacts
UpdateShell
EnergyAngles
EnergyDihedrals
PrintEnergies
PrintEnergies_Parallel
CopyPotential
CheckOverlaps
EnergyBox
EnergyBoxPolymer

*rotamers.c

AlloRotamer
ReadRotamers
Rot2PolymerScratch
AddDefaultCB
SubstituteDefaultCB
Rot2Polymer

* stempering.c
FreeStempering
STempering
DoSTempering
OrderTemperatures
PrintAverageEnergies
PrintHistogram
Restart
FindClosestT
PrintStatistics


* montegrappa.c

main
Welcome


GRAPPINO =================================================================


FILES

...cape...
grappino.c 
???
???

* grappino.c
...cape...
main
Parse
AppendPotentialComments
AlloAtoms
GrappinOWelcome

* ???.c 
 ...cape...


MHISTOGRAM ===============================================================
