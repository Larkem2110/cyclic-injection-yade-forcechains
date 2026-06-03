#!/usr/bin/env yade
# -*- coding: utf-8 -*-

# ============================================================================
# Cyclic fluid injection into a Pre-Sheared gouge layer in YADE DEM
# ============================================================================
# Purpose
# -------
# This script continues from the final state of the direct shear simulation
# saved as ``laststate.yade``. It reloads the sheared granular fault gouge,
# stabilizes the normal stress, relaxes the shear stress to a prescribed fraction
# of the steady state shear stress, and then applies cyclic pore-pressure injection
# using PeriodicFlowEngine.
#
# Main stages
# -----------
# 1. Load the final sheared state from the direct shear test.
# 2. Reconfigure recorders and the time step for the injection simulation.
# 3. Maintain the target normal stress using a servo-controlled top wall.
# 4. Relax the shear stress to a prescribed fraction of its initial value.
# 5. Activate the fluid engine and impose cyclic pressure boundary conditions.
#
# Notes
# -----
# - This script assumes that the direct shear script already defined the YADE
#   engines and labels: ``flowEng``, ``saveSolid``, ``saveFluid``, ``recData``,
#   ``recMark``, ``servo``, and ``newton``.
# - The top and bottom walls are body IDs 0 and 1, respectively.
# - The normal stress is maintained at sigmaN during injection.
# - The shear stress is first relaxed to rTau × tau_ss.
# ============================================================================

from yade import pack, plot, export
from math import pi, copysign


# ============================================================================
# 1. Load previous direct-shear simulation state
# ============================================================================

prevSim = "laststate"
O.load(prevSim + ".yade")

# Target normal stress used by the normal servo-control [Pa].
# This value could also be computed automatically from the loaded state

sigmaN = 5.0e6

# Target shear stress ratio. The shear stress is relaxed to rTau × the initial
# shear stress measured after loading the direct-shear state (tau_ss).
rTau = 0.8

output = prevSim + "_Cyclic_injection"


# ============================================================================
# 2. Redefine time step and recorders
# ============================================================================

# ---------------------------------------------------------------------------
# Time stepper
# ---------------------------------------------------------------------------
# The automatic time stepper is disabled and a fixed time step is imposed. This
# is useful when a controlled time scale is needed for the cyclic injection.
timer = utils.typedEngine("GlobalStiffnessTimeStepper")
timer.active = False
timer.dead = True
O.dt = 0.5 * utils.PWaveTimeStep()

# ---------------------------------------------------------------------------
# Solid VTK recorder
# ---------------------------------------------------------------------------
saveSolid.dead = True
saveSolid.iterPeriod = int(1)
saveSolid.fileName = output + "."
saveSolid.skipNondynamic = 1
saveSolid.recorders = ["spheres", "colors", "velocity", "bstresses"]

# ---------------------------------------------------------------------------
# Fluid VTK recorder
# ---------------------------------------------------------------------------
def saveFlowVTK():
    """Save fluid-pressure and flow fields from the PeriodicFlowEngine."""
    flowEng.saveVtk(folder=output + "_vtkFluid")


saveFluid.dead = True
saveFluid.iterPeriod = int(1)


# ============================================================================
# 3. Data recorder
# ============================================================================

inputP = 0.0


def dataRecorder():
    """
    Record macroscopic mechanical and hydraulic quantities.

    Recorded quantities include wall stresses, contact stress tensor components,
    wall displacement, sample height, volume, porosity, contact index, imposed
    pore pressure, kinetic energy, and unbalanced force.
    """
    global inputP

    h = 0.0
    vol = 0.0
    vol_s = 0.0
    nb_s = 0.0

    h = O.bodies[0].state.pos[1] - O.bodies[1].state.pos[1]
    vol = h * O.cell.hSize[0, 0] * O.cell.hSize[2, 2]
    contactStress = getStress(vol)

    for o in O.bodies:
        if isinstance(o.shape, Sphere) and o.shape.color[0] != 1:
            nb_s += 1
            vol_s += 4.0 * pi / 3.0 * (o.shape.radius) ** 3

    n = 1.0 - vol_s / vol

    nbFrictCont = 0.0
    for i in O.interactions:
        if i.isReal and i.phys.cohesionBroken:
            nbFrictCont += 1

    plot.addData(
        iter=O.iter,
        stress_upWall0=abs(O.forces.f(0)[0] / (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])),
        stress_upWall1=abs(O.forces.f(0)[1] / (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])),
        stress_upWall2=abs(O.forces.f(0)[2] / (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])),
        contactStress00=contactStress[0, 0],
        contactStress01=contactStress[0, 1],
        contactStress02=contactStress[0, 2],
        contactStress10=contactStress[1, 0],
        contactStress11=contactStress[1, 1],
        contactStress12=contactStress[1, 2],
        contactStress20=contactStress[2, 0],
        contactStress21=contactStress[2, 1],
        contactStress22=contactStress[2, 2],
        xW=O.bodies[0].state.pos[0],
        height=h,
        volume=vol,
        porosity=n,
        k=2.0 * nbFrictCont / nb_s,
        p=inputP,
        Ek=kineticEnergy(),
        unbF=unbalancedForce(),
    )


# ============================================================================
# 4. Marker-particle recorder
# ============================================================================

markerSpheres = []
for o in O.bodies:
    if o.shape.color[2] == 0:
        markerSpheres.append(o)

print("nb or markers=", len(markerSpheres))


def markerRecorder():
    """Record marker-particle positions and velocities during injection."""
    global intOnFracPlane

    inFile = open(output + "_markerPosAndVel_" + str(O.iter), "a")
    for s in markerSpheres:
        inFile.write(
            str(s.state.pos[0]) + "\t"
            + str(s.state.pos[1]) + "\t"
            + str(s.state.pos[2]) + "\t"
            + str(s.state.vel[0]) + "\t"
            + str(s.state.vel[1]) + "\t"
            + str(s.state.vel[2]) + "\n"
        )
    inFile.close()


recMark.dead = True
recMark.iterPeriod = 1


# ============================================================================
# 5. Normal-stress servo-control of the top wall
# ============================================================================

# One step is performed with frozen bodies to recompute the wall forces in the
# loaded state. This ensures that the servo starts from an updated force state.
for o in O.bodies:
    o.dynamic = False

saveSolid.dead = False
recMark.dead = False
servo.dead = True
O.run(1, 1)
saveSolid.dead = True
recMark.dead = True
servo.dead = False
O.save(output + "_" + str(O.iter) + ".yade")

for o in O.bodies:
    if isinstance(o.shape, Sphere):
        o.dynamic = True

# Normal stiffness of the contacts involving the top wall.
nStiff = 0.0
for i in O.interactions.withBody(O.bodies[0].id):
    nStiff += i.phys.kn

print("normal stiffness=", nStiff)

initFnPlaten = O.forces.f(0)[1]
initSigmaN = initFnPlaten / (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])
print("normal stress (platen) =", initSigmaN)


def servo():
    """Maintain the target normal stress sigmaN on the top wall."""
    fnDesired = sigmaN * (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])

    nBoundaryVel = copysign(
        min(0.1, abs(0.35 * (O.forces.f(0)[1] - fnDesired) / nStiff / O.dt)),
        O.forces.f(0)[1] - fnDesired,
    )

    O.bodies[0].state.vel[1] = nBoundaryVel


# ============================================================================
# 6. Mechanical stabilization before shear-stress relaxation
# ============================================================================

print("stabilizing | iter=", O.iter)

newton.damping = 0.3
O.bodies[0].state.vel[0] = 0

recData.dead = False
recData.iterPeriod = 1

O.run(int(1e3), 1)

saveSolid.dead = False
recMark.dead = False
O.run(1, 1)
saveSolid.dead = True
recMark.dead = True

plot.saveDataTxt(output)
O.save(output + "_" + str(O.iter) + ".yade")

currentSN = O.forces.f(0)[1] / (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])
unbF = unbalancedForce()
print("\ncurrent Normal stress =", currentSN, " | unbF=", unbF)


# ============================================================================
# 7. Shear-stress control
# ============================================================================

# Shear stiffness of the contacts involving the top wall.
sStiff = 0.0
for i in O.interactions.withBody(O.bodies[0].id):
    sStiff += i.phys.ks

print("shear stiffness=", sStiff)

initFsPlaten = O.forces.f(0)[0]
initTau = initFsPlaten / (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])
print("shear stress (platen) =", initTau)

maxSBVel = 0.1
gain = 0.35


def servoShear():
    """
    Relax and then maintain the shear stress at rTau × the initial value (tau_ss).

    The velocity of the top wall in the x-direction is controlled using the
    difference between the current shear force and the desired shear force.
    """
    global maxSBVel

    fsDesired = rTau * initTau * (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])

    sBoundaryVel = copysign(
        min(maxSBVel, abs(gain * (O.forces.f(0)[0] - fsDesired) / sStiff / O.dt)),
        O.forces.f(0)[0] - fsDesired,
    )

    O.bodies[0].state.vel[0] = sBoundaryVel


# Insert the shear servo into the engine list.
O.engines = O.engines[:5] + [PyRunner(command="servoShear()", iterPeriod=1, label="servoShear")] + O.engines[5:]

print("shear control now | iter=", O.iter)

while True:
    O.run(100, 1)

    currentTau = O.forces.f(0)[0] / (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])
    targetTau = rTau * initTau

    if (unbalancedForce() < 0.001) and ((abs(abs(currentTau) - abs(targetTau)) / abs(targetTau)) < 0.001):
        print("stress state reached| iter=", O.iter, " | shear stress (platen) =", currentTau)
        print("stabilizing")

        O.run(int(1e3), 1)
        plot.saveDataTxt(output)

        saveSolid.dead = False
        recMark.dead = False
        O.run(1, 1)
        saveSolid.dead = True
        recMark.dead = True

        O.save(output + "_" + str(O.iter) + ".yade")
        break

currentSN = O.forces.f(0)[1] / (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])
unbF = unbalancedForce()
print("\ncurrent Normal stress =", currentSN, " | unbF=", unbF)

recMark.dead = True


# ============================================================================
# 8. Cyclic fluid injection
# ============================================================================

print("FLUID NOW | iter=", O.iter)

# Servo parameters are slightly relaxed during the fluid-injection phase
# to ensure the stability of the wall at post-failure stage.
maxSBVel = 0.15
gain = 0.15

# Activate fluid solver and impose pressure at the ymin boundary.
flowEng.isActivated = 1
flowEng.useSolver = 3
flowEng.bndCondIsPressure = [0, 0, 1, 0, 0, 0]
flowEng.bndCondValue = [0, 0, 0, 0, 0, 0]

# Cyclic pressure protocol, (from Pmin to Pmax and back to Pmin gives one cycle)
deltaP = 1e5       # pressure increment [Pa] = 0.1 MPa
inputP = deltaP
Pmax = 1.1e6       # maximum imposed pressure [Pa] = 1.1 MPa; this corresponds to
                   # the gouge critical pressure previously identified from a
                   # monotonic injection test.
Pmin = 1e5         # minimum imposed pressure [Pa] = 0.1 MPa
increasing = True

# Fluid properties.
# In YADE, a positive fluidBulkModulus uses the compressible-flow scheme,
# whereas a negative value uses the incompressible scheme.
flowEng.fluidBulkModulus = 2.2e9
flowEng.permeabilityFactor = 1
flowEng.viscosity = 1

newton.damping = 0.0

iterInit = O.iter
injectIter = O.iter
recIter = O.iter

# Number of DEM iterations between two pressure steps.
deltaIter = 3000

# Total duration of the cyclic injection stage.
iterMax = 500 * deltaIter

print(" Fixed time step O.dt = ", O.dt)

while True:
    O.run(1, 1)

    # -----------------------------------------------------------------------
    # Update the imposed injection pressure
    # -----------------------------------------------------------------------
    if O.iter >= int(injectIter + deltaIter):
        injectIter = O.iter

        if increasing:
            inputP += deltaP
            if inputP >= Pmax:
                inputP = Pmax
                increasing = False
        else:
            inputP -= deltaP
            if inputP <= Pmin:
                inputP = Pmin
                increasing = True
                # Set increasing = False here if only one full cycle is desired.

        flowEng.bndCondValue = [0.0, 0.0, inputP, 0.0, 0.0, 0.0]
        flowEng.updateBCs()

        print("updateBCs! inputP=", inputP, " | increasing=", increasing)

    # -----------------------------------------------------------------------
    # Save simulation data and VTK outputs
    # -----------------------------------------------------------------------
    if O.iter >= int(recIter + int(deltaIter / 2.0)):
        sn = O.forces.f(0)[1] / (O.cell.hSize[0, 0] * O.cell.hSize[2, 2])
        unbF = unbalancedForce()
        print("normal stress =", sn, " | unbF =", unbF)

        recIter = O.iter

        saveSolid.dead = False
        saveFluid.dead = False
        O.run(1, 1)
        saveSolid.dead = True
        saveFluid.dead = True

        plot.saveDataTxt(output)
        O.save(output + "_" + str(O.iter) + ".yade")
        print("saving data!")

    # -----------------------------------------------------------------------
    # End of cyclic injection simulation
    # -----------------------------------------------------------------------
    if O.iter >= int(iterInit + iterMax):
        print("iter=", O.iter, " -> END!")
        print(" Final time step O.dt = ", O.dt)
        plot.saveDataTxt(output)
        O.save(output + "_" + str(O.iter) + ".yade")
        break
