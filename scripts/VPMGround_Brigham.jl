###This is a copy/paste of the usefull parts of windcraft
#------------ Activate Enviornment -----------------------------------------
import Pkg 
projectpath = joinpath(@__DIR__, "..") * "/"
Pkg.activate(projectpath)
# ------------ MODULES ---------------------------------------------------------
# Load simulation engine
import FLOWUnsteady
# reload("FLOWUnsteady")
uns = FLOWUnsteady
vlm = uns.vlm
gt = uns.gt

using PyPlot
using Statistics

# ------------ GLOBAL VARIABLES ------------------------------------------------
# Default path where to save data
# extdrive_path = "/media/edoalvar/MyExtDrive/simulationdata7/"
extdrive_path = projectpath;
# extdrive_path = "temps/"

#Print What's Happening as we go
v_lvl = 0;
verbose = true;
gt.verbalize("Defining parameters: ", v_lvl, verbose)
save_folder = "202220224"
save_path = joinpath(projectpath,"data",save_folder)
run_name = "altitude1D_00"

# ------------ GEOMETRIC PARAMETERS (meters) -------------------------------
# ----- Main Wing (wing system)----- #
#Set Global Origin and Axis
vehicleorigin = [0.0; 0.0; 0.0]
vehicleaxis = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

xfoil = false;
includerotors = true;
inlinerotors = true;
numrotors = 1;
rotor_file = "DJI-II.csv";  #same as single rotor case.

rotorpitch       = 22;      #degrees
numbladeelements = 15 #todo how many?  #Check this in paper
data_path = uns.def_data_path;

# --- Rotors --- #
# rotorspacing    = 2.4
# if includerotors == true
    toprotorxpos    = -1.4
    bottomrotorxpos = toprotorxpos #-1.35

    #if we have inline rotors
    # if numrotors == 1
        n_rotors    = 1
        rotorpos1   = [toprotorxpos; 0.0; 0.0]
        rotorposs   = rotorpos1
        rotorccw    = [1]
    # elseif numrotors == 2
    #     n_rotors = 2
    #     rotorpos1   = [toprotorxpos; 0.0; rotorspacing/2.0]
    #     rotorpos2   = [toprotorxpos; 0.0; -rotorspacing/2.0]
    #     rotorposs   = [rotorpos1 rotorpos2]
    #     if counterrotate == true
    #         rotorccw        = [1,2]
    #     else
    #         rotorccw        = [1,1]
    #     end
    # else
    #     gt.verbalize("number specified rotors: $numrotors.  4 being defined...", v_lvl, verbose)
    #     n_rotors        = 4
    #     rotorpos1  = [toprotorxpos; 0.5*rotorspacing; rotorspacing/2.0]
    #     rotorpos2  = [toprotorxpos; 0.5*rotorspacing; -rotorspacing/2.0]
    #     rotorpos3  = [toprotorxpos; -0.5*rotorspacing; rotorspacing/2.0]
    #     rotorpos4  = [toprotorxpos; -0.5*rotorspacing; -rotorspacing/2.0]
    #     rotorposs  = [rotorpos1 rotorpos2 rotorpos3 rotorpos4]
    #     rotorccw        = [1,2,1,2]
    # end

# end #if we want rotors
"""
angle theta output in degrees
"""
function thetastar(tstar; theta0=0.0, thetan=360.0, omegan=[4.0/9.0;20.0/27.0;4.0/9.0], tn=[0.0;0.5;1.0])
    #find index of max tn less than tstar and set to index variable m.
    if tstar==tn[end]
        m = length(tn)
    elseif tstar==tn[1]
        m = 1
    else
        m = maximum(findall(x->x<tstar,tn))
    end
    # m = searchsortednearest(tn,tstar)
    #loop from 1 to m and sum areas.
    areas = 0.0
    den = 0.0
    for i = 2:length(tn)
        if i <= m
            areas += (tn[i]-tn[i-1])*(omegan[i]+omegan[i-1])/2.0
        end
        den += (tn[i]-tn[i-1])*(omegan[i]+omegan[i-1])/2.0
    end

    #add integral of final partial area up to tstar
    if m < length(tn) #if we haven't already integrated over everything
        integral = ( tstar - tn[m] ) * ( omegan[m] + (omegan[m+1]-omegan[m])/2.0 * (tstar-tn[m])/(tn[m+1]-tn[m]) )
    else
        integral = 0.0
    end

    #solve for scale factor, which
    sf = thetan / den
    #add theta0 (constant from integration)
    theta = sf * (areas + integral) + theta0

    return theta
end

function omegastar(tstar; h=1e-8, optargs...)
    return (thetastar(tstar+h; optargs...) - thetastar(tstar; optargs...))/h
end

"""
Transform time of thetastar and make it periodic to extend the angle output
outside the range [0, 360]
"""
function thetastar_periodic(t, nrevs; optargs...)
    tstar = (t*nrevs)%1    # Convert general time to time of one revolution
    return 360*floor(t*nrevs) + thetastar(tstar; optargs...)
end

function omegastar_periodic(tstar, nrevs; h=1e-8, optargs...)
    return (thetastar_periodic(tstar+h, nrevs; optargs...) - thetastar_periodic(tstar, nrevs; optargs...))/h
end

############################################################################
# ASSEMBLY: This block actually builds the geometry using the above parameters
############################################################################

gt.verbalize("Generating geometry...", v_lvl, verbose)

# --- Start Whole System Assembly --- #
# Create system assembly
gt.verbalize("Initializing System Assembly...", v_lvl, verbose)
system = vlm.WingSystem()

gt.verbalize("Initializing VLM System...", v_lvl, verbose)
# System to solve through the VLM solver
vlm_system = vlm.WingSystem()

#System to add rotor to
mainwingsystem = vlm.WingSystem()
_fun_fun(X,t) = [0.0, 0.0, -1e-12]
vlm.setVinf(mainwingsystem, Vinf_fun)

#make dummy rotor in case there aren't any
rotors = vlm.Rotor[]
rotor_systems = ()

#if we want rotors, generate them and add them to the wing system here
gt.verbalize("Generating Rotors...", v_lvl, verbose)
# --- Generate ROTORS --- #
R, B = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]]
##more rotor definition stuff:
n = 200.0 #target RPS
Vinf = 40.0 #freestream velocity
vind = sqrt( Vinf^2 + (n*0.84)^2 ) #velocity at 70% blade span
# Simulation parameters
J = Vinf/(n*2.0*R)                      # Advance ratio Vinf/(nD)
rho = 1.225                         # (kg/m^3) air density
mu = 1.81e-5                        # (kg/ms) air dynamic viscosity
ReD07 = rho*2.0*R*0.7*vind/mu            # Diameter-based Reynolds at 70% span
ReD = ReD07/0.7                     # Diameter-based Reynolds

# Generates base rotors (one on each rotation orientation)
props = vlm.Rotor[]
gt.verbalize("Generating first propeller...", v_lvl, verbose)
@time push!(props, uns.generate_rotor(rotor_file; pitch=rotorpitch,
n=numbladeelements, CW=true, ReD=ReD,
verbose=verbose, xfoil=xfoil,
data_path=data_path,
# plot_disc=plot_disc,
v_lvl=v_lvl+2))

gt.verbalize("Generating second propeller...", v_lvl, verbose)
@time push!(props, vlm.Rotor(!props[1].CW, props[1].r,
                        props[1].chord, props[1].theta,
                        props[1].LE_x, props[1].LE_z,
                        props[1].B, props[1].airfoils))
@time vlm.initialize(props[2], props[1].m)

# --- Assemble Rotors --- #
gt.verbalize("Adding Rotors to Main Wing...", v_lvl, verbose)
rotors = vlm.Rotor[]
for i = 1:n_rotors
    copy_prop = props[rotorccw[i]]
    this_prop = deepcopy(copy_prop) # Alternates rotation orientation
    this_O = rotorposs[:,i]
    vlm.setcoordsystem(this_prop, this_O, vehicleaxis; user=true)
    vlm.setVinf(this_prop, Vinf_fun)
    # Rotates props to be tip to tip #?what does this do?
    # vlm.rotate(this_prop, (-1)^(!CW_w) * init_ori_prop)

    # Adds the original polars that don't get copied in deepcopy
    this_prop.airfoils = copy_prop.airfoils
    this_prop._polars = copy_prop._polars
    this_prop._polarroot = copy_prop._polarroot
    this_prop._polartip = copy_prop._polartip

    push!(rotors, this_prop)
end

for (i, rotor) in enumerate(rotors)
    vlm.setVinf(rotor,Vinf_fun)
    rotor._wingsystem.Vinf = Vinf_fun
    vlm.addwing(mainwingsystem, "rotor$i", rotor)
end

println("Sherlock!\n\trotors[1]._wingsystem.Vinf = $(rotors[1]._wingsystem.Vinf)")

vlm.setVinf(mainwingsystem,Vinf_fun)

# --- Define rotor system --- #
gt.verbalize("Creating Rotor System...", v_lvl, verbose)
# Rotors grouped by systems of the same RPM
rotor_systems = (rotors,)

 # ------------ GEOMETRY OUTPUTS------------------------------------------------

#Rotate to be in correct starting pos.
O = zeros(3)
Oaxis = [1.0 0 0; 0 1 0; 0 0 1]
vlm.setcoordsystem(system, O, Oaxis)

gt.verbalize("Creating Wake System...", v_lvl, verbose)
# Wake-shedding system (`vlm_system`+`rotors`)
wake_system = vlm.WingSystem()
vlm.addwing(wake_system, "SolveVLM", vlm_system)
if includerotors == true
    for (i, rotor) in enumerate(rotors)
        vlm.addwing(wake_system, "Rotor$i", rotor)
    end
end

gt.verbalize("Completed Geometry Generation and Assembly...", v_lvl, verbose)

#Add geometry to vehicle system
vlm.setVinf(system,Vinf_fun)
vehicle = uns.QVLMVehicle(   system;
# tilting_systems = tilting_systems,
rotor_systems   = rotor_systems,
# vlm_system      = vlm_system,
wake_system     = wake_system,
)


function generate_maneuver_windcraft_kinematic(nrevs;
    disp_plot=false,
    includetail=false,
    includewing=false,
    includecontrols=false,
    includerotors=false,
    optargs...)

    omegamean = pi/180 * mean(omegastar.(range(0, 1, length=361); optargs...))

    ############################################################################
    # AIRCRAFT VELOCITY
    ############################################################################
    """
    Receives a nondimensional time between 0 and 1, and returns the
    vector of velocity of the vehicle at that instant.
    """
    function Vvehicle(t)

        theta = thetastar_periodic(t, nrevs; optargs...)*pi/180
        omega = omegastar_periodic(t, nrevs; optargs...)*pi/180
        scaling = omega/omegamean               # Scales velocity to a maximum value of 1

        Vcomp = [0.0, cos(theta), -sin(theta)]  # Counter-clockwise rotation
        return scaling*Vcomp
    end

    ############################################################################
    # AIRCRAFT ANGLES
    ############################################################################
    """
    Receives a nondimensional time between 0 and 1, and returns the angle
    (in degrees) of the vehicle.
    Returns: (angle_aircraft_x, angle_aircraft_y, angle_aircraft_z)
    """
    function anglevehicle(t)

        angle = [0.0, 0.0, thetastar_periodic(abs(t), nrevs; optargs...)]

        if t < 0.0
        return -angle
        else
        return angle
        end

    end
end


RPMref          = 1600;
t_per_rev       = 60/RPMref;      #60s / rotations per minute -> s/rotation
nrevs           = 40;
ttot            = nrevs*t_per_rev
t               = 0:t_per_rev:ttot;
maneuver = generate_maneuver_windcraft_kinematic(nrevs);
Vmean           = 2*pi*R/(nrevs*t_per_rev) # (m/s) mean velocity along a full circle
Vref            = Vmean    #determine from paper?
angle = 0;
# angle           = (anglevehicle(t),); #It either hates this or above, 321
# anglevehicle    = 0;
theta0          = 0.0
thetan          = 360.0
omegan          = [4.0/9.0, 20.0/27.0, 4.0/9.0]
tn              = [0.0, 0.5, 1.0]

my_RPM_function(t) = 1.0
RPM = (my_RPM_function,)
angle = ()

Vvehicle(t) = zeros(3)
angleofvehicle(t) = zeros(3)

maneuver = uns.KinematicManeuver(angle, RPM, Vvehicle, angleofvehicle)
#angle, RPM, velocity(t), 

revinit         = 0.25              # Part of revolution where to start the simulation
# tinit           = revinit*t_per_rev
tinit           = 0;
angle1 = maneuver.anglevehicle(tinit/ttot)
angle2 = maneuver.anglevehicle(tinit/ttot + 1e-12)
Winit = pi/180 * (angle2-angle1)/(ttot*1e-12)
Vinit = Vref*maneuver.Vvehicle(tinit/ttot)

# maneuver = generate_maneuver_windcraft_kinematic(nrevs;
#         disp_plot        = true,
#         includetail      = false,
#         includewing      = false,
#         includecontrols  = false,
#         includerotors    = true,
#         theta0           = 0,
#         thetan           = thetan,
#         omegan           = omegan,
#         tn               = tn)



simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
    Vinit=Vinit, Winit=Winit, t=tinit)

# Visualize maneuver
gt.verbalize("STEPPING THROUGH MANEUVER", v_lvl, verbose)

n_steps = 36*40;

strn = uns.run_simulation(simulation,n_steps;  #10 degrees for 40 rev
                          Vinf=Vinf_fun,
                          save_path = save_path,
                          run_name = run_name,
                          verbose = verbose,
                          sigma_vlm_solver=-1,
                          sigma_vlm_surf=R/numbladeelements*1.3,
                          sigma_rotor_surf=R/numbladeelements*1.3,
                          sigmafactor_vpm=1.3,
                          sigmafactor_vpmonvlm=1.0
                          )
                          #this instead of below, find necessary parameters
# strn = uns.visualize_kinematics(simulation, nsteps, save_path;
#                                 run_name=run_name,
#                                 save_vtk_optsargs=save_vtk_optsargs,
#                                 prompt=prompt, verbose=verbose, v_lvl=v_lvl,
#                                 paraview=false
# )

# Save ground
for (i, ground) in enumerate(grounds)
    gt.save(ground, run_name*"_Ground$i"; path=save_path)
    strn *= run_name*"_Ground$i.vtk;"
end

# Call paraview
if paraview
    run(`paraview --data="$save_path/$strn"`)
end

#--------------Visualize Geometry of rotor----------------------
# gt.verbalize("GEOMETRY VISUALIZATION", v_lvl, verbose)

#     # Generate the Geometry in order to visualize
#     (system, vlm_system, rotors,
#     tilting_systems, rotor_systems,
#     wake_system) = generate_geometry_windcraft(; xfoil=xfoil,
#                 data_path=uns.def_data_path,
#                 run_name=run_name,
#                 v_lvl=v_lvl+1,
#                 verbose=verbose, optargs...)

#     # Set up dummy values of RPM and freestream
#     _fun(x,t) = [1,0,0]
    vlm.set_fun(system, _fun)
#     for rotor in rotors
#     vlm.setRPM(rotor, 6000)
#     end

#     # Save visualization of geometry
#     gt.create_path(save_path, prompt)
#     strn = vlm.save(system, run_name; save_horseshoes=false, path=save_path)

#     # Call Paraview
#     if paraview
#         run(`paraview --data="$save_path/$strn"`)
#     end

