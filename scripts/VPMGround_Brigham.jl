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

# ------------ GLOBAL VARIABLES ------------------------------------------------
# Default path where to save data
# extdrive_path = "/media/edoalvar/MyExtDrive/simulationdata7/"
extdrive_path = projectpath;
# extdrive_path = "temps/"

#Print What's Happening as we go
gt.verbalize("Defining parameters: ", v_lvl, verbose)

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

# --- Rotors --- #
# rotorspacing    = 2.4
if includerotors == true
    toprotorxpos    = -1.4
    bottomrotorxpos = toprotorxpos #-1.35

    #if we have inline rotors
    if inlinerotors == true
        if numrotors == 1
            n_rotors    = 1
            rotorpos1   = [toprotorxpos; 0.0; 0.0]
            rotorposs   = rotorpos1
            rotorccw    = [1]
        elseif numrotors == 2
            n_rotors = 2
            rotorpos1   = [toprotorxpos; 0.0; rotorspacing/2.0]
            rotorpos2   = [toprotorxpos; 0.0; -rotorspacing/2.0]
            rotorposs   = [rotorpos1 rotorpos2]
            if counterrotate == true
                rotorccw        = [1,2]
            else
                rotorccw        = [1,1]
            end
        else
            gt.verbalize("number specified rotors: $numrotors.  4 being defined...", v_lvl, verbose)
            n_rotors        = 4
            rotorpos1  = [toprotorxpos; 0.5*rotorspacing; rotorspacing/2.0]
            rotorpos2  = [toprotorxpos; 0.5*rotorspacing; -rotorspacing/2.0]
            rotorpos3  = [toprotorxpos; -0.5*rotorspacing; rotorspacing/2.0]
            rotorpos4  = [toprotorxpos; -0.5*rotorspacing; -rotorspacing/2.0]
            rotorposs  = [rotorpos1 rotorpos2 rotorpos3 rotorpos4]
            rotorccw        = [1,2,1,2]
        end
    #if we have offset rotors
    else
        n_rotors        = 8
        toprotorxpos    = -1.4
        bottomrotorxpos = -1.35
        rotorspacing    = 2.4
        # Position of rotors on main wing numbering starts on top right (lookingforward from windcraft) and reads counterclockwise
        rotorpos1 = [toprotorxpos; 1.5*rotorspacing; rotorspacing/2.0]
        rotorpos2 = [toprotorxpos; 0.5*rotorspacing; rotorspacing/2.0]
        rotorpos3 = [toprotorxpos; -0.5*rotorspacing; rotorspacing/2.0]
        rotorpos4 = [toprotorxpos; -1.5*rotorspacing; rotorspacing/2.0]
        rotorpos5 = [bottomrotorxpos; -1.5*rotorspacing; -rotorspacing/2.0]
        rotorpos6 = [bottomrotorxpos; -0.5*rotorspacing; -rotorspacing/2.0]
        rotorpos7 = [bottomrotorxpos; 0.5*rotorspacing; -rotorspacing/2.0]
        rotorpos8 = [bottomrotorxpos; 1.5*rotorspacing; -rotorspacing/2.0]
        rotorposs = [rotorpos1 rotorpos2 rotorpos3 rotorpos4 rotorpos5 rotorpos6 rotorpos7 rotorpos8]
        rotorccw  = [1,1,2,2,1,1,2,2]
    end


end #if we want rotors

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

#make dummy rotor in case there aren't any
rotors = vlm.Rotor[]
rotor_systems = ()

#if we want rotors, generate them and add them to the wing system here
if includerotors == true
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
        vlm.addwing(mainwingsystem, "rotor$i", rotor)
    end

    # --- Define rotor system --- #
    gt.verbalize("Creating Rotor System...", v_lvl, verbose)
    # Rotors grouped by systems of the same RPM
    rotor_systems = (rotors,)
end #if adding rotors

 # ------------ GEOMETRY OUTPUTS------------------------------------------------

#Rotate to be in correct starting pos.
O = zeros(3)
if circlepath == true
    Oaxis = [0.0 -1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]
else
    Oaxis = [1.0 0 0; 0 1 0; 0 0 1]
end
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
vehicle = uns.QVLMVehicle(   system;
# tilting_systems = tilting_systems,
rotor_systems   = rotor_systems,
# vlm_system      = vlm_system,
wake_system     = wake_system,
)


# Visualize maneuver
gt.verbalize("STEPPING THROUGH MANEUVER", v_lvl, verbose)

strn = uns.run_simulation(sim::Simulation,n_steps = 36*40;  #10 degrees for 40 rev
                          save_path = save_path,
                          run_name = run_name,
                          verbose = verbose,
                          optargs ...)  
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
#     Vinf(x,t) = [1,0,0]
#     vlm.setVinf(system, Vinf)
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