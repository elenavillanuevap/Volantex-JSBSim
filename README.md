# Volantex-JSBSim

JSBSim flight dynamics model of the Volantex Ranger 2400 UAV (ex2400_trim).

## Repository structure

Volantex-JSBSim/
├── aircraft/
│   └── ex2400_trim/
│       ├── ex2400_trim.xml         # Main JSBSim aircraft definition: geometry/metrics,
│       │                           # mass and inertia, ground reactions, propulsion links,
│       │                           # flight control system (pitch/roll/yaw/throttle channels),
│       │                           # and the full aerodynamic model (lift, drag, side force,
│       │                           # roll, pitch and yaw coefficients and derivatives)
│       ├── ex2400_trim_init.xml    # Initial conditions used to start the trim search
│       │                           # (airspeed, altitude, position, attitude)
│       ├── GBSap.xml               # Autopilot definition referenced by ex2400_trim.xml
│       └── GNCUtilities_TEST.xml   # System file (GNC utilities) referenced by ex2400_trim.xml
│        
├── engine/
│   ├── engine_ex2400.xml           # Engine definition used by the aircraft's propulsion block
│   └── prop_ex2400.xml             # Propeller/thruster definition used by the engine
└── scripts/
    ├── trim/
    │   └── ex2400_trim_trim.xml    # Trim run: finds the steady, level-flight equilibrium
    │   
    ├── reset/
    │   ├── ex2400_ailzero.xml      # Resets the aileron to zero deflection after trim
    │   ├── ex2400_elevzero.xml     # Resets the elevator to zero deflection after trim
    │   └── ex2400_rudzero.xml      # Resets the rudder to zero deflection after trim
    ├── doublet/
    │   ├── ex2400_trim_aildoublet2.xml   # Doublet excitation on the aileron
    │   ├── ex2400_trim_elevdoublet2.xml  # Doublet excitation on the elevator
    │   └── ex2400_trim_ruddoublet2.xml   # Doublet excitation on the rudder
    └── chirp/
        ├── chirp_ail_1.xml ... chirp_ail_4.xml    # 4 chirp excitation configurations for aileron
        ├── chirp_elev_1.xml ... chirp_elev_4.xml  # 4 chirp excitation configurations for elevator
        └── chirp_rud_1.xml ... chirp_rud_4.xml    # 4 chirp excitation configurations for rudder
        
