# GBS3_trimmed2

This set corresponds to:

- `aircraft/GBS3_trimmed2/GBS3_trimmed2.xml`
- `aircraft/GBS3_trimmed2/GBS3_trimmed2_init.xml`
- `scripts/GBS3_trimmed2_script.xml`
- `engine/Zenoah_G-26A_GBS3_trimmed2.xml`

## Changes with respect to the original

### Inertia and mass

- `Ixz`: `0.0951 -> 0.0951`
- `CG.x`: `14.4881 in -> 19.9020 in`

### Aerodynamics

Only these two drag terms change:

- `CDo`: `0.1000 -> 0.053876`
- `Drag_induced`: `0.0877 -> 0.014905`

### Propulsion and FCS

- Engine: `Zenoah_G-26A -> Zenoah_G-26A_GBS3_trimmed2`
- Power: `1789.68 W -> 121.70 W`
- The throttle channel 

### Init

`GBS3_trimmed2_init.xml`:

- `vt = 20 kt`
- `altitude = 1000 ft`
- `phi = -0.11 deg`
- `theta = 7.81 deg`
- `alpha = 7.81 deg`

### Script

`GBS3_trimmed2_script.xml`:

- total duration: `5 s`
- `dt = 0.01`
- keeps the seeded controls until `t < 1.0 s`
- performs an early stabilizing trim at `t = 0.05 s`
- performs the full trim at `t = 1.0 s`

Control:

- `aileron = -0.01`
- `rudder = 0.00`
- `pitch-trim = -0.04`
- `throttle = 1.000000`

## Summary

 The changes with respect to the original model are:

- `CG`
- `CDo`
- `Drag_induced`
- engine power
- throttle channel 
- trim-specific init and script

