# GBS3_trimmed3

This set corresponds to:

- `aircraft/GBS3_trimmed3/GBS3_trimmed3.xml`
- `aircraft/GBS3_trimmed3/GBS3_trimmed3_init.xml`
- `scripts/GBS3_trimmed3_script.xml`
- `engine/Zenoah_G-26A_GBS3_trimmed3.xml`

## Changes with respect to the original

### Inertia and mass

- `Ixz`: `0.0951 -> 0.0951`
- `CG.x`: `14.4881 in -> 19.9020 in`

### Aerodynamics

These terms change:

- `CDo`: `0.1000 -> 0.053824`
- `Drag_induced`: `0.0877 -> 0.014949`
- `Cldr`: `-0.1070 -> -0.0961`

### Propulsion and FCS

- Engine: `Zenoah_G-26A -> Zenoah_G-26A_GBS3_trimmed3`
- Power: `1789.68 W -> 121.75 W`
- The throttle channel 

### Init

`GBS3_trimmed3_init.xml`:

- `vt = 20 kt`
- `altitude = 1000 ft`
- `phi = 0.0 deg`
- `theta = 9.566408 deg`
- `alpha = 9.566408 deg`

### Script

`GBS3_trimmed3_script.xml`:

- total duration: `5 s`
- `dt = 0.01`
- keeps the seeded controls until `t < 1.0 s`
- performs an early stabilizing trim at `t = 0.05 s`
- performs the full trim at `t = 1.0 s`

Control:

- `aileron = 0.0`
- `rudder = 0.0`
- `pitch-trim = -0.001348`
- `throttle = 0.494078`

## Summary

`GBS3_trimmed3` in aerodynamics changes `CDo` and `Drag_induced`, and it also needs a change in `Cldr`.

