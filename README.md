# mitransient-TAL-converter

This script converts [mitransient](https://github.com/diegoroyo/mitransient) NLOS simulation files to 
[TAL](https://github.com/diegoroyo/tal) compatible HDF5 files, allowing its easy integration for performing
NLOS reconstructions. 

While TAL provides tools for simulating NLOS scenes, these are limited to planar relay walls, uniform scanning patterns,
and hidden scenes manually defined through YAML files. On the other hand, mitransient allows defining more complex scenes,
with generic relay surfaces, hidden geometry, sensors and illumination, being also able to perform simulations on
scenes created with Blender, which are easier to edit than text configuration files. However, it does not provide NLOS
reconstruction tools.

With this script you can freely design a NLOS scene compatible with mitransient, and then perform the desired 
reconstructions with tal.

## Usage

    python mitransient-TAL-converter.py <xml scene file>

