# mitransient-tal-converter

This script converts [mitransient](https://github.com/diegoroyo/mitransient) NLOS simulation files to 
[tal](https://github.com/diegoroyo/tal) compatible HDF5 files, allowing its easy integration for performing
NLOS reconstructions. 

While tal provides tools for simulating NLOS scenes, these are limited to planar relay walls, uniform scanning patterns,
and hidden scenes manually defined through YAML files. On the other hand, mitransient allows defining more complex scenes,
with generic relay surfaces, hidden geometry, sensors and illumination, being also able to perform simulations on
scenes created with Blender, which are easier to edit than text configuration files. However, it does not provide NLOS
reconstruction tools.

With this script you can freely design a NLOS scene compatible with mitransient, and then perform the desired 
reconstructions with tal.

## Usage

To use the converter, simply create a mitransient compatible XML file, detailing the characteristics of the desired NLOS
scene, then run the script as:

    python mitransient-tal-converter.py <xml scene file> -o <output file> -v <Mitsuba variant>

This will perform the simulation of the NLOS capture process, extract the scene configuration from the XML and save 
everything in a tal compatible HDF5 file. By using the option ```--dryrun```, the simulation will not be performed, and 
instead only the scene parameters will be saved to the HDF5 file.

Also note that if you use ```rbg``` or ```spectral``` Mitsuba variants for rendering, the different color channels will
be summed into a single one, as required for tal. When using ```polarized``` variants, the complete Stokes vector will
be saved, however tal does not provide any reconstruction algorithm that leverages polarization, and will only use 
the intensity component.
