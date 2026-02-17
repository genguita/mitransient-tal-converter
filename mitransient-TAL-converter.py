import sys
import os
import argparse
import h5py
import numpy as np


def main(args):
    import mitsuba as mi
    import drjit as dr
    mi.set_variant(args.variant)
    import mitransient as mitr

    scene = mi.load_file(args.scene_file)
    integrator = scene.integrator()
    sensor = scene.sensors()[0]
    film = sensor.film()
    emitter = scene.emitters()[0]
    scene_params = mi.traverse(scene)

    # Fill TAL capture data info
    TAL_dict = dict()

    # Transient sequence parameters
    T = film.temporal_bins
    TAL_dict['t_start'] = film.start_opl
    TAL_dict['delta_t'] = film.bin_width_opl
    TAL_dict['t_account_first_and_last_bounces'] = integrator.account_first_and_last_bounces

    # Data format parameters
    scan_size = film.size().numpy()
    is_exhaustive = integrator.capture_type == 3
    if is_exhaustive:
        H_shape = (T, film.laser_scan_width, film.laser_scan_height, scan_size[0], scan_size[1])
    else:
        H_shape = (T, scan_size[0], scan_size[1])

    TAL_dict['H_format'] = 2 if is_exhaustive else 1 # T_Lx_Ly_Sx_Sy if exhaustive, T_Sx_Sy if single of confocal
    TAL_dict['sensor_grid_format'] = 2 # X_Y_3 (3D points)
    TAL_dict['laser_grid_format'] = 2

    # Sensor and laser positions and point grids
    # TODO: get laser and sensor, intersect with the scene to obtain these values
    TAL_dict['sensor_xyz'] = sensor.m_to_world.translation().numpy().flatten()
    TAL_dict['sensor_grid_xyz'] = 0
    TAL_dict['sensor_grid_normals'] = 0
    TAL_dict['laser_xyz'] = mi.traverse(emitter)['to_world'].translation().numpy().flatten()
    TAL_dict['laser_grid_xyz'] = 0
    TAL_dict['laser_grid_normals'] = 0

    transient_data = np.zeros(H_shape)
    if not args.dryrun:
        _, transient_data = mi.render(scene)
        # TODO: reshape into TAL's format

    TAL_dict['H'] = transient_data

    # Write TAL compatible data to a HDF5 file
    out_file = h5py.File(args.output_file, 'w')
    # TODO: full hdf5 file
    out_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="mitransient-TAL-converter")
    parser.add_argument("scene_file", type=str,
                        help="XML file describing the scene to be rendered with mitransient")
    parser.add_argument("-v", "--variant", type=str, default="llvm_mono",
                        help="Mitsuba variant used for rendering (default = llvm_mono)")
    parser.add_argument("-d", "--dryrun", action="store_true",
                        help="If set, the script does not render the scene, but fills the rest of the HDF5 file")
    parser.add_argument("-o", "--output_file", type=str, default="./output.hdf5",
                        help="Path to save the output HDF5 file (default = ./output.hdf5)")
    args = parser.parse_args(sys.argv[1:])
    main(args)

