import sys
import time
import argparse
import numpy as np
import h5py


def get_sensor_laser_intersections(scene, scan_size, laser_scan_size, force_equal_grids, is_single_capture):
    """
        Computes sensor and laser intersections with the scene to fill the sensor/laser_grid position and normals
        required for tal
    """
    import mitsuba as mi
    import drjit as dr

    # Get sensor and emitter
    sensor = scene.sensors()[0]
    emitter = scene.emitters()[0]

    # Sample rays from the sensor
    scan_x, scan_y = dr.meshgrid(
        dr.linspace(mi.Float, 0.0, 1.0, scan_size[0], endpoint=False),
        dr.linspace(mi.Float, 0.0, 1.0, scan_size[1], endpoint=False))
    points: mi.Point2f = mi.Point2f(scan_x, scan_y)
    sensor_rays, _ = sensor.sample_ray(mi.Float(0.0), mi.Float(0.0), points, mi.Point2f(0.0), mi.Bool(True))

    # Intersect with the scene to obtain the scanned points
    si_sensor = scene.ray_intersect(sensor_rays, ray_flags=mi.RayFlags.All, coherent=mi.Bool(True))

    si_laser = None
    if is_single_capture:
        trafo = emitter.world_transform()
        laser_dir: mi.Vector3f = trafo @ mi.Vector3f(0, 0, 1)

        laser_ray = mi.Ray3f(mi.Point3f(trafo.translation()), laser_dir)
        si_laser = scene.ray_intersect(laser_ray)

        assert dr.all(si_laser.is_valid()), \
            'The emitter is not pointing at the scene!'
    elif not force_equal_grids:
        # Else, sample rays from the emitter to obtain illuminated points
        # Dummy emitter, with custom FOV
        dummy_emitter = mi.load_dict({
            "type": "projector",
            "fov": scene.integrator().illumination_scan_fov,
            "to_world": mi.ScalarTransform4f(emitter.world_transform().matrix.to_numpy()[..., 0]),
        })

        # Generate a ray for each of the illuminated points
        laser_x, laser_y = dr.meshgrid(
            dr.linspace(mi.Float, 0.0, 1.0, laser_scan_size[0], endpoint=False),
            dr.linspace(mi.Float, 0.0, 1.0, laser_scan_size[1], endpoint=False),
        )
        points: mi.Point2f = mi.Point2f(laser_x, laser_y)
        laser_rays, _ = dummy_emitter.sample_ray(mi.Float(0.0), mi.Float(0.0), mi.Point2f(0.0), points, mi.Bool(True))

        # Intersect with the scene to obtain the illuminated points
        si_laser = scene.ray_intersect(laser_rays, ray_flags=mi.RayFlags.All, coherent=mi.Bool(True))
    else:
        # If scanning is equal, scanned and illuminated points are the same
        si_laser = si_sensor

    # Transform point positions and normals to tal compatible format
    sensor_grid_xyz = dr.reshape(mi.TensorXf, si_sensor.p, (scan_size[0], scan_size[1], 3)).numpy()
    sensor_grid_normals = dr.reshape(mi.TensorXf, si_sensor.n, (scan_size[0], scan_size[1], 3)).numpy()
    sensor_valid_mask = dr.reshape(mi.TensorXb, si_sensor.is_valid(), (scan_size[0], scan_size[1])).numpy()
    laser_grid_xyz = dr.reshape(mi.TensorXf, si_laser.p, (laser_scan_size[0], laser_scan_size[1], 3)).numpy()
    laser_grid_normals = dr.reshape(mi.TensorXf, si_laser.n, (laser_scan_size[0], laser_scan_size[1], 3)).numpy()
    laser_valid_mask = dr.reshape(mi.TensorXb, si_laser.is_valid(), (laser_scan_size[0], laser_scan_size[1])).numpy()

    return sensor_grid_xyz, sensor_grid_normals, sensor_valid_mask, laser_grid_xyz, laser_grid_normals, laser_valid_mask


def main(args):
    import mitsuba as mi
    import drjit as dr
    mi.set_variant(args.variant)
    import mitransient as mitr

    is_polarized = 'polarized' in args.variant
    if is_polarized:
        print("NOTE: tal does not implement any reconstruction algorithms that take polarization into account. "
              "If you want to use the output file for reconstruction, you must discard everything but the "
              "intensity component of the Stokes vector.")

    # Load the scene
    scene = mi.load_file(args.scene_file)
    integrator = scene.integrator()
    sensor = scene.sensors()[0]
    film = sensor.film()
    emitter = scene.emitters()[0]

    # Fill tal capture data info
    tal_dict = dict()

    # Transient sequence parameters
    T = film.temporal_bins
    tal_dict['t_start'] = film.start_opl
    tal_dict['delta_t'] = film.bin_width_opl
    tal_dict['t_accounts_first_and_last_bounces'] = integrator.account_first_and_last_bounces

    # Data format parameters
    scan_size = film.size().numpy()
    force_equal_scan = integrator.force_equal_grids
    is_single = integrator.capture_type == 1
    if is_single:
        laser_scan_size = (1, 1)
    elif force_equal_scan:
        laser_scan_size = scan_size
    else:
        laser_scan_size = (film.laser_scan_width if film.laser_scan_width > 1 else 1,
                           film.laser_scan_height if film.laser_scan_height > 1 else 1)
    is_exhaustive = integrator.capture_type == 3
    if is_exhaustive:
        H_shape = (T, laser_scan_size[0], laser_scan_size[1], scan_size[0], scan_size[1])
    else:
        H_shape = (T, scan_size[0], scan_size[1])

    tal_dict['H_format'] = 2 if is_exhaustive else 1 # T_Lx_Ly_Sx_Sy if exhaustive, T_Sx_Sy if single of confocal
    tal_dict['sensor_grid_format'] = 2 # X_Y_3 (3D points)
    tal_dict['laser_grid_format'] = 2

    # Sensor and laser positions and point grids
    sensor_grid_xyz, sensor_grid_normals, sensor_grid_valid, laser_grid_xyz, laser_grid_normals, laser_grid_valid = (
        get_sensor_laser_intersections(scene, scan_size, laser_scan_size, force_equal_scan, is_single))

    tal_dict['sensor_xyz'] = sensor.m_to_world.translation().numpy().flatten()
    tal_dict['sensor_grid_xyz'] = sensor_grid_xyz
    tal_dict['sensor_grid_normals'] = sensor_grid_normals
    tal_dict['laser_xyz'] = mi.traverse(emitter)['to_world'].translation().numpy().flatten()
    tal_dict['laser_grid_xyz'] = laser_grid_xyz
    tal_dict['laser_grid_normals'] = laser_grid_normals

    # Render the NLOS scene
    transient_data = np.zeros(H_shape)
    if not args.dryrun:
        dr.print("Rendering the NLOS scene...")
        start = time.time()
        _, transient_data = mi.render(scene)
        transient_data = np.array(transient_data)
        print(f"Rendering done, took {time.time() - start:.3f} seconds")

        # Reshape to match tal's H format
        transient_data = np.moveaxis(transient_data, -2, 0) # Time dimension should be first
        if not is_polarized:
            # If simulated than one channel (RGB), sum them, except for polarized variants
            transient_data = np.sum(transient_data, axis=-1)

        if is_exhaustive:
            # Swap sensor scan and laser illumination dimensions
            transient_data = transient_data.swapaxes(1, 3).swapaxes(2, 4)

            # Zero-out invalid laser and sensor scans (laser or sensor ray does not intersect the scene)
            transient_data[:, ~laser_grid_valid, :, :] = 0.0
            transient_data[:, :, :, ~sensor_grid_valid] = 0.0
        else:
            # Zero-out invalid sensor scans (sensor ray does not intersect the scene)
            transient_data[:,  ~sensor_grid_valid] = 0.0

    tal_dict['H'] = transient_data

    # Write tal compatible data to a HDF5 file
    out_file = h5py.File(args.output_file, 'w')
    for key, value in tal_dict.items():
        out_file[key] = value
    out_file.close()
    print(f"Saved tal compatible HDF5 file to {args.output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="mitransient-tal-converter")
    parser.add_argument("scene_file", type=str,
                        help="XML file describing the scene to be rendered with mitransient")
    parser.add_argument("-v", "--variant", type=str, default="llvm_ad_mono",
                        help="Mitsuba variant used for rendering (default = llvm_ad_mono)")
    parser.add_argument("-d", "--dryrun", action="store_true",
                        help="If set, the script does not render the scene, but fills the rest of the HDF5 file")
    parser.add_argument("-o", "--output_file", type=str, default="./output.hdf5",
                        help="Path to save the output HDF5 file (default = ./output.hdf5)")
    args = parser.parse_args(sys.argv[1:])
    main(args)

