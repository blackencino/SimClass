import { is_safe_divide, safe_divide_or } from "./common";
import { MAX_NEIGHBORS, Neighborhood } from "./neighborhood";

//------------------------------------------------------------------------------
function iisph_pseudo_ap_density_stars_and_pseudo_diagonals(
    particle_count: number,
    dt: number,
    target_density: number,
    density_stars: Float32Array,
    pseudo_diagonals: Float32Array,

    densities: Float32Array,

    fluid_volumes: Float32Array,
    fluid_velocities: Float32Array,
    fluid_neighborhood: Neighborhood,

    solid_volumes: Float32Array,
    solid_velocities: Float32Array,
    solid_neighborhood: Neighborhood
) {
    for (
        let particle_index = 0;
        particle_index < particle_count;
        ++particle_index
    ) {
        const fluid_nbhd_count = fluid_neighborhood.counts[particle_index];
        const solid_nbhd_count = solid_neighborhood.counts[particle_index];
        if (fluid_nbhd_count + solid_nbhd_count < 1) {
            density_stars[particle_index] = densities[particle_index];
            pseudo_diagonals[particle_index] = 0.0;
            return;
        }

        const self_mass = fluid_volumes[particle_index] * target_density;
        const self_density = densities[particle_index];
        const self_velocity_x = fluid_velocities[particle_index * 2];
        const self_velocity_y = fluid_velocities[particle_index * 2 + 1];

        let sum_m_gradw_x = 0.0;
        let sum_m_gradw_y = 0.0;
        let sum_m_gradw_dot_gradw = 0.0;

        let divergence = 0.0;

        // Loop over fluids
        if (fluid_nbhd_count > 0) {
            for (let j = 0; j < fluid_nbhd_count; ++j) {
                const neighbor_particle_index =
                    fluid_neighborhood.indices[
                        particle_index * MAX_NEIGHBORS + j
                    ];
                const grad_w_x =
                    fluid_neighborhood.kernel_gradients[
                        (particle_index * MAX_NEIGHBORS + j) * 2
                    ];
                const grad_w_y =
                    fluid_neighborhood.kernel_gradients[
                        (particle_index * MAX_NEIGHBORS + j) * 2 + 1
                    ];
                const neighbor_mass =
                    fluid_volumes[neighbor_particle_index] * target_density;

                sum_m_gradw_x += neighbor_mass * grad_w_x;
                sum_m_gradw_y += neighbor_mass * grad_w_y;
                sum_m_gradw_dot_gradw +=
                    neighbor_mass * (grad_w_x * grad_w_x + grad_w_y * grad_w_y);

                const delta_vel_x =
                    self_velocity_x -
                    fluid_velocities[neighbor_particle_index * 2];
                const delta_vel_y =
                    self_velocity_y -
                    fluid_velocities[neighbor_particle_index * 2 + 1];
                divergence +=
                    neighbor_mass *
                    (delta_vel_x * grad_w_x + delta_vel_y * grad_w_y);
            }
        }

        // Loop over solids
        // Note that the sum of dot products is omitted, this is due to the fact that the
        // solids are not moved in this simulation phase by the fluids.
        if (solid_nbhd_count > 0) {
            for (let j = 0; j < solid_nbhd_count; ++j) {
                const neighbor_particle_index =
                    solid_neighborhood.indices[
                        particle_index * MAX_NEIGHBORS + j
                    ];
                const grad_w_x =
                    solid_neighborhood.kernel_gradients[
                        (particle_index * MAX_NEIGHBORS + j) * 2
                    ];
                const grad_w_y =
                    solid_neighborhood.kernel_gradients[
                        (particle_index * MAX_NEIGHBORS + j) * 2 + 1
                    ];
                const neighbor_mass =
                    solid_volumes[neighbor_particle_index] * target_density;

                sum_m_gradw_x += neighbor_mass * grad_w_x;
                sum_m_gradw_y += neighbor_mass * grad_w_y;

                const delta_vel_x =
                    self_velocity_x -
                    solid_velocities[neighbor_particle_index * 2];
                const delta_vel_y =
                    self_velocity_y -
                    solid_velocities[neighbor_particle_index * 2 + 1];
                divergence +=
                    neighbor_mass *
                    (delta_vel_x * grad_w_x + delta_vel_y * grad_w_y);
            }
        }

        const self_density_sqr = self_density * self_density;

        sum_m_gradw_dot_gradw *= self_mass;
        if (
            is_safe_divide(sum_m_gradw_x, self_density) &&
            is_safe_divide(sum_m_gradw_y, self_density) &&
            is_safe_divide(sum_m_gradw_dot_gradw, self_density_sqr)
        ) {
            sum_m_gradw_x /= self_density;
            sum_m_gradw_y /= self_density;
            sum_m_gradw_dot_gradw /= self_density_sqr;
            const aii =
                sum_m_gradw_x * sum_m_gradw_x +
                sum_m_gradw_y * sum_m_gradw_y +
                sum_m_gradw_dot_gradw;
            pseudo_diagonals[particle_index] = -aii;
        } else {
            pseudo_diagonals[particle_index] = 0.0;
        }

        density_stars[particle_index] = self_density + dt * divergence;
    }
}

//------------------------------------------------------------------------------
function iisph_pseudo_ap_compute_pseudo_pressure_displacements_partial(
    particle_count: number,
    boundary: boolean,
    target_density: number,

    pseudo_pressure_displacements: Float32Array,

    pseudo_pressures: Float32Array,
    densities: Float32Array,

    neighbor_volumes: Float32Array,
    neighborhood: Neighborhood
): void {
    for (
        let particle_index = 0;
        particle_index < particle_count;
        ++particle_index
    ) {
        const nbhd_count = neighborhood.counts[particle_index];
        if (nbhd_count < 1) {
            if (!boundary) {
                pseudo_pressure_displacements[2 * particle_index] = 0.0;
                pseudo_pressure_displacements[2 * particle_index + 1] = 0.0;
            }
            return;
        }

        const self_pseudo_pressure = pseudo_pressures[particle_index];
        const self_density = densities[particle_index];
        const p_over_rho_sqr = safe_divide_or(
            self_pseudo_pressure,
            self_density * self_density,
            0.0
        );

        let pseudo_pressure_displacement_x = 0.0;
        let pseudo_pressure_displacement_y = 0.0;
        if (boundary) {
            const coeff = p_over_rho_sqr;
            for (let j = 0; j < nbhd_count; ++j) {
                const neighbor_particle_index =
                    neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
                const neighbor_mass =
                    target_density * neighbor_volumes[neighbor_particle_index];
                const grad_w_x =
                    neighborhood.kernel_gradients[
                        (particle_index * MAX_NEIGHBORS + j) * 2
                    ];
                const grad_w_y =
                    neighborhood.kernel_gradients[
                        (particle_index * MAX_NEIGHBORS + j) * 2 + 1
                    ];

                pseudo_pressure_displacement_x +=
                    -neighbor_mass * coeff * grad_w_x;
                pseudo_pressure_displacement_y +=
                    -neighbor_mass * coeff * grad_w_y;
            }
        } else {
            for (let j = 0; j < nbhd_count; ++j) {
                const neighbor_particle_index =
                    neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
                const neighbor_mass =
                    target_density * neighbor_volumes[neighbor_particle_index];
                const grad_w_x =
                    neighborhood.kernel_gradients[
                        (particle_index * MAX_NEIGHBORS + j) * 2
                    ];
                const grad_w_y =
                    neighborhood.kernel_gradients[
                        (particle_index * MAX_NEIGHBORS + j) * 2 + 1
                    ];

                const neighbor_pseudo_pressure =
                    pseudo_pressures[neighbor_particle_index];
                const neighbor_density = densities[neighbor_particle_index];
                const coeff =
                    p_over_rho_sqr +
                    safe_divide_or(
                        neighbor_pseudo_pressure,
                        neighbor_density * neighbor_density,
                        0.0
                    );

                pseudo_pressure_displacement_x +=
                    -neighbor_mass * coeff * grad_w_x;
                pseudo_pressure_displacement_y +=
                    -neighbor_mass * coeff * grad_w_y;
            }
        }

        if (boundary) {
            pseudo_pressure_displacements[
                particle_index * 2
            ] += pseudo_pressure_displacement_x;
            pseudo_pressure_displacements[
                particle_index * 2 + 1
            ] += pseudo_pressure_displacement_y;
        } else {
            pseudo_pressure_displacements[
                particle_index * 2
            ] = pseudo_pressure_displacement_x;
            pseudo_pressure_displacements[
                particle_index * 2 + 1
            ] = pseudo_pressure_displacement_y;
        }
    }
}

function iisph_pseudo_ap_compute_pseudo_pressure_displacements(
    particle_count: number,
    target_density: number,

    pseudo_pressure_displacements: Float32Array,

    pseudo_pressures: Float32Array,
    densities: Float32Array,

    fluid_volumes: Float32Array,
    fluid_neighborhood: Neighborhood,
    solid_volumes: Float32Array,
    solid_neighborhood: Neighborhood
) {
    iisph_pseudo_ap_compute_pseudo_pressure_displacements_partial(
        particle_count,
        false,
        target_density,
        pseudo_pressure_displacements,
        pseudo_pressures,
        densities,
        fluid_volumes,
        fluid_neighborhood
    );

    iisph_pseudo_ap_compute_pseudo_pressure_displacements_partial(
        particle_count,
        true,
        target_density,
        pseudo_pressure_displacements,
        pseudo_pressures,
        densities,
        solid_volumes,
        solid_neighborhood
    );
}