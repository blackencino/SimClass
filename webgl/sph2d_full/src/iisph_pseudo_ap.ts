import { is_safe_divide, safe_divide_or, Config, sqr } from "./common";
import { State, Fluid_temp_data, Solid_temp_data } from "./state";
import {
    MAX_NEIGHBORS,
    Neighborhood,
    compute_all_neighborhoods,
    compute_all_neighborhood_kernels,
} from "./neighborhood";
import {
    compute_densities,
    integrate_velocities_in_place,
    fill_float_array,
    User_forces_function,
} from "./sph_common";

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
    for (let particle_index = 0; particle_index < particle_count; ++particle_index) {
        const fluid_nbhd_count = fluid_neighborhood.counts[particle_index];
        const solid_nbhd_count = solid_neighborhood.counts[particle_index];
        if (fluid_nbhd_count + solid_nbhd_count < 1) {
            density_stars[particle_index] = densities[particle_index];
            pseudo_diagonals[particle_index] = 0.0;
            continue;
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
        for (let j = 0; j < fluid_nbhd_count; ++j) {
            const neighbor_particle_index =
                fluid_neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
            const grad_w_x =
                fluid_neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2];
            const grad_w_y =
                fluid_neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2 + 1];
            const neighbor_mass = fluid_volumes[neighbor_particle_index] * target_density;

            sum_m_gradw_x += neighbor_mass * grad_w_x;
            sum_m_gradw_y += neighbor_mass * grad_w_y;
            sum_m_gradw_dot_gradw += neighbor_mass * (grad_w_x * grad_w_x + grad_w_y * grad_w_y);

            const delta_vel_x = self_velocity_x - fluid_velocities[neighbor_particle_index * 2];
            const delta_vel_y = self_velocity_y - fluid_velocities[neighbor_particle_index * 2 + 1];
            divergence += neighbor_mass * (delta_vel_x * grad_w_x + delta_vel_y * grad_w_y);
        }

        // Loop over solids
        // Note that the sum of dot products is omitted, this is due to the fact that the
        // solids are not moved in this simulation phase by the fluids.
        for (let j = 0; j < solid_nbhd_count; ++j) {
            const neighbor_particle_index =
                solid_neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
            const grad_w_x =
                solid_neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2];
            const grad_w_y =
                solid_neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2 + 1];
            const neighbor_mass = solid_volumes[neighbor_particle_index] * target_density;

            sum_m_gradw_x += neighbor_mass * grad_w_x;
            sum_m_gradw_y += neighbor_mass * grad_w_y;

            const delta_vel_x = self_velocity_x - solid_velocities[neighbor_particle_index * 2];
            const delta_vel_y = self_velocity_y - solid_velocities[neighbor_particle_index * 2 + 1];
            divergence += neighbor_mass * (delta_vel_x * grad_w_x + delta_vel_y * grad_w_y);
        }

        const numer = -((sum_m_gradw_x * sum_m_gradw_x + sum_m_gradw_y * sum_m_gradw_y) +
                        self_mass * sum_m_gradw_dot_gradw);
        const denom = self_density * self_density;

        pseudo_diagonals[particle_index] = safe_divide_or(numer, denom, 0.0);
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
    for (let particle_index = 0; particle_index < particle_count; ++particle_index) {
        const nbhd_count = neighborhood.counts[particle_index];
        if (nbhd_count < 1) {
            if (!boundary) {
                pseudo_pressure_displacements[2 * particle_index] = 0.0;
                pseudo_pressure_displacements[2 * particle_index + 1] = 0.0;
            }
            continue;
        }

        const self_pseudo_pressure = pseudo_pressures[particle_index];
        const self_density = densities[particle_index];
        const p_over_rho_sqr = safe_divide_or(
            self_pseudo_pressure,
            sqr(self_density),
            0.0
        );

        let pseudo_pressure_displacement_x = 0.0;
        let pseudo_pressure_displacement_y = 0.0;
        if (boundary) {
            const coeff = p_over_rho_sqr;
            for (let j = 0; j < nbhd_count; ++j) {
                const neighbor_particle_index =
                    neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
                const neighbor_mass = target_density * neighbor_volumes[neighbor_particle_index];
                const grad_w_x =
                    neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2];
                const grad_w_y =
                    neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2 + 1];

                pseudo_pressure_displacement_x += -neighbor_mass * coeff * grad_w_x;
                pseudo_pressure_displacement_y += -neighbor_mass * coeff * grad_w_y;
            }
        } else {
            for (let j = 0; j < nbhd_count; ++j) {
                const neighbor_particle_index =
                    neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
                const neighbor_mass = target_density * neighbor_volumes[neighbor_particle_index];
                const grad_w_x =
                    neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2];
                const grad_w_y =
                    neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2 + 1];

                const neighbor_pseudo_pressure = pseudo_pressures[neighbor_particle_index];
                const neighbor_density = densities[neighbor_particle_index];
                const coeff =
                    p_over_rho_sqr +
                    safe_divide_or(
                        neighbor_pseudo_pressure,
                        sqr(neighbor_density),
                        0.0
                    );

                pseudo_pressure_displacement_x += -neighbor_mass * coeff * grad_w_x;
                pseudo_pressure_displacement_y += -neighbor_mass * coeff * grad_w_y;
            }
        }

        if (boundary) {
            pseudo_pressure_displacements[particle_index * 2] += pseudo_pressure_displacement_x;
            pseudo_pressure_displacements[particle_index * 2 + 1] += pseudo_pressure_displacement_y;
        } else {
            pseudo_pressure_displacements[particle_index * 2] = pseudo_pressure_displacement_x;
            pseudo_pressure_displacements[particle_index * 2 + 1] = pseudo_pressure_displacement_y;
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

function iisph_pseudo_ap_iterate_pseudo_pressures_in_place(
    particle_count: number,
    target_density: number,
    omega: number,
    pseudo_pressures: Float32Array,
    density_stars: Float32Array,
    pseudo_pressure_displacements: Float32Array,
    pseudo_diagonals: Float32Array,
    fluid_volumes: Float32Array,
    fluid_neighborhood: Neighborhood,
    solid_volumes: Float32Array,
    solid_neighborhood: Neighborhood
): Array<number> {
    if (particle_count < 1) {
        return [0.0, 0.0];
    }

    // Also computes the average density variation
    let error_sum = 0.0;
    let error_max = 0.0;

    for (let particle_index = 0; particle_index < particle_count; ++particle_index) {
        const fluid_nbhd_count = fluid_neighborhood.counts[particle_index];
        const solid_nbhd_count = solid_neighborhood.counts[particle_index];
        if (fluid_nbhd_count < 1 && solid_nbhd_count < 1) {
            pseudo_pressures[particle_index] = 0.0;
            continue;
        }

        // This is checking to see whether a sentinel value was set during
        // the computation of the pseudo_diagonal. We'll check proper
        // proportional safe divides below.
        const pseudo_diagonal = pseudo_diagonals[particle_index];
        if (pseudo_diagonal == 0.0) {
            pseudo_pressures[particle_index] = 0.0;
            continue;
        }

        let Api = 0.0;
        const self_pseudo_pressure_displacement_x =
            pseudo_pressure_displacements[particle_index * 2];
        const self_pseudo_pressure_displacement_y =
            pseudo_pressure_displacements[particle_index * 2 + 1];

        for (let j = 0; j < fluid_nbhd_count; ++j) {
            const neighbor_particle_index =
                fluid_neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
            const neighbor_mass = target_density * fluid_volumes[neighbor_particle_index];
            const grad_w_x =
                fluid_neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2];
            const grad_w_y =
                fluid_neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2 + 1];

            const neighbor_pseudo_pressure_displacement_x =
                pseudo_pressure_displacements[neighbor_particle_index * 2];
            const neighbor_pseudo_pressure_displacement_y =
                pseudo_pressure_displacements[neighbor_particle_index * 2 + 1];

            const dp_x =
                self_pseudo_pressure_displacement_x - neighbor_pseudo_pressure_displacement_x;
            const dp_y =
                self_pseudo_pressure_displacement_y - neighbor_pseudo_pressure_displacement_y;

            Api += neighbor_mass * (dp_x * grad_w_x + dp_y * grad_w_y);
        }

        for (let j = 0; j < solid_nbhd_count; ++j) {
            const neighbor_particle_index =
                solid_neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
            const neighbor_mass = target_density * solid_volumes[neighbor_particle_index];
            const grad_w_x =
                solid_neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2];
            const grad_w_y =
                solid_neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2 + 1];

            Api +=
                neighbor_mass *
                (self_pseudo_pressure_displacement_x * grad_w_x +
                    self_pseudo_pressure_displacement_y * grad_w_y);
        }

        const source = target_density - density_stars[particle_index];
        const error = Math.max(0.0, safe_divide_or(Api - source, target_density, 0.0));
        error_max = Math.max(error, error_max);
        error_sum += error;

        const old_pseudo_pressure = pseudo_pressures[particle_index];
        const numer = omega * (source - Api);
        pseudo_pressures[particle_index] = Math.max(
            0.0,
            old_pseudo_pressure + safe_divide_or(numer, pseudo_diagonal, 0.0)
        );
    }

    return [error_sum / particle_count, error_max];
}

//------------------------------------------------------------------------------
// Remember that this one is special! Don't replace it with the common version!!
function iisph_pseudo_ap_integrate_velocities_and_positions_in_place(
    particle_count: number,
    dt: number,
    velocities: Float32Array,
    positions: Float32Array,
    pseudo_pressure_displacements: Float32Array
): void {
    for (let i = 0; i < particle_count; ++i) {
        const velocity_star_x = velocities[i*2];
        const velocity_star_y = velocities[i*2+1];

        const old_position_x = positions[i*2];
        const old_position_y = positions[i*2+1];

        const displacement_x = pseudo_pressure_displacements[i*2];
        const displacement_y = pseudo_pressure_displacements[i*2+1];

        positions[i*2] = old_position_x + dt * velocity_star_x + displacement_x;
        positions[i*2+1] = old_position_y + dt * velocity_star_y + displacement_y;
        if (is_safe_divide(velocity_star_x, dt) &&
            is_safe_divide(velocity_star_y, dt)) {
            velocities[i*2] = velocity_star_x + displacement_x / dt;
            velocities[i*2+1] = velocity_star_y + displacement_y / dt;
        }
    }
}

export function iisph_pseudo_ap_sub_step(
    global_time: number,
    dt: number,
    config: Config,
    fluid_state: State,
    fluid_temp_data: Fluid_temp_data,
    solid_state: State,
    solid_temp_data: Solid_temp_data,
    user_forces: User_forces_function
): void {
    const particle_count = fluid_state.positions.length / 2;

    compute_all_neighborhoods(config, fluid_state, fluid_temp_data, solid_state, solid_temp_data);
    compute_all_neighborhood_kernels(config, fluid_temp_data.self_neighborhood);
    compute_all_neighborhood_kernels(config, fluid_temp_data.other_neighborhood);
    user_forces(
        global_time,
        dt,
        config,
        fluid_state,
        fluid_temp_data,
        solid_state,
        solid_temp_data
    );

    compute_densities(
        particle_count,
        config.params.support,
        config.params.target_density,
        fluid_temp_data.densities,
        fluid_temp_data.volumes,
        fluid_temp_data.self_neighborhood,
        solid_temp_data.volumes,
        fluid_temp_data.other_neighborhood
    );

    integrate_velocities_in_place(
        particle_count,
        dt,
        config.params.target_density,
        fluid_state.velocities,
        fluid_temp_data.volumes,
        fluid_temp_data.external_forces
    );

    iisph_pseudo_ap_density_stars_and_pseudo_diagonals(
        particle_count,
        dt,
        config.params.target_density,
        fluid_temp_data.density_stars,
        fluid_temp_data.pseudo_diagonals,
        fluid_temp_data.densities,
        fluid_temp_data.volumes,
        fluid_state.velocities,
        fluid_temp_data.self_neighborhood,
        solid_temp_data.volumes,
        solid_state.velocities,
        fluid_temp_data.other_neighborhood
    );

    fill_float_array(particle_count, 0.0, fluid_temp_data.pseudo_pressures);

    let iter = 0;
    for (; iter < config.params.iisph.max_pressure_iterations; ++iter) {
        iisph_pseudo_ap_compute_pseudo_pressure_displacements(
            particle_count,
            config.params.target_density,
            fluid_temp_data.pseudo_pressure_displacements,
            fluid_temp_data.pseudo_pressures,
            fluid_temp_data.densities,
            fluid_temp_data.volumes,
            fluid_temp_data.self_neighborhood,
            solid_temp_data.volumes,
            fluid_temp_data.other_neighborhood
        );

        const [error_average, error_max] = iisph_pseudo_ap_iterate_pseudo_pressures_in_place(
            particle_count,
            config.params.target_density,
            config.params.iisph.omega,
            fluid_temp_data.pseudo_pressures,
            fluid_temp_data.density_stars,
            fluid_temp_data.pseudo_pressure_displacements,
            fluid_temp_data.pseudo_diagonals,
            fluid_temp_data.volumes,
            fluid_temp_data.self_neighborhood,
            solid_temp_data.volumes,
            fluid_temp_data.other_neighborhood
        );

        // console.log("\t Pseudo Pressure iter, error avg, max: ",
        //            iter,
        //            error_average,
        //            error_max);

        if (
            iter > 1 &&
            error_average <= config.params.iisph.error_average_threshold &&
            error_max <= config.params.iisph.error_max_threshold
        ) {
            break;
        }
    }
    //console.log("IISPH PSEUDO AP pressure iters: ", iter);

    iisph_pseudo_ap_compute_pseudo_pressure_displacements(
        particle_count,
        config.params.target_density,
        fluid_temp_data.pseudo_pressure_displacements,
        fluid_temp_data.pseudo_pressures,
        fluid_temp_data.densities,
        fluid_temp_data.volumes,
        fluid_temp_data.self_neighborhood,
        solid_temp_data.volumes,
        fluid_temp_data.other_neighborhood
    );

    iisph_pseudo_ap_integrate_velocities_and_positions_in_place(
        particle_count,
        dt,
        fluid_state.velocities,
        fluid_state.positions,
        fluid_temp_data.pseudo_pressure_displacements
    );
}