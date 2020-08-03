import { is_safe_divide, safe_divide_or, kernel_sph_w, kernel_sph_gradw, Config } from "./common";
import { MAX_NEIGHBORS, Neighborhood } from "./neighborhood";
import { State, Fluid_temp_data, Solid_temp_data } from "./state";

//------------------------------------------------------------------------------
export function fill_float_array(count: number, value: number, values: Float32Array) {
    for (let i = 0; i < count; ++i) {
        values[i] = value;
    }
}

//------------------------------------------------------------------------------
// Velocities with only external forces
export function integrate_velocities_in_place(
    particle_count: number,
    dt: number,
    target_density: number,
    velocities: Float32Array,
    volumes: Float32Array,
    forces: Float32Array
): void {
    for (let particle_index = 0; particle_index < particle_count; ++particle_index) {
        const denom = target_density * volumes[particle_index];
        const numer_x = dt * forces[particle_index * 2];
        const numer_y = dt * forces[particle_index * 2 + 1];
        if (is_safe_divide(numer_x, denom) && is_safe_divide(numer_y, denom)) {
            velocities[particle_index * 2] += numer_x / denom;
            velocities[particle_index * 2 + 1] += numer_y / denom;
        }
    }
}

//------------------------------------------------------------------------------
export function compute_volumes(
    particle_count: number,
    support: number,
    volumes: Float32Array,
    neighbor_counts: Uint8Array,
    neighbor_kernels: Float32Array
): void {
    const muchness_init = kernel_sph_w(0.0, support);
    for (let particle_index = 0; particle_index < particle_count; ++particle_index) {
        const nbhd_count = neighbor_counts[particle_index];
        if (nbhd_count < 1) {
            volumes[particle_index] = 1.0 / muchness_init;
            continue;
        }

        let muchness = muchness_init;
        for (let j = 0; j < nbhd_count; ++j) {
            muchness += neighbor_kernels[particle_index * MAX_NEIGHBORS + j];
        }

        // CJH HACK numer should be 1.0f
        volumes[particle_index] = 1.7 / muchness;
    }
}

//------------------------------------------------------------------------------
export function compute_constant_volumes(
    particle_count: number,
    target_density: number,
    mass_per_particle: number,
    volumes: Float32Array
) {
    // mass = volume * density,
    // volume = mass / density
    fill_float_array(particle_count, mass_per_particle / target_density, volumes);
}

//------------------------------------------------------------------------------
// Densities
function compute_densities_partial(
    particle_count: number,
    boundary: boolean,
    support: number,
    target_density: number,
    densities: Float32Array,
    self_volumes: Float32Array,
    neighbor_volumes: Float32Array,
    neighborhood: Neighborhood
): void {
    const w0 = kernel_sph_w(0.0, support);
    for (let particle_index = 0; particle_index < particle_count; ++particle_index) {
        let density = 0.0;
        if (!boundary) {
            density = self_volumes[particle_index] * target_density * w0;
        }
        const nbhd_count = neighborhood.counts[particle_index];
        if (nbhd_count < 1) {
            if (!boundary) {
                densities[particle_index] = density;
            }
            continue;
        }

        for (let j = 0; j < nbhd_count; ++j) {
            const neighbor_particle_index =
                neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
            const mass = target_density * neighbor_volumes[neighbor_particle_index];
            density += neighborhood.kernels[particle_index * MAX_NEIGHBORS + j] * mass;
        }

        if (boundary) {
            densities[particle_index] += density;
        } else {
            densities[particle_index] = density;
        }
    }
}

export function compute_densities(
    particle_count: number,
    support: number,
    target_density: number,
    densities: Float32Array,
    fluid_volumes: Float32Array,
    fluid_neighborhood: Neighborhood,
    solid_volumes: Float32Array,
    solid_neighborhood: Neighborhood
): void {
    compute_densities_partial(
        particle_count,
        false,
        support,
        target_density,
        densities,
        fluid_volumes,
        fluid_volumes,
        fluid_neighborhood
    );

    compute_densities_partial(
        particle_count,
        true,
        support,
        target_density,
        densities,
        fluid_volumes,
        solid_volumes,
        solid_neighborhood
    );
}

//------------------------------------------------------------------------------
export function cfl_maximum_time_step(
    particle_count: number,
    support: number,
    factor: number,
    velocities: Float32Array
) {
    let max_squared_vel = 0.0;
    for (let i = 0; i < particle_count; ++i) {
        const vx = velocities[i * 2];
        const vy = velocities[i * 2 + 1];
        max_squared_vel = Math.max(max_squared_vel, vx * vx + vy * vy);
    }

    const numer = factor * support;
    const denom = Math.sqrt(max_squared_vel);

    return safe_divide_or(numer, denom, 1.0);
}

//------------------------------------------------------------------------------
export type User_forces_function = (
    global_time: number,
    dt: number,
    config: Config,
    fluid_state: State,
    fluid_temp_data: Fluid_temp_data,
    solid_state: State,
    solid_temp_data: Solid_temp_data
) => void;

