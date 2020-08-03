import { MAX_NEIGHBORS, Neighborhood } from "./neighborhood";

//------------------------------------------------------------------------------
function compute_divergences_partial(
    particle_count: number,
    boundary: boolean,
    target_density: number,
    divergences: Float32Array,
    self_velocities: Float32Array,
    neighbor_volumes: Float32Array,
    neighbor_velocities: Float32Array,
    neighborhood: Neighborhood
): void {
    for (let particle_index = 0; particle_index < particle_count; ++particle_index) {
        const nbhd_count = neighborhood.counts[particle_index];
        if (nbhd_count < 1) {
            if (!boundary) {
                divergences[particle_index] = 0.0;
            }
            return;
        }

        const self_velocity_x = self_velocities[particle_index * 2];
        const self_velocity_y = self_velocities[particle_index * 2 + 1];

        let divergence = 0.0;

        for (let j = 0; j < nbhd_count; ++j) {
            const neighbor_particle_index =
                neighborhood.indices[particle_index * MAX_NEIGHBORS + j];
            const neighbor_mass = target_density * neighbor_volumes[neighbor_particle_index];
            const grad_w_x =
                neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2];
            const grad_w_y =
                neighborhood.kernel_gradients[(particle_index * MAX_NEIGHBORS + j) * 2 + 1];
            const neighbor_velocity_x = neighbor_velocities[neighbor_particle_index * 2];
            const neighbor_velocity_y = neighbor_velocities[neighbor_particle_index * 2 + 1];

            const dv_x = self_velocity_x - neighbor_velocity_x;
            const dv_y = self_velocity_y - neighbor_velocity_y;

            divergence += neighbor_mass * (dv_x * grad_w_x + dv_y * grad_w_y);
        }

        if (!boundary) {
            divergences[particle_index] = divergence;
        } else {
            divergences[particle_index] += divergence;
        }
    }
}

export function compute_divergences(
    particle_count: number,
    target_density: number,
    divergences: Float32Array,
    fluid_volumes: Float32Array,
    fluid_velocities: Float32Array,
    fluid_neighborhood: Neighborhood,
    solid_volumes: Float32Array,
    solid_velocities: Float32Array,
    solid_neighborhood: Neighborhood
): void {
    compute_divergences_partial(
        particle_count,
        false,
        target_density,
        divergences,
        fluid_velocities,
        fluid_volumes,
        fluid_velocities,
        fluid_neighborhood
    );

    compute_divergences_partial(
        particle_count,
        true,
        target_density,
        divergences,
        fluid_velocities,
        solid_volumes,
        solid_velocities,
        solid_neighborhood
    );
}


//------------------------------------------------------------------------------
export function predict_density_stars_from_divergences(
    particle_count: number,
    dt: number,
    density_stars: Float32Array,
    densities: Float32Array,
    divergences: Float32Array
): void {
    for (let particle_index = 0; particle_index < particle_count; ++particle_index) {
        density_stars[particle_index] = Math.max(
            0.0,
            densities[particle_index] + dt * divergences[particle_index]
        );
    }
}

//------------------------------------------------------------------------------
export function predict_density_stars(
    particle_count: number,
    dt: number,
    target_density: number,
    density_stars: Float32Array,
    densities: Float32Array,
    fluid_volumes: Float32Array,
    fluid_velocities: Float32Array,
    fluid_neighborhood: Neighborhood,
    solid_volumes: Float32Array,
    solid_velocities: Float32Array,
    solid_neighborhood: Neighborhood
): void {
    compute_divergences(
        particle_count,
        target_density,
        density_stars,
        fluid_volumes,
        fluid_velocities,
        fluid_neighborhood,
        solid_volumes,
        solid_velocities,
        solid_neighborhood
    );

    predict_density_stars_from_divergences(
        particle_count,
        dt,
        density_stars,
        densities,
        density_stars
    );
}

