import seedrandom = require("seedrandom");

import { emit_solid_box, compute_grid_coords, morton, Parameters, Config } from "./common";
import {
    MAX_NEIGHBORS,
    Neighborhood,
    Block,
    compute_sorted_index_pairs,
    compute_block_map,
    compute_neighborhood_parts,
    compute_self_neighborhoods,
    compute_other_neighborhoods,
    compute_all_neighborhood_kernels,
    compute_all_neighborhoods,
    visit_neighbor_pairs,
} from "./neighborhood";

import { State, Fluid_temp_data, Solid_temp_data, create_simulation_state_of_size } from "./state";
import {
    compute_volumes,
    compute_constant_volumes,
    fill_float_array,
    cfl_maximum_time_step,
    User_forces_function,
    compute_densities,
} from "./sph_common";
import { iisph_pseudo_ap_sub_step } from "./iisph_pseudo_ap";

//import "es6-shim";

//------------------------------------------------------------------------------
function simple_simulation_step(
    state: State,
    delta_time: number,
    radius: number,
    width: number,
    height: number
): void {
    const POINT_COUNT = state.positions.length / 2;
    // X-coordinates
    for (let i = 0; i < 2 * POINT_COUNT; i += 2) {
        state.positions[i] += delta_time * state.velocities[i];
        if (state.positions[i] - radius < 0) {
            // Move point back onto rectangle and make sure it is moving
            // in positive direction
            state.positions[i] = radius - (state.positions[i] - radius);
            state.velocities[i] = Math.abs(state.velocities[i]);
        } else if (state.positions[i] + radius > width) {
            // Move point back onto rectangle and make sure it is moving
            // in negative direction
            state.positions[i] = width - (state.positions[i] + radius - width) - radius;
            state.velocities[i] = -Math.abs(state.velocities[i]);
        }
    }
    // Y-coordinates
    for (let i = 1; i < 2 * POINT_COUNT; i += 2) {
        state.positions[i] += delta_time * state.velocities[i];
        if (state.positions[i] - radius < 0) {
            // Move point back onto rectangle and make sure it is moving
            // in positive direction
            state.positions[i] = radius - (state.positions[i] - radius);
            state.velocities[i] = Math.abs(state.velocities[i]);
        } else if (state.positions[i] + radius > height) {
            // Move point back onto rectangle and make sure it is moving
            // in negative direction
            state.positions[i] = height - (state.positions[i] + radius - height) - radius;
            state.velocities[i] = -Math.abs(state.velocities[i]);
        }
    }
}

//------------------------------------------------------------------------------
function random_box_initial_state(
    count: number,
    width: number,
    height: number,
    seed: string,
    speed_gain: number
): State {
    const state = create_simulation_state_of_size(count);
    const my_rng = seedrandom(seed);
    // const my_rng = () => 0.5;
    for (let i = 0; i < count; ++i) {
        state.ids[i] = i;

        // Random point in the width/height box
        state.positions[2 * i] = width * my_rng();
        state.positions[2 * i + 1] = height * my_rng();

        const speed = speed_gain * (1.0 + 3.0 * my_rng());
        const angle = 2.0 * Math.PI * my_rng();
        state.velocities[2 * i] = speed * Math.cos(angle);
        state.velocities[2 * i + 1] = speed * Math.sin(angle);
    }
    for (let i = 0; i < 3 * count; ++i) {
        state.colors[i] = my_rng();
    }
    return state;
}

//------------------------------------------------------------------------------
function world_walls_initial_solid_state(
    support: number,
    width: number,
    height: number,
    seed: string
): State {
    const H = support;
    const Hhalf = H / 2;
    const R = 0.99 * Hhalf;
    const border_size = H * 3;

    const positions_array: Array<number> = [];

    // left wall
    emit_solid_box(positions_array, [-border_size, -border_size], [0.0, height + border_size], R);

    // right wall
    emit_solid_box(
        positions_array,
        [width, -border_size],
        [width + border_size, height + border_size],
        R
    );

    // floor wall
    emit_solid_box(
        positions_array,
        [-0.5 * border_size, -border_size],
        [width + 0.5 * border_size, 0.0],
        R
    );

    // ceiling wall
    emit_solid_box(
        positions_array,
        [-0.5 * border_size, height],
        [width + 0.5 * border_size, height + border_size],
        R
    );

    const tiny_block = (unit: Array<number> | Float32Array, roff: number): void => {
        const r = 0.04 + Math.abs(roff);
        emit_solid_box(
            positions_array,
            [width * (unit[0] - r), height * (unit[1] - r)],
            [width * (unit[0] + r), height * (unit[1] + r)],
            R
        );
    };

    const my_rng = seedrandom(seed);
    // const my_rng = () => 0.5;
    const my_rng_dist = () => 0.06 * (2 * my_rng() - 1);

    const N = 4;
    for (let by = 0; by < N; ++by) {
        const byf = (by + 1) / (N + 1);

        for (let bx = 0; bx < N; ++bx) {
            const bxf = (bx + 1) / (N + 1);
            tiny_block([bxf + my_rng_dist(), byf + my_rng_dist()], my_rng_dist());
        }
    }

    const positions = new Float32Array(positions_array);
    const count = positions.length / 2;
    const ids = new Int32Array(count);
    const velocities = new Float32Array(count * 2);
    const colors = new Float32Array(count * 3);
    for (let i = 0; i < count; ++i) {
        ids[i] = i;

        velocities[2 * i] = 0.0;
        velocities[2 * i + 1] = 0.0;

        colors[3 * i] = 0.784;
        colors[3 * i + 1] = 0.608;
        colors[3 * i + 2] = 0.078;
    }

    return new State(ids, positions, velocities, colors);
}

//------------------------------------------------------------------------------
function dam_break_initial_state(
    support: number,
    width: number,
    height: number,
    solid_state: State,
    solid_temp_data: Solid_temp_data
): State {
    const H = support;
    const Hhalf = H * 0.5;
    const R = 0.99 * Hhalf;
    const extent = 2.0 * support;

    const box_min = [0.0, 0.0];
    const box_max = [width * 3.0 / 4.0, (3.25 * height) / 4.0];

    // const max_count = 100 + estimate_solid_box_emission_count(box_min, box_max,
    // R); let state = create_simulation_state(max_count);

    // emit all the points at first.
    const positions_array: Array<number> = [];
    emit_solid_box(positions_array, box_min, box_max, R);
    const positions = new Float32Array(positions_array);

    const emitted_count = positions_array.length / 2;
    let grid_coords = new Int32Array(emitted_count * 2);
    compute_grid_coords(emitted_count, extent, grid_coords, positions);

    let sorted_index_pairs = new Uint32Array(emitted_count * 2);
    compute_sorted_index_pairs(emitted_count, sorted_index_pairs, grid_coords);

    let block_map = new Map<number, Block>();
    compute_block_map(emitted_count, block_map, sorted_index_pairs);

    let killed = new Uint8Array(emitted_count);
    let killed_count = 0;
    killed.fill(0);
    let visitor = (
        index_a: number,
        index_b: number,
        pos_a_x: number,
        pos_a_y: number,
        pos_b_x: number,
        pos_b_y: number,
        dx: number,
        dy: number,
        distance: number
    ): void => {
        if (killed[index_a] == 0) {
            killed[index_a] = 1;
            ++killed_count;
        }
    };

    visit_neighbor_pairs(
        positions,
        grid_coords,
        sorted_index_pairs,
        block_map,

        solid_state.positions,
        solid_temp_data.sorted_index_pairs,
        solid_temp_data.block_map,

        2.0 * R,

        visitor
    );

    const final_emit_count = emitted_count - killed_count;
    let output_state = create_simulation_state_of_size(final_emit_count);
    let emit_i = 0;
    for (let i = 0; i < emitted_count; ++i) {
        if (!killed[i]) {
            output_state.ids[emit_i] = emit_i;
            output_state.positions[emit_i * 2] = positions_array[i * 2];
            output_state.positions[emit_i * 2 + 1] = positions_array[i * 2 + 1];
            output_state.velocities[emit_i * 2] = 0.0;
            output_state.velocities[emit_i * 2 + 1] = 0.0;
            output_state.colors[emit_i * 3] = 0.1;
            output_state.colors[emit_i * 3 + 1] = 0.2;
            output_state.colors[emit_i * 3 + 2] = 0.9;
            ++emit_i;
        }
    }

    return output_state;
}

//------------------------------------------------------------------------------
function color_from_density(
    count: number,
    target_density: number,
    colors: Float32Array,
    densities: Float32Array
): void {
    for (let i = 0; i < count; ++i) {
        const norm_d = densities[i] / (1.1 * target_density);
        let c = 1.0 - norm_d;
        c = c < 0.0 ? 0.0 : c > 1.0 ? 1.0 : c;
        colors[i * 3] = c;
        colors[i * 3 + 1] = c;
        colors[i * 3 + 2] = 1;
    }
}

function color_from_neighbor_count(
    count: number,
    colors: Float32Array,
    neighbor_counts: Uint8Array
): void {
    for (let i = 0; i < count; ++i) {
        const j = neighbor_counts[i];
        if (j == 0) {
            colors[i * 3] = 1;
            colors[i * 3 + 1] = 0;
            colors[i * 3 + 2] = 0;
        } else if (j == 1) {
            colors[i * 3] = 0;
            colors[i * 3 + 1] = 1;
            colors[i * 3 + 2] = 0;
        } else if (j == 2) {
            colors[i * 3] = 0;
            colors[i * 3 + 1] = 0;
            colors[i * 3 + 2] = 1;
        } else if (j == 3) {
            colors[i * 3] = 1;
            colors[i * 3 + 1] = 1;
            colors[i * 3 + 2] = 0;
        } else if (j == 4) {
            colors[i * 3] = 1;
            colors[i * 3 + 1] = 0;
            colors[i * 3 + 2] = 1;
        } else if (j == 5) {
            colors[i * 3] = 0;
            colors[i * 3 + 1] = 1;
            colors[i * 3 + 2] = 1;
        } else {
            colors[i * 3] = 1;
            colors[i * 3 + 1] = 1;
            colors[i * 3 + 2] = 1;
        }

        // const d = 5 * (neighbor_counts[i] / MAX_NEIGHBORS);
        // colors[3 * i] = d;
        // colors[3 * i + 1] = d;
        // colors[3 * i + 2] = 1.0;
    }
}

//------------------------------------------------------------------------------
// Every one of the regular SPH algorithms, excluding PCISPH, has the ability
// to be called with an adaptive timestep. Getting this to work just right is
// difficult, as many of the algorithms are sensitive to very small timesteps
// I don't know why this is the case, and is worthy of investigation, but
// meanwhile... this std adaptive timestep keeps the timesteps in some discrete
// number of sub steps so that the dt is never too small.
function std_adaptive_time_step(
    sub_step: (
        global_time: number,
        dt: number,
        config: Config,
        fluid_state: State,
        fluid_temp_data: Fluid_temp_data,
        solid_state: State,
        solid_temp_data: Solid_temp_data,
        user_forces: User_forces_function
    ) => void,
    min_sub_step_denom: number,
    max_sub_step_denom: number,
    cfl_factor: number,
    global_time: number,
    config: Config,
    fluid_state: State,
    fluid_temp_data: Fluid_temp_data,
    solid_state: State,
    solid_temp_data: Solid_temp_data,
    user_forces: User_forces_function
): void {
    if (min_sub_step_denom % max_sub_step_denom != 0) {
        throw "INVALID SUB_STEP DENOMS. Min denom must be perfectly divisible by max denom";
    }

    const particle_count = fluid_state.positions.length / 2;
    compute_constant_volumes(
        particle_count,
        config.params.target_density,
        config.mass_per_particle,
        fluid_temp_data.volumes
    );

    const min_sub_time_step = config.params.seconds_per_step / min_sub_step_denom;
    const max_sub_time_step = config.params.seconds_per_step / max_sub_step_denom;
    const sub_step_divisor = min_sub_step_denom / max_sub_step_denom;

    let remaining_sub_steps = min_sub_step_denom;
    let time = global_time;

    const clamp = (a: number, lo: number, hi: number): number => {
        return a < lo ? lo : a > hi ? hi : a;
    };

    while (remaining_sub_steps > 0) {
        const cfl_step = cfl_maximum_time_step(
            particle_count,
            config.params.support,
            cfl_factor,
            fluid_state.velocities
        );

        let num_sub_steps = Math.ceil(cfl_step / min_sub_time_step);
        num_sub_steps = clamp(num_sub_steps, 1, sub_step_divisor);
        num_sub_steps = Math.min(num_sub_steps, remaining_sub_steps);

        const sub_time_step = num_sub_steps * min_sub_time_step;

        //const sub_time_step = 1.0 / 30.0 / 30.0;

        sub_step(
            time,
            sub_time_step,
            config,
            fluid_state,
            fluid_temp_data,
            solid_state,
            solid_temp_data,
            user_forces
        );

        remaining_sub_steps -= num_sub_steps;
        time += sub_time_step;
    }
}

function gravity_forces(
    global_time: number,
    dt: number,
    config: Config,
    fluid_state: State,
    fluid_temp_data: Fluid_temp_data,
    solid_state: State,
    solid_temp_data: Solid_temp_data
): void {
    const count = fluid_state.positions.length / 2;
    for (let i = 0; i < count; ++i) {
        fluid_temp_data.external_forces[2 * i + 0] = 0.0;
        fluid_temp_data.external_forces[2 * i + 1] =
            -config.mass_per_particle * config.params.gravity;
    }
}

//------------------------------------------------------------------------------
export class Simple_simulation {
    config: Config;
    solid_state: State;
    solid_temp_data: Solid_temp_data;
    fluid_state: State;
    fluid_temp_data: Fluid_temp_data;
    accumulated_time: number;

    constructor(params: Parameters) {
        this.config = new Config(params);
        const HHalf = 0.5 * this.config.params.support;

        this.solid_state = world_walls_initial_solid_state(
            this.config.params.support,
            this.config.params.width,
            this.config.params.height,
            "world_walls_initial_solid_state"
        );

        const solid_count = this.solid_state.positions.length / 2;
        this.solid_temp_data = new Solid_temp_data(solid_count);
        compute_neighborhood_parts(this.config, this.solid_state, this.solid_temp_data);
        compute_self_neighborhoods(
            this.config,
            this.solid_temp_data.self_neighborhood,
            this.solid_state,
            this.solid_temp_data
        );
        compute_all_neighborhood_kernels(this.config, this.solid_temp_data.self_neighborhood);
        compute_volumes(
            solid_count,
            this.config.params.support,
            this.solid_temp_data.volumes,
            this.solid_temp_data.self_neighborhood.counts,
            this.solid_temp_data.self_neighborhood.kernels
        );

        compute_constant_volumes(
            solid_count,
            this.config.params.target_density,
            this.config.mass_per_particle,
            this.solid_temp_data.volumes
        );

        color_from_neighbor_count(
            solid_count,
            this.solid_state.colors,
            this.solid_temp_data.self_neighborhood.counts
        );

        this.fluid_state = dam_break_initial_state(
            this.config.params.support,
            this.config.params.width,
            this.config.params.height,
            this.solid_state,
            this.solid_temp_data
        );
        // this.fluid_state = random_box_initial_state(
        //     1000,
        //     this.config.params.width,
        //     this.config.params.height,
        //     "random_box_initial_state",
        //     30.0
        // );

        const fluid_count = this.fluid_state.positions.length / 2;

        // Make a fluid temp state
        this.fluid_temp_data = new Fluid_temp_data(fluid_count);
        compute_all_neighborhoods(
            this.config,
            this.fluid_state,
            this.fluid_temp_data,
            this.solid_state,
            this.solid_temp_data
        );
        compute_all_neighborhood_kernels(this.config, this.fluid_temp_data.self_neighborhood);
        compute_all_neighborhood_kernels(this.config, this.fluid_temp_data.other_neighborhood);

        compute_densities(
            fluid_count,
            this.config.params.support,
            this.config.params.target_density,
            this.fluid_temp_data.densities,
            this.fluid_temp_data.volumes,
            this.fluid_temp_data.self_neighborhood,
            this.solid_temp_data.volumes,
            this.fluid_temp_data.other_neighborhood
        );

        color_from_density(
            fluid_count,
            this.config.params.target_density,
            this.fluid_state.colors,
            this.fluid_temp_data.densities
        );

        this.accumulated_time = 0;
    }

    reset(): void {
        this.fluid_state = dam_break_initial_state(
            this.config.params.support,
            this.config.params.width,
            this.config.params.height,
            this.solid_state,
            this.solid_temp_data
        );
        // this.fluid_state = random_box_initial_state(
        //     1000,
        //     this.config.params.width,
        //     this.config.params.height,
        //     "random_box_initial_state",
        //     30.0
        // );

        const fluid_count = this.fluid_state.positions.length / 2;

        // Make a fluid temp state
        this.fluid_temp_data = new Fluid_temp_data(fluid_count);
        compute_all_neighborhoods(
            this.config,
            this.fluid_state,
            this.fluid_temp_data,
            this.solid_state,
            this.solid_temp_data
        );
        compute_all_neighborhood_kernels(this.config, this.fluid_temp_data.self_neighborhood);
        compute_all_neighborhood_kernels(this.config, this.fluid_temp_data.other_neighborhood);
        compute_densities(
            fluid_count,
            this.config.params.support,
            this.config.params.target_density,
            this.fluid_temp_data.densities,
            this.fluid_temp_data.volumes,
            this.fluid_temp_data.self_neighborhood,
            this.solid_temp_data.volumes,
            this.fluid_temp_data.other_neighborhood
        );

        color_from_density(
            fluid_count,
            this.config.params.target_density,
            this.fluid_state.colors,
            this.fluid_temp_data.densities
        );
        this.accumulated_time = 0;
    }

    step(delta_time: number): void {
        const count = this.fluid_state.positions.length / 2;

        std_adaptive_time_step(
            iisph_pseudo_ap_sub_step,
            60,
            4,
            0.2,
            this.accumulated_time,
            this.config,
            this.fluid_state,
            this.fluid_temp_data,
            this.solid_state,
            this.solid_temp_data,
            gravity_forces
        );

        const solid_count = this.solid_state.positions.length / 2;
        fill_float_array(count * 3, 0.25, this.fluid_state.colors);
        fill_float_array(solid_count * 3, 0.2, this.solid_state.colors);

        for (let i = 0; i < count; i += 1) {

            let fluid_nbhd_count = this.fluid_temp_data.self_neighborhood.counts[i];
            for (let j = 0; j < fluid_nbhd_count; ++j) {
                let nbr_i = this.fluid_temp_data.self_neighborhood.indices[i * MAX_NEIGHBORS + j];

                let kernel = this.fluid_temp_data.self_neighborhood.kernels[i * MAX_NEIGHBORS + j];
                let volume = this.fluid_temp_data.volumes[nbr_i];
                let k = kernel * volume;

                this.fluid_state.colors[nbr_i * 3] += 0.5 * k;
                this.fluid_state.colors[nbr_i * 3 + 1] += 1.0 * k;
                this.fluid_state.colors[nbr_i * 3 + 2] += 0.5 * k;
            }

            let solid_nbhd_count = this.fluid_temp_data.other_neighborhood.counts[i];
            for (let j = 0; j < solid_nbhd_count; ++j) {
                let nbr_i = this.fluid_temp_data.other_neighborhood.indices[i * MAX_NEIGHBORS + j];

                let kernel = this.fluid_temp_data.other_neighborhood.kernels[i * MAX_NEIGHBORS + j];
                let volume = this.solid_temp_data.volumes[nbr_i];
                let k = kernel * volume;
                k *= 5;

                this.solid_state.colors[nbr_i * 3] = 0.5;
                this.solid_state.colors[nbr_i * 3 + 1] += 0.5 * k;
                this.solid_state.colors[nbr_i * 3 + 2] += 1.0 * k;
            }
        }
    }
}