import seedrandom = require("seedrandom");

import { emit_solid_box, compute_grid_coords, morton,
Parameters, Config } from "./common";
import { MAX_NEIGHBORS, Neighborhood, Block,
    compute_sorted_index_pairs,
    compute_block_map,
    compute_neighborhood_parts,
    compute_self_neighborhoods,
    compute_other_neighborhoods,
    compute_all_neighborhood_kernels,
    compute_all_neighborhoods,
    visit_neighbor_pairs } from "./neighborhood";

import { State, Fluid_temp_data, Solid_temp_data, create_simulation_state_of_size,
    random_box_initial_state } from "./state";

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
            state.positions[i] =
                width - (state.positions[i] + radius - width) - radius;
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
            state.positions[i] =
                height - (state.positions[i] + radius - height) - radius;
            state.velocities[i] = -Math.abs(state.velocities[i]);
        }
    }
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
    emit_solid_box(
        positions_array,
        [-border_size, -border_size],
        [0.0, height + border_size],
        R
    );

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

    const tiny_block = (
        unit: Array<number> | Float32Array,
        roff: number
    ): void => {
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
            tiny_block(
                [bxf + my_rng_dist(), byf + my_rng_dist()],
                my_rng_dist()
            );
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
    const box_max = [width / 2.0, (3.25 * height) / 4.0];

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
function color_from_neighbor_count(
    count: number,
    colors: Float32Array,
    neighbor_counts: Uint8Array
): void {
    for (let i = 0; i < count; ++i) {
        const d = 2.5 * (neighbor_counts[i] / MAX_NEIGHBORS);
        colors[3 * i] = d;
        colors[3 * i + 1] = d;
        colors[3 * i + 2] = 1.0;
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

    constructor(
        params: Parameters
    ) {
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
        compute_neighborhood_parts(this.config,
                                   this.solid_state,
                                   this.solid_temp_data);
        compute_self_neighborhoods(this.config,
                                   this.solid_temp_data.self_neighborhood,
                                   this.solid_state,
                                   this.solid_temp_data);
        compute_all_neighborhood_kernels(this.config,
                                         this.solid_temp_data.self_neighborhood);


        // this.fluid_state = dam_break_initial_state(
        //     {
        //         support: support,
        //         width: width,
        //         height: height,
        //     },
        //     this.solid_state,
        //     this.solid_temp_data
        // );
        this.fluid_state = random_box_initial_state(1000,
            this.config.params.width,
            this.config.params.height,
            "random_box_initial_state",
            30.0
        );

        const fluid_count = this.fluid_state.positions.length / 2;

        // Make a fluid temp state
        this.fluid_temp_data = new Fluid_temp_data(fluid_count * 5);
        compute_all_neighborhoods(this.config,
                                  this.fluid_state,
                                  this.fluid_temp_data,
                                  this.solid_state,
                                  this.solid_temp_data);
        compute_all_neighborhood_kernels(this.config, this.fluid_temp_data.self_neighborhood);
        compute_all_neighborhood_kernels(this.config, this.fluid_temp_data.other_neighborhood);
        this.accumulated_time = 0;
    }

    step(delta_time: number): void {
        simple_simulation_step(
            this.fluid_state,
            delta_time,
            this.config.draw_radius,
            this.config.params.width,
            this.config.params.height
        );

        const count = this.fluid_state.positions.length / 2;
        compute_all_neighborhoods(this.config,
                                  this.fluid_state,
                                  this.fluid_temp_data,
                                  this.solid_state,
                                  this.solid_temp_data);
        compute_all_neighborhood_kernels(this.config, this.fluid_temp_data.self_neighborhood);
        compute_all_neighborhood_kernels(this.config, this.fluid_temp_data.other_neighborhood);

        color_from_neighbor_count(
            count,
            this.fluid_state.colors,
            this.fluid_temp_data.self_neighborhood.counts
        );
    }
}
