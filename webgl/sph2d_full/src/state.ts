import seedrandom = require("seedrandom");
import { Config } from "./common";
import {
    Block,
    Neighborhood
} from "./neighborhood";

//------------------------------------------------------------------------------
export class State {
    constructor(
        public ids: Int32Array,
        public positions: Float32Array,
        public velocities: Float32Array,
        public colors: Float32Array
    ) {}
}

//------------------------------------------------------------------------------
export function create_simulation_state_of_size(count: number): State {
    return new State(
        new Int32Array(count),
        new Float32Array(2 * count),
        new Float32Array(2 * count),
        new Float32Array(3 * count)
    );
}

export class Base_temp_data {
    grid_coords: Int32Array;
    sorted_index_pairs: Uint32Array;
    block_map: Map<number, Block>;
    self_neighborhood: Neighborhood;
    volumes: Float32Array;

    constructor(public readonly max_count: number) {
        this.grid_coords = new Int32Array(max_count * 2);
        this.sorted_index_pairs = new Uint32Array(max_count * 2);
        this.block_map = new Map<number, Block>();
        this.self_neighborhood = new Neighborhood(max_count);
        this.volumes = new Float32Array(max_count);
    }
}

export class Fluid_temp_data extends Base_temp_data {
    other_neighborhood: Neighborhood;
    constructor(max_count: number) {
        super(max_count);
        this.other_neighborhood = new Neighborhood(max_count);
    }
}

export class Solid_temp_data extends Base_temp_data {
    constructor(max_count: number) {
        super(max_count);
    }
}

//------------------------------------------------------------------------------
export function random_box_initial_state(
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