import { Config } from "./common";
import { Block, Neighborhood } from "./neighborhood";

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

    densities: Float32Array;
    density_stars: Float32Array;
    external_forces: Float32Array;
    pseudo_pressures: Float32Array;
    pseudo_pressure_displacements: Float32Array;
    pseudo_diagonals: Float32Array;

    constructor(max_count: number) {
        super(max_count);
        this.other_neighborhood = new Neighborhood(max_count);
        this.densities = new Float32Array(max_count);
        this.density_stars = new Float32Array(max_count);
        this.external_forces = new Float32Array(max_count * 2);
        this.pseudo_pressures = new Float32Array(max_count);
        this.pseudo_pressure_displacements = new Float32Array(max_count * 2);
        this.pseudo_diagonals = new Float32Array(max_count);
    }
}

export class Solid_temp_data extends Base_temp_data {
    constructor(max_count: number) {
        super(max_count);
    }
}