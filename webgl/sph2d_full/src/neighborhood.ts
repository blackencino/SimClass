import {
    morton,
    compute_grid_coords,
    kernel_sph_w,
    kernel_sph_gradw,
    Config,
} from "./common";

export const MAX_NEIGHBORS: number = 32;

export class Neighborhood {
    counts: Uint8Array;
    indices: Uint32Array;
    distances: Float32Array;
    vectors_to: Float32Array;
    kernels: Float32Array;
    kernel_gradients: Float32Array;

    constructor(public readonly max_count: number) {
        this.counts = new Uint8Array(max_count);
        this.indices = new Uint32Array(max_count * MAX_NEIGHBORS);
        this.distances = new Float32Array(max_count * MAX_NEIGHBORS);
        this.vectors_to = new Float32Array(max_count * 2 * MAX_NEIGHBORS);
        this.kernels = new Float32Array(max_count * MAX_NEIGHBORS);
        this.kernel_gradients = new Float32Array(max_count * 2 * MAX_NEIGHBORS);
    }
}

export class Block {
    constructor(public begin: number, public end: number) {}
}

//------------------------------------------------------------------------------
function create_index_pairs(
    count: number,
    index_pairs: Uint32Array,
    grid_coords: Int32Array
): void {
    for (let i = 0; i < count; ++i) {
        index_pairs[i * 2] = morton(grid_coords[i * 2], grid_coords[i * 2 + 1]);
        index_pairs[i * 2 + 1] = i;
    }
}

//------------------------------------------------------------------------------
function index_pairs_less_than_elements(
    a0: number,
    a1: number,
    b0: number,
    b1: number
): boolean {
    if (a0 < b0) {
        return true;
    } else if (a0 > b0) {
        return false;
    } else if (a1 < b1) {
        return true;
    } else {
        return false;
    }
    // return a0 < b0 ? true : a0 > b0 ? false : a1 < b1;
}

function index_pairs_less_than_pivot(
    index_pairs: Uint32Array,
    a: number,
    pivot0: number,
    pivot1: number
): boolean {
    const a0 = index_pairs[a * 2];
    const a1 = index_pairs[a * 2 + 1];
    return index_pairs_less_than_elements(a0, a1, pivot0, pivot1);
}

function index_pairs_greater_than_pivot(
    index_pairs: Uint32Array,
    a: number,
    pivot0: number,
    pivot1: number
): boolean {
    const a0 = index_pairs[a * 2];
    const a1 = index_pairs[a * 2 + 1];
    return index_pairs_less_than_elements(pivot0, pivot1, a0, a1);
}

function index_pairs_less_than(
    index_pairs: Uint32Array,
    a: number,
    b: number
): boolean {
    const a0 = index_pairs[a * 2];
    const a1 = index_pairs[a * 2 + 1];
    const b0 = index_pairs[b * 2];
    const b1 = index_pairs[b * 2 + 1];
    return index_pairs_less_than_elements(a0, a1, b0, b1);
}

function index_pairs_swap(
    index_pairs: Uint32Array,
    a: number,
    b: number
): void {
    let tmp0 = index_pairs[a * 2];
    let tmp1 = index_pairs[a * 2 + 1];
    index_pairs[a * 2] = index_pairs[b * 2];
    index_pairs[a * 2 + 1] = index_pairs[b * 2 + 1];
    index_pairs[b * 2] = tmp0;
    index_pairs[b * 2 + 1] = tmp1;
}

function index_pairs_partition(
    index_pairs: Uint32Array,
    lo: number,
    hi: number
): number {
    let mid = Math.floor((lo + hi) / 2);

    if (index_pairs_less_than(index_pairs, mid, lo)) {
        index_pairs_swap(index_pairs, lo, mid);
    }
    if (index_pairs_less_than(index_pairs, hi, lo)) {
        index_pairs_swap(index_pairs, lo, hi);
    }
    if (index_pairs_less_than(index_pairs, mid, hi)) {
        index_pairs_swap(index_pairs, mid, hi);
    }

    const pivot0 = index_pairs[hi * 2];
    const pivot1 = index_pairs[hi * 2 + 1];

    let i = lo - 1;
    let j = hi + 1;
    while (true) {
        do {
            i = i + 1;
        } while (index_pairs_less_than_pivot(index_pairs, i, pivot0, pivot1));

        do {
            j = j - 1;
        } while (
            index_pairs_greater_than_pivot(index_pairs, j, pivot0, pivot1)
        );

        if (i >= j) {
            return j;
        }

        index_pairs_swap(index_pairs, i, j);
    }
}

// The c standard library implementation of qsort uses 4 as the cutoff for the switch
// to insertion sort.
const SORT_INSERTION_THRESH: number = 4;

function index_pairs_insertionsort(
    index_pairs: Uint32Array,
    lo: number,
    hi: number
): void {
    for (let i = lo + 1; i <= hi; ++i) {
        let pivot0 = index_pairs[i * 2];
        let pivot1 = index_pairs[i * 2 + 1];
        let j: number;
        for (
            j = i - 1;
            j >= lo &&
            index_pairs_greater_than_pivot(index_pairs, j, pivot0, pivot1);
            --j
        ) {
            index_pairs[(j + 1) * 2] = index_pairs[j * 2];
            index_pairs[(j + 1) * 2 + 1] = index_pairs[j * 2 + 1];
        }

        index_pairs[(j + 1) * 2] = pivot0;
        index_pairs[(j + 1) * 2 + 1] = pivot1;
    }
}

function index_pairs_quicksort(
    index_pairs: Uint32Array,
    lo: number,
    hi: number
): void {
    const count = hi - lo + 1;
    if (count < 2) {
        return;
    } else if (count < SORT_INSERTION_THRESH) {
        index_pairs_insertionsort(index_pairs, lo, hi);
    } else {
        let p = index_pairs_partition(index_pairs, lo, hi);
        index_pairs_quicksort(index_pairs, lo, p);
        index_pairs_quicksort(index_pairs, p + 1, hi);
    }
}

export function compute_sorted_index_pairs(
    count: number,
    sorted_index_pairs: Uint32Array,
    grid_coords: Int32Array
): void {
    create_index_pairs(count, sorted_index_pairs, grid_coords);
    index_pairs_quicksort(sorted_index_pairs, 0, count - 1);
}

//------------------------------------------------------------------------------
export function compute_block_map(
    count: number,
    map: Map<number, Block>,
    sorted_index_pairs: Uint32Array
): void {
    map.clear();

    if (count < 1) {
        return;
    }

    // Loop over the sorted index pairs. Each time the z-index changes
    // from one sorted index pair to the next, begin a new block.
    let prev_z_index = sorted_index_pairs[0];
    let block = new Block(0, count);
    for (let i = 1; i < count; ++i) {
        const z_index = sorted_index_pairs[i * 2];
        if (prev_z_index != z_index) {
            // i is the beginning of a new block.

            // End the old block and add it to the map.
            block.end = i;
            const block_count = block.end - block.begin + 1;
            if (block_count <= 0 || block_count > count) {
                throw "Bad block count: " + block_count.toString();
            }
            map.set(prev_z_index, block);

            // Make a new block
            block = new Block(i, count);

            // Remember this prev z index.
            prev_z_index = z_index;
        }
    }

    // Add the last (or possible only) block.
    const block_count = block.end - block.begin + 1;
    if (block_count <= 0 || block_count > count) {
        throw "Bad block count: " + block_count.toString();
    }
    map.set(prev_z_index, block);
}

class Neighbor {
    constructor(
        public index: number,
        public distance: number,
        public vector_to_x: number,
        public vector_to_y: number
    ) {}
}

function neighbor_compare(a: Neighbor, b: Neighbor): number {
    if (a.distance < b.distance) {
        return -1;
    } else if (a.distance > b.distance) {
        return 1;
    } else if (a.index < b.index) {
        return -1;
    } else if (a.index > b.index) {
        return 1;
    } else {
        return 0;
    }
}

export function compute_neighborhoods(
    count: number,
    max_distance: number,

    neighbor_counts: Uint8Array,
    neighbor_indices: Uint32Array,
    neighbor_distances: Float32Array,
    neighbor_vectors_to: Float32Array,

    self_positions: Float32Array,
    self_grid_coords: Int32Array,

    other_positions: Float32Array,
    other_sorted_index_pairs: Uint32Array,
    other_block_map: Map<number, Block>
): void {
    let neighbors: Array<Neighbor> = [];

    for (let i = 0; i < count; ++i) {
        const pos_x = self_positions[i * 2];
        const pos_y = self_positions[i * 2 + 1];

        const grid_coord_i = self_grid_coords[i * 2];
        const grid_coord_j = self_grid_coords[i * 2 + 1];

        neighbors.length = 0;
        for (let j = grid_coord_j - 1; j <= grid_coord_j + 1; ++j) {
            for (let ii = grid_coord_i - 1; ii <= grid_coord_i + 1; ++ii) {
                const other_z_index = morton(ii, j);
                const found_block = other_block_map.get(other_z_index);
                if (found_block === undefined) {
                    continue;
                }

                // sip == sorted_index_pair
                for (
                    let sip = found_block.begin;
                    sip != found_block.end;
                    ++sip
                ) {
                    const other_i = other_sorted_index_pairs[sip * 2 + 1];

                    // Don't count our self if the positions array is the same.
                    if (other_positions === self_positions && other_i == i) {
                        continue;
                    }

                    const vector_to_x = other_positions[other_i * 2] - pos_x;
                    const vector_to_y =
                        other_positions[other_i * 2 + 1] - pos_y;
                    const distance = Math.hypot(vector_to_x, vector_to_y);
                    if (distance >= max_distance) {
                        continue;
                    }

                    neighbors.push(
                        new Neighbor(
                            other_i,
                            distance,
                            vector_to_x,
                            vector_to_y
                        )
                    );
                }
            }
        }
        const nbhd_count = Math.min(neighbors.length, MAX_NEIGHBORS);
        // if (nbhd_count < neighbors.length) {
        //     console.log("Warning, overfull neighborhood");
        // }

        neighbor_counts[i] = nbhd_count;
        neighbors.sort(neighbor_compare);
        for (let neighbor_i = 0; neighbor_i < nbhd_count; ++neighbor_i) {
            const neighbor = neighbors[neighbor_i];
            neighbor_indices[i * MAX_NEIGHBORS + neighbor_i] = neighbor.index;
            neighbor_distances[i * MAX_NEIGHBORS + neighbor_i] =
                neighbor.distance;
            neighbor_vectors_to[(i * MAX_NEIGHBORS + neighbor_i) * 2] =
                neighbor.vector_to_x;
            neighbor_vectors_to[(i * MAX_NEIGHBORS + neighbor_i) * 2 + 1] =
                neighbor.vector_to_y;
        }
    }
}

//------------------------------------------------------------------------------
export function compute_neighbor_kernels(
    particle_count: number,
    support: number,
    neighbor_kernels: Float32Array,
    neighbor_counts: Uint8Array,
    neighbor_distances: Float32Array
): void {
    for (
        let particle_index = 0;
        particle_index < particle_count;
        ++particle_index
    ) {
        const neighbor_count = neighbor_counts[particle_index];
        for (let j = 0; j < neighbor_count; ++j) {
            neighbor_kernels[particle_index * MAX_NEIGHBORS + j] = kernel_sph_w(
                neighbor_distances[particle_index * MAX_NEIGHBORS + j],
                support
            );
        }
    }
}

export function compute_neighbor_kernel_gradients(
    particle_count: number,
    support: number,
    neighbor_kernel_gradients: Float32Array,
    neighbor_counts: Uint8Array,
    neighbor_vectors_to: Float32Array
) {
    for (
        let particle_index = 0;
        particle_index < particle_count;
        ++particle_index
    ) {
        const neighbor_count = neighbor_counts[particle_index];
        for (let j = 0; j < neighbor_count; ++j) {
            const gradw = kernel_sph_gradw(
                neighbor_vectors_to[(particle_index * MAX_NEIGHBORS + j) * 2],
                neighbor_vectors_to[
                    (particle_index * MAX_NEIGHBORS + j) * 2 + 1
                ],
                support
            );
            neighbor_kernel_gradients[
                (particle_index * MAX_NEIGHBORS + j) * 2
            ] = gradw[0];
            neighbor_kernel_gradients[
                (particle_index * MAX_NEIGHBORS + j) * 2 + 1
            ] = gradw[1];
        }
    }
}

//------------------------------------------------------------------------------
// Call a visitor for each pair of particles within a certain distance
export function visit_neighbor_pairs(
    positions_a: Float32Array,
    grid_coords_a: Int32Array,
    sorted_index_pairs_a: Uint32Array,
    block_map_a: Map<number, Block>,

    positions_b: Float32Array,
    sorted_index_pairs_b: Uint32Array,
    block_map_b: Map<number, Block>,

    max_distance: number,

    visitor: (
        index_a: number,
        index_b: number,
        pos_a_x: number,
        pos_a_y: number,
        pos_b_x: number,
        pos_b_y: number,
        dx: number,
        dy: number,
        distance: number
    ) => void
): void {
    // Loop over the blocks of A
    block_map_a.forEach((block_a: Block): void => {
        const block_begin_index_a = sorted_index_pairs_a[block_a.begin * 2 + 1];
        const grid_coord_a_x = grid_coords_a[block_begin_index_a * 2];
        const grid_coord_a_y = grid_coords_a[block_begin_index_a * 2 + 1];

        for (
            let grid_coord_j = grid_coord_a_y - 1;
            grid_coord_j <= grid_coord_a_y + 1;
            ++grid_coord_j
        ) {
            for (
                let grid_coord_i = grid_coord_a_x - 1;
                grid_coord_i <= grid_coord_a_x + 1;
                ++grid_coord_i
            ) {
                const z_index_b = morton(grid_coord_i, grid_coord_j);

                const block_b = block_map_b.get(z_index_b);
                if (block_b === undefined) {
                    continue;
                }

                for (
                    let block_sip_a = block_a.begin;
                    block_sip_a < block_a.end;
                    ++block_sip_a
                ) {
                    const index_a = sorted_index_pairs_a[block_sip_a * 2 + 1];
                    const pos_a_x = positions_a[index_a * 2];
                    const pos_a_y = positions_a[index_a * 2 + 1];

                    for (
                        let block_sip_b = block_b.begin;
                        block_sip_b < block_b.end;
                        ++block_sip_b
                    ) {
                        const index_b =
                            sorted_index_pairs_b[block_sip_b * 2 + 1];
                        // Don't call visitor on self.
                        if (positions_a === positions_b && index_a == index_b) {
                            continue;
                        }

                        const pos_b_x = positions_b[index_b * 2];
                        const pos_b_y = positions_b[index_b * 2 + 1];

                        const dx = pos_b_x - pos_a_x;
                        const dy = pos_b_y - pos_a_y;

                        const distance = Math.hypot(dx, dy);

                        if (distance < max_distance) {
                            visitor(
                                index_a,
                                index_b,
                                pos_a_x,
                                pos_a_y,
                                pos_b_x,
                                pos_b_y,
                                dx,
                                dy,
                                distance
                            );
                        }
                    }
                }
            }
        }
    });
}

//------------------------------------------------------------------------------
interface Expected_state {
    positions: Float32Array;
}

interface Expected_temp_data {
    grid_coords: Int32Array;
    sorted_index_pairs: Uint32Array;
    block_map: Map<number, Block>;
    self_neighborhood: Neighborhood;
}

//------------------------------------------------------------------------------
export function compute_neighborhood_parts(
    config: Config,
    state: Expected_state,
    temp_data: Expected_temp_data
): void {
    const count = state.positions.length / 2;
    const cell_size = config.params.support * 2.0;
    compute_grid_coords(
        count,
        cell_size,
        temp_data.grid_coords,
        state.positions
    );
    compute_sorted_index_pairs(
        count,
        temp_data.sorted_index_pairs,
        temp_data.grid_coords
    );
    compute_block_map(count, temp_data.block_map, temp_data.sorted_index_pairs);
}

export function compute_self_neighborhoods(
    config: Config,
    neighborhood: Neighborhood,
    state: Expected_state,
    temp_data: Expected_temp_data
): void {
    const count = state.positions.length / 2;
    const cell_size = config.params.support * 2.0;
    compute_neighborhoods(
        count,
        cell_size,
        neighborhood.counts,
        neighborhood.indices,
        neighborhood.distances,
        neighborhood.vectors_to,

        // self
        state.positions,
        temp_data.grid_coords,

        // self
        state.positions,
        temp_data.sorted_index_pairs,
        temp_data.block_map
    );
}

export function compute_other_neighborhoods(
    config: Config,
    neighborhood: Neighborhood,
    state: Expected_state,
    temp_data: Expected_temp_data,
    other_state: Expected_state,
    other_temp_data: Expected_temp_data
): void {
    const count = state.positions.length / 2;
    const cell_size = config.params.support * 2.0;
    compute_neighborhoods(
        count,
        cell_size,
        neighborhood.counts,
        neighborhood.indices,
        neighborhood.distances,
        neighborhood.vectors_to,

        // self
        state.positions,
        temp_data.grid_coords,

        // self
        other_state.positions,
        other_temp_data.sorted_index_pairs,
        other_temp_data.block_map
    );
}

export function compute_all_neighborhood_kernels(
    config: Config,
    neighborhood: Neighborhood
): void {
    const count = neighborhood.counts.length;
    compute_neighbor_kernels(
        count,
        config.params.support,
        neighborhood.kernels,
        neighborhood.counts,
        neighborhood.distances
    );

    compute_neighbor_kernel_gradients(
        count,
        config.params.support,
        neighborhood.kernel_gradients,
        neighborhood.counts,
        neighborhood.vectors_to
    );
}

//------------------------------------------------------------------------------
interface Expected_fluid_temp_data extends Expected_temp_data {
    other_neighborhood: Neighborhood;
}

export function compute_all_neighborhoods(
    config: Config,
    fluid_state: Expected_state,
    fluid_temp_data: Expected_fluid_temp_data,
    solid_state: Expected_state,
    solid_temp_data: Expected_temp_data
): void {
    compute_neighborhood_parts(config, fluid_state, fluid_temp_data);
    compute_self_neighborhoods(
        config,
        fluid_temp_data.self_neighborhood,
        fluid_state,
        fluid_temp_data
    );
    compute_other_neighborhoods(
        config,
        fluid_temp_data.other_neighborhood,
        fluid_state,
        fluid_temp_data,
        solid_state,
        solid_temp_data
    );
}
