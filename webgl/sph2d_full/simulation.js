"use strict";

import "./seedrandom.js";

//------------------------------------------------------------------------------
class State {
    constructor(ids, positions, velocities, colors) {
        this.ids = ids;
        this.positions = positions;
        this.velocities = velocities;
        this.colors = colors;
    }
}

//------------------------------------------------------------------------------
function create_simulation_state(count) {
    return new State(
        new Uint32Array(count),
        new Float32Array(2 * count),
        new Float32Array(2 * count),
        new Float32Array(3 * count)
    );
}

//------------------------------------------------------------------------------
function random_box_initial_state(
    state,
    {
        width = 512,
        height = 512,
        seed = "random_box_initial_state",
        speed_gain = 30.0,
    } = {}
) {
    const POINT_COUNT = state.positions.length / 2;
    const my_rng = new Math.seedrandom(seed);
    for (let i = 0; i < POINT_COUNT; ++i) {
        state.ids[i] = i;

        // Random point in the width/height box
        state.positions[2 * i] = width * my_rng();
        state.positions[2 * i + 1] = height * my_rng();

        const speed = speed_gain * (1.0 + 3.0 * my_rng());
        const angle = 2.0 * Math.PI * my_rng();
        state.velocities[2 * i] = speed * Math.cos(angle);
        state.velocities[2 * i + 1] = speed * Math.sin(angle);
    }
    for (let i = 0; i < 3 * POINT_COUNT; ++i) {
        state.colors[i] = my_rng();
    }
}

const MAX_NEIGHBORS = 32;

class Neighborhood {
    constructor(max_count) {
        this.max_count = max_count;

        this.counts = new Uint8Array(max_count);
        this.indices = new Uint32Array(max_count * MAX_NEIGHBORS);
        this.distances = new Float32Array(max_count * MAX_NEIGHBORS);
        this.vectors_to = new Float32Array(max_count * 2 * MAX_NEIGHBORS);
        this.kernels = new Float32Array(max_count * MAX_NEIGHBORS);
        this.kernel_gradients = new Float32Array(max_count * 2 * MAX_NEIGHBORS);
    }
}

class Fluid_temp_data {
    constructor(max_count) {
        this.max_count = max_count;

        this.grid_coords = new Int32Array(max_count * 2);
        this.index_pairs = new Uint32Array(max_count * 2);
        this.block_map = new Map();

        this.fluid_neighborhood = new Neighborhood(max_count, 2);
        this.solid_neighborhood = new Neighborhood(max_count, 2);
    }
}

class Solid_temp_data {
    constructor(max_count) {
        this.max_count = max_count;

        this.grid_coords = new Int32Array(max_count * 2);
        this.index_pairs = new Uint32Array(max_count * 2);
        this.block_map = new Map();
    }
}

//------------------------------------------------------------------------------
function compute_grid_coords(
    count,
    dimension,
    cell_size,
    grid_coords,
    positions
) {
    for (let i = 0; i < count * dimension; ++i) {
        grid_coords[i] = Math.floor(positions[i] / cell_size);
    }
}

//------------------------------------------------------------------------------
// Morton lookup tables.
// Based on
// http://graphics.stanford.edu/~seander/bithacks.html#InterleaveTableLookup
var X = [0, 1],
    Y = [0, 2];
for (let i = 4; i < 0xffff; i <<= 2) {
    for (let j = 0, l = X.length; j < l; ++j) {
        X.push(X[j] | i);
        Y.push((X[j] | i) << 1);
    }
}

// Only works for 24 bit input numbers (up to 16777215).
function morton(x, y) {
    return (
        (Y[y & 0xff] | X[x & 0xff]) +
        (Y[(y >> 8) & 0xff] | X[(x >> 8) & 0xff]) * 0x10000 +
        (Y[(y >> 16) & 0xff] | X[(x >> 16) & 0xff]) * 0x100000000
    );
}

//------------------------------------------------------------------------------
function create_index_pairs(count, dimension, index_pairs, grid_coords) {
    if (dimension != 2) {
        throw "We only support 2-dimensional z-indices at the moment";
    }

    for (let i = 0; i < count; ++i) {
        index_pairs[i * 2] = morton(grid_coords[i * 2], grid_coords[i * 2 + 1]);
        index_pairs[i * 2 + 1] = i;
    }
}

//------------------------------------------------------------------------------
function index_pairs_less_than_elements(a0, a1, b0, b1) {
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

function index_pairs_less_than_pivot(index_pairs, a, pivot0, pivot1) {
    const a0 = index_pairs[a * 2];
    const a1 = index_pairs[a * 2 + 1];
    // console.log("index_pairs_less_than_pivot", a0, a1, pivot0, pivot1);
    return index_pairs_less_than_elements(a0, a1, pivot0, pivot1);
}

function index_pairs_greater_than_pivot(index_pairs, a, pivot0, pivot1) {
    const a0 = index_pairs[a * 2];
    const a1 = index_pairs[a * 2 + 1];
    return index_pairs_less_than_elements(pivot0, pivot1, a0, a1);
}

function index_pairs_less_than(index_pairs, a, b) {
    const a0 = index_pairs[a * 2];
    const a1 = index_pairs[a * 2 + 1];
    const b0 = index_pairs[b * 2];
    const b1 = index_pairs[b * 2 + 1];
    return index_pairs_less_than_elements(a0, a1, b0, b1);
}

function index_pairs_swap(index_pairs, a, b) {
    // console.log("swap: ", a, ", ", b);
    let tmp0 = index_pairs[a * 2];
    let tmp1 = index_pairs[a * 2 + 1];
    index_pairs[a * 2] = index_pairs[b * 2];
    index_pairs[a * 2 + 1] = index_pairs[b * 2 + 1];
    index_pairs[b * 2] = tmp0;
    index_pairs[b * 2 + 1] = tmp1;
}

function index_pairs_partition(index_pairs, lo, hi) {
    // let mid = Math.floor((lo + hi) / 2);

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

    // const pivot0 = index_pairs[mid * 2];
    // const pivot1 = index_pairs[mid * 2 + 1];
    const pivot0 = index_pairs[hi * 2];
    const pivot1 = index_pairs[hi * 2 + 1];

    // console.log("lo, hi: ", lo, ", ", hi, ", pivot: ", pivot0, ", ", pivot1);

    let i = lo - 1;
    let j = hi + 1;
    while (true) {
        do {
            i = i + 1;
            // console.log("i: ", i);
        } while (index_pairs_less_than_pivot(index_pairs, i, pivot0, pivot1));

        do {
            j = j - 1;
            // console.log("j: ", j);
        } while (
            index_pairs_greater_than_pivot(index_pairs, j, pivot0, pivot1)
        );

        if (i >= j) {
            return j;
        }

        index_pairs_swap(index_pairs, i, j);
    }
}

function index_pairs_quicksort(index_pairs, lo, hi) {
    if (lo < hi) {
        let p = index_pairs_partition(index_pairs, lo, hi);
        index_pairs_quicksort(index_pairs, lo, p);
        index_pairs_quicksort(index_pairs, p + 1, hi);
    }
}

function sort_index_pairs(count, index_pairs) {
    index_pairs_quicksort(index_pairs, 0, count - 1);
    // for (let i = 0; i < count; ++i) {
    //     console.log(
    //         "[" +
    //             index_pairs[i * 2].toString() +
    //             ", " +
    //             index_pairs[i * 2 + 1].toString() +
    //             "]"
    //     );
    // }
}

//------------------------------------------------------------------------------
class Block {
    constructor(begin, end) {
        this.begin = begin;
        this.end = end;
    }
}

//------------------------------------------------------------------------------
function compute_block_map(count, map, sorted_index_pairs) {
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
    constructor(index, distance, vector_to_x, vector_to_y) {
        this.index = index;
        this.distance = distance;
        this.vector_to_x = vector_to_x;
        this.vector_to_y = vector_to_y;
    }
}

function neighbor_compare(a, b) {
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

function create_neighborhoods(
    count,
    max_distance,

    neighbor_counts,
    neighbor_indices,
    neighbor_distances,
    neighbor_vectors_to,

    self_positions,
    self_grid_coords,

    other_positions,
    other_sorted_index_pairs,
    other_block_map
) {
    let neighbors = [];

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
        if (nbhd_count < neighbors.length) {
            console.log("Warning, overfull neighborhood at position: ", i);
        }

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
function simple_simulation_step(
    state,
    { delta_time = 1.0, radius = 1.0, width = 512, height = 512 } = {}
) {
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
function estimate_solid_box_emission_count(box_min, box_max, radius) {
    const size = [box_max[0] - box_min[0], box_max[1] - box_min[1]];
    const num_x = Math.ceil(size[0] / (2 * radius));
    const num_y = Math.ceil(size[1] / (Math.sqrt(3.0) * radius));
    return num_x * num_y;
}

//------------------------------------------------------------------------------
function emit_solid_box(positions, box_min, box_max, radius) {
    let push_x = false;
    const dx = 2 * radius;
    const dy = Math.sqrt(3.0) * radius;

    const start_point = [box_min[0] + radius, box_min[1] + radius];
    const point = [start_point[0], start_point[1]];
    const stop_max = [box_max[0] - radius, box_max[1] - radius];
    for (; point[1] <= stop_max[1]; point[1] += dy) {
        // Do a row.
        point[0] = push_x ? start_point[0] + radius : start_point[0];
        for (; point[0] <= stop_max[0]; point[0] += dx) {
            positions.push(point[0]);
            positions.push(point[1]);
        }
        push_x = !push_x;
    }
}

//------------------------------------------------------------------------------
function world_walls_initial_solid_state({
    support: H,
    width = 512,
    height = 512,
    seed = "world_walls_initial_solid_state",
}) {
    const Hhalf = H / 2;
    const R = 0.99 * Hhalf;
    const border_size = H * 3;

    const positions_array = [];

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

    const tiny_block = (unit, roff) => {
        const r = 0.04 + Math.abs(roff);
        emit_solid_box(
            positions_array,
            [width * (unit[0] - r), height * (unit[1] - r)],
            [width * (unit[0] + r), height * (unit[1] + r)],
            R
        );
    };

    const my_rng = new Math.seedrandom(seed);
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
    const ids = new Uint32Array(count);
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
// Call a visitor for each pair of particles within a certain distance
function visit_neighbor_pairs(
    positions_a,
    grid_coords_a,
    sorted_index_pairs_a,
    block_map_a,

    positions_b,
    sorted_index_pairs_b,
    block_map_b,

    max_distance,

    visitor
) {
    // Loop over the blocks of A
    for (let block_a of block_map_a.values()) {
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
    }
}

//------------------------------------------------------------------------------
function dam_break_initial_state(
    { support = 0.01, width = 512, height = 512 },
    solid_state,
    solid_temp_data
) {
    const H = support;
    const Hhalf = H * 0.5;
    const R = 0.99 * Hhalf;
    const extent = 2.0 * support;

    const box_min = new Float32Array([0.0, 0.0]);
    const box_max = new Float32Array([width / 2.0, (3.25 * height) / 4.0]);

    // const max_count = 100 + estimate_solid_box_emission_count(box_min, box_max,
    // R); let state = create_simulation_state(max_count);

    // emit all the points at first.
    const positions_array = [];
    emit_solid_box(positions_array, box_min, box_max, R);

    const emitted_count = positions_array.length / 2;
    let grid_coords = new Int32Array(emitted_count * 2);
    compute_grid_coords(emitted_count, 2, extent, grid_coords, positions_array);

    let index_pairs = new Uint32Array(emitted_count * 2);
    create_index_pairs(emitted_count, 2, index_pairs, grid_coords);
    sort_index_pairs(emitted_count, index_pairs);

    let block_map = new Map();
    compute_block_map(emitted_count, block_map, index_pairs);

    let killed = new Uint8Array(emitted_count);
    let killed_count = 0;
    killed.fill(0);
    let visitor = (
        index_a,
        index_b,
        pos_a_x,
        pos_a_y,
        pos_b_x,
        pos_b_y,
        dx,
        dy,
        distance
    ) => {
        if (killed[index_a] == 0) {
            killed[index_a] = 1;
            ++killed_count;
        }
    };

    visit_neighbor_pairs(
        positions_array,
        grid_coords,
        index_pairs,
        block_map,

        solid_state.positions,
        solid_temp_data.index_pairs,
        solid_temp_data.block_map,

        2.0 * R,

        visitor
    );

    const final_emit_count = emitted_count - killed_count;
    let output_state = create_simulation_state(final_emit_count);
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
function color_from_neighbor_count(count, colors, neighbor_counts) {
    for (let i = 0; i < count; ++i) {
        const d = 2.5 * (neighbor_counts[i] / MAX_NEIGHBORS);
        colors[3 * i] = d;
        colors[3 * i + 1] = d;
        colors[3 * i + 2] = 1.0;
    }
}

//------------------------------------------------------------------------------
export class Simple_simulation {
    constructor({ count = 1000, support = 0.01, width = 512, height = 512 }) {
        this.support = support;
        const HHalf = 0.5 * support;
        this.radius = 0.99 * HHalf;
        this.extent = 2.0 * support;
        this.width = width;
        this.height = height;

        // Make fluid state
        // this.fluid_state = create_simulation_state(count);
        // random_box_initial_state(this.fluid_state, {
        //     width: width,
        //     height: height,
        // });
        // this.fluid_state =
        // Make solid state
        this.solid_state = world_walls_initial_solid_state({
            support: support,
            width: width,
            height: height,
        });

        const solid_count = this.solid_state.positions.length / 2;

        this.solid_temp_data = new Solid_temp_data(count * 5);

        compute_grid_coords(
            solid_count,
            2,
            this.extent,
            this.solid_temp_data.grid_coords,
            this.solid_state.positions
        );
        create_index_pairs(
            solid_count,
            2,
            this.solid_temp_data.index_pairs,
            this.solid_temp_data.grid_coords
        );
        sort_index_pairs(solid_count, this.solid_temp_data.index_pairs);

        compute_block_map(
            solid_count,
            this.solid_temp_data.block_map,
            this.solid_temp_data.index_pairs
        );

        // this.fluid_state = dam_break_initial_state(
        //     {
        //         support: support,
        //         width: width,
        //         height: height,
        //     },
        //     this.solid_state,
        //     this.solid_temp_data
        // );
        this.fluid_state = create_simulation_state(1000);
        random_box_initial_state(this.fluid_state, {
            width: width,
            height: height,
        });

        let fluid_count = this.fluid_state.positions.length / 2;

        this.accumulated_time = 0.0;

        // Make a fluid temp state
        this.fluid_temp_data = new Fluid_temp_data(fluid_count * 5, 2);
    }

    step(delta_time) {
        simple_simulation_step(this.fluid_state, {
            delta_time: delta_time,
            radius: this.radius,
            width: this.width,
            height: this.height,
        });

        const count = this.fluid_state.positions.length / 2;

        compute_grid_coords(
            count,
            2,
            this.extent,
            this.fluid_temp_data.grid_coords,
            this.fluid_state.positions
        );
        create_index_pairs(
            count,
            2,
            this.fluid_temp_data.index_pairs,
            this.fluid_temp_data.grid_coords
        );
        sort_index_pairs(count, this.fluid_temp_data.index_pairs);

        compute_block_map(
            count,
            this.fluid_temp_data.block_map,
            this.fluid_temp_data.index_pairs
        );

        create_neighborhoods(
            count,
            this.extent,

            this.fluid_temp_data.fluid_neighborhood.counts,
            this.fluid_temp_data.fluid_neighborhood.indices,
            this.fluid_temp_data.fluid_neighborhood.distances,
            this.fluid_temp_data.fluid_neighborhood.vectors_to,

            this.fluid_state.positions,
            this.fluid_temp_data.grid_coords,
            this.fluid_state.positions,
            this.fluid_temp_data.index_pairs,
            this.fluid_temp_data.block_map
        );

        color_from_neighbor_count(
            count,
            this.fluid_state.colors,
            this.fluid_temp_data.fluid_neighborhood.counts
        );
    }
}