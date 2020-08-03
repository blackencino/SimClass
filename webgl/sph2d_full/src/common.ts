


//------------------------------------------------------------------------------
export function is_safe_divide(numer: number, denom: number): boolean {
    if (denom == 0) {
        return false;
    } else if (Math.abs(denom) >= 1) {
        return true;
    } else {
        const DENOM_EPSILON = 1.0e-20; // Number.MIN_VALUE
        const mr = Math.abs(denom) / DENOM_EPSILON;
        if (mr > Math.abs(numer)) {
            return true;
        } else {
            return false;
        }
    }
}

export function safe_divide_or(
    numer: number,
    denom: number,
    instead: number
): number {
    if (denom == 0) {
        return instead;
    } else if (Math.abs(denom) >= 1) {
        return numer / denom;
    } else {
        const DENOM_EPSILON = 1.0e-20; // Number.MIN_VALUE
        const mr = Math.abs(denom) / DENOM_EPSILON;
        if (mr > Math.abs(numer)) {
            return numer / denom;
        } else {
            return instead;
        }
    }
}

//------------------------------------------------------------------------------
export function compute_grid_coords(
    count: number,
    cell_size: number,
    grid_coords: Int32Array,
    positions: Float32Array
): void {
    for (let i = 0; i < count * 2; ++i) {
        grid_coords[i] = Math.floor(positions[i] / cell_size);
    }
}

//------------------------------------------------------------------------------
// Morton lookup tables.
// Based on
// http://graphics.stanford.edu/~seander/bithacks.html#InterleaveTableLookup
const X: Array<number> = [0, 1];
const Y: Array<number> = [0, 2];
for (let i = 4; i < 0xffff; i <<= 2) {
    for (let j = 0, l = X.length; j < l; ++j) {
        X.push(X[j] | i);
        Y.push((X[j] | i) << 1);
    }
}

// Only works for 24 bit input numbers (up to 16777215).

// clamping this to 32 bits.
// export function morton(x: number, y: number): number {
//     return (
//         (Y[y & 0xff] | X[x & 0xff]) +
//         (Y[(y >> 8) & 0xff] | X[(x >> 8) & 0xff]) * 0x10000 +
//         (Y[(y >> 16) & 0xff] | X[(x >> 16) & 0xff]) * 0x100000000
//     );
// }
const _clamped_index = new Uint32Array(2);
export function morton(x: number, y: number): number {
    _clamped_index[0] = (
        (Y[y & 0xff] | X[x & 0xff]) +
        (Y[(y >> 8) & 0xff] | X[(x >> 8) & 0xff]) * 0x10000 //+
        //(Y[(y >> 16) & 0xff] | X[(x >> 16) & 0xff]) * 0x100000000
    );
    _clamped_index[1] = 0;
    return _clamped_index[0];
}

export function sqr(x: number): number {
    return x * x;
}

//------------------------------------------------------------------------------
function kernel_sph_wa(q: number): number {
    return 1.0 + q * q * (q * (3.0 / 4.0) - 3.0 / 2.0);
}

function kernel_sph_wb(q_in: number): number {
    const q = 2.0 - q_in;
    return (q * q * q) / 4.0;
}

export function kernel_sph_w(r: number, h: number): number {
    const q = r / h;
    const k = h * h * h * Math.PI;
    if (q <= 1) {
        return kernel_sph_wa(q) / k;
    } else if (q <= 2) {
        return kernel_sph_wb(q) / k;
    } else {
        return 0.0;
    }
}

function kernel_sph_dwa(r: number, h: number) {
    return (-3.0 * r) / (h * h) + ((9.0 / 4.0) * r * r) / (h * h * h);
}

function kernel_sph_dwb(r: number, h: number) {
    return (-3.0 * sqr(2.0 - r / h)) / (4.0 * h);
}

function kernel_sph_dw(r: number, h: number) {
    const q = r / h;
    const k = h * h * h * Math.PI;
    if (q <= 1) {
        return kernel_sph_dwa(r, h) / k;
    } else if (q <= 2) {
        return kernel_sph_dwb(r, h) / k;
    } else {
        return 0.0;
    }
}

export function kernel_sph_gradw(
    dp_x: number,
    dp_y: number,
    h: number
): Array<number> {
    const r = Math.hypot(dp_x, dp_y);
    if (is_safe_divide(dp_x, r) && is_safe_divide(dp_y, r)) {
        const dw = kernel_sph_dw(r, h);
        return [(dp_x / r) * dw, (dp_y / r) * dw];
    } else {
        return [0.0, 0.0];
    }
}


//------------------------------------------------------------------------------
export function estimate_solid_box_emission_count(
    box_min: Array<number> | Float32Array,
    box_max: Array<number> | Float32Array,
    radius: number
): number {
    const size = [box_max[0] - box_min[0], box_max[1] - box_min[1]];
    const num_x = Math.ceil(size[0] / (2 * radius));
    const num_y = Math.ceil(size[1] / (Math.sqrt(3.0) * radius));
    return num_x * num_y;
}

//------------------------------------------------------------------------------
export function emit_solid_box(
    positions: Array<number>,
    box_min: Array<number> | Float32Array,
    box_max: Array<number> | Float32Array,
    radius: number
): void {
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

export type Parameters = {
    seed: number;
    seconds_per_step: number;
    sub_steps: number;
    width: number;
    height: number;
    support: number;
    gravity: number;
    viscosity: number;
    target_density: number;
    iisph: {
        max_pressure_iterations: number;
        error_average_threshold: number;
        error_max_threshold: number;
        omega: number;
    };
    num_batch_steps: number;
}

export function default_parameters(): Parameters {
    return {
        seed: 1,
        seconds_per_step: 1.0/30.0,
        sub_steps: 30,
        width: 1.0,
        height: 1.0,
        support: 0.025,
        gravity: 9.81,
        viscosity: 0.01,
        target_density: 1000.0,
        iisph: {
            max_pressure_iterations: 30,
            error_average_threshold: 0.005,
            error_max_threshold: 0.05,
            omega: 0.5
        },
        num_batch_steps: 100
    };
}

export class Config {
    draw_radius: number;
    target_concentration: number;
    mass_per_particle: number;
    constructor(public readonly params: Parameters) {
        const H = params.support;
        const Hhalf = H * 0.5;
        const R = 0.99 * Hhalf;
        this.draw_radius = R;

        let muchness0 = kernel_sph_w(0.0, H);

        let positions: Array<number> = [];
        emit_solid_box(positions,
                       [-15.0 * R, -15.0*R],
                       [15.0*R, 15.0*R],
                       R);
        const emitted_count = positions.length / 2;

        // Find the emitted point closest to zero
        let best_point = [positions[0], positions[1]];
        let best_r2 = Math.hypot(best_point[0], best_point[1]);
        let best_index = 0;
        for (let i = 1; i < emitted_count; ++i) {
            let point = [positions[i*2], positions[i*2+1]];
            const r2 = Math.hypot(point[0], point[1]);
            if (r2 < best_r2) {
                best_point = point;
                best_r2 = r2;
                best_index = i;
            }
        }
        console.log("Best index: ", best_index);
        console.log("Best r2: ", best_r2);

        let num_neighbors = 0;
        for (let i = 0; i < emitted_count; ++i) {
            if (i == best_index) {
                continue;
            }

            const point_x = positions[i*2];
            const point_y = positions[i*2+1];

            const r = Math.hypot(point_x - best_point[0],
                                 point_y - best_point[1]);

            const w = kernel_sph_w(r, H);
            if (w > 0.0) {
                muchness0 += w;
                ++num_neighbors;
            }
        }

        this.target_concentration = muchness0;
        this.mass_per_particle = this.params.target_density / muchness0;
        console.log("radius: ", this.draw_radius);
        console.log("mass: ", this.mass_per_particle);
        console.log("support: ", this.params.support);
        console.log("target concentration: ", this.target_concentration);
        console.log("num_neighbors: ", num_neighbors);
        console.log("best_point: ", best_point[0], best_point[1]);
    }
};
