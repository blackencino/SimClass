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
import { Config } from "./common";

//------------------------------------------------------------------------------
// needed for a water wheel.
// solid state has a rest position
// solid state has a translation
// solid state has a rotation
// at the beginning of each step we...
// ... compute the transformed positions
// ... compute the transformed velocities
// ... compute the neighborhood parts
// ... compute the solid self neighborhood
//
// to compute rest position we need a...
// ... annular inner radius
// ... annular outer radius
// ... spoke radius
// ... num inner spokes
// ... num outer spokes
// then we compute a bounding box,
//

class Spoke_shape {
    axis_x: number;
    axis_y: number;

    constructor(
        public readonly inner_radius: number,
        public readonly outer_radius: number,
        public readonly thickness: number,
        angle: number
    ) {
        this.axis_x = Math.cos(angle);
        this.axis_y = Math.sin(angle);
    }

    test_point(point_x: number, point_y: number): boolean {
        const proj = point_x * this.axis_x + point_y * this.axis_y;
        if (proj < this.inner_radius || proj > this.outer_radius) {
            return false;
        }

        // get the projected point
        const proj_pt_x = this.axis_x * proj;
        const proj_pt_y = this.axis_y * proj;

        // get the separating vector
        const dp_x = point_x - proj_pt_x;
        const dp_y = point_y - proj_pt_y;

        if (Math.hypot(dp_x, dp_y) <= 0.5 * this.thickness) {
            return true;
        } else {
            return false;
        }
    }
}

class Ring_shape {
    constructor(public readonly inner_radius: number, public readonly outer_radius: number) {}

    test_point(point_x: number, point_y: number): boolean {
        const r = Math.hypot(point_x, point_y);
        if (r >= this.inner_radius && r <= this.outer_radius) {
            return true;
        } else {
            return false;
        }
    }
}

class Water_wheel_shape {
    shapes: Array<Spoke_shape | Ring_shape>;

    constructor(
        wheel_radius: number,
        outer_spoke_radius: number,
        wheel_thickness: number,
        inner_spoke_thickness: number,
        outer_spoke_thickness: number,
        num_inner_spokes: number,
        num_outer_spokes: number
    ) {
        this.shapes = [];
        this.shapes.push(new Ring_shape(wheel_radius - wheel_thickness, wheel_radius));

        for (let i = 0; i < num_inner_spokes; ++i) {
            const angle = (Math.PI * 2.0 * i) / num_inner_spokes;
            this.shapes.push(
                new Spoke_shape(
                    0.0,
                    wheel_radius - 0.5 * wheel_thickness,
                    inner_spoke_thickness,
                    angle
                )
            );
        }

        for (let i = 0; i < num_outer_spokes; ++i) {
            const angle = (Math.PI * 2.0 * i) / num_outer_spokes;
            this.shapes.push(
                new Spoke_shape(
                    wheel_radius - 0.5 * wheel_thickness,
                    outer_spoke_radius,
                    outer_spoke_thickness,
                    angle
                )
            );
        }
    }

    test_point(point_x: number, point_y: number): boolean {
        return this.shapes.some((shape: Spoke_shape | Ring_shape): boolean => {
            return shape.test_point(point_x, point_y);
        });
    }
}

//------------------------------------------------------------------------------
function water_wheel_rest_positions(
    wheel_radius: number,
    outer_spoke_radius: number,
    wheel_thickness: number,
    inner_spoke_thickness: number,
    outer_spoke_thickness: number,
    num_inner_spokes: number,
    num_outer_spokes: number,
    point_radius: number
): Float32Array {
    const wheel = new Water_wheel_shape(
        wheel_radius,
        outer_spoke_radius,
        wheel_thickness,
        inner_spoke_thickness,
        outer_spoke_thickness,
        num_inner_spokes,
        num_outer_spokes
    );

    let push_x = false;
    const dx = 2 * point_radius;
    const dy = Math.sqrt(3.0) * point_radius;

    const max_r = outer_spoke_radius + point_radius;
    const box_min = [-max_r, -max_r];
    const box_max = [max_r, max_r];

    let positions: Array<number> = [];

    const start_point = [box_min[0] + point_radius, box_min[1] + point_radius];
    const point = [start_point[0], start_point[1]];
    const stop_max = [box_max[0] - point_radius, box_max[1] - point_radius];
    for (; point[1] <= stop_max[1]; point[1] += dy) {
        // Do a row.
        point[0] = push_x ? start_point[0] + point_radius : start_point[0];
        for (; point[0] <= stop_max[0]; point[0] += dx) {
            if (wheel.test_point(point[0], point[1])) {
                positions.push(point[0]);
                positions.push(point[1]);
            }
        }
        push_x = !push_x;
    }

    console.log("Rest positions count:", positions.length / 2);

    return new Float32Array(positions);
}

//------------------------------------------------------------------------------
function compute_moving_water_wheel_state(
    count: number,
    angle: number,
    angular_velocity: number,
    origin_x: number,
    origin_y: number,
    state: State,
    rest_positions: Float32Array
) {
    const cos_ang = Math.cos(angle);
    const sin_ang = Math.sin(angle);

    const axis_0_x = cos_ang;
    const axis_0_y = sin_ang;
    const axis_1_x = -sin_ang;
    const axis_1_y = cos_ang;

    for (let i = 0; i < count; ++i) {
        const rest_pos_x = rest_positions[i * 2];
        const rest_pos_y = rest_positions[i * 2 + 1];

        const pos_x = origin_x + rest_pos_x * axis_0_x + rest_pos_y * axis_1_x;
        const pos_y = origin_y + rest_pos_x * axis_0_y + rest_pos_y * axis_1_y;

        const r = Math.hypot(rest_pos_x * rest_pos_y);

        const vel_x = angular_velocity * axis_1_x * r;
        const vel_y = angular_velocity * axis_1_y * r;

        state.positions[i * 2] = pos_x;
        state.positions[i * 2 + 1] = pos_y;
        state.velocities[i * 2] = vel_x;
        state.velocities[i * 2 + 1] = vel_y;

        // state.positions[i*2] = origin_x + rest_pos_x;
        // state.positions[i*2+1] = origin_y + rest_pos_y;
        // state.velocities[i*2] = 0.0;
        // state.velocities[i*2+1] = 0.0;
    }
}

//------------------------------------------------------------------------------
function compute_border(
    width: number,
    height: number,
    border: number,
    point_radius: number
): Float32Array {
    let push_x = false;
    const dx = 2 * point_radius;
    const dy = Math.sqrt(3.0) * point_radius;

    const box_min = [-border, -border];
    const box_max = [width + border, height + border];

    let positions: Array<number> = [];

    const start_point = [box_min[0] + point_radius, box_min[1] + point_radius];
    const point = [start_point[0], start_point[1]];
    const stop_max = [box_max[0] - point_radius, box_max[1] - point_radius];
    for (; point[1] <= stop_max[1]; point[1] += dy) {
        // Do a row.
        point[0] = push_x ? start_point[0] + point_radius : start_point[0];
        for (; point[0] <= stop_max[0]; point[0] += dx) {
            if (!(point[0] > 0 && point[0] < width && point[1] > 0 && point[1] < height)) {
                positions.push(point[0]);
                positions.push(point[1]);
            }
        }
        push_x = !push_x;
    }

    return new Float32Array(positions);
}

//------------------------------------------------------------------------------
export class Water_wheel {
    state: State;
    temp_data: Solid_temp_data;
    wheel_rest_positions: Float32Array;

    constructor(
        public readonly config: Config,
        wheel_radius: number,
        outer_spoke_radius: number,
        wheel_thickness: number,
        inner_spoke_thickness: number,
        outer_spoke_thickness: number,
        num_inner_spokes: number,
        num_outer_spokes: number,

        angle: number,
        angular_velocity: number,
        public readonly origin_x: number,
        public readonly origin_y: number
    ) {
        this.wheel_rest_positions = water_wheel_rest_positions(
            wheel_radius,
            outer_spoke_radius,
            wheel_thickness,
            inner_spoke_thickness,
            outer_spoke_thickness,
            num_inner_spokes,
            num_outer_spokes,
            this.config.draw_radius
        );
        const wheel_count = this.wheel_rest_positions.length / 2;
        const walls = compute_border(
            this.config.params.width,
            this.config.params.height,
            6 * this.config.draw_radius,
            this.config.draw_radius
        );
        const walls_count = walls.length / 2;
        const total_count = wheel_count + walls_count;
        this.temp_data = new Solid_temp_data(total_count);
        this.state = create_simulation_state_of_size(total_count);
        for (let i = 0; i < total_count; ++i) {
            this.state.ids[i] = i;
            this.state.colors[3 * i] = 0.9;
            this.state.colors[3 * i + 1] = 0.7;
            this.state.colors[3 * i + 2] = 0.2;
            if (i >= wheel_count) {
                this.state.positions[i * 2] = walls[(i - wheel_count) * 2];
                this.state.positions[i * 2 + 1] = walls[(i - wheel_count) * 2 + 1];
            }
            this.state.velocities[i * 2] = 0.0;
            this.state.velocities[i * 2 + 1] = 0.0;
        }

        this.set_pose(angle, angular_velocity);
        compute_volumes(
            total_count,
            this.config.params.support,
            this.temp_data.volumes,
            this.temp_data.self_neighborhood.counts,
            this.temp_data.self_neighborhood.kernels
        );
    }

    set_pose(angle: number, angular_velocity: number): void {
        const wheel_count = this.wheel_rest_positions.length / 2;
        compute_moving_water_wheel_state(
            wheel_count,
            angle,
            angular_velocity,
            this.origin_x,
            this.origin_y,
            this.state,
            this.wheel_rest_positions
        );

        compute_neighborhood_parts(this.config, this.state, this.temp_data);
        compute_self_neighborhoods(
            this.config,
            this.temp_data.self_neighborhood,
            this.state,
            this.temp_data
        );
        compute_all_neighborhood_kernels(this.config, this.temp_data.self_neighborhood);
        compute_volumes(
            wheel_count,
            this.config.params.support,
            this.temp_data.volumes,
            this.temp_data.self_neighborhood.counts,
            this.temp_data.self_neighborhood.kernels
        );
    }
}