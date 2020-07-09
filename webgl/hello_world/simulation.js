
"use strict";

//------------------------------------------------------------------------------
function create_simulation_state(count) {
    return {
        ids: new Uint32Array(count),
        positions: new Float32Array(2 * count),
        velocities: new Float32Array(2 * count),
        colors: new Float32Array(3 * count),
    };
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
    let i;
    for (i = 0; i < POINT_COUNT; ++i) {
        state.ids[i] = i;

        // Random point in the width/height box
        state.positions[2 * i] = width * my_rng();
        state.positions[2 * i + 1] = height * my_rng();

        const speed = speed_gain * (1.0 + 3.0 * my_rng());
        const angle = 2.0 * Math.PI * my_rng();
        state.velocities[2 * i] = speed * Math.cos(angle);
        state.velocities[2 * i + 1] = speed * Math.sin(angle);
    }
    for (i = 0; i < 3 * POINT_COUNT; ++i) {
        state.colors[i] = my_rng();
    }
}
//------------------------------------------------------------------------------
function simple_simulation_step(
    state,
    { delta_time = 1.0, radius = 1.0, width = 512, height = 512 } = {}
) {
    const POINT_COUNT = state.positions.length / 2;
    let i;
    // X-coordinates
    for (i = 0; i < 2 * POINT_COUNT; i += 2) {
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
    for (i = 1; i < 2 * POINT_COUNT; i += 2) {
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
    let bx;
    let by;
    for (by = 0; by < N; ++by) {
        const byf = (by + 1) / (N + 1);

        for (bx = 0; bx < N; ++bx) {
            const bxf = (bx + 1) / (N + 1);
            tiny_block([bxf + my_rng_dist(), byf + my_rng_dist()], my_rng_dist());
        }
    }

    const positions = new Float32Array(positions_array);
    const count = positions.length / 2;
    const ids = new Uint32Array(count);
    const velocities = new Float32Array(count * 2);
    const colors = new Float32Array(count * 3);
    let i;
    for (i = 0; i < count; ++i) {
        ids[i] = i;

        velocities[2 * i] = 0.0;
        velocities[2 * i + 1] = 0.0;

        colors[3 * i] = 0.784;
        colors[3 * i + 1] = 0.608;
        colors[3 * i + 2] = 0.078;
    }

    return {
        ids: ids,
        positions: positions,
        velocities: velocities,
        colors: colors,
    };
}
//------------------------------------------------------------------------------
class Simple_simulation {
    constructor({ count = 1000, support = 0.01, width = 512, height = 512 }) {
        this.support = support;
        const HHalf = 0.5 * support;
        this.radius = 0.99 * HHalf;
        this.extent = 2.0 * support;
        this.width = width;
        this.height = height;

        // Make fluid state
        this.state = create_simulation_state(count);
        random_box_initial_state(this.state, { width: width, height: height });

        // Make solid state
        this.solid_state = world_walls_initial_solid_state({
            support: support,
            width: width,
            height: height,
        });

        this.accumulated_time = 0.0;
    }

    step(delta_time) {
        simple_simulation_step(this.state, {
            delta_time: delta_time,
            radius: this.radius,
            width: this.width,
            height: this.height,
        });
    }
}

module.exports = {
    Simple_simulation
};

