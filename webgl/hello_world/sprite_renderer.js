
"use strict";

//------------------------------------------------------------------------------
function get_gl_context_bundle(canvas) {
    const options = {
        // no need for alpha channel or depth buffer in this program
        alpha: false,
        depth: false,
    };
    const gl =
        canvas.getContext("webgl", options) ||
        canvas.getContext("experimental-webgl", options);
    if (!gl) {
        throw "Browser does not support WebGL";
    }

    const instancing_ext = gl.getExtension("ANGLE_instanced_arrays");
    if (!instancing_ext) {
        throw "Browser does not support ANGLE_instanced_arrays";
    }

    const vertex_array_objects_ext = gl.getExtension("OES_vertex_array_object");
    if (!vertex_array_objects_ext) {
        throw "Browser does not support OES_vertex_array_object";
    }

    return {
        gl: gl,
        instancing: instancing_ext,
        vertex_array_objects: vertex_array_objects_ext,
    };
}

//------------------------------------------------------------------------------
function create_program(gl, vertex_shader_source, fragment_shader_source) {
    const vertex_shader = gl.createShader(gl.VERTEX_SHADER);
    gl.shaderSource(vertex_shader, vertex_shader_source);
    gl.compileShader(vertex_shader);
    if (!gl.getShaderParameter(vertex_shader, gl.COMPILE_STATUS)) {
        throw "Error in vertex shader:  " + gl.getShaderInfoLog(vertex_shader);
    }

    const fragment_shader = gl.createShader(gl.FRAGMENT_SHADER);
    gl.shaderSource(fragment_shader, fragment_shader_source);
    gl.compileShader(fragment_shader);
    if (!gl.getShaderParameter(fragment_shader, gl.COMPILE_STATUS)) {
        throw "Error in fragment shader:  " +
            gl.getShaderInfoLog(fragment_shader);
    }

    const program = gl.createProgram();
    gl.attachShader(program, vertex_shader);
    gl.attachShader(program, fragment_shader);
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        throw "Link error in program:  " + gl.getProgramInfoLog(program);
    }
    return program;
}

//------------------------------------------------------------------------------
function init_program_info(gl) {
    const vertex_source = `
        attribute vec2 a_quad_corner;
        attribute vec2 a_center;
        attribute vec3 a_color;
        varying vec3 v_color;
        varying vec2 v_point_coord;
        uniform float u_point_size;

        uniform mat4 u_projection;
        uniform mat4 u_modelview;
        void main() {
            v_point_coord = a_quad_corner;
            vec2 pos2d = a_center + 0.5 * u_point_size * a_quad_corner;
            gl_Position = u_projection * u_modelview * vec4(pos2d, 0, 1);
            v_color = a_color;
        }`;

    const fragment_source = `
        precision mediump float;
        varying vec3 v_color;
        varying vec2 v_point_coord;
        void main() {
            float distance_from_center = length(v_point_coord);
            if (distance_from_center >= 1.0) {
                discard;
            }
            gl_FragColor = vec4(v_color, 1.0);
        }`;

    const program = create_program(
        gl,
        vertex_source,
        fragment_source);

    return {
        program: program,
        attribute_locations: {
            a_quad_corner: gl.getAttribLocation(program, "a_quad_corner"),
            a_center: gl.getAttribLocation(program, "a_center"),
            a_color: gl.getAttribLocation(program, "a_color")
        },
        uniform_locations: {
            u_projection: gl.getUniformLocation(program, "u_projection"),
            u_modelview: gl.getUniformLocation(program, "u_modelview"),
            u_point_size: gl.getUniformLocation(program, "u_point_size")
        }
    };
}
//------------------------------------------------------------------------------
class Sprite_renderer {
    constructor(gl_context_bundle, program_info, state) {
        this.gl_context_bundle = gl_context_bundle;
        const gl = this.gl_context_bundle.gl;
        const gl_inst = this.gl_context_bundle.instancing;
        const gl_vao = this.gl_context_bundle.vertex_array_objects;

        this.program_info = program_info;
        this.buffers = {
            quad_corners: gl.createBuffer(),
            centers: gl.createBuffer(),
            colors: gl.createBuffer(),
        };

        this.point_count = state.positions.length / 2;

        const attrs = this.program_info.attribute_locations;
        const unis = this.program_info.uniform_locations;

        this.vao = gl_vao.createVertexArrayOES();
        gl_vao.bindVertexArrayOES(this.vao);

        const bind = ({
            buffer,
            data,
            draw_type = gl.STATIC_DRAW,
            attribute,
            dimension,
            divisor,
        }) => {
            gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
            gl.bufferData(gl.ARRAY_BUFFER, data, draw_type);
            gl.vertexAttribPointer(attribute, dimension, gl.FLOAT, false, 0, 0);
            gl_inst.vertexAttribDivisorANGLE(attribute, divisor);
            gl.enableVertexAttribArray(attribute);
        };

        // Make quad corners
        const quad_corners = new Float32Array([
            -1.0,
            -1.0,
            1.0,
            -1.0,
            1.0,
            1.0,
            -1.0,
            1.0,
        ]);
        bind({
            buffer: this.buffers.quad_corners,
            data: quad_corners,
            attribute: attrs.a_quad_corner,
            dimension: 2,
            divisor: 0,
        });

        bind({
            buffer: this.buffers.centers,
            data: state.positions,
            draw_type: gl.STREAM_DRAW,
            attribute: attrs.a_center,
            dimension: 2,
            divisor: 1,
        });

        bind({
            buffer: this.buffers.colors,
            data: state.colors,
            attribute: attrs.a_color,
            dimension: 3,
            divisor: 1,
        });
    }

    update_buffers(state) {
        //assert(state.positions.length / 2 === this.point_count);
        const gl = this.gl_context_bundle.gl;
        gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.centers);
        gl.bufferData(gl.ARRAY_BUFFER, state.positions, gl.STREAM_DRAW);
    }

    render({
        do_array_colors = true,
        point_size = 10.0,
        projection,
        modelview,
        fixed_color = [1, 0, 0],
    } = {}) {
        const gl = this.gl_context_bundle.gl;
        const gl_inst = this.gl_context_bundle.instancing;
        const gl_vao = this.gl_context_bundle.vertex_array_objects;
        const attrs = this.program_info.attribute_locations;
        const unis = this.program_info.uniform_locations;

        gl.useProgram(this.program_info.program);
        gl_vao.bindVertexArrayOES(this.vao);

        if (do_array_colors) {
            gl.enableVertexAttribArray(attrs.a_color);
        } else {
            gl.disableVertexAttribArray(attrs.a_color);
            gl.vertexAttrib3f(
                attrs.a_color,
                fixed_color[0],
                fixed_color[1],
                fixed_color[2]
            );
        }

        // Other uniforms
        gl.uniform1f(unis.u_point_size, point_size);
        gl.uniformMatrix4fv(unis.u_projection, false, projection);
        gl.uniformMatrix4fv(unis.u_modelview, false, modelview);

        // Instancing!
        gl_inst.drawArraysInstancedANGLE(
            gl.TRIANGLE_FAN,
            0,
            4,
            this.point_count
        );
    }
}
//------------------------------------------------------------------------------
class Simple_simulation_renderer {
    constructor(gl_context_bundle, simple_simulation) {
        this.gl_context_bundle = gl_context_bundle;
        this.program_info = init_program_info(this.gl_context_bundle.gl);
        this.simulation = simple_simulation;

        this.fluid_renderer = new Sprite_renderer(
            this.gl_context_bundle,
            this.program_info,
            this.simulation.state
        );
        this.fluid_modelview = glMatrix.mat4.create();
        glMatrix.mat4.identity(this.fluid_modelview);

        this.solid_renderer = new Sprite_renderer(
            this.gl_context_bundle,
            this.program_info,
            this.simulation.solid_state
        );
        this.solid_modelview = glMatrix.mat4.create();
        glMatrix.mat4.identity(this.solid_modelview);
    }

    update_buffers(simulation) {
        this.fluid_renderer.update_buffers(simulation.state);
        //this.solid_renderer.update_buffers(simulation.solid_state);
    }

    render({
        do_array_colors = true,
        point_size_gain = 1.0,
        fluid_fixed_color = [0, 0, 1],
        solid_fixed_color = [1, 0, 0],
        width = 512,
        height = 512,
    } = {}) {
        const gl = this.gl_context_bundle.gl;
        gl.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);

        const projection = glMatrix.mat4.create();
        const border = 50.0;
        glMatrix.mat4.ortho(projection, -border, width + border, -border, height + border, -1.0, 1.0);

        this.solid_renderer.render({
            do_array_colors: do_array_colors,
            point_size: 2.0 * point_size_gain * this.simulation.radius,
            projection: projection,
            modelview: this.solid_modelview,
            fixed_color: solid_fixed_color,
        });

        this.fluid_renderer.render({
            do_array_colors: do_array_colors,
            point_size: 2.0 * point_size_gain * this.simulation.radius,
            projection: projection,
            modelview: this.solid_modelview,
            fixed_color: fluid_fixed_color,
        });
    }
}

module.exports = {
    get_gl_context_bundle,
    Sprite_renderer,
    Simple_simulation_renderer
};
