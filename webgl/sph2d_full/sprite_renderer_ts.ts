//import glMatrix = require("gl-matrix");
//import from './gl-matrix.js'

export class Gl_context_bundle {
    constructor(
        public readonly gl: WebGLRenderingContext,
        public readonly instancing: ANGLE_instanced_arrays,
        public readonly vertex_array_objects: OES_vertex_array_object
    ) {}
}

//------------------------------------------------------------------------------
export function get_gl_context_bundle(
    canvas: HTMLCanvasElement
): Gl_context_bundle {
    const options = {
        // no need for alpha channel or depth buffer in this program
        alpha: false,
        depth: false,
    };
    const gl =
        (canvas.getContext("webgl", options) as WebGLRenderingContext) ||
        (canvas.getContext(
            "experimental-webgl",
            options
        ) as WebGLRenderingContext);
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

    return new Gl_context_bundle(gl, instancing_ext, vertex_array_objects_ext);
}

//------------------------------------------------------------------------------
function create_program(
    gl: WebGLRenderingContext,
    vertex_shader_source: string,
    fragment_shader_source: string
): WebGLProgram {
    const vertex_shader = gl.createShader(gl.VERTEX_SHADER);
    if (!vertex_shader) {
        throw "Could not create vertex shader";
    }
    gl.shaderSource(vertex_shader, vertex_shader_source);
    gl.compileShader(vertex_shader);
    if (!gl.getShaderParameter(vertex_shader, gl.COMPILE_STATUS)) {
        throw "Error in vertex shader:  " + gl.getShaderInfoLog(vertex_shader);
    }

    const fragment_shader = gl.createShader(gl.FRAGMENT_SHADER);
    if (!fragment_shader) {
        throw "Could not create fragment shader";
    }
    gl.shaderSource(fragment_shader, fragment_shader_source);
    gl.compileShader(fragment_shader);
    if (!gl.getShaderParameter(fragment_shader, gl.COMPILE_STATUS)) {
        throw "Error in fragment shader:  " +
            gl.getShaderInfoLog(fragment_shader);
    }

    const program = gl.createProgram();
    if (!program) {
        throw "Could not create GL program";
    }
    gl.attachShader(program, vertex_shader);
    gl.attachShader(program, fragment_shader);
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        throw "Link error in program:  " + gl.getProgramInfoLog(program);
    }
    return program;
}

class Program_info {
    program: WebGLProgram;
    attribute_locations: {
        a_quad_corner: GLint;
        a_center: GLint;
        a_color: GLint;
    };
    uniform_locations: {
        u_projection: WebGLUniformLocation | null;
        u_modelview: WebGLUniformLocation | null;
        u_point_size: WebGLUniformLocation | null;
    };

    constructor(gl: WebGLRenderingContext) {
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

        this.program = create_program(gl, vertex_source, fragment_source);
        this.attribute_locations = {
            a_quad_corner: gl.getAttribLocation(this.program, "a_quad_corner"),
            a_center: gl.getAttribLocation(this.program, "a_center"),
            a_color: gl.getAttribLocation(this.program, "a_color"),
        };
        this.uniform_locations = {
            u_projection: gl.getUniformLocation(this.program, "u_projection"),
            u_modelview: gl.getUniformLocation(this.program, "u_modelview"),
            u_point_size: gl.getUniformLocation(this.program, "u_point_size"),
        };
    }
}
//------------------------------------------------------------------------------
export class Sprite_renderer {
    buffers: {
        quad_corners: WebGLBuffer;
        centers: WebGLBuffer;
        colors: WebGLBuffer;
    };

    point_count: number;

    vao: WebGLVertexArrayObjectOES;

    constructor(
        public readonly gl_context_bundle: Gl_context_bundle,
        public readonly program_info: Program_info,
        positions: Float32Array,
        colors: Float32Array
    ) {
        this.gl_context_bundle = gl_context_bundle;
        const gl = this.gl_context_bundle.gl;
        const gl_inst = this.gl_context_bundle.instancing;
        const gl_vao = this.gl_context_bundle.vertex_array_objects;

        this.program_info = program_info;
        const quad_corners_buffer = gl.createBuffer();
        if (!quad_corners_buffer) {
            throw "Could not create quad_corners buffer";
        }

        const centers_buffer = gl.createBuffer();
        if (!centers_buffer) {
            throw "Could not create centers buffer";
        }

        const colors_buffer = gl.createBuffer();
        if (!colors_buffer) {
            throw "Could not create colors_buffer";
        }

        this.buffers = {
            quad_corners: quad_corners_buffer,
            centers: centers_buffer,
            colors: colors_buffer,
        };

        this.point_count = positions.length / 2;

        const attrs = this.program_info.attribute_locations;
        const unis = this.program_info.uniform_locations;

        const vao = gl_vao.createVertexArrayOES();
        if (!vao) {
            throw "Could not create vertex array object";
        }
        this.vao = vao;

        gl_vao.bindVertexArrayOES(this.vao);

        const bind = (
            buffer: WebGLBuffer,
            data: Float32Array,
            draw_type: number,
            attribute: number,
            dimension: number,
            divisor: number
        ) => {
            gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
            gl.bufferData(gl.ARRAY_BUFFER, data, draw_type);
            gl.vertexAttribPointer(attribute, dimension, gl.FLOAT, false, 0, 0);
            gl_inst.vertexAttribDivisorANGLE(attribute, divisor);
            gl.enableVertexAttribArray(attribute);
        };

        // Make quad corners
        const quad_corners_data = new Float32Array([
            -1.0,
            -1.0,
            1.0,
            -1.0,
            1.0,
            1.0,
            -1.0,
            1.0,
        ]);
        bind(
            this.buffers.quad_corners,
            quad_corners_data,
            gl.STATIC_DRAW,
            attrs.a_quad_corner,
            2,
            0
        );

        bind(
            this.buffers.centers,
            positions,
            gl.STREAM_DRAW,
            attrs.a_center,
            2,
            1
        );

        bind(this.buffers.colors, colors, gl.STREAM_DRAW, attrs.a_color, 3, 1);
    }

    // Right now only positions are updated - we'll change that
    update_buffers(
        positions: Float32Array,
        colors: Float32Array | null = null
    ) {
        const gl = this.gl_context_bundle.gl;

        const input_count = positions.length / 2;
        if (input_count != this.point_count) {
            gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.centers);
            gl.bufferData(gl.ARRAY_BUFFER, positions, gl.STREAM_DRAW);

            if (colors) {
                if (colors.length < input_count * 3) {
                    throw "Invalid colors array length";
                }
                gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.colors);
                gl.bufferData(gl.ARRAY_BUFFER, colors, gl.STREAM_DRAW);
            }
            this.point_count = input_count;
        } else {
            gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.centers);
            gl.bufferSubData(gl.ARRAY_BUFFER, 0, positions);

            if (colors) {
                if (colors.length < input_count * 3) {
                    throw "Invalid colors array length";
                }
                gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.colors);
                gl.bufferSubData(gl.ARRAY_BUFFER, 0, colors);
            }
        }
    }

    render(
        do_array_colors: boolean,
        point_size: number,
        projection: glMatrix.mat4,
        modelview: glMatrix.mat4,
        fixed_color: Array<number> = [1, 0, 0]
    ) {
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
interface Physical_state {
    positions: Float32Array;
    colors: Float32Array;
}

//------------------------------------------------------------------------------
export class Simple_simulation_renderer {
    program_info: Program_info;
    solid_renderer: Sprite_renderer;
    solid_modelview: glMatrix.mat4;
    fluid_renderer: Sprite_renderer;
    fluid_modelview: glMatrix.mat4;

    constructor(
        public readonly gl_context_bundle: Gl_context_bundle,
        solid_state: Physical_state,
        fluid_state: Physical_state,
        public readonly radius: number
    ) {
        this.gl_context_bundle = gl_context_bundle;
        this.program_info = new Program_info(this.gl_context_bundle.gl);

        this.solid_renderer = new Sprite_renderer(
            this.gl_context_bundle,
            this.program_info,
            solid_state.positions,
            solid_state.colors
        );
        this.solid_modelview = glMatrix.mat4.create();
        glMatrix.mat4.identity(this.solid_modelview);

        this.fluid_renderer = new Sprite_renderer(
            this.gl_context_bundle,
            this.program_info,
            fluid_state.positions,
            fluid_state.colors
        );
        this.fluid_modelview = glMatrix.mat4.create();
        glMatrix.mat4.identity(this.fluid_modelview);
    }

    update_buffers(fluid_state: Physical_state) {
        this.fluid_renderer.update_buffers(
            fluid_state.positions,
            fluid_state.colors
        );
    }

    render(
        do_array_colors: boolean,
        point_size_gain: number,
        solid_fixed_color: Array<number>,
        fluid_fixed_color: Array<number>,
        width: number,
        height: number
    ) {
        const gl = this.gl_context_bundle.gl;
        gl.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);

        const projection = glMatrix.mat4.create();
        const border = 6 * this.radius;
        glMatrix.mat4.ortho(
            projection,
            -border,
            width + border,
            -border,
            height + border,
            -1.0,
            1.0
        );

        this.solid_renderer.render(
            do_array_colors,
            2.0 * this.radius * point_size_gain,
            projection,
            this.solid_modelview,
            solid_fixed_color
        );

        this.fluid_renderer.render(
            do_array_colors,
            2.0 * this.radius * point_size_gain,
            projection,
            this.fluid_modelview,
            fluid_fixed_color
        );
    }
}