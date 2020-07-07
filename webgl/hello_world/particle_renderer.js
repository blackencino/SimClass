//------------------------------------------------------------------------------
function Particle_renderer(width, height, density) {
    this.width = width;
    this.height = height;

    this.uniforms = {};

    this.canvas = document.querySelector("#glcanvas");
    this.gl = this.canvas.getContext("webgl");
    if (!this.gl) {
        alert(
            "Unable to initialize WebGL. Your browser or machine may not support it."
        );
        return;
    }

    this.instancing_ext = this.gl.getExtension('ANGLE_instanced_arrays');
    if (!this.instancing_ext) {
        alert(
              "Unable to initialize ANGLE_instanced_arrays extension."
              );
        return;
    }

    this.program_info = init_program_info(this.gl);
    this.buffers = init_buffers(this.gl);
}

//------------------------------------------------------------------------------
Particle_renderer.prototype.render = function() {
    draw_scene(
        this.gl,
        this.program_info,
        this.buffers,
        this.uniforms
    );
};

//------------------------------------------------------------------------------
function init_program_info(gl) {
    const vertex_source = `
        vec2 rotate(vec2 v, float a) {
            float s = sin(a);
            float c = cos(a);
            mat2 m = mat2(c, -s, s, c);
            return m * v;
        }

        attribute vec2 quad_corner;
        attribute vec2 center;
        attribute float angle;
        attribute vec4 rgba;

        uniform mat4x4 modelview;
        uniform mat4x4 projection;

        uniform float radius;

        varying medp vec4 v_rgba;
        varying highp vec2 v_ndc;

        void main() {
            v_ndc = quad_corner;
            v_rgba = rgba;

            vec2 pos2d = center + radius * rotate(quad_corner, -angle);
            gl_Position = projection * modelview * vec4(pos2d, 0, 1);
        }
    `;

    const fragment_source = `
        precision highp float;

        varying medp vec4 v_rgba;
        varying highp vec2 v_ndc;

        void main() {
            float r = length(v_ndc);
            if (r > 1) {
                discard;
            }

            // Make a little white triangle on right
            float x_edge = mix(0.3, 0.0, v_ndc.x);
            float tri_alpha = step(0.5, v_ndc.x) *
                              step(-x_edge, v_ndc.y) *
                              (1 - step(x_edge, v_ndc.y));

            if (tri_alpha == 0 && r > 0.75) {
                discard;
            }

            gl_FragColor = mix(v_rgba, vec4(1, 1, 1, 1), tri_alpha);
        }
    `;

    const program = init_shader_program(gl, vertex_source, fragment_source);
    const program_info = {
        program: program,
        attrib_locations: {
             quad_corner: gl.getAttribLocation(program, "quad_corner"),
             center: gl.getAttribLocation(program, "center"),
             angle: gl.getAttribLocation(program, "angle"),
             rgba: gl.getAttribLocation(program, "rgba")
        },
        uniform_locations: {
            modelview: gl.getUniformLocation(program, "modelview"),
            projection: gl.getUniformLocation(program, "projection"),
            radius: gl.getUniformLocation(program, "radius")
        }
    };

    return program_info;
}

//------------------------------------------------------------------------------
function init_shader_program(gl, vertex_source, fragment_source) {
    const vertex_shader = load_shader(gl, gl.VERTEX_SHADER, vertex_source);
    const fragment_shader = load_shader(
        gl,
        gl.FRAGMENT_SHADER,
        fragment_source
    );
    const program = gl.createProgram();
    gl.attachShader(program, vertex_shader);
    gl.attachShader(program, fragment_shader);
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        alert(
            "Unable to initialize the shader program: " +
                gl.getProgramInfoLog(program)
        );
        return null;
    }
    return program;
}

//------------------------------------------------------------------------------
function load_shader(gl, type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        alert(
            "An error occurred compiling the shaders: " +
                gl.getShaderInfoLog(shader)
        );
        gl.deleteShader(shader);
        return null;
    }
    return shader;
}

//------------------------------------------------------------------------------
function init_buffers(gl, quad_corners, centers, angles, rgbas) {
    const quad_corner_buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, quad_corner_buffer);
    gl.bufferData(gl.ARRAY_BUFFER, quad_corners, gl.STATIC_DRAW);

    gl.bindBuffer(gl.ARRAY_BUFFER,

    return {
        position: position_buffer
    };
}

//------------------------------------------------------------------------------
function draw_scene(gl, program_info, buffers, textures, uniforms) {
    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    // bind positions
    {
        const num_components = 2;
        const type = gl.FLOAT;
        const normalize = false;
        const stride = 0;
        const offset = 0;
        gl.bindBuffer(gl.ARRAY_BUFFER, buffers.position);
        gl.vertexAttribPointer(
            program_info.attrib_locations.position,
            num_components,
            type,
            normalize,
            stride,
            offset
        );
        gl.enableVertexAttribArray(program_info.attrib_locations.position);
    }

    // Uniforms
    gl.useProgram(program_info.program);

    // Draw the slab
    {
        const offset = 0;
        const vertex_count = 4;
        gl.drawArrays(gl.TRIANGLE_STRIP, offset, vertex_count);
    }
}

module.exports = {
    Particle_renderer
};