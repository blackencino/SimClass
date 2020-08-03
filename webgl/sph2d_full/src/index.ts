import {
    Gl_context_bundle,
    get_gl_context_bundle,
    Simple_simulation_renderer,
} from "./sprite_renderer";
import { Simple_simulation } from "./simulation";
import { default_parameters, Parameters, Config } from "./common";

//------------------------------------------------------------------------------
// Initialize the program.
// This function is called after the page has been loaded.
function init() {
    var canvas: HTMLCanvasElement;
    var gl_context_bundle: Gl_context_bundle;
    var simulation: Simple_simulation;
    var simulation_renderer: Simple_simulation_renderer;
    try {
        var canvas_or_null = document.getElementById("webglcanvas") as HTMLCanvasElement;
        if (!canvas_or_null) {
            throw "Cannot get webglcanvas element";
        }
        canvas = canvas_or_null;

        gl_context_bundle = get_gl_context_bundle(canvas);
    } catch (e) {
        const canvas_holder = document.getElementById("canvas-holder");
        if (canvas_holder) {
            canvas_holder.innerHTML = "<p>Sorry, could not get a WebGL graphics context.</p>";
        }
        return;
    }

    const params = default_parameters();
    params.support = 0.025;
    params.width = 1.0;
    params.height = 1.0;
    simulation = new Simple_simulation(params);

    try {
        simulation_renderer = new Simple_simulation_renderer(
            gl_context_bundle,
            simulation.solid_state,
            simulation.fluid_state,
            simulation.config.draw_radius
        );
    } catch (e) {
        const canvas_holder = document.getElementById("canvas-holder");
        if (canvas_holder) {
            canvas_holder.innerHTML =
                "<p>Sorry, could not initialize the WebGL graphics context:" + e + "</p>";
        }
        return;
    }

    var animate_checkbox = document.getElementById("animateCheckbox") as HTMLInputElement;
    var color_checkbox = document.getElementById("colorCheckbox") as HTMLInputElement;
    var size_choice = document.getElementById("sizeChoice") as HTMLSelectElement;
    var animating = true;

    animate_checkbox.checked = true;
    color_checkbox.checked = true;
    size_choice.value = "1.0";

    var prev_time = performance.now();
    var do_frame = (new_time: number): void => {
        if (new_time > prev_time) {
            const delta_time = new_time - prev_time;
            prev_time = new_time;
            const delta_time_seconds = Math.min(3.0, 0.001 * delta_time);
            simulation.step(delta_time_seconds);
            simulation_renderer.update_buffers(simulation.fluid_state);
        }

        simulation_renderer.render(
            Boolean(color_checkbox.checked),
            Number(size_choice.value),
            [1, 0, 0],
            [0, 0, 1],
            params.width,
            params.height
        );

        if (Boolean(animate_checkbox.checked)) {
            requestAnimationFrame(do_frame);
        }
    };

    var regular_change = () => {
        if (!animate_checkbox.checked) {
            prev_time = performance.now();
            do_frame(prev_time);
        }
    };

    color_checkbox.onchange = regular_change;
    size_choice.onchange = regular_change;

    animate_checkbox.onchange = () => {
        if (animate_checkbox.checked) {
            prev_time = performance.now();
            do_frame(prev_time);
        }
    };

    var step_button = document.getElementById("stepButton");
    if (step_button) {
        step_button.onclick = () => {
            simulation.step(0.0);
            //simulation_renderer.update_buffers(simulation.fluid_state);
            simulation_renderer.fluid_renderer.update_buffers(
                simulation.fluid_state.positions,
                simulation.fluid_state.colors
            );
            // simulation_renderer.solid_renderer.update_buffers(
            //     simulation.solid_state.positions,
            //     simulation.solid_state.colors
            // );

            simulation_renderer.render(
                Boolean(color_checkbox.checked),
                Number(size_choice.value),
                [1, 0, 0],
                [0, 0, 1],
                params.width,
                params.height
            );
        };
    }

    var reset_button = document.getElementById("resetButton");
    if (reset_button) {
        reset_button.onclick = () => {
            simulation.reset();
            simulation_renderer.update_buffers(simulation.fluid_state);

            simulation_renderer.render(
                Boolean(color_checkbox.checked),
                Number(size_choice.value),
                [1, 0, 0],
                [0, 0, 1],
                params.width,
                params.height
            );
        };
    }

    prev_time = performance.now();
    do_frame(prev_time);
}

init();