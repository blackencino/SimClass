import { Particle_renderer } from "./particle_renderer.js";

function main() {
    const glcanvas = document.querySelector("#glcanvas");
    // glcanvas.style = "margin-left: -50%;";
    const particle_renderer = new Particle_renderer(...);
    particle_renderer.render();
}