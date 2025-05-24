/*!
 * Raindrop Canvas Effect
 * 
 * © Mandala Logics – https://operalimina.co.uk
 * GitHub: https://github.com/Mandala-Logics/raindrop-canvas
 * Buy Me a Coffee: https://www.buymeacoffee.com/operalimina
 * 
 * This code is free to use for non-commercial purposes.
 * If you use it, please credit operalimina.co.uk.
 */

import { RaindropCanvas } from './raindrops.js';
import { DropPolygon } from './raindrops.js'; 

let context;
let canvas;
let dropTimer;
let trailTimer;
let sys;
let mainInner;
let rect;
let contain_div;

function ResizeCanvas() {

    rect = contain_div.getBoundingClientRect();

    canvas.width = rect.width;
    canvas.height = rect.height;

    if (sys) sys.OnCanvasResize(canvas); // all the raindrop system to know how big the canvas is now, for removal of cells outside.
}

function DrawRaindrops(polygons) {
context.fillStyle = 'black';
context.fillRect(0, 0, canvas.width, canvas.height);

context.globalCompositeOperation = 'destination-out';
context.fillStyle = 'rgba(0, 0, 0, 1)';

polygons.forEach(polygon => {
    polygon.FillSmooth(context);
});

context.globalCompositeOperation = 'source-over';
}

window.onload = () => {

    // Create the color shift background layer
    let div = document.createElement("div");
    div.classList.add("color-shift-layer");
    document.body.appendChild(div);

    contain_div = document.createElement("div");
    contain_div.classList.add("mask-container");
    document.body.appendChild(contain_div);

    // Setup canvas and context
    canvas = document.createElement("canvas");
    canvas.classList.add("mask-layer");
    contain_div.appendChild(canvas);
    context = canvas.getContext("2d");

    // Create the color shift background layer
    div = document.createElement("div");
    div.classList.add("blur-layer");
    document.body.appendChild(div);

    // Initialize system
    sys = new RaindropCanvas(DrawRaindrops);

    sys.ApplySettings({
    maxCellCount: 200,
    smoothLoopLimit: 2e3,
    canvasMargin : 20
    });

    // Set canvas size and bind resize handler
    mainInner = document.querySelector('.md-main__inner.md-grid');
    ResizeCanvas();
    window.addEventListener("resize", ResizeCanvas);

    context.fillStyle = 'black';
    context.fillRect(0, 0, canvas.width, canvas.height);

    // Start drop spawners
    dropTimer = window.setInterval(() => {
        sys.SpawnDrop(Math.random() * 10 + 10,
            Math.random() * rect.width + rect.left,
            Math.random() * rect.height + rect.top); 
    }, 1000);

    trailTimer = window.setInterval(() => {
    sys.ShedDrops();
    }, 500);

    // Start animation loop
    sys.Start();

};