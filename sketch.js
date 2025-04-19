// --- Stable Fluids (Stam, 1999) Implementation ---
// This code implements the core steps of Jos Stam's "Stable Fluids" algorithm:
// add_source, diffuse, advect, and project. See: https://graphics.cs.cmu.edu/nsp/course/15-464/Spring11/papers/StamFluidforGames.pdf

let N = 128; // Grid size (N x N), as in Stam's paper
let size = (N + 2) * (N + 2); // Include boundary cells
let dt = 0.1; // Time step
let diff = 0.0001; // Diffusion rate
let visc = 0.0001; // Viscosity
let dens, dens_prev; // Density fields
let u, v, u_prev, v_prev; // Velocity fields
let colorModeEnabled = false;

// Convert 2D indices to 1D for array access (Stam's IX function)
function IX(x, y) {
    return x + (N + 2) * y;
}

// --- p5.js Setup ---
function setup() {
    createCanvas(512, 512);
    initializeFields();
    setupColorToggle();
}

// --- p5.js Draw Loop ---
function draw() {
    background(0);
    handleMouseInput();
    stepFluid();
    renderDens(dens);
    resetPrevFields();
}

// --- Initialize all simulation fields to zero ---
function initializeFields() {
    dens = new Array(size).fill(0);
    dens_prev = new Array(size).fill(0);
    u = new Array(size).fill(0);
    v = new Array(size).fill(0);
    u_prev = new Array(size).fill(0);
    v_prev = new Array(size).fill(0);
}

// --- Attach event listener for color toggle UI ---
function setupColorToggle() {
    const colorToggle = document.getElementById('colorToggle');
    if (colorToggle) {
        colorToggle.addEventListener('change', function() {
            colorModeEnabled = colorToggle.checked;
        });
    }
}

// --- Handle mouse input to inject density and velocity ---
function handleMouseInput() {
    if (mouseIsPressed && mouseX >= 0 && mouseX < width && mouseY >= 0 && mouseY < height) {
        let x = floor(map(mouseX, 0, width, 1, N));
        let y = floor(map(mouseY, 0, height, 1, N));
        injectDensityAndVelocity(x, y);
    }
}

// --- Inject density and velocity at a given grid location ---
function injectDensityAndVelocity(x, y) {
    for (let i = -2; i <= 2; i++) {
        for (let j = -2; j <= 2; j++) {
            let idx = IX(x + i, y + j);
            if (idx > 0 && idx < size) {
                dens[idx] += 100; // Inject density (see add_source in Stam)
                u[idx] += (movedX || 0) * 2; // Inject velocity
                v[idx] += (movedY || 0) * 2;
            }
        }
    }
}

// --- Advance the fluid simulation by one time step ---
function stepFluid() {
    vel_step(N, u, v, u_prev, v_prev, visc, dt);
    dens_step(N, dens, dens_prev, u, v, diff, dt);
}

// --- Reset previous fields to zero for the next frame ---
function resetPrevFields() {
    for (let i = 0; i < size; i++) {
        u_prev[i] = 0;
        v_prev[i] = 0;
        dens_prev[i] = 0;
    }
}

// --- Render the density field as grayscale or color ---
function renderDens(d) {
    let cellW = width / N;
    let cellH = height / N;
    for (let i = 1; i <= N; i++) {
        for (let j = 1; j <= N; j++) {
            let val = constrain(d[IX(i, j)], 0, 255);
            if (colorModeEnabled) {
                colorMode(HSB, 255);
                fill(val, 200, 255);
                colorMode(RGB, 255);
            } else {
                fill(val);
            }
            noStroke();
            rect((i - 1) * cellW, (j - 1) * cellH, cellW, cellH);
        }
    }
}

// --- Fluid Solver Functions (see Stam, 1999) ---
// --- Density step: add source, diffuse, and advect density ---
function dens_step(N, x, x0, u, v, diff, dt) {
    add_source(N, x, x0, dt); // Add density sources
    [x0, x] = [x, x0];
    diffuse(N, 0, x, x0, diff, dt); // Diffuse density
    [x0, x] = [x, x0];
    advect(N, 0, x, x0, u, v, dt); // Advect density by velocity field
}

// --- Velocity step: add source, diffuse, project, advect, project ---
function vel_step(N, u, v, u0, v0, visc, dt) {
    add_source(N, u, u0, dt); // Add velocity sources
    add_source(N, v, v0, dt);
    [u0, u] = [u, u0];
    diffuse(N, 1, u, u0, visc, dt); // Diffuse velocity u
    [v0, v] = [v, v0];
    diffuse(N, 2, v, v0, visc, dt); // Diffuse velocity v
    project(N, u, v, u0, v0); // Project to make velocity field mass conserving
    [u0, u] = [u, u0];
    [v0, v] = [v, v0];
    advect(N, 1, u, u0, u0, v0, dt); // Advect velocity u
    advect(N, 2, v, v0, u0, v0, dt); // Advect velocity v
    project(N, u, v, u0, v0); // Project again
}

// --- Add source to field x from s (see add_source in Stam) ---
function add_source(N, x, s, dt) {
    for (let i = 0; i < (N + 2) * (N + 2); i++) x[i] += dt * s[i];
}

// --- Set boundary conditions (see set_bnd in Stam) ---
function set_bnd(N, b, x) {
    for (let i = 1; i <= N; i++) {
        x[IX(0, i)]     = b === 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = b === 1 ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)]     = b === 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = b === 2 ? -x[IX(i, N)] : x[IX(i, N)];
    }
    x[IX(0, 0)]         = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)]     = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)]     = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

// --- Linear solver for diffusion and projection (see lin_solve in Stam) ---
function lin_solve(N, b, x, x0, a, c) {
    for (let k = 0; k < 20; k++) {
        for (let i = 1; i <= N; i++) {
            for (let j = 1; j <= N; j++) {
                x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
            }
        }
        set_bnd(N, b, x);
    }
}

// --- Diffuse field x0 into x (see diffuse in Stam) ---
function diffuse(N, b, x, x0, diff, dt) {
    let a = dt * diff * N * N;
    lin_solve(N, b, x, x0, a, 1 + 4 * a);
}

// --- Advect field d0 into d by velocity field (see advect in Stam) ---
function advect(N, b, d, d0, u, v, dt) {
    let dt0 = dt * N;
    for (let i = 1; i <= N; i++) {
        for (let j = 1; j <= N; j++) {
            let x = i - dt0 * u[IX(i, j)];
            let y = j - dt0 * v[IX(i, j)];
            if (x < 0.5) x = 0.5;
            if (x > N + 0.5) x = N + 0.5;
            let i0 = floor(x);
            let i1 = i0 + 1;
            if (y < 0.5) y = 0.5;
            if (y > N + 0.5) y = N + 0.5;
            let j0 = floor(y);
            let j1 = j0 + 1;
            let s1 = x - i0;
            let s0 = 1 - s1;
            let t1 = y - j0;
            let t0 = 1 - t1;
            // Bilinear interpolation (see eq. 12 in Stam)
            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(N, b, d);
}

// --- Project velocity field to be mass conserving (see project in Stam) ---
function project(N, u, v, p, div) {
    for (let i = 1; i <= N; i++) {
        for (let j = 1; j <= N; j++) {
            div[IX(i, j)] = -0.5 * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(N, 0, div);
    set_bnd(N, 0, p);
    lin_solve(N, 0, p, div, 1, 4);
    for (let i = 1; i <= N; i++) {
        for (let j = 1; j <= N; j++) {
            u[IX(i, j)] -= 0.5 * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
            v[IX(i, j)] -= 0.5 * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
        }
    }
    set_bnd(N, 1, u);
    set_bnd(N, 2, v);
}

// Export core functions for unit testing
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        IX,
        add_source,
        lin_solve,
        diffuse,
        advect,
        set_bnd,
        project,
        dens_step,
        vel_step
    };
}