// --- UI Helper Functions ---
// Handles color toggle UI
function setupColorToggle() {
    const colorToggle = document.getElementById('colorToggle');
    if (colorToggle) {
        colorToggle.addEventListener('change', function() {
            colorModeEnabled = colorToggle.checked;
        });
    }
}

// Handles mouse input for injecting density and velocity
function handleMouseInput() {
    if (mouseIsPressed && mouseX >= 0 && mouseX < width && mouseY >= 0 && mouseY < height) {
        let x = floor(map(mouseX, 0, width, 1, N));
        let y = floor(map(mouseY, 0, height, 1, N));
        injectDensityAndVelocity(x, y);
    }
}

// Injects density and velocity at a given grid location
function injectDensityAndVelocity(x, y) {
    for (let i = -2; i <= 2; i++) {
        for (let j = -2; j <= 2; j++) {
            let idx = IX(x + i, y + j);
            if (idx > 0 && idx < size) {
                dens[idx] += 100;
                u[idx] += (movedX || 0) * 2;
                v[idx] += (movedY || 0) * 2;
            }
        }
    }
}

// Resets previous fields to zero for the next frame
function resetPrevFields() {
    for (let i = 0; i < size; i++) {
        u_prev[i] = 0;
        v_prev[i] = 0;
        dens_prev[i] = 0;
    }
}
