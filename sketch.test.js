const { IX, add_source, lin_solve, set_bnd } = require('./sketch');

describe('Stable Fluids Core Functions', () => {
  test('IX maps 2D indices to 1D array index correctly', () => {
    expect(IX(1, 1)).toBe(1 + (128 + 2) * 1);
    expect(IX(0, 0)).toBe(0);
    expect(IX(128, 128)).toBe(128 + (128 + 2) * 128);
  });

  test('add_source adds scaled source to field', () => {
    const N = 4;
    const size = (N + 2) * (N + 2);
    let x = new Array(size).fill(0);
    let s = new Array(size).fill(2);
    add_source(N, x, s, 0.5);
    for (let i = 0; i < size; i++) {
      expect(x[i]).toBeCloseTo(1);
    }
  });

  test('lin_solve solves a simple linear system', () => {
    const N = 4;
    const size = (N + 2) * (N + 2);
    let x = new Array(size).fill(0);
    let x0 = new Array(size).fill(1);
    set_bnd(N, 0, x); // Set boundary conditions before solving
    lin_solve(N, 0, x, x0, 1, 4);
    set_bnd(N, 0, x); // Set boundary conditions after solving
    // The solution should be stable and not NaN or Infinity
    for (let i = 0; i < size; i++) {
      expect(Number.isFinite(x[i])).toBe(true);
    }
  });
});
