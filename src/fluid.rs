use super::{ITER, N};

const NSIZE: usize = (N * N) as usize;

/**
 * flat index
 */
pub fn ix(i: i64, j: i64) -> usize {
    (i + j * N) as usize
}

/**
 * set bounds by average neighbors
 */
fn set_bnd(b: i64, x: &mut [f64; NSIZE]) {
    for i in 1..(N - 1) {
        x[ix(i, 0)] = if b == 2 { -x[ix(i, 1)] } else { x[ix(i, 1)] };
        x[ix(i, N - 1)] = if b == 2 {
            -x[ix(i, N - 2)]
        } else {
            x[ix(i, N - 2)]
        };
    }
    for j in 1..(N - 1) {
        x[ix(0, j)] = if b == 1 { -x[ix(1, j)] } else { x[ix(1, j)] };
        x[ix(N - 1, j)] = if b == 1 {
            -x[ix(N - 2, j)]
        } else {
            x[ix(N - 2, j)]
        };
    }

    x[ix(0, 0)] = 0.5f64 * (x[ix(1, 0)] + x[ix(0, 1)]);
    x[ix(0, N - 1)] = 0.5f64 * (x[ix(1, N - 1)] + x[ix(0, N - 2)]);
    x[ix(N - 1, 0)] = 0.5f64 * (x[ix(N - 2, 0)] + x[ix(N - 1, 1)]);
    x[ix(N - 1, N - 1)] = 0.5f64 * (x[ix(N - 2, N - 1)] + x[ix(N - 1, N - 2)]);
}

/**
 * linear solving
 */
fn lin_solve(b: i64, x: &mut [f64; NSIZE], x0: &mut [f64; NSIZE], a: f64, c: f64) {
    let c_recip = 1. / c;
    for _ in 0..ITER {
        for j in 1..(N - 1) {
            for i in 1..(N - 1) {
                x[ix(i, j)] = (x0[ix(i, j)]
                    + a * (x[ix(i + 1, j)] + x[ix(i - 1, j)] + x[ix(i, j + 1)] + x[ix(i, j - 1)]))
                    * c_recip;
            }
        }
        set_bnd(b, x);
    }
}

/**
 * diffuse density and velocity
 */
fn diffuse(b: i64, x: &mut [f64; NSIZE], x0: &mut [f64; NSIZE], diff: f64, dt: f64) {
    let a = dt * diff * i64::pow(N - 2, 2) as f64;
    lin_solve(b, x, x0, a, 1. + 4. * a);
}

/**
 * project density and velocity
 */
fn project(
    vx: &mut [f64; NSIZE],
    vy: &mut [f64; NSIZE],
    p: &mut [f64; NSIZE],
    div: &mut [f64; NSIZE],
) {
    for j in 1..(N - 1) {
        for i in 1..(N - 1) {
            div[ix(i, j)] = -0.5f64
                * (vx[ix(i + 1, j)] - vx[ix(i - 1, j)] + vy[ix(i, j + 1)] - vy[ix(i, j - 1)])
                / N as f64;
            p[ix(i, j)] = 0.;
        }
    }

    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1., 4.);

    for j in 1..(N - 1) {
        for i in 1..(N - 1) {
            vx[ix(i, j)] -= 0.5f64 * (p[ix(i + 1, j)] - p[ix(i - 1, j)]) * N as f64;
            vy[ix(i, j)] -= 0.5f64 * (p[ix(i, j + 1)] - p[ix(i, j - 1)]) * N as f64;
        }
    }

    set_bnd(1, vx);
    set_bnd(2, vy);
}

/**
 * short for advect x velocity
 */
fn advect_x(b: i64, d: &mut [f64; NSIZE], vx: &mut [f64; NSIZE], vy: &mut [f64; NSIZE], dt: f64) {
    let (mut i0, mut i1, mut j0, mut j1): (f64, f64, f64, f64);

    let dtx = dt * (N - 2) as f64;
    let dty = dt * (N - 2) as f64;

    let (mut s0, mut s1, mut t0, mut t1): (f64, f64, f64, f64);
    let (mut tmp1, mut tmp2, mut x, mut y): (f64, f64, f64, f64);

    let nfloat = N as f64;
    let (mut ifloat, mut jfloat): (f64, f64);

    for j in 1..(N - 1) {
        jfloat = j as f64;
        for i in 1..(N - 1) {
            ifloat = i as f64;

            tmp1 = dtx * vx[ix(i, j)];
            tmp2 = dty * vy[ix(i, j)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            x = if x < 0.5f64 {
                0.5f64
            } else if x > nfloat + 0.5f64 {
                nfloat + 0.5f64
            } else {
                x
            };
            i0 = f64::floor(x);
            i1 = i0 + 1.0f64;
            y = if y < 0.5f64 {
                0.5f64
            } else if y > nfloat + 0.5f64 {
                nfloat + 0.5f64
            } else {
                y
            };
            j0 = f64::floor(y);
            j1 = j0 + 1.0f64;

            s1 = x - i0;
            s0 = 1.0f64 - s1;
            t1 = y - j0;
            t0 = 1.0f64 - t1;

            let i0i = i0 as i64;
            let i1i = i1 as i64;
            let j0i = j0 as i64;
            let j1i = j1 as i64;

            d[ix(i, j)] = s0 * (t0 * vx[ix(i0i, j0i)] + t1 * vx[ix(i0i, j1i)])
                + s1 * (t0 * vx[ix(i1i, j0i)] + t1 * vx[ix(i1i, j1i)]);
        }
    }
    set_bnd(b, d);
}

/**
 * short for advect y velocity
 */
fn advect_y(b: i64, d: &mut [f64; NSIZE], vx: &mut [f64; NSIZE], vy: &mut [f64; NSIZE], dt: f64) {
    let (mut i0, mut i1, mut j0, mut j1): (f64, f64, f64, f64);

    let dtx = dt * (N - 2) as f64;
    let dty = dt * (N - 2) as f64;

    let (mut s0, mut s1, mut t0, mut t1): (f64, f64, f64, f64);
    let (mut tmp1, mut tmp2, mut x, mut y): (f64, f64, f64, f64);

    let nfloat = N as f64;
    let (mut ifloat, mut jfloat): (f64, f64);

    for j in 1..(N - 1) {
        jfloat = j as f64;
        for i in 1..(N - 1) {
            ifloat = i as f64;

            tmp1 = dtx * vx[ix(i, j)];
            tmp2 = dty * vy[ix(i, j)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            x = if x < 0.5f64 {
                0.5f64
            } else if x > nfloat + 0.5f64 {
                nfloat + 0.5f64
            } else {
                x
            };
            i0 = f64::floor(x);
            i1 = i0 + 1.0f64;
            y = if y < 0.5f64 {
                0.5f64
            } else if y > nfloat + 0.5f64 {
                nfloat + 0.5f64
            } else {
                y
            };
            j0 = f64::floor(y);
            j1 = j0 + 1.0f64;

            s1 = x - i0;
            s0 = 1.0f64 - s1;
            t1 = y - j0;
            t0 = 1.0f64 - t1;

            let i0i = i0 as i64;
            let i1i = i1 as i64;
            let j0i = j0 as i64;
            let j1i = j1 as i64;

            d[ix(i, j)] = s0 * (t0 * vy[ix(i0i, j0i)] + t1 * vy[ix(i0i, j1i)])
                + s1 * (t0 * vy[ix(i1i, j0i)] + t1 * vy[ix(i1i, j1i)]);
        }
    }
    set_bnd(b, d);
}

/**
 * advect density and maybe velocity
 */
fn advect(
    b: i64,
    d: &mut [f64; NSIZE],
    d0: &mut [f64; NSIZE],
    vx: &mut [f64; NSIZE],
    vy: &mut [f64; NSIZE],
    dt: f64,
) {
    let (mut i0, mut i1, mut j0, mut j1): (f64, f64, f64, f64);

    let dtx = dt * (N - 2) as f64;
    let dty = dt * (N - 2) as f64;

    let (mut s0, mut s1, mut t0, mut t1): (f64, f64, f64, f64);
    let (mut tmp1, mut tmp2, mut x, mut y): (f64, f64, f64, f64);

    let nfloat = N as f64;
    let (mut ifloat, mut jfloat): (f64, f64);

    for j in 1..(N - 1) {
        jfloat = j as f64;
        for i in 1..(N - 1) {
            ifloat = i as f64;

            tmp1 = dtx * vx[ix(i, j)];
            tmp2 = dty * vy[ix(i, j)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            x = if x < 0.5f64 {
                0.5f64
            } else if x > nfloat + 0.5f64 {
                nfloat + 0.5f64
            } else {
                x
            };
            i0 = f64::floor(x);
            i1 = i0 + 1.0f64;
            y = if y < 0.5f64 {
                0.5f64
            } else if y > nfloat + 0.5f64 {
                nfloat + 0.5f64
            } else {
                y
            };
            j0 = f64::floor(y);
            j1 = j0 + 1.0f64;

            s1 = x - i0;
            s0 = 1.0f64 - s1;
            t1 = y - j0;
            t0 = 1.0f64 - t1;

            let i0i = i0 as i64;
            let i1i = i1 as i64;
            let j0i = j0 as i64;
            let j1i = j1 as i64;

            d[ix(i, j)] = s0 * (t0 * d0[ix(i0i, j0i)] + t1 * d0[ix(i0i, j1i)])
                + s1 * (t0 * d0[ix(i1i, j0i)] + t1 * d0[ix(i1i, j1i)]);
        }
    }
    set_bnd(b, d);
}

/**
 * Fluid Cube Map
 */
pub struct FluidCube {
    pub size: i64,
    pub dt: f64,
    pub diff: f64,
    pub visc: f64,

    pub s: [f64; NSIZE],
    pub density: [f64; NSIZE],

    pub vx: [f64; NSIZE],
    pub vy: [f64; NSIZE],

    pub vx0: [f64; NSIZE],
    pub vy0: [f64; NSIZE],
}

impl FluidCube {
    /**
     * new Fluid Cube
     */
    pub fn new(dt: f64, diff: f64, visc: f64) -> FluidCube {
        FluidCube {
            size: N,
            dt,
            diff,
            visc,

            s: [0.; NSIZE],
            density: [0.; NSIZE],

            vx: [0.; NSIZE],
            vy: [0.; NSIZE],

            vx0: [0.; NSIZE],
            vy0: [0.; NSIZE],
        }
    }

    /**
     * add some density
     */
    pub fn add_density(&mut self, x: i64, y: i64, amount: f64) {
        let index = ix(x, y);
        self.density[index] += amount;
    }

    /**
     * add some velocity
     */
    pub fn add_velocity(&mut self, x: i64, y: i64, amount_x: f64, amount_y: f64) {
        let index = ix(x, y);
        self.vx[index] += amount_x;
        self.vy[index] += amount_y;
    }

    /**
     * step forward in time
     */
    pub fn step(&mut self) {
        let visc = self.visc;
        let diff = self.diff;
        let dt = self.dt;
        let vx = &mut self.vx;
        let vy = &mut self.vy;
        let vx0 = &mut self.vx0;
        let vy0 = &mut self.vy0;
        let s = &mut self.s;
        let density = &mut self.density;

        diffuse(1, vx0, vx, visc, dt);
        diffuse(1, vy0, vy, visc, dt);

        project(vx0, vy0, vx, vy);

        advect_x(1, vx, vx0, vy0, dt);
        advect_y(2, vy, vx0, vy0, dt);

        project(vx, vy, vx0, vy0);

        diffuse(0, s, density, diff, dt);
        advect(0, density, s, vx, vy, dt);
    }

    /**
     * fade density and keep it between 0. and 255.
     */
    pub fn fade_d(&mut self) {
        let f = 0.5;
        for i in 0..(N * N) {
            let d = self.density[i as usize];
            self.density[i as usize] = if d - f < 0. {
                0.
            } else if d - f > 255. {
                255.
            } else {
                d - f
            };
        }
    }
}
