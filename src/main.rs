mod fluid;

extern crate find_folder;
extern crate fps_counter;
extern crate noise;
extern crate piston_window;

use fluid::*;
use fps_counter::*;
use noise::*;
use piston_window::*;
use rand::Rng;

const N: i64 = 128; // number of cells
const ITER: i64 = 16; // times to iter linear solving
const SCALE: i64 = 2; // draw scale
const PI: f64 = std::f64::consts::PI; // well...

fn main() {
    let mut window: PistonWindow =
        WindowSettings::new("Fluid Simulation", [(N * SCALE) as u32, (N * SCALE) as u32])
            .exit_on_esc(true)
            .resizable(false)
            .build()
            .unwrap_or_else(|e| panic!("Failed to build PistonWindow: {}", e));
    let noise = Perlin::new();
    let mut rng = rand::thread_rng();
    let mut fluid = FluidCube::new(0.1, 0., 1e-7);
    let mut fps_counter = FPSCounter::new();
    let mut frames: usize;

    let assets = find_folder::Search::ParentsThenKids(3, 3)
        .for_folder("assets")
        .unwrap();
    let mut glyphs = window
        .load_font(assets.join("RobotoMono-Thin.ttf"))
        .unwrap();

    let mut t_offset = 0.;
    let c_x = N / 2;
    let c_y = N / 2;

    // main loop (drawing performance issues...)
    window.set_lazy(false);
    window.set_bench_mode(true);
    while let Some(e) = window.next() {
        frames = fps_counter.tick(); // get current fps
        t_offset += 1e-2; // offset for perlin noise

        // adding density to fluid cube map
        for i in -1..2 {
            for j in -1..2 {
                let a = rng.gen_range(50., 150.);
                fluid.add_density(c_x + i, c_y + j, a);
            }
        }
        // adding velocity to fluid cube map
        for _i in 0..2 {
            let angle = noise.get([t_offset, 0.]) * 2. * PI;
            let a_x = angle.cos();
            let a_y = angle.sin();
            fluid.add_velocity(c_x, c_y, a_x, a_y);
        }

        fluid.step(); // step forward in time
        fluid.fade_d(); // fade density

        // render density as black and white rectangles
        window.draw_2d(&e, |c, g, device| {
            let transform = c.transform.trans(10., 20.);
            clear([0.; 4], g);
            // for all cubes
            for i in 0..N {
                for j in 0..N {
                    let x = (i * SCALE) as f64;
                    let y = (j * SCALE) as f64;
                    let d = (fluid.density[ix(i, j)] / 255.) as f32;
                    let scale = SCALE as f64;
                    if d > 0. {
                        // would comment 100
                        let sq = rectangle::square(x, y, scale);
                        rectangle([d, d, d, 1.], sq, c.transform, g);
                    }
                }
            }
            // render fps on window as text
            text::Text::new_color([1.; 4], 7)
                .draw(
                    &mut format!("fps : {}", frames),
                    &mut glyphs,
                    &c.draw_state,
                    transform,
                    g,
                )
                .unwrap();
            glyphs.factory.encoder.flush(device); // update glyphs before rendering
        });
    }
}
