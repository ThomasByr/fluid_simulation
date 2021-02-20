# Fluid Simulation using Rust

1. [Why this subject ?](#why-this-subject-)
2. [From C to Processing to Rust](#from-c-to-processing-to-rust)
3. [Dependencies](#dependencies)
4. [Changelog](#changelog)

> Inspired by Fluid Simulation for Dummies, [see Mike Ash paper](https://mikeash.com/pyblog/fluid-simulation-for-dummies.html)
## Why this subject ?
Well it appears that fluid simulation is off topic for french *classes préparatoires aux grandes écoles*. So why not try it ? Furthermore, being borderline beyond the scope does not mean this subject is not susceptible to be found in written exams. Besides that, Navier-Stokes equations are not unfamiliar and linear solving is something I should be aware of. So there it is. Not very accurate simulation but simulation nonetheless. Tweak some constants in [main.rs](src/main.rs) and run project with ``cargo run --release``. This project is under the GNU-GPL v3 public license.

## From C to Processing to Rust
What have changed ?? Not much. I thought it would be fun to give Rust a try. Many prople have noticed a great increase in performances. The first idea came as mentionned from Mike Ash paper, which has inspired Daniel Shiffman aka *The Coding Train* and creator of the Processing language. Minor changes have occured : C needed the user to manually allocate and free memory, Java and Rust don't. There is no need to work with pointers in Java but we do need mutable references in Rust if we want to dynamically change values in arrays. The advect function also now comes in 3 different forms because Rust does not allow more than one mutable reference borrow at a time (for instance I could not call a function of twice the same mutable reference). Speaking about performances, that damn line
```rs
rectangle([d, d, d, 1.], sq, c.transform, g);
```
is costing me a hundred frames on my potato laptop. I so have about 20 to 30 and I manage to double it when I switch device. But why is drawing a thousand rectangles so slow ?

## Dependencies
The project will need RobotoMono Thin font which can be found inside of the assets folder alongside with its usage license. Inside of [Cargo.toml](Cargo.toml) are soft dependencies wich will be automatically met by cargo when first building the project. So there is **Piston Window** used for the rendering and based on OpenGL and Glutin Window, **noise** which implements 2 to 4 dimentionnal Perlin Noise, **rand** for standard randomness, **fps counter** which is a small structure which coounts frames and **find colder** to find folder the font is located in. The whole thing runs with Rust 2018.

## Changelog
* v0.1.0
  * initial release and commit
  * why is Piston so slow ?
  * does not draw black rectangles anymore
  * some comments have been added under torture
