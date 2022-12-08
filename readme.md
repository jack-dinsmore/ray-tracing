# Compton Scattering of Accretion Disk Emission off Black Hole Coronas in Curved Spacetime

Jack T. Dinsmore. November 20, 2022

Finished Dec 8, 2022

Physics 360 Final Project. Stanford University.

## Description

The purpose of this project is to generate images of a simple black hole system consisting of an accretion disk and a Compton-scattering corona. Images are generated for the X-ray and optical bands by numerically integrating the equations of motion of light from the observer to the source under the laws of General Relativity. This technique is called ray-tracing.

For a description of the projects and a summary, see the **writeup** directory. For my final presentation, see **presentation**. Generated data is contained in the **data** directory and the code required to generate image from that data is in **imager.**

Finally, the code to actually do the ray-tracing is in the **raytracer** directory. It is flexible enough to use any relativistic metric, though Minkowsky, Schwarzschild, Kerr, and Morris-Thorne are already implemented. Likewise, other sources of light may be added to the code by modifying the **source.rs** and **observer.rs** files. The code is written in Rust &mdash; a relativity recent programming language which is blazing fast and safer than any other alternative of the same performance (i.e., C and C++). The Rust documentation is phenomenal and you should definitely try the online [Rust book](https://doc.rust-lang.org/book/) if you would like to learn it. Experience with another type-safe language like C, C++, or Java is helpful, since Rust has a steep learning curve.



## Usage

To run the code, in the **raytracer** directory run
```
cargo run --release
```
To generate images, in the **imager** directory run 
```
python image.py
```

## Changelog

### Dec 7
Initial Commit. Original code moved to **old** directory. Added new code, data, and summary.