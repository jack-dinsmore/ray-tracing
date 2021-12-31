use crate::{Christoffels, Photon, State, Vec3};
use crate::skybox::Skybox;
use bmp::Pixel;

const NUM_THREADS : usize = 6;

pub type Color = (u8, u8, u8);

pub trait View<C> {
    fn run(&mut self, max_iterations: u32, dtau: f64) -> bool;
    fn render(self, skybox: Skybox) -> ColorArray;
}

pub struct ViewLine<C> {
    width: u32,
    height: u32,
    photons: Vec<Photon<C>>,
}

pub struct ViewPlane<C> {
    width: u32,
    height: u32,
    pub photons: Vec<Photon<C>>,
}

pub struct ViewCube<C> {
    length: u32,
    planes: [ViewPlane<C>; 6],
}

pub struct ColorArray {
    pixels: Vec<Vec<Color>>,
    width: u32,
    height: u32,
}


fn cross (a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1],
     a[2] * b[0] - a[0] * b[2],
     a[0] * b[1] - a[1] * b[0]]
}

fn opposite(v: Vec3) -> Vec3 {
    [-v[0], -v[1], -v[2]]
}


impl ColorArray {
    fn new(pixels: Vec<Vec<Color>>, width: u32, height: u32) -> ColorArray {
        ColorArray {pixels, width, height}
    }

    pub fn save(&self, filename: &str) -> std::io::Result<()> {
        let mut img = bmp::Image::new(self.width, self.height);

        for (x, y) in img.coordinates() {
            let color = self.pixels[x as usize][y as usize];
            img.set_pixel(x, y, bmp::px!(color.0, color.1, color.2));
        }
        img.save(filename)
    }
}

impl<C> ViewPlane<C> where C: Christoffels {
    pub fn new<T>(pos: Vec3, normal: Vec3, up: Vec3, width: u32, height: u32) -> ViewPlane<T>
        where T: Christoffels {
        // Pos is in spherical coords; normal, up are cartesian (local)
        let right = cross(normal, up);
        let mut photons = Vec::new();
        for px in 0..width {
            for py in 0..height {
                let right_amount = match width > 1 {
                    true => (px as f64 - (width as f64 - 1.0) / 2.0) / (width as f64 - 1.0),
                    false => 0.0,
                };
                let up_amount = match height > 1 {
                    true => (py as f64 - (height as f64 - 1.0) / 2.0) / (height as f64 - 1.0),
                    false => 0.0,
                };
                let cart_vel = [
                    up[0] * up_amount + right[0] * right_amount + normal[0],
                    up[1] * up_amount + right[1] * right_amount + normal[1],
                    up[2] * up_amount + right[2] * right_amount + normal[2],
                ];

                let st = pos[1].sin();
                let ct = pos[1].cos();
                let sp = pos[2].sin();
                let cp = pos[2].cos();
                let mut vel = [ st * (cart_vel[0] * cp + cart_vel[1] * sp) + cart_vel[2] * ct,
                           (ct * (cart_vel[0] * cp + cart_vel[1] * sp) - cart_vel[2] * st) / pos[0],
                           (cart_vel[1] * cp - cart_vel[0] * sp) / (pos[0] * st)
                ]; //dr / dt, dtheta / dt, dphi / dt
                let norm = f64::sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
                vel[0] /= norm;
                vel[1] /= norm;
                vel[2] /= norm;

                let pos4 = [0.0, pos[0], pos[1], pos[2]];
                let metric = T::get_metric(pos4);

                let ldotv = metric[(0, 1)] * vel[0] + metric[(0, 2)] * vel[1] + metric[(0, 3)] * vel[2];
                let mut vgv = 0.0;
                for i in 1..4 {
                    for j in 1..4 {
                        vgv += vel[i-1] * vel[j-1] * metric[(i, j)];
                    }
                }
                let time = (-ldotv + (ldotv * ldotv - metric[(0,0)]*vgv).sqrt()) / metric[(0,0)];

                photons.push(Photon::<T>::new(
                    pos4,
                    [time, vel[0], vel[1], vel[2]],
                ));
            }
        }
        ViewPlane { width, height, photons }
    }

    fn run_slice<T>(slice: &mut [Photon<T>], max_iterations: u32, dtau: f64) -> bool
        where T: Christoffels + std::fmt::Debug {
        
        for iteration in 0..max_iterations {
            println!("{:?}", slice[0]);
            let mut ended = true;
            for photon in slice.iter_mut() {
                match photon.update(iteration, dtau) {
                    State::Running => ended = false,
                    State::Dead => (),
                    State::Escape => ()
                };
            }
            if ended {
                return true;
            }
        }
        println!("Simulation did not succeed");
        false
    }
}

impl<C> View<C> for ViewPlane<C> where C: Christoffels + Send + std::fmt::Debug {
    fn run(&mut self, max_iterations: u32, dtau: f64) -> bool{
        let slice_length = self.photons.len() / NUM_THREADS + 1;
        let mut chunks = self.photons.chunks_mut(slice_length);
        println!("{:?}", chunks);
        
        if NUM_THREADS > 1 {
            // Threaded code
            crossbeam::scope(|scope| {
                let mut handles = Vec::new();
                handles.reserve(NUM_THREADS);

                for chunk in &mut chunks {
                    handles.push(scope.spawn(move |_| {
                        ViewPlane::<C>::run_slice(chunk, max_iterations, dtau)
                    }));
                }

                let mut result = true;
                for handle in handles {
                    result = result && handle.join().expect("Thread has panicked");
                }
                result
            })
            .expect("A child thread panicked")
        }
        else {
            // Unthreaded code
            let mut result = true;
            for chunk in &mut chunks {
                result = result && ViewPlane::<C>::run_slice(chunk, max_iterations, dtau);
            }
            result
        }
    }

    fn render(mut self, skybox: Skybox) -> ColorArray {
        let mut pixels = Vec::new();
        self.photons.reverse();
        for _ in 0..self.width {
            let mut line = Vec::new();
            for _ in 0..self.height {
                line.push(skybox.call(self.photons.pop().expect("Photon list was too short")));
            }
            pixels.push(line);
        }
        ColorArray::new(pixels, self.width, self.height)
    }
}

/*impl<C> ViewCube<C> where C: Christoffels {
    pub fn new<T>(pos: Vec3, top: Vec3, forward: Vec3, length: u32) -> ViewCube<T>
        where T: Christoffels {
        let right = cross(forward, top);
        let left = opposite(right);
        let bottom = opposite(top);
        let back = opposite(forward);
        ViewCube { length,
            planes: [ViewPlane::<T>::new(pos, top, forward, length, length),
                     ViewPlane::<T>::new(pos, bottom, forward, length, length),
                     ViewPlane::<T>::new(pos, left, forward, length, length),
                     ViewPlane::<T>::new(pos, right, forward, length, length),
                     ViewPlane::<T>::new(pos, forward, top, length, length),
                     ViewPlane::<T>::new(pos, back, top, length, length)]}
    }
}

impl<C> View<C> for ViewCube<C> where C: Christoffels {
    fn update(&mut self, iteration: u32, dtau: f64) -> bool {
        let mut ended = true;
        for plane in &mut self.planes {
            ended = ended && plane.update(iteration, dtau);
        }
        ended
    }

    fn get_output(self) -> ViewOutput<C> {
        let mut output = Vec::new();
        for plane in self.planes {
            output.append(&mut plane.get_output());
        }
        output
    }
}*/