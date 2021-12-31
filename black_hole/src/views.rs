use crate::{Christoffels, Photon, State, Vec3, Vec4};
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
    normal: Vec3,
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


fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1],
     a[2] * b[0] - a[0] * b[2],
     a[0] * b[1] - a[1] * b[0]]
}

fn add3 (a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn mul3 (a: Vec3, f: f64) -> Vec3 {
    [a[0] * f, a[1] * f, a[2] * f]
}

fn opposite(v: Vec3) -> Vec3 {
    [-v[0], -v[1], -v[2]]
}

fn length(v: Vec3) -> f64 {
    (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt()
}

fn get_initial_state<T>(pos: Vec3, cart_vel: Vec3) -> (Vec4, Vec4)
    where T : Christoffels {
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
    let time = (-ldotv - (ldotv * ldotv - metric[(0,0)]*vgv).sqrt()) / metric[(0,0)];
    (pos4, [time, vel[0], vel[1], vel[2]])
}

fn rotate_pos_vel(pos: Vec4, vel: Vec4, normal: Vec3, psi: f64) -> (Vec3, Vec3) {
    let cart_pos = [pos[1] * pos[2].sin() * pos[3].cos(),
                    pos[1] * pos[2].sin() * pos[3].sin(),
                    pos[1] * pos[2].cos()];
    let dot = (cart_pos[0] * normal[0] + cart_pos[1] * normal[1] + cart_pos[2] * normal[2])
            / (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    if dot <= 0.0 {
        return rotate_pos_vel(pos, vel, opposite(normal), -psi);
    }

    let scale_pos = [cart_pos[0] / dot, cart_pos[1] / dot, cart_pos[2] / dot];
    let x = add3(scale_pos, opposite(normal));
    let y = mul3(cross(normal, x), 1.0 / length(normal));

    let offset = add3(mul3(x, psi.cos()), mul3(y, psi.sin()));
    let new_pos = add3(normal, offset);
    
    //[pos[0], 1.0, theta, phi]
    (mul3(new_pos, pos[1].abs() / length(new_pos)), [0.0, 0.0, 0.0])// TO DO: Vel is not rotated!!!
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


impl<C> ViewLine<C> where C: Christoffels {
    pub fn new<T>(pos: Vec3, normal: Vec3, up: Vec3, width: u32, height: u32) -> ViewLine<T>
        where T: Christoffels {
        // Pos is in spherical coords; normal, up are cartesian (local)
        let right = cross(normal, up);
        let right = mul3(right, length(up) / length(right) * (width as f64) / (height as f64));

        let mut photons = Vec::new();
        let new_length = (f64::sqrt((width * width + height * height) as f64) / 2.0) as u32 + 2;
        for px in 0..new_length {
            let right_amount = match new_length > 1 {
                true => px as f64 / (new_length as f64 - 1.0),
                false => 0.0,
            };
            let right_amount = right_amount * new_length as f64 / width as f64;

            let cart_vel = [
                right[0] * right_amount + normal[0],
                right[1] * right_amount + normal[1],
                right[2] * right_amount + normal[2],
            ];

            let (pos4, vel4) = get_initial_state::<T>(pos, cart_vel);

            photons.push(Photon::<T>::new(pos4, vel4));
        }
        ViewLine { width, height, normal, photons }
    }

    fn run_slice<T>(slice: &mut [Photon<T>], max_iterations: u32, dtau: f64) -> bool
        where T: Christoffels + std::fmt::Debug {
        
        for iteration in 0..max_iterations {
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

impl<C> View<C> for ViewLine<C> where C: Christoffels + Send + std::fmt::Debug {
    fn run(&mut self, max_iterations: u32, dtau: f64) -> bool{
        let slice_length = self.photons.len() / NUM_THREADS + 1;
        let mut chunks = self.photons.chunks_mut(slice_length);
        
        if NUM_THREADS > 1 {
            // Threaded code
            crossbeam::scope(|scope| {
                let mut handles = Vec::new();
                handles.reserve(NUM_THREADS);

                for chunk in &mut chunks {
                    handles.push(scope.spawn(move |_| {
                        Self::run_slice(chunk, max_iterations, dtau)
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
                result = result && Self::run_slice(chunk, max_iterations, dtau);
            }
            result
        }
    }

    fn render(self, skybox: Skybox) -> ColorArray {
        let mut pixels = Vec::new();

        for px in 0..self.width {
            let mut line = Vec::new();
            for py in 0..self.height {
                let rpx = px as f64 - (self.width as f64 - 1.0) / 2.0;
                let rpy = py as f64 - (self.height as f64 - 1.0) / 2.0;
                let psi = f64::atan2(rpy, rpx);
                let radial_length = f64::sqrt(rpx * rpx + rpy * rpy);
                let fraction = radial_length - (radial_length as usize) as f64;
                let radial_length = radial_length as usize;
                let photon_below = &self.photons[radial_length];
                let photon_above = &self.photons[radial_length + 1];

                if let State::Dead = photon_below.state {
                    line.push((0, 0, 0));
                    continue;
                }
                if let State::Dead = photon_above.state {
                    line.push((0, 0, 0));
                    continue;
                }
                
                let (below_pos, below_vel) = rotate_pos_vel(photon_below.pos, photon_below.vel, self.normal, psi);
                let (above_pos, above_vel) = rotate_pos_vel(photon_above.pos, photon_above.vel, self.normal, psi);

                let cart_pos = add3(mul3(below_pos, 1.0 - fraction), mul3(above_pos, fraction));
                let cart_vel = add3(mul3(below_vel, 1.0 - fraction), mul3(above_vel, fraction));
                let pos_time = photon_below.pos[0] * (1.0 - fraction) + photon_above.pos[0] * fraction;
                let vel_time = photon_below.vel[0] * (1.0 - fraction) + photon_above.vel[0] * fraction;

                line.push(skybox.call(
                    [pos_time,
                     length(cart_pos),
                     f64::acos(cart_pos[2] / length(cart_pos)),
                     f64::atan2(cart_pos[1], cart_pos[0])],
                    [0.0, 0.0, 0.0, 0.0]));
            }
            pixels.push(line);
        }
        ColorArray::new(pixels, self.width, self.height)
    }
}


impl<C> ViewPlane<C> where C: Christoffels {
    pub fn new<T>(pos: Vec3, normal: Vec3, up: Vec3, width: u32, height: u32) -> ViewPlane<T>
        where T: Christoffels {
        // Pos is in spherical coords; normal, up are cartesian (local)
        let right = cross(normal, up);
        let right = mul3(right, length(up) / length(right) * (width as f64) / (height as f64));
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

                let (pos4, vel4) = get_initial_state::<T>(pos, cart_vel);

                photons.push(Photon::<T>::new(pos4, vel4));
            }
        }
        ViewPlane { width, height, photons }
    }

    fn run_slice<T>(slice: &mut [Photon<T>], max_iterations: u32, dtau: f64) -> bool
        where T: Christoffels + std::fmt::Debug {
        
        for iteration in 0..max_iterations {
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
        
        if NUM_THREADS > 1 {
            // Threaded code
            crossbeam::scope(|scope| {
                let mut handles = Vec::new();
                handles.reserve(NUM_THREADS);

                for chunk in &mut chunks {
                    handles.push(scope.spawn(move |_| {
                        Self::run_slice(chunk, max_iterations, dtau)
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
                result = result && Self::run_slice(chunk, max_iterations, dtau);
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
                let photon = self.photons.pop().expect("Photon list was too short");
                if let State::Dead = photon.state {
                    line.push((0, 0, 0));
                    continue;
                }
                line.push(skybox.call(photon.pos, photon.vel));
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