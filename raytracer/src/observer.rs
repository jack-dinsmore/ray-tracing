use std::mem::MaybeUninit;

use crate::util::*;
use crate::metrics::{Metric, State};
use crate::source::{AccretionDisk, RESCALE_FOR_XRAY};

const SIZE: usize = 32; // 32
pub const PHOTON_BATCH_SIZE: usize = 0x100;
pub const CORONA_INTERACTION: usize = 0x08;
const NUM_COVERS_PER: usize = 10; // 10
const NUM_SLICES: usize = 7; // 7
pub const IS_KERR: bool = true; // Determines whether to cut disk
const LENGTH_SCALE: f64 = 0.9; // Decrease in size ratio for each rectangle of the image

pub const WIDTH: usize = 16 * SIZE;
pub const HEIGHT: usize = 9 * SIZE;
const MAX_COUNTS: usize = WIDTH * HEIGHT * NUM_SLICES * NUM_COVERS_PER;
const OFFSET_CHECK: usize = 0x10;
const TEMP_RECORD: usize = 8;
const EPSILON: f64 = 1e-3; // Addition to position to prevent singularities for very close to event horizon.
const POLE_PROTECTION: f64 = 3.0; // Controls the delta t factor near theta = 0
const IMAGE_WIDTH: f64 = 1.8; // 1.8
const FALLOFF_SIG: f64 = 1000.0;
const PIXEL_WIDTH: f64 = IMAGE_WIDTH / WIDTH as f64;

pub trait Observer {
    fn update(&mut self, photon_data: &Vec<PhotonData>) -> u32;
    fn next_photons(&mut self) -> [Option<Photon>; PHOTON_BATCH_SIZE];
}

pub struct Simple<M: Metric> {
    pos: Vec4, // t, r, theta, phi
    look: Vec3,
    up: Vec3,
    right: Vec3,
    num_counts: usize,
    metric: M,
}

impl<M: Metric> Simple<M> {
    pub fn new(pos: Vec3, look: Vec3, metric: M) -> Self {
        let right = normalize(cross(look, [0.0, 0.0, 1.0]));
        let up = normalize(cross(right, look));
        Self {
            pos: [0.0, pos[0], pos[1], pos[2]],
            look,
            up,
            right,
            num_counts: 0,
            metric,
        }
    }

    fn theta_phi_to_vel(&self, pos: Vec4, dir: (f64, f64)) -> Vec4 {
        let v3 = normalize(add3(self.look, add3(mul3(self.up, dir.0), mul3(self.right, dir.1))));
        
        // Use the reverse velocity
        let v3 = cart_to_spher_vel([pos[1], pos[2], pos[3]], opposite(v3));
        get_vel_from_metric(v3, &self.metric.get_metric(self.pos))
    }

    fn get_theta_phi_from_count(&self, mut count: usize) -> Option<((f64, f64), (usize, usize))> {
        let mut width = WIDTH;
        let mut height = HEIGHT;
        let mut index = 0;
        count /= NUM_COVERS_PER;
        loop {
            // Loop through the number of trials in each pixel
            if count >= width * height {
                count -= width * height;
                width = (width as f64 * LENGTH_SCALE) as usize;
                height = (height as f64 * LENGTH_SCALE) as usize;
                index += 1;
                if index == NUM_SLICES || height <= 2 {
                    return None;
                }
                continue;
            }
            let i = count / width;
            let j = count % width;
            return Some(((
                PIXEL_WIDTH * (i as f64 - (height / 2) as f64 + 0.001 * rand::random::<f64>()),
                PIXEL_WIDTH * (j as f64 - (width / 2) as f64 + 0.001 * rand::random::<f64>())
            ), (i + (HEIGHT - height) / 2, j + (WIDTH - width) / 2)));
        }
    }
}

impl<M: Metric> Observer for Simple<M> {
    fn update(&mut self, _photon_data: &Vec<PhotonData>) -> u32 {
        ((self.num_counts * 100) / MAX_COUNTS) as u32
    }

    fn next_photons(&mut self) -> [Option<Photon>; PHOTON_BATCH_SIZE] {
        let mut array: [MaybeUninit<Option<Photon>>; PHOTON_BATCH_SIZE] = unsafe { MaybeUninit::uninit().assume_init() };

        for (i, element) in array.iter_mut().enumerate() {
            *element = match self.get_theta_phi_from_count(self.num_counts + i) {
                Some((theta_phi, pixel)) => MaybeUninit::new(
                    Some(Photon::new(
                        self.pos,
                        pixel,
                        self.theta_phi_to_vel(self.pos, theta_phi),
                    ))
                ),
                None => MaybeUninit::new(None)
            };
        }

        self.num_counts += PHOTON_BATCH_SIZE;
        unsafe { std::mem::transmute::<_, [Option<Photon>; PHOTON_BATCH_SIZE]>(array) }
    }
}

pub struct Photon {
    pos: Vec4,
    vel: Vec4,
    pixel: (usize, usize),
    temps: [(f64, f64); TEMP_RECORD], // Temp, luminosity
    temp_index: usize,
    depth: f64,
    compton_scatter: Option<f64>, // Compton shift
}

#[derive(Clone)]
pub struct PhotonData {
    pub optical_color: (f64, f64, f64),
    pub xray_color: (f64, f64, f64),
    pub pixel: (usize, usize),
}

impl PhotonData {
    fn temp_to_color(temp: f64) -> (f64, f64, f64) {
        const TEMP_TO_COLOR: [(u8, u8, u8); 111] = [
            (255, 56, 0),
            (255, 71, 0),
            (255, 83, 0),
            (255, 93, 0),
            (255, 101, 0),
            (255, 109, 0),
            (255, 115, 0),
            (255, 121, 0),
            (255, 126, 0),
            (255, 131, 0),
            (255, 138, 18),
            (255, 142, 33),
            (255, 147, 44),
            (255, 152, 54),
            (255, 157, 63),
            (255, 161, 72),
            (255, 165, 79),
            (255, 169, 87),
            (255, 173, 94),
            (255, 177, 101),
            (255, 180, 107),
            (255, 184, 114),
            (255, 187, 120),
            (255, 190, 126),
            (255, 193, 132),
            (255, 196, 137),
            (255, 199, 143),
            (255, 201, 148),
            (255, 204, 153),
            (255, 206, 159),
            (255, 209, 163),
            (255, 211, 168),
            (255, 213, 173),
            (255, 215, 177),
            (255, 217, 182),
            (255, 219, 186),
            (255, 221, 190),
            (255, 223, 194),
            (255, 225, 198),
            (255, 227, 202),
            (255, 228, 206),
            (255, 230, 210),
            (255, 232, 213),
            (255, 233, 217),
            (255, 235, 220),
            (255, 236, 224),
            (255, 238, 227),
            (255, 239, 230),
            (255, 240, 233),
            (255, 242, 236),
            (255, 243, 239),
            (255, 244, 242),
            (255, 245, 245),
            (255, 246, 247),
            (255, 248, 251),
            (255, 249, 253),
            (254, 249, 255),
            (252, 247, 255),
            (249, 246, 255),
            (247, 245, 255),
            (245, 243, 255),
            (243, 242, 255),
            (240, 241, 255),
            (239, 240, 255),
            (237, 239, 255),
            (235, 238, 255),
            (233, 237, 255),
            (231, 236, 255),
            (230, 235, 255),
            (228, 234, 255),
            (227, 233, 255),
            (225, 232, 255),
            (224, 231, 255),
            (222, 230, 255),
            (221, 230, 255),
            (220, 229, 255),
            (218, 229, 255),
            (217, 227, 255),
            (216, 227, 255),
            (215, 226, 255),
            (214, 225, 255),
            (212, 225, 255),
            (211, 224, 255),
            (210, 223, 255),
            (209, 223, 255),
            (208, 222, 255),
            (207, 221, 255),
            (207, 221, 255),
            (206, 220, 255),
            (205, 220, 255),
            (207, 218, 255),
            (207, 218, 255),
            (206, 217, 255),
            (205, 217, 255),
            (204, 216, 255),
            (204, 216, 255),
            (203, 215, 255),
            (202, 215, 255),
            (202, 214, 255),
            (201, 214, 255),
            (200, 213, 255),
            (200, 213, 255),
            (199, 212, 255),
            (198, 212, 255),
            (198, 212, 255),
            (197, 211, 255),
            (197, 211, 255),
            (197, 210, 255),
            (196, 210, 255),
            (195, 210, 255),
            (195, 209, 255)
        ];

        let float_index = temp / 100.0 - 10.0;
        if float_index < 0.0 {
            let falloff = (-float_index * float_index / FALLOFF_SIG).exp();
            return (
                TEMP_TO_COLOR[0].0 as f64 / 255.0 * falloff,
                TEMP_TO_COLOR[0].1 as f64 / 255.0 * falloff,
                TEMP_TO_COLOR[0].2 as f64 / 255.0 * falloff
            );
        } else if float_index > TEMP_TO_COLOR.len() as f64 - 1.0 {
            let falloff = (-(float_index - TEMP_TO_COLOR.len() as f64 - 1.0).powi(2) / FALLOFF_SIG).exp();
            return (
                TEMP_TO_COLOR[TEMP_TO_COLOR.len()- 1].0 as f64 / 255.0 * falloff,
                TEMP_TO_COLOR[TEMP_TO_COLOR.len()- 1].1 as f64 / 255.0 * falloff,
                TEMP_TO_COLOR[TEMP_TO_COLOR.len()- 1].2 as f64 / 255.0 * falloff
            );
        }
        let int_index = float_index as usize;
        let frac = float_index - int_index as f64;
        let out = (
            (TEMP_TO_COLOR[int_index + 1].0 as f64 * frac + TEMP_TO_COLOR[int_index].0 as f64 * (1.0 - frac)) / 255.0,
            (TEMP_TO_COLOR[int_index + 1].1 as f64 * frac + TEMP_TO_COLOR[int_index].1 as f64 * (1.0 - frac)) / 255.0,
            (TEMP_TO_COLOR[int_index + 1].2 as f64 * frac + TEMP_TO_COLOR[int_index].2 as f64 * (1.0 - frac)) / 255.0,
        );
        out
    }
}

impl Photon {
    pub fn new(pos: Vec4, pixel: (usize, usize), vel: Vec4) -> Self {
        Self {
            pos,
            vel,
            pixel,
            temps: [(0.0, 0.0); TEMP_RECORD],
            temp_index: 0,
            depth: 1.0,
            compton_scatter: None,
        }
    }

    fn get_data(&self) -> PhotonData {
        let mut optical_color = (0.0, 0.0, 0.0);
        let mut xray_color = (0.0, 0.0, 0.0);
        for temp_index in 0..self.temp_index {
            let (temp, lum) = self.temps[temp_index];
            let temp_color = PhotonData::temp_to_color(temp / 8.62e-5); // Convert to kelvin
            match self.compton_scatter {
                Some(lum_shift) => {
                    xray_color = (
                        xray_color.0 + temp_color.0 * lum * lum_shift / RESCALE_FOR_XRAY,
                        xray_color.1 + temp_color.1 * lum * lum_shift / RESCALE_FOR_XRAY,
                        xray_color.2 + temp_color.2 * lum * lum_shift / RESCALE_FOR_XRAY
                    );
                }
                None => {
                    optical_color = (
                        optical_color.0 + temp_color.0 * lum,
                        optical_color.1 + temp_color.1 * lum,
                        optical_color.2 + temp_color.2 * lum
                    );
                }
            }
        }

        PhotonData {
            optical_color,
            xray_color,
            pixel: self.pixel,
        }
    }

    pub fn run<M: Metric>(mut self, max_iterations: usize, dtau: f64, metric: &M, source: &AccretionDisk, rng: &fastrand::Rng) -> PhotonData {
        let mut iteration = 0;
        loop {
            let angle_dist = std::f64::consts::PI / 2.0 - (std::f64::consts::PI / 2.0 - self.pos[2]).abs();
            let use_dtau = dtau * f64::min(
                self.pos[1] - metric.get_horizon() + EPSILON,
                POLE_PROTECTION * angle_dist + EPSILON);
            let (dp1, dv1) = metric.christoffel(self.pos).accel(self.vel);
            let old_above = self.pos[2] < std::f64::consts::PI / 2.0;
            // Use subtraction because we're back-propagating.
            self.pos = sub4(self.pos, mul4(dp1, use_dtau));
            self.vel = sub4(self.vel, mul4(dv1, use_dtau));

            let new_above = self.pos[2] < std::f64::consts::PI / 2.0;
            if old_above ^ new_above {
                // Crossed the disk
                let redshift = (-metric.get_metric(self.pos)[0]).sqrt();
                let (temp, lum, depth_inc) = source.disk_collision(self.pos, self.vel, redshift);
                self.temps[self.temp_index] = (temp, lum * self.depth);
                self.depth *= depth_inc;
                self.temp_index += 1;
                if (self.temp_index) >= TEMP_RECORD {
                    break;
                }
            }

            if iteration % CORONA_INTERACTION == 0 && self.compton_scatter.is_none() {
                let collision_prob = source.corona_prob(self.pos);
                if rng.f64() < collision_prob * (dp1[0].abs() * use_dtau * CORONA_INTERACTION as f64) {
                    // Compton interacted!
                    let (energy_factor, new_vel) = source.corona_collide(self.pos, self.vel, &metric.get_metric(self.pos), &rng);
                    // Get rid of data accumulated so far
                    self.temp_index = 0;
                    self.vel = new_vel;
                    self.compton_scatter = Some(energy_factor);
                    self.depth = 1.0;
                }
            }

            let offset_check = if self.pos[1] > 3.0 {
                OFFSET_CHECK
            } else {
                1
            };
            if iteration % offset_check == 0 {
                let metric = metric.get_metric(self.pos);
                let linear_offset = dot4(self.vel, matvecmul(&metric, self.vel));
                let quad_offset = dot4(self.vel, matvecmul(&matsquare(&metric), self.vel));
                self.vel = add4(self.vel, mul4(matvecmul(&metric, self.vel), -0.5 * linear_offset / quad_offset));
                
                let new_offset = dot4(self.vel, matvecmul(&metric, self.vel));
                self.vel[0] = (-new_offset / metric[0] + self.vel[0] * self.vel[0]).sqrt();

                // #[cfg(debug_assertions)]
                // {
                //     if self.pos[0].is_nan() {
                //         break;
                //     }
                //     if linear_offset > 0.1 {
                //         println!("WARNING: large offset of {} detected, {}", linear_offset, new_offset);
                //     }
                // }
            }

            // Check location
            match metric.get_state(self.pos) {
                State::Dead => {
                    break;
                },
                State::Escape => {
                    break
                },
                State::Running => (),
            };
            if iteration >= max_iterations {
                break;
            }

            iteration += 1;
        }

        // Convert to photon data
        self.get_data()
    }
}