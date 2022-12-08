use crate::util::*;

// Fix R_S at 1.
const ALPHA: f64 = 500.0;
const MASS: f64 = 2.4e12;
const REDSHIFT: f64 = 3.5;
const CORONA_DENSITY_SCALE: f64 = 5.0e-26;
const CORONA_PROTON_GAMMA_MINUS_ONE: f64 = 0.1;

const CORONA_DENSITY_HEIGHT: f64 = 2.25 / CORONA_PROTON_GAMMA_MINUS_ONE;
const CORONA_ELECTRON_GAMMA: f64 = 1836.0 * CORONA_PROTON_GAMMA_MINUS_ONE;
pub const RESCALE_FOR_XRAY: f64 = CORONA_ELECTRON_GAMMA * CORONA_ELECTRON_GAMMA;

#[derive(Clone)]
pub struct AccretionDisk {
    temp_scale: f64,
    ang_vel_at_horizon: f64,
    tau_scale: f64,
    lum_scale: f64,
    corona_scale: f64,
}

impl AccretionDisk {
    pub fn corona() -> Self {
        Self {
            temp_scale: 3200.0 * MASS.powf(-0.25), // eV
            ang_vel_at_horizon: 0.5f64.sqrt(),
            tau_scale: 3.7e-3 * ALPHA * MASS.powf(-1.0/8.0) * (2.0 * std::f64::consts::PI).sqrt(),
            lum_scale: 1.14e26 / MASS,
            corona_scale: 1.0,
        }
    }
    pub fn flat() -> Self {
        Self {
            temp_scale: 3200.0 * MASS.powf(-0.25), // eV
            ang_vel_at_horizon: 0.0,
            tau_scale: 3.7e-3 * ALPHA * MASS.powf(-1.0/8.0) * (2.0 * std::f64::consts::PI).sqrt(),
            lum_scale: 1.14e26 / MASS,
            corona_scale: 0.0,
        }
    }
    pub fn thick() -> Self {
        Self {
            temp_scale: 3200.0 * MASS.powf(-0.25), // eV
            ang_vel_at_horizon: 0.5f64.sqrt(),
            tau_scale: 3.7e-3 * ALPHA * MASS.powf(-1.0/8.0) * (2.0 * std::f64::consts::PI).sqrt() * 100.0,
            lum_scale: 1.14e26 / MASS,
            corona_scale: 0.0,
        }
    }
    pub fn thin() -> Self {
        Self {
            temp_scale: 3200.0 * MASS.powf(-0.25), // eV
            ang_vel_at_horizon: 0.5f64.sqrt(),
            tau_scale: 3.7e-3 * ALPHA * MASS.powf(-1.0/8.0) * (2.0 * std::f64::consts::PI).sqrt(),
            lum_scale: 1.14e26 / MASS,
            corona_scale: 0.0,
        }
    }

    // Get the temperature, luminosity, and depth value after this point.
    pub fn disk_collision(&self, pos: Vec4, vel: Vec4, grav_redshift: f64) -> (f64, f64, f64) {
        let disk_vel = self.ang_vel_at_horizon / pos[1].sqrt();
        let disk_vel = [-pos[3].sin() * disk_vel, pos[3].cos() * disk_vel, 0.0];
        let light_3vel = spher_to_cart_vel(
            [pos[1], pos[2], pos[3]],
            [vel[1] / vel[0], vel[2] / vel[0], vel[3] / vel[0]]
        );
        let vel_redshift = (1.0 + dot3(disk_vel, light_3vel) / (1.0 - dot3(disk_vel, light_3vel))).sqrt();
        let depth = self.tau_scale * (pos[1]).powf(1.25) / light_3vel[2].abs();
            
        let temp = self.temp_scale / pos[1].sqrt() * grav_redshift * vel_redshift / REDSHIFT;
        let lum = self.lum_scale / (pos[1] * pos[1]);
        (temp, lum, (-depth).exp())
    }

    // Get the probability of collision per distance unit
    pub fn corona_prob(&self, pos: Vec4) -> f64 {
        let scale_factor = CORONA_DENSITY_HEIGHT / pos[1];
        let corona_density = self.corona_scale * CORONA_DENSITY_SCALE * (1.0 + scale_factor + scale_factor * scale_factor * 0.5);
        2.16e8 * MASS * corona_density // But also energy dependence`
    }

    // Get the energy shift and the new velocity.
    pub fn corona_collide(&self, pos: Vec4, _vel: Vec4, metric: &Matrix4, rng: &fastrand::Rng) -> (f64, Vec4) {
        // let theta = cdf_theta.asin();
        // let phi = cdf_phi * std::f64::consts::PI * 2.0;
        // let new_vel_cart = [
        //     theta.sin() * phi.cos(),
        //     theta.sin() * phi.sin(),
        //     theta.cos()
        // ];

        let mut new_vel_cart = [0.0, 0.0, 0.0];
        const NUM_RANDS: usize = 5;
        for i in 0..3 {
            for _ in 0..NUM_RANDS {
                new_vel_cart[i] += rng.f64() - 0.5;
            }
        }

        let energy = RESCALE_FOR_XRAY;
        let vel3 = cart_to_spher_vel([pos[1],pos[2],pos[3]], normalize(new_vel_cart));
        let vel4 = get_vel_from_metric(vel3, metric);
        (energy, vel4)
    }
}

/// Find the root of a function using Newton's method
pub fn find_root(f: impl Fn(f64)->f64, df: impl Fn(f64)->f64, start: f64) -> f64 {
    let mut x = start;
    let mut error = f(x);
    loop {
        x -= error / df(x);
        error = f(x);
        if error.abs() < 0.001 {
            return x;
        }
    }
}