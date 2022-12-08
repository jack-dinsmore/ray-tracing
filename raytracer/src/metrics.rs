use crate::util::{Vec4, Matrix4};

// https://arxiv.org/pdf/0904.4184.pdf

const SPACETIME_EDGE: f64 = 100.0;

pub enum State {
    Dead,
    Escape,
    Running
}

pub struct Christoffel {
    pub numbers: [f64; 64],
}

pub trait Metric: Copy + Clone + Send + 'static {
    fn get_state(&self, pos: Vec4) -> State;
    fn christoffel(&self, pos: Vec4) -> Christoffel;
    fn get_metric(&self, pos: Vec4) -> Matrix4;
    fn get_horizon(&self) -> f64;
}

#[derive(Debug, Copy, Clone)]
pub struct Minkowski { radius: f64 }
impl Minkowski {
    pub fn new(radius: f64) -> Self { Self { radius } }
}
#[derive(Debug, Copy, Clone)]
pub struct Schwarzschild {
}
impl Schwarzschild {
    pub fn new() -> Self { Self {} }
}
#[derive(Debug, Copy, Clone)]
pub struct MorrisThorne {
    b0: f64
}
impl MorrisThorne {
    pub fn new(b0: f64) -> Self { Self { b0 } }
}
#[derive(Debug, Copy, Clone)]
pub struct Kerr {
    a: f64,
    horizon: f64
}
impl Kerr {
    pub fn new(a: f64) -> Self {
        Self {
            a,
            horizon: 1.0 / 2.0 + (1.0 - a * a).sqrt() / 2.0
        }
    }
}

impl Christoffel {
    fn new(numbers: [f64; 64]) -> Christoffel {
        Christoffel {numbers}
    }

    fn get(&self, mu: usize, alpha: usize, beta: usize) -> f64 {
        match alpha > beta {
            true => self.numbers[mu * 16 + beta * 4 + alpha],
            false => self.numbers[mu * 16 + alpha * 4 + beta]
        }
    }

    pub fn accel(&self, vel: Vec4) -> (Vec4, Vec4) {
        let mut sum = [0.0; 4];
        for mu in 0..4 {
            for alpha in 0..4 {
                for beta in 0..4 {
                    sum[mu] -= self.get(mu, alpha, beta) * vel[alpha] * vel[beta];
                }
            }
        }
        (vel, sum)
    }
}

impl Metric for Minkowski {
    fn get_horizon(&self) -> f64 { 0.0 }
    fn get_state(&self, pos: Vec4) -> State {
        if pos[0].abs() > 1e6 || pos[1].abs() < self.radius {
            return State::Dead;
        }
        if pos[1].abs() > SPACETIME_EDGE {
            return State::Escape;
        }
        State::Running
    }

    fn christoffel(&self, pos: Vec4) -> Christoffel {
        let st = pos[2].sin().abs();
        let ct = pos[2].cos();
        Christoffel::new([
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, -pos[1], 0.0,
            0.0, 0.0, 0.0, -pos[1] * st * st,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0/pos[1], 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -st * ct,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0/pos[1],
            0.0, 0.0, 0.0, ct / st,
            0.0, 0.0, 0.0, 0.0,
        ])
    }

    fn get_metric(&self, pos: Vec4) -> Matrix4 {
        let st = pos[2].sin().abs();
        let r2 = pos[1] * pos[1];
        [
            -1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, r2, 0.0,
            0.0, 0.0, 0.0, r2 * st * st
        ]
    }
}

impl Metric for Schwarzschild {
    fn get_horizon(&self) -> f64 { 1.0 }
    fn get_state(&self, pos: Vec4) -> State {
        if pos[0].abs() > 1e6 {
            return State::Dead;
        }
        if pos[1].abs() < 1.0 {
            return State::Dead;
        }
        if pos[1].abs() > SPACETIME_EDGE {
            return State::Escape;
        }
        State::Running
    }

    fn christoffel(&self, pos: Vec4) -> Christoffel {
        let st = pos[2].sin().abs();
        let ct = pos[2].cos();
        let diff = pos[1] - 1.0;
        Christoffel::new([
            0.0, 1.0 / (2.0 * pos[1] * diff), 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,

            1.0 * diff / (2.0 * pos[1] * pos[1] * pos[1]), 0.0, 0.0, 0.0,
            0.0, -1.0 / (2.0 * pos[1] * diff), 0.0, 0.0,
            0.0, 0.0, -diff, 0.0,
            0.0, 0.0, 0.0, -diff * st * st,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0/pos[1], 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -st * ct,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0/pos[1],
            0.0, 0.0, 0.0, ct / st,
            0.0, 0.0, 0.0, 0.0,
        ])
    }

    fn get_metric(&self, pos: Vec4) -> Matrix4 {
        let st = pos[2].sin().abs();
        let r2 = pos[1] * pos[1];
        let f = 1.0 - 1.0 / pos[1];
        [
            -f, 0.0, 0.0, 0.0,
            0.0, 1.0 / f, 0.0, 0.0,
            0.0, 0.0, r2, 0.0,
            0.0, 0.0, 0.0, r2 * st * st
        ]
    }
}

impl Metric for MorrisThorne {
    fn get_horizon(&self) -> f64 { self.b0 }

    fn get_state(&self, pos: Vec4) -> State {
        if pos[0].abs() > 1e6 {
            return State::Dead;
        }
        if pos[1].abs() > SPACETIME_EDGE {
            return State::Escape;
        }
        State::Running
    }

    fn christoffel(&self, pos: Vec4) -> Christoffel {
        let st = pos[2].sin().abs();
        let ct = pos[2].cos();
        let quad = pos[1] / (self.b0 * self.b0 + pos[1] * pos[1]);
        Christoffel::new([
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, -pos[1], 0.0,
            0.0, 0.0, 0.0, -pos[1] * st * st,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, quad, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -st * ct,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, quad,
            0.0, 0.0, 0.0, ct / st,
            0.0, 0.0, 0.0, 0.0,
        ])
    }

    fn get_metric(&self, pos: Vec4) -> Matrix4 {
        let st = pos[2].sin().abs();
        let fake_r2 = self.b0 * self.b0 + pos[1] * pos[1];
        [
            -1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, fake_r2, 0.0,
            0.0, 0.0, 0.0, fake_r2 * st * st
        ]
    }
}

impl Metric for Kerr {
    fn get_horizon(&self) -> f64 { 1.0 }

    fn get_state(&self, pos: Vec4) -> State {
        if pos[0].abs() > 1e6 {
            return State::Dead;
        }
        if pos[1].abs() < self.horizon {
            return State::Dead;
        }
        if pos[1].abs() > SPACETIME_EDGE {
            return State::Escape;
        }
        State::Running
    }

    fn christoffel(&self, pos: Vec4) -> Christoffel {
        let st = pos[2].sin().abs();
        let ct = pos[2].cos();
        let r2 = pos[1] * pos[1];
        let a2 = self.a * self.a;
        let sigma = r2 + a2 * ct * ct;
        let sigma2 = sigma * sigma;
        let delta = r2 - 1.0 * pos[1] + a2;
        let adiff = r2 - a2 * ct * ct;
        let A = (r2 + a2) * (r2 + a2) - a2 * delta * st * st;
        Christoffel::new([
            0.0,
            1.0 * (r2 + a2) * adiff / (2.0 * sigma2 * delta),
            -1.0 * a2 * pos[1] * ct * st / sigma2,
            0.0,

            0.0,
            0.0,
            0.0,
            1.0 * self.a * st * st * (a2 * ct * ct * (a2 - r2) - r2 * (a2 + 3.0 * r2)) / (2.0 * sigma2 * delta),

            0.0,
            0.0,
            1.0 * self.a * a2 * pos[1] * st * st * st * ct / sigma2,
            0.0,

            0.0,
            0.0,
            0.0,
            0.0,

            // r
            1.0 * delta * adiff / (2.0 * sigma2 * sigma),
            0.0,
            0.0,
            -delta * 1.0 * self.a * st * st * adiff / (2.0 * sigma * sigma2),

            0.0,
            (2.0 * pos[1] * a2 * st * st - 1.0 * adiff) / (2.0 * sigma * delta),
            -a2 * st * ct / sigma,
            0.0,

            0.0,
            0.0,
            -pos[1] * delta / sigma,
            0.0,

            0.0,
            0.0,
            0.0,
            delta * st * st / (2.0 * sigma2 * sigma) * (-2.0 * pos[1] * sigma2 + 1.0 * a2 * st * st * adiff),

            // theta
            -1.0 * a2 * pos[1] * ct * st / (sigma * sigma2),
            0.0,
            0.0,
            1.0 * self.a * pos[1] * (r2 + a2) * st * ct / (sigma * sigma2),

            0.0,
            a2 * ct * st / (sigma * delta),
            pos[1] / sigma,
            0.0,

            0.0,
            0.0,
            - a2 * st * ct / sigma,
            0.0,

            0.0,
            0.0,
            0.0,
            -st * ct / (sigma * sigma2) * (A * sigma + (r2 + a2) * 1.0 * a2 * pos[1] * st * st),

            // phi
            0.0,
            1.0 * self.a * adiff / (2.0 * sigma2 * delta),
            -1.0 * self.a * pos[1] * ct / sigma2,
            0.0,

            0.0,
            0.0,
            0.0,
            (2.0 * pos[1] * sigma2 + 1.0 * (a2 * a2 * st * st * ct * ct - r2 * (sigma + r2 + a2))) / (2.0 * sigma2 * delta),

            0.0,
            0.0,
            0.0,
            ct / st / sigma2 * (sigma2 + 1.0 * a2 * pos[1] * st * st),

            0.0,
            0.0,
            0.0,
            0.0,
        ])
    }

    fn get_metric(&self, pos: Vec4) -> Matrix4 {
        let st = pos[2].sin().abs();
        let ct = pos[2].cos().abs();
        let r2 = pos[1] * pos[1];
        let a2 = self.a * self.a;
        let sigma = r2 + a2 * ct * ct;
        let delta = r2 - 1.0 * pos[1] + a2;
        [
            -(1.0 - 1.0 * pos[1] / sigma), 0.0, 0.0, 1.0 * self.a * pos[1] * st * st / sigma,
            0.0, sigma / delta, 0.0, 0.0,
            0.0, 0.0, sigma, 0.0,
            1.0 * self.a * pos[1] * st * st / sigma, 0.0, 0.0, (r2 + a2 + 1.0*pos[1] * a2*st*st / sigma) * st * st
        ]
    }
}