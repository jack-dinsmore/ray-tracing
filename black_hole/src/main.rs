use std::marker::PhantomData;
use nalgebra::{Matrix4, Vector4};
mod metrics;
mod skybox;
mod views;
use views::View;


const PI : f64 = 3.141592653589793238462643383;

type Vec4 = [f64; 4];
type Vec3 = [f64; 3];


fn add (a: Vec4, b: Vec4) -> Vec4 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]
}

fn mul (a: Vec4, f: f64) -> Vec4 {
    [a[0] * f, a[1] * f, a[2] * f, a[3] * f]
}

#[derive(Debug)]
pub enum State {
    Running,
    Dead,
    Escape,
}

#[derive(Debug)]
pub struct Photon<C> {
    pos: Vec4,
    vel: Vec4,
    state: State,
    phantom: PhantomData<C>,
}

pub struct Christoffel {
    pub numbers: [f64; 64],
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

    fn accel(&self, vel: Vec4) -> (Vec4, Vec4) {
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

pub trait Christoffels {
    fn get_state(pos: Vec4) -> State;
    fn christoffel(pos: Vec4) -> Christoffel;
    fn get_metric(pos: Vec4) -> Matrix4<f64>;
}

impl<C> Photon<C> where C: Christoffels  {
    pub fn new<T>(pos: Vec4, vel: Vec4) -> Photon<T> {
        Photon {pos, vel, state: State::Running, phantom: PhantomData::<T> }
    }

    pub fn update(&mut self, iteration: u32, dtau: f64) -> &State {
        // Don't execute non-running photons
        match &self.state {
            State::Running => (),
            State::Escape => return &self.state,
            State::Dead => return &self.state,
        }

        let (dp1, dv1) = C::christoffel(self.pos).accel(self.vel);
        self.pos = add(self.pos, mul(dp1, dtau));
        self.vel = add(self.vel, mul(dv1, dtau));

        /*// Runge Kutta 4
        let (dp1, dv1) = C::christoffel(self.pos).accel(self.vel);
        let (dp2, dv2) = C::christoffel(add(self.pos, mul(dp1, dtau / 2.0)))
            .accel(add(self.vel, mul(dv1, dtau / 2.0)));
        let (dp3, dv3) = C::christoffel(add(self.pos, mul(dp2, dtau / 2.0)))
            .accel(add(self.vel, mul(dv2, dtau / 2.0)));
        let (dp4, dv4) = C::christoffel(add(self.pos, mul(dp3, dtau)))
            .accel(add(self.vel, mul(dv3, dtau)));

        self.pos = add(self.pos, add(mul(add(dp1, dp4), dtau/6.0), mul(add(dp2, dp3), dtau/3.0)));
        self.vel = add(self.vel, add(mul(add(dv1, dv4), dtau/6.0), mul(add(dv2, dv3), dtau/3.0)));*/

        // Check velocity
        if iteration % 0x100 == 0 {
            let metric = C::get_metric(self.pos);
            let mut vel = Vector4::new(self.vel[0], self.vel[1], self.vel[2], self.vel[3]);
            let offset = -0.5 * Vector4::dot(&vel, &(metric * vel))
                / Vector4::dot(&vel, &(metric * metric * vel));
            vel += offset * metric * vel;
            self.vel = [vel[0], vel[1], vel[2], vel[3]];
        }
        self.state = C::get_state(self.pos);
        &self.state
    }
}

fn main() {
    let pos = [1.0, PI / 2.0, 0.0];// r, theta, phi

    let mut plane: views::ViewLine<metrics::Minkowski3> = views::ViewLine::<metrics::Minkowski3>::new(
        pos, [-2.0, 0.0, 0.0], [0.0, 0.0, 1.0], 640, 480);
    
    let dtau = 1e-3;

    plane.run(1_000_000, dtau);

    let skybox = skybox::Skybox::blank();
    let image = plane.render(skybox);
    if let Err(_) = image.save("out.bmp") {
        println!("Image saving threw an error.");
    }
    println!("Success");
}