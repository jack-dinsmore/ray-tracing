use nalgebra::Matrix4;
use crate::{Christoffels, State, Vec4, Christoffel};

// https://arxiv.org/pdf/0904.4184.pdf

const R_S : f64 = 0.2;

#[derive(Debug)]
pub struct Minkowski {}
#[derive(Debug)]
pub struct Minkowski3 {}
#[derive(Debug)]
pub struct Schwarzschild {}
#[derive(Debug)]
pub struct Schwarzschild3 {}

impl Christoffels for Minkowski {
    fn get_state(pos: Vec4) -> State {
        if pos[0].abs() > 1e6 {
            return State::Dead;
        }
        if pos[1].abs() > 20.0 {
            return State::Escape;
        }
        State::Running
    }

    fn christoffel(pos: Vec4) -> Christoffel {
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

    fn get_metric(pos: Vec4) -> Matrix4<f64> {
        let st = pos[2].sin().abs();
        let r2 = pos[1] * pos[1];
        Matrix4::new(
            -1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, r2, 0.0,
            0.0, 0.0, 0.0, r2 * st * st)
    }
}

impl Christoffels for Minkowski3 {
    fn get_state(pos: Vec4) -> State {
        if pos[0].abs() > 1e6 {
            return State::Dead;
        }
        if pos[1].abs() > 20.0 {
            return State::Escape;
        }
        State::Running
    }

    fn christoffel(pos: Vec4) -> Christoffel {
        Christoffel::new([
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, -pos[1], 0.0,
            0.0, 0.0, 0.0, -pos[1],

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0/pos[1],
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
        ])
    }

    fn get_metric(pos: Vec4) -> Matrix4<f64> {
        let r2 = pos[1] * pos[1];
        Matrix4::new(
            -1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, r2, 0.0,
            0.0, 0.0, 0.0, r2)
    }
}

impl Christoffels for Schwarzschild {
    fn get_state(pos: Vec4) -> State {
        if pos[0].abs() > 1e6 {
            return State::Dead;
        }
        if pos[1].abs() < R_S {
            return State::Dead;
        }
        if pos[1].abs() > 20.0 {
            return State::Escape;
        }
        State::Running
    }

    fn christoffel(pos: Vec4) -> Christoffel {
        let st = pos[2].sin().abs();
        let ct = pos[2].cos();
        let diff = pos[1] - R_S;
        Christoffel::new([
            0.0, R_S / (2.0 * pos[1] * diff), 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,

            R_S * diff / (2.0 * pos[1] * pos[1] * pos[1]), 0.0, 0.0, 0.0,
            0.0, -R_S / (2.0 * pos[1] * diff), 0.0, 0.0,
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

    fn get_metric(pos: Vec4) -> Matrix4<f64> {
        let st = pos[2].sin().abs();
        let r2 = pos[1] * pos[1];
        let f = 1.0 - R_S / pos[1];
        Matrix4::new(
            -f, 0.0, 0.0, 0.0,
            0.0, 1.0 / f, 0.0, 0.0,
            0.0, 0.0, r2, 0.0,
            0.0, 0.0, 0.0, r2 * st * st)
    }
}

impl Christoffels for Schwarzschild3 {
    fn get_state(pos: Vec4) -> State {
        if pos[0].abs() > 1e6 {
            return State::Dead;
        }
        if pos[1].abs() < R_S {
            return State::Dead;
        }
        if pos[1].abs() > 20.0 {
            return State::Escape;
        }
        State::Running
    }

    fn christoffel(pos: Vec4) -> Christoffel {
        let diff = pos[1] - R_S;
        Christoffel::new([
            0.0, R_S / (2.0 * pos[1] * diff), 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,

            R_S * diff / (2.0 * pos[1] * pos[1] * pos[1]), 0.0, 0.0, 0.0,
            0.0, -R_S / (2.0 * pos[1] * diff), 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -diff,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,

            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0/pos[1],
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
        ])
    }

    fn get_metric(pos: Vec4) -> Matrix4<f64> {
        let r2 = pos[1] * pos[1];
        let f = 1.0 - R_S / pos[1];
        Matrix4::new(
            -f, 0.0, 0.0, 0.0,
            0.0, 1.0 / f, 0.0, 0.0,
            0.0, 0.0, r2, 0.0,
            0.0, 0.0, 0.0, r2)
    }
}