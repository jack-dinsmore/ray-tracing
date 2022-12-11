#![allow(dead_code)]
#![allow(non_snake_case)]

mod metrics;
mod observer;
mod util;
mod engine;
mod source;

use observer::Simple;
use engine::Engine;

const THETA: f64 = 3.14 / 2.0 - 0.3f64;
const START_POS: [f64; 3] = [10.0, THETA, 0.0];
const MAX_ITER: usize = 1_000_000;

fn flat() {
    let metric = metrics::Minkowski::new(0.0);
    let source = source::AccretionDisk::flat();

    let observer = Simple::new(START_POS, [-THETA.sin(), 0.0, -THETA.cos()], metric);
    let mut engine = Engine::new(observer, "flat".to_owned());
    engine.run(MAX_ITER, 5e-3, metric, source);
}

fn mink() {
    let metric = metrics::Minkowski::new(1.0);
    let source = source::AccretionDisk::thick();

    let observer = Simple::new(START_POS, [-THETA.sin(), 0.0, -THETA.cos()], metric);
    let mut engine = Engine::new(observer, "minkowski".to_owned());
    engine.run(MAX_ITER, 5e-3, metric, source);
}

fn thick() {
    let metric = metrics::Schwarzschild::new();
    let source = source::AccretionDisk::thick();

    
    let observer = Simple::new(START_POS, [-THETA.sin(), 0.0, -THETA.cos()], metric);
    let mut engine = Engine::new(observer, "thick".to_owned());
    let photon_count = engine.run(MAX_ITER, 1e-3, metric, source);
    println!("{} photons run successfully", photon_count);
}

fn thin() {
    let metric = metrics::Schwarzschild::new();
    let source = source::AccretionDisk::thin();
    
    let observer = Simple::new(START_POS, [-THETA.sin(), 0.0, -THETA.cos()], metric);
    let mut engine = Engine::new(observer, "thin".to_owned());
    let photon_count = engine.run(MAX_ITER, 2e-3, metric, source);
    println!("{} photons run successfully", photon_count);
}

fn sch() {
    let metric = metrics::Schwarzschild::new();
    let source = source::AccretionDisk::corona();

    let observer = Simple::new(START_POS, [-THETA.sin(), 0.0, -THETA.cos()], metric);
    let mut engine = Engine::new(observer, "schwarzschild".to_owned());
    let photon_count = engine.run(MAX_ITER, 2e-3, metric, source);
    println!("{} photons run successfully", photon_count);
}

fn kerr() {
    let metric = metrics::Kerr::new (0.8);
    let source = source::AccretionDisk::corona();

    let observer = Simple::new(START_POS, [-THETA.sin(), 0.0, -THETA.cos()], metric);
    let mut engine = Engine::new(observer, "kerr".to_owned());
    let photon_count = engine.run(MAX_ITER, 2e-3, metric, source);
    println!("{} photons run successfully", photon_count);
}

fn main() {
    // flat();
    // mink();
    // thick();
    // thin();
    // sch();
    kerr();
}