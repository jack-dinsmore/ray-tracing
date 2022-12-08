use crate::{Vec4};
use crate::views::Color;
use perlin_noise::PerlinNoise;

const MAX_L: i64 = 10;
const GAS_COLOR: Color = (255, 89, 203);

pub struct Skybox {
    color_fn: Box<dyn Fn(Vec4, Vec4) -> Color>,
}

impl Skybox {
    pub fn demo() -> Skybox {
        Skybox {color_fn: Box::new(|pos: Vec4, _vel: Vec4| -> Color {
            let brightness = (1.0 - pos[2].cos().abs()).powf(4.0);
            let red = f64::sin(4.0 * pos[3]).powi(2);
            let green = f64::cos(4.0 * pos[3]).powi(2);
            ((brightness * red * 255.0) as u8, (brightness * green * 255.0) as u8, 0)
        })}
    }
    pub fn black() -> Skybox {
        Skybox {color_fn: Box::new(|_pos: Vec4, _vel: Vec4| -> Color {
            (0, 0, 0)
        })}
    }

    pub fn gas() -> Skybox {
        let perlin = PerlinNoise::new();
        Skybox {color_fn: Box::new(move |pos: Vec4, _vel: Vec4| -> Color {
            let amp = perlin.get3d([f64::abs(pos[2].sin() * pos[3].cos()), f64::abs(pos[2].sin() * pos[3].sin()), f64::abs(pos[2].cos())]);
            let amp = f64::max(0.0, f64::min(1.0, (amp-0.5) * 3.0 + 0.5));
            ((amp * (GAS_COLOR.0 as f64)) as u8, (amp * (GAS_COLOR.1 as f64)) as u8, (amp * (GAS_COLOR.2 as f64)) as u8)
        })}
    }

    pub fn call(&self, pos: Vec4, vel: Vec4) -> Color {
        (self.color_fn)(pos, vel)
    }
}