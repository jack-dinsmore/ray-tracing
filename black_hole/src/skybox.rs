use crate::{Vec4, PI};
use crate::views::Color;

pub struct Skybox {
    color_fn: fn (Vec4, Vec4) -> Color,
}

impl Skybox {
    pub fn blank() -> Skybox {
        Skybox {color_fn: |pos: Vec4, _vel: Vec4| -> Color {
            let brightness = (1.0 - pos[2].cos().abs()).powf(4.0);
            let red = f64::sin(4.0 * pos[3]).powi(2);
            let green = f64::cos(4.0 * pos[3]).powi(2);
            ((brightness * red * 255.0) as u8, (brightness * green * 255.0) as u8, 0)
        }}
    }

    pub fn call(&self, pos: Vec4, vel: Vec4) -> Color{
        (self.color_fn)(pos, vel)
    }
}