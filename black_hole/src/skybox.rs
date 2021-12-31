use crate::{Vec4, Photon};
use crate::views::Color;

pub struct Skybox {
    color_fn: fn (Vec4, Vec4) -> Color,
}

impl Skybox {
    pub fn blank() -> Skybox {
        Skybox {color_fn: |pos: Vec4, _vel: Vec4| -> Color {
            let brightness = pos[2].sin();
            let red = pos[3].sin().powi(2);
            let green = pos[3].cos().powi(2);
            ((brightness * red * 255.0) as u8, (brightness * green * 255.0) as u8, 0)
        }}
    }

    pub fn call<C>(&self, photon: Photon<C>) -> Color{
        (self.color_fn)(photon.pos, photon.vel)
    }
}