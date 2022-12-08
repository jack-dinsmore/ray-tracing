pub use crate::metrics::Metric;

pub type Vec4 = [f64; 4];
pub type Vec3 = [f64; 3];
pub type Matrix4 = [f64; 16];

fn get_elem(mat: &Matrix4, elem: (usize, usize)) -> f64{
    mat[elem.0 + 4 * elem.1]
}

pub fn add4(a: Vec4, b: Vec4) -> Vec4 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]
}

pub fn sub4(a: Vec4, b: Vec4) -> Vec4 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]]
}

pub fn mul4(a: Vec4, f: f64) -> Vec4 {
    [a[0] * f, a[1] * f, a[2] * f, a[3] * f]
}

pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1],
     a[2] * b[0] - a[0] * b[2],
     a[0] * b[1] - a[1] * b[0]]
}

pub fn add3 (a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub fn mul3 (a: Vec3, f: f64) -> Vec3 {
    [a[0] * f, a[1] * f, a[2] * f]
}

pub fn opposite(v: Vec3) -> Vec3 {
    [-v[0], -v[1], -v[2]]
}

pub fn length(v: Vec3) -> f64 {
    (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt()
}

pub fn dot4(a: Vec4, b: Vec4) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
}

pub fn dot3(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn normalize(v: Vec3) -> Vec3 {
    mul3(v, 1.0 / length(v))
}

pub fn matsquare(m: &Matrix4) -> Matrix4 {
    let mut new_mat = [0.0; 16];
    for i in 0..4 {
        for j in 0..4 {
            new_mat[i * 4 + j] = get_elem(m, (i,0)) * get_elem(m, (0,j)) + get_elem(m, (i,1)) * get_elem(m, (1,j)) + get_elem(m, (i,2)) * get_elem(m, (2,j)) + get_elem(m, (i,3)) * get_elem(m, (3,j));
        }
    }
    new_mat
}

pub fn matvecmul(m: &Matrix4, v: Vec4) -> Vec4 {
    [
        get_elem(m, (0, 0)) * v[0] + get_elem(m, (0, 1)) * v[1] + get_elem(m, (0, 2)) * v[2] + get_elem(m, (0, 3)) * v[3],
        get_elem(m, (1, 0)) * v[0] + get_elem(m, (1, 1)) * v[1] + get_elem(m, (1, 2)) * v[2] + get_elem(m, (1, 3)) * v[3],
        get_elem(m, (2, 0)) * v[0] + get_elem(m, (2, 1)) * v[1] + get_elem(m, (2, 2)) * v[2] + get_elem(m, (2, 3)) * v[3],
        get_elem(m, (3, 0)) * v[0] + get_elem(m, (3, 1)) * v[1] + get_elem(m, (3, 2)) * v[2] + get_elem(m, (3, 3)) * v[3]
    ]
}

pub fn cart_to_spher_vel(pos: Vec3, cart_vel: Vec3) -> Vec3 {
    let st = pos[1].sin();
    let ct = pos[1].cos();
    let sp = pos[2].sin();
    let cp = pos[2].cos();
    [
        st * (cart_vel[0] * cp + cart_vel[1] * sp) + cart_vel[2] * ct,
        (ct * (cart_vel[0] * cp + cart_vel[1] * sp) - cart_vel[2] * st) / pos[0],
        (cart_vel[1] * cp - cart_vel[0] * sp) / (pos[0] * st)
    ]
}

pub fn spher_to_cart_vel(pos: Vec3, spher_vel: Vec3) -> Vec3 {
    let st = pos[1].sin();
    let ct = pos[1].cos();
    let sp = pos[2].sin();
    let cp = pos[2].cos();
    [
        spher_vel[0] * st * cp + pos[0] * spher_vel[1] * ct * cp - pos[0] * spher_vel[2] * st * sp,
        spher_vel[0] * st * sp + pos[0] * spher_vel[1] * ct * sp + pos[0] * spher_vel[2] * st * cp,
        spher_vel[0] * ct - pos[0] * spher_vel[1] * st,
    ]
}

/*pub fn get_initial_state(pos: Vec3, cart_vel: Vec3) -> (Vec4, Vec4) {
    let vel = cart_to_spher_vel(pos, cart_vel);

    let pos4 = [0.0, pos[0], pos[1], pos[2]];
    let metric = CHRISTOFFELS.get_metric(pos4);

    let ldotv = get_elem(&metric, (0, 1)) * vel[0] + get_elem(&metric, (0, 2)) * vel[1] + get_elem(&metric, (0, 3)) * vel[2];
    let mut vgv = 0.0;
    for i in 1..4 {
        for j in 1..4 {
            vgv += vel[i-1] * vel[j-1] * get_elem(&metric, (i, j));
        }
    }
    let time = (-ldotv - (ldotv * ldotv - get_elem(&metric, (0,0))*vgv).sqrt()) / get_elem(&metric, (0,0));
    (pos4, [time, vel[0], vel[1], vel[2]])
}*/

/*pub fn rotate_pos_vel(pos: Vec4, vel: Vec4, normal: Vec3, psi: f64) -> (Vec3, Vec3) {
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
}*/

pub fn get_vel_from_metric(vel: Vec3, g: &Matrix4) -> Vec4 {
    let mut v4 = [0.0, vel[0], vel[1], vel[2]];
    let spatial_norm = dot4(v4, matvecmul(&g, v4));
    v4[0] = (-spatial_norm / g[0]).sqrt();
    v4
}