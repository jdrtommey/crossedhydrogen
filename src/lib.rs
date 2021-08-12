/// This small module computes the energy of hydrogen (& positronium) energy levels in arbitrary
/// crossed electric and magnetic fields.
/// The energies are given in atomic units and are for the |n,n1,n2,m> basis.
/// The function here takes electic and magnetic fields in SI units and returns an energy in atomic
/// units of hartrees.
///
/// The function assumes that physically allowed values of the quantum numbers of the states are
/// given.
use std::fmt;

// define new type to impliment neccessary operations on the [f64;3] input type for
// the fields.
type Vecc = [f64; 3];

pub trait Operations<Rhs = Self> {
    type Output;

    fn add(self, rhs: Rhs) -> Self::Output;
    fn subtract(self, rhs: Rhs) -> Self::Output;
    fn multiply(self, rhs: Rhs) -> Self::Output;
    fn cross(self, rhs: Rhs) -> Self::Output;
    fn dot(self, rhs: Rhs) -> f64;
    fn modulus(self) -> f64;
    fn dot_angle(self, rhs: Rhs) -> f64;
    fn float_multiply(self, value: f64) -> Self::Output;
}
impl Operations for Vecc {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        [self[0] + other[0], self[1] + other[1], self[2] + other[2]]
    }
    fn subtract(self, other: Self) -> Self {
        [self[0] - other[0], self[1] - other[1], self[2] - other[2]]
    }
    fn multiply(self, other: Self) -> Self {
        [self[0] * other[0], self[1] * other[1], self[2] * other[2]]
    }
    fn cross(self, other: Self) -> Self {
        [
            self[1] * other[2] - self[2] * other[1],
            self[2] * other[0] - self[0] * other[2],
            self[0] * other[1] - self[1] * other[0],
        ]
    }
    fn dot(self, other: Self) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
    }
    fn modulus(self) -> f64 {
        f64::sqrt(self[0] * self[0] + self[1] * self[1] + self[2] * self[2])
    }
    fn dot_angle(self, other: Self) -> f64 {
        (self.dot(other) / (self.modulus() * other.modulus())).acos()
    }
    fn float_multiply(self, value: f64) -> Self {
        [self[0] * value, self[1] * value, self[2] * value]
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct HydrogenError {}

impl fmt::Display for HydrogenError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Invalid input parameters")
    }
}

fn check_inputs(n: i64, n1: i64, n2: i64, m: i64) -> bool {
    if n == n1 + n2 + i64::abs(m) + 1 {
        true
    } else {
        false
    }
}

// Constants to convert the SI input vectors into atomic unit vectors.
const SI_TO_ATOMIC_ELEC: f64 = 1.0 / 5.14220674763e-11;
const SI_TO_ATOMIC_MAG: f64 = 1.0 / 2.35051756758e5;

//second order function which returns the energy of hydrogen in crossed electric and magnetic
//fields up to second order.
pub fn hydrogen_energy(
    electric_field: [f64; 3],
    magnetic_field: [f64; 3],
    n: i64,
    n1: i64,
    n2: i64,
    m: i64,
) -> Result<f64, HydrogenError> {
    match check_inputs(n, n1, n2, m) {
        true => {}
        false => return Err(HydrogenError {}),
    }

    if electric_field == [0.0, 0.0, 0.0] && magnetic_field == [0.0, 0.0, 0.0] {
        return Ok(-1.0 / (2.0 * (n * n) as f64));
    }
    let atomic_electric_field = electric_field.float_multiply(SI_TO_ATOMIC_ELEC);
    let atomic_magnetic_field = magnetic_field.float_multiply(SI_TO_ATOMIC_MAG);
    let n_dash = (m + n2 - n1) as f64 / 2.0;
    let n_double_dash = (m - n2 + n1) as f64 / 2.0;
    let (omega1, omega2) =
        get_normalised_vectors(atomic_electric_field, atomic_magnetic_field, n as f64);
    let (alpha1, alpha2) = get_angles(magnetic_field, omega1, omega2);

    let n_float = n as f64;
    let n_dash_float = n_dash as f64;
    let n_double_dash_float = n_double_dash as f64;

    let e1 = _energy_term_1(n as f64, omega1, omega2, n_dash_float, n_double_dash_float);
    let e2 = _energy_term_2(
        n_float,
        n_dash_float,
        n_double_dash_float,
        atomic_electric_field,
        alpha1,
        alpha2,
    );
    let e3 = _energy_term_3(
        n_float,
        n_dash_float,
        n_double_dash_float,
        alpha1,
        alpha2,
        atomic_magnetic_field,
    );
    let e4 = _energy_term_4(n_float, n_dash_float, n_double_dash_float, alpha1, alpha2);

    Ok(e1 + e2 + e3 + e4)
}
fn _energy_term_1(
    n: f64,
    omega1: [f64; 3],
    omega2: [f64; 3],
    n_dash: f64,
    n_double_dash: f64,
) -> f64 {
    -1.0 / (2.0 * n.powi(2) as f64) + omega1.modulus() * n_dash + omega2.modulus() * n_double_dash
}
//second term in the energy
fn _energy_term_2(
    n: f64,
    n_dash: f64,
    n_double_dash: f64,
    electric_field: [f64; 3],
    alpha1: f64,
    alpha2: f64,
) -> f64 {
    (-1.0 / 16.0)
        * n.powi(4)
        * electric_field.dot(electric_field)
        * (17.0 * n.powi(2) + 19.0
            - 12.0
                * (n_dash.powi(2)
                    + n_dash * n_double_dash * f64::cos(alpha1 + alpha2)
                    + n_double_dash.powi(2)))
}
//third term in the energy
fn _energy_term_3(
    n: f64,
    n_dash: f64,
    n_double_dash: f64,
    alpha1: f64,
    alpha2: f64,
    magnetic_field: [f64; 3],
) -> f64 {
    let coef = (1.0 / 48.0) * n.powi(2) * magnetic_field.dot(magnetic_field);
    let term =
        7.0 * n.powi(2) + 5.0 + 4.0 * n_dash * n_double_dash * f64::sin(alpha1) * f64::sin(alpha2);
    coef * term
}
fn _energy_term_4(n: f64, n_dash: f64, n_double_dash: f64, alpha1: f64, alpha2: f64) -> f64 {
    let t1 = (n.powi(2) - 1.0) * ((f64::cos(alpha1)).powi(2) + (f64::cos(alpha2)).powi(2));
    let t2 = -12.0
        * (n_dash.powi(2) * (f64::cos(alpha1)).powi(2)
            - n_dash * n_double_dash * f64::cos(alpha1) * (f64::cos(alpha2))
            + n_double_dash.powi(2) * (f64::cos(alpha2)).powi(2));
    t1 + t2
}
// The input electric and magnetic field vectors need to be converted in to the correct form and in
// atomic units.
fn get_normalised_vectors(
    electric_field: [f64; 3],
    magnetic_field: [f64; 3],
    n: f64,
) -> (Vecc, Vecc) {
    let omega1 = magnetic_field
        .subtract(electric_field.float_multiply(3.0 * n))
        .float_multiply(0.5);
    let omega2 = magnetic_field
        .add(electric_field.float_multiply(3.0 * n))
        .float_multiply(0.5);
    return (omega1, omega2);
}
fn get_angles(magnetic_field: [f64; 3], omega1: [f64; 3], omega2: [f64; 3]) -> (f64, f64) {
    (
        magnetic_field.dot_angle(omega1),
        magnetic_field.dot_angle(omega2),
    )
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_add() {
        let veca = [1.0, 1.0, 1.0];
        let vecb = [2.0, 2.0, 2.0];
        assert_eq!(veca.add(vecb), [3.0, 3.0, 3.0]);
        let veca = [1.0, 1.0, 1.0];
        let vecb = [2.0, 5.0, -2.0];
        assert_eq!(veca.add(vecb), [3.0, 6.0, -1.0]);
    }
    #[test]
    fn test_subtract() {
        let veca = [1.0, 1.0, 1.0];
        let vecb = [2.0, 2.0, 2.0];
        assert_eq!(veca.subtract(vecb), [-1.0, -1.0, -1.0]);
        let veca = [1.0, 1.0, 1.0];
        let vecb = [2.0, 5.0, -2.0];
        assert_eq!(veca.subtract(vecb), [-1.0, -4.0, 3.0]);
    }
    #[test]
    fn test_multiply() {
        let veca = [1.0, 1.0, 1.0];
        let vecb = [2.0, 2.0, 2.0];
        assert_eq!(veca.multiply(vecb), [2.0, 2.0, 2.0]);
        let veca = [1.0, 3.0, 1.0];
        let vecb = [2.0, 5.0, -2.0];
        assert_eq!(veca.multiply(vecb), [2.0, 15.0, -2.0]);
    }
    #[test]
    fn test_cross() {
        let veca = [1.0, 2.0, 2.0];
        let vecb = [2.0, 1.0, 2.0];
        assert_eq!(veca.cross(vecb), [2.0, 2.0, -3.0]);
        let vecb = [2.0, 1.0, 2.0];
        let veca = [1.0, 2.0, 2.0];
        assert_eq!(vecb.cross(veca), [-2.0, -2.0, 3.0]);
        let veca = [1.0, 0.0, 0.0];
        let vecb = [0.0, 1.0, 0.0];
        assert_eq!(veca.cross(vecb), [0.0, 0.0, 1.0]);
    }
    #[test]
    fn test_dot() {
        let veca = [1.0, 2.0, 2.0];
        let vecb = [2.0, 1.0, 2.0];
        assert_eq!(veca.dot(vecb), 8.0);
    }
    #[test]
    fn test_modulus() {
        let veca = [1.0, 2.0, 2.0];
        assert_eq!(veca.modulus(), 3.0);
    }
    #[test]
    fn test_dot_angle() {
        let veca = [1.0, 0.0, 0.0];
        let vecb = [0.0, 1.0, 0.0];
        assert_eq!(veca.dot_angle(vecb), std::f64::consts::PI / 2.0);
        let veca = [1.0, 0.0, 0.0];
        let vecb = [1.0, 0.0, 0.0];
        assert_eq!(veca.dot_angle(vecb), 0.0);
    }
    #[test]
    fn test_float_multiply() {
        let veca = [1.0, 0.0, 2.0];
        let float = 2.0;
        assert_eq!(veca.float_multiply(float), [2.0, 0.0, 4.0]);
    }
    #[test]
    fn test_normalised_vectors() {
        let vec_elec = [1.0, 0.0, 0.0];
        let vec_mag = [0.0, 1.0, 0.0];
        let n = 20;
        let (omega1, omega2) = get_normalised_vectors(vec_elec, vec_mag, n as f64);
        assert_eq!(omega1, [-30.0, 0.5, 0.0]);
        assert_eq!(omega2, [30.0, 0.5, 0.0]);
        let n = 50;
        let (omega1, omega2) = get_normalised_vectors(vec_elec, vec_mag, n as f64);
        assert_eq!(omega1, [-75.0, 0.5, 0.0]);
        assert_eq!(omega2, [75.0, 0.5, 0.0]);
        let vec_elec = [1.0, 2.0, 1.0];
        let vec_mag = [2.0, 1.0, 0.0];
        let n = 10;
        let (omega1, omega2) = get_normalised_vectors(vec_elec, vec_mag, n as f64);
        assert_eq!(omega1, [1.0, 0.5, 0.0].add([-15.0, -30.0, -15.0]));
        assert_eq!(omega2, [1.0, 0.5, 0.0].add([15.0, 30.0, 15.0]));
    }
    //test the individual energy terms.
    #[test]
    fn test_energies() {
        let n = 6;
        let m = 0.0;
        let n2 = 10.0;
        let n1 = 6.0;
        let n_dash = (m + n2 - n1) / 2.0;
        let n_double_dash = (m - n2 + n1) / 2.0;
        let electric = [4.0, 6.0, -2.0];
        let magnetic = [2.0, 4.0, 6.0];
        let (omega1, omega2) = get_normalised_vectors(electric, magnetic, n as f64);
        assert_eq!(omega1, [-35.0, -52.0, 21.0]);
        assert_eq!(omega2, [37.0, 56.0, -15.0]);
        assert_eq!(omega1.modulus(), f64::sqrt(4370.0));
        assert_eq!(omega2.modulus(), f64::sqrt(4730.0));
        let (alpha1, alpha2) = get_angles(magnetic, omega1, omega2);
        let alpha1_calc = f64::acos(-2.0 * f64::sqrt(15295.0) / 805.0);
        let alpha2_calc = f64::acos(52.0 * f64::sqrt(16555.0) / 16555.0);
        assert_eq!(alpha1, alpha1_calc);
        assert_eq!(alpha2, alpha2_calc);

        assert_eq!(electric.dot(electric), 56.0);
        let energy_1 = _energy_term_1(n as f64, omega1, omega2, n_dash, n_double_dash);
        assert_eq!(
            energy_1,
            (-1.0 / 72.0) + 2.0 * f64::sqrt(4370.0) - 2.0 * f64::sqrt(4730.0)
        );
        let energy2 = _energy_term_2(n as f64, n_dash, n_double_dash, electric, alpha1, alpha2);
        let energy2_ana = (-1.0 * f64::powi(6.0, 4) * 56.0 / 16.0)
            * (17.0 * 36.0 + 19.0 - 12.0 * (8.0 - 4.0 * f64::cos(alpha1_calc + alpha2_calc)));
        assert_eq!(energy2, energy2_ana);
        let energy3 = _energy_term_3(n as f64, n_dash, n_double_dash, alpha1, alpha2, magnetic);
        let energy3_ana =
            (36.0 * 56.0 / 48.0) * (7.0 * 36.0 + 5.0 - 16.0 * f64::sin(alpha1) * f64::sin(alpha2));
        assert_eq!(energy3, energy3_ana);

        let energy4 = _energy_term_4(n as f64, n_dash, n_double_dash, alpha1, alpha2);
        let energy4_ana = 35.0 * (f64::cos(alpha1).powi(2) + f64::cos(alpha2).powi(2))
            - 12.0
                * (4.0 * f64::cos(alpha1).powi(2)
                    + 4.0 * f64::cos(alpha1) * f64::cos(alpha2)
                    + 4.0 * f64::cos(alpha2).powi(2));
        assert_eq!(energy4, energy4_ana);
    }

    #[test]
    fn test_inputs() {
        let n = 10;
        let n1 = 0;
        let n2 = 9;
        let m = 0;

        let res = hydrogen_energy([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], n, n1, n2, m);
        assert_eq!(res, Ok(-1.0 / (2.0 * n as f64 * n as f64)));
    }
    #[test]
    fn test_inputs_wrong() {
        let n = 10;
        let n1 = 0;
        let n2 = 9;
        let m = 1;

        let res = hydrogen_energy([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], n, n1, n2, m);
        assert_eq!(res, Err(HydrogenError {}));
    }
}
