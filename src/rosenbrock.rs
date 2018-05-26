// Copyright 2018 Stefan Kroboth
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
// http://opensource.org/licenses/MIT>, at your option. This file may not be
// copied, modified, or distributed except according to those terms.

//! # Rosenbrock function

use num::{Float, FromPrimitive};

/// 2D Rosenbrock test function
///
/// Parameters are usually: `a = 1` and `b = 100`
pub fn rosenbrock<T: Float + FromPrimitive>(param: &[T], a: T, b: T) -> T {
    if let &[x, y] = param {
        (a - x).powi(2) + b * (y - x.powi(2)).powi(2)
    } else {
        panic!("rosenbrock function only accepts 2 parameters.");
    }
}

/// Derivative of 2D Rosenbrock function
pub fn rosenbrock_derivative<T: Float + FromPrimitive>(param: &[T], a: T, b: T) -> Vec<T> {
    let num2 = T::from_f64(2.0).unwrap();
    let num4 = T::from_f64(4.0).unwrap();
    if let &[x, y] = param {
        let mut out = Vec::with_capacity(2);
        out.push(-num2 * a + num4 * b * x.powi(3) - num4 * b * x * y + num2 * x);
        out.push(num2 * b * (y - x.powi(2)));
        out
    } else {
        panic!("rosenbrock function only accepts 2 parameters.");
    }
}

/// Hessian of 2D Rosenbrock function
pub fn rosenbrock_hessian<T: Float + FromPrimitive>(param: &[T], _a: T, b: T) -> Vec<T> {
    let num2 = T::from_f64(2.0).unwrap();
    let num4 = T::from_f64(4.0).unwrap();
    let num12 = T::from_f64(12.0).unwrap();
    if let &[x, y] = param {
        let mut out = Vec::with_capacity(4);
        // d/dxdx
        out.push(num12 * b * x.powi(2) - num4 * b * y + num2);
        // d/dxdy
        out.push(-num4 * b * x);
        // d/dydx
        out.push(-num4 * b * x);
        // d/dydy
        out.push(num2 * b);
        out
    } else {
        panic!("rosenbrock_hessian only accepts 2 parameters.");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std;

    #[test]
    fn test_rosenbrock_optimum() {
        assert!(rosenbrock(&[1.0_f32, 1.0_f32], 1.0, 100.0).abs() < std::f32::EPSILON);
        assert!(rosenbrock(&[1.0, 1.0], 1.0, 100.0).abs() < std::f64::EPSILON);
    }

    #[test]
    fn test_rosenbrock_derivative() {
        let res: Vec<f64> = rosenbrock_derivative(&[1.0, 1.0], 1.0, 100.0);
        for elem in &res {
            assert!((elem - 0.0).abs() < std::f64::EPSILON);
        }
        let res: Vec<f32> = rosenbrock_derivative(&[1.0_f32, 1.0_f32], 1.0_f32, 100.0_f32);
        for elem in &res {
            assert!((elem - 0.0).abs() < std::f32::EPSILON);
        }
    }

    #[test]
    #[should_panic]
    fn test_rosenbrock_not_2d() {
        rosenbrock(&[1.0, 1.0, 1.0], 1.0, 100.0);
    }
}
