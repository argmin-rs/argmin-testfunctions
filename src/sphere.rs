// Copyright 2018 Stefan Kroboth
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
// http://opensource.org/licenses/MIT>, at your option. This file may not be
// copied, modified, or distributed except according to those terms.

//! # Sphere function

use num::{Float, FromPrimitive};
use std::iter::Sum;

/// Sphere test function
pub fn sphere<T: Float + FromPrimitive + Sum>(param: &[T]) -> T {
    let num2 = T::from_f64(2.0).unwrap();
    param.iter().map(|x| x.powf(num2)).sum()
}

/// Derivative of sphere test function
pub fn sphere_derivative<T: Float + FromPrimitive>(param: &[T]) -> Vec<T> {
    let num2 = T::from_f64(2.0).unwrap();
    param.iter().map(|x| num2 * *x).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std;

    #[test]
    fn test_sphere_optimum() {
        assert!(sphere(&[0.0_f32, 0.0_f32]).abs() < std::f32::EPSILON);
        assert!(sphere(&[0.0_f64, 0.0_f64]).abs() < std::f64::EPSILON);
    }
}
