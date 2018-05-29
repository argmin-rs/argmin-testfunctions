// Copyright 2018 Stefan Kroboth
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
// http://opensource.org/licenses/MIT>, at your option. This file may not be
// copied, modified, or distributed except according to those terms.

//! # Test functions for optimization algorithms

extern crate num;

mod ackley;
mod rastrigin;
mod rosenbrock;
mod sphere;
mod zero;

pub use ackley::*;
pub use rastrigin::*;
pub use rosenbrock::*;
pub use sphere::*;
pub use zero::*;
