extern crate std;

use super::{
    fq::{FROBENIUS_COEFF_FQ6_C1, Fq},
    fq3::{ Fq3Extension },
};

use crate::{
    generics::fields::{
        fp6_as_2_over_3::{
            Fp6,
            Fp6Extension,
        },
    },
};

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Fq6Extension();

impl Fp6Extension for Fq6Extension {
    type Fp3P = Fq3Extension;
    const FROBENIUS_COEFFICIENTS_C1: [Fq; 6] = FROBENIUS_COEFF_FQ6_C1;
}

pub type Fq6 = Fp6<Fq6Extension>;

#[test]
fn fq6_field_tests() {
    use ff::PrimeField;

    crate::tests::field::random_field_tests::<Fq6>();
    crate::tests::field::random_frobenius_tests::<Fq6, _>(super::fq::Fq::char(), 13);
}