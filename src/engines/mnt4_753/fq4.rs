extern crate std;

use super::{
    fq::{Fq, FROBENIUS_COEFF_FQ4_C1},
    fq2::{ Fq2Extension },
};

use crate::{
    generics::fields::{
        fp4_as_2_over_2::{
            Fp4,
            Fp4Extension,
        },
    },
};

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Fq4Extension();

impl Fp4Extension for Fq4Extension {
    type Fp2P = Fq2Extension;

    const FROBENIUS_COEFFICIENTS_C1: [Fq; 4] = FROBENIUS_COEFF_FQ4_C1;
}

pub type Fq4 = Fp4<Fq4Extension>;

#[cfg(test)]


#[test]
fn fq4_field_tests() {
    use ff::PrimeField;

    crate::tests::field::random_field_tests::<Fq4>();
    crate::tests::field::random_frobenius_tests::<Fq4, _>(super::fq::Fq::char(), 13);
}
