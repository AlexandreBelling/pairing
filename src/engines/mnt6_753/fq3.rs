use super::{
    fq::{
        Fq,
        FQ3_NQR_T,
        FQ3_T_MINUS_1,
        NON_RESIDUE,
        FROBENIUS_COEFF_FQ3_C1,
        FROBENIUS_COEFF_FQ3_C2
    },
};

use crate::generics::fields::fp3::{ Fp3, Fp3Extension };

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Fq3Extension();

impl Fp3Extension for Fq3Extension {
    type Fp = Fq;
    const FROBENIUS_COEFFICIENTS_C1: [Fq; 3] = FROBENIUS_COEFF_FQ3_C1;
    const FROBENIUS_COEFFICIENTS_C2: [Fq; 3] = FROBENIUS_COEFF_FQ3_C2;
    const NON_RESIDUE: Fq = NON_RESIDUE;
    const QUADRATIC_NONRESIDUE_TO_T: (Fq, Fq, Fq) = FQ3_NQR_T;
    const T_MINUS_1_OVER_2: &'static [u64] = &FQ3_T_MINUS_1;
}

pub type Fq3 = Fp3<Fq3Extension>;

#[test]
fn fq3_field_tests() {
    use ff::PrimeField;

    crate::tests::field::random_field_tests::<Fq3>();
    crate::tests::field::random_sqrt_tests::<Fq3>();
    crate::tests::field::random_frobenius_tests::<Fq3, _>(super::fq::Fq::char(), 13);
}
