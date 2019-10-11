use crate::{
    generics::fields::fp2::{ 
        Fp2, 
        Fp2Extension 
    },
};

use super::fq::{
    Fq, 
    FROBENIUS_COEFF_FQ2_C1,
    NON_RESIDUE,
};

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Fq2Extension();

impl Fp2Extension for Fq2Extension {
    type Fp = Fq;

    const FROBENIUS_COEFFICIENTS: [Fq; 2] = FROBENIUS_COEFF_FQ2_C1;
    const NON_RESIDUE: Fq = NON_RESIDUE;
}

pub type Fq2 = Fp2<Fq2Extension>;

impl Fp2<Fq2Extension> {
    #[inline(always)]
    pub fn mul_by_nonresidue(&mut self) {
        ::std::mem::swap(&mut self.c0, &mut self.c1);
        Fq2Extension::mul_by_nonresidue(&mut self.c0);
    }
}

#[cfg(test)]
use crate::{
    rand::{
        Rand,
        SeedableRng,
        XorShiftRng
    },
    ff::{
        Field,
        SqrtField,
    },
};

#[test]
fn fq2_frobenius_map_test() {
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    for _ in 1..100 {
        let a = Fq2::rand(&mut rng);
        let mut b = a;
        for i in 1..3 {
            println!("{}", i);
            let mut c = a;
            c.frobenius_map(i);
            b.frobenius_map(1);
            assert_eq!(b, c);
        }
        assert_eq!(a, b);
    }
}

#[test]
fn test_fq2_legendre() {
    use ff::LegendreSymbol::*;

    assert_eq!(Zero, Fq2::zero().legendre());
    // i^2 = -1
    let mut m1 = Fq2::one();
    m1.negate();
    assert_eq!(QuadraticResidue, m1.legendre());
    m1.mul_by_nonresidue();
    assert_eq!(QuadraticNonResidue, m1.legendre());
}

#[test]
fn fq2_field_tests() {
    use ff::PrimeField;

    crate::tests::field::random_field_tests::<Fq2>();
    crate::tests::field::random_sqrt_tests::<Fq2>();
    crate::tests::field::random_frobenius_tests::<Fq2, _>(super::fq::Fq::char(), 13);
}
