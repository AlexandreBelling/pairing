#![allow(unused_imports)]

use ff::{
    Field, 
};

use crate::BitIterator;
use rand::{Rand, Rng};
use std::{
    cmp::Ordering,
    fmt::Debug,
    mem::swap,
};

use super::fp2::{
    Fp2Extension,
    Fp2,
};

pub trait Fp4Extension: 'static + Copy + Debug + Eq + Sync + Send { 
    type Fp2P: Fp2Extension;

    const FROBENIUS_COEFFICIENTS_C1: [<Self::Fp2P as Fp2Extension>::Fp; 4];

    const NON_RESIDUE: <Self::Fp2P as Fp2Extension>::Fp;

    #[inline(always)]
    fn mul_by_nonresidue(x: &mut Fp2<Self::Fp2P>) {
        swap(&mut x.c0, &mut x.c1);
        Self::Fp2P::mul_by_nonresidue(&mut x.c0);
    }
}

/// An element of Fp2, represented by c0 + c1 * u.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Fp4<P: Fp4Extension> {
    pub c0: Fp2<P::Fp2P>,
    pub c1: Fp2<P::Fp2P>,
}

impl<P: Fp4Extension> ::std::fmt::Display for Fp4<P> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fp4({} + {} * w)", self.c0, self.c1)
    }
}

impl<P: Fp4Extension> Rand for Fp4<P> {
    fn rand<R: Rng>(rng: &mut R) -> Self {
        Self {
            c0: rng.gen(),
            c1: rng.gen(),
        }
    }
}

impl<P: Fp4Extension> Fp4<P> {

    #[inline(always)]
    pub fn conjugate(&mut self) {
        self.c1.negate();
    }

    // When the Fp4 element is known to be an r-th root of 
    // unity, we can use this function instead of squaring
    // TODO: Implement NAF
    #[inline(always)]
    pub fn cyclotomic_exp<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one();
        let mut found_one = false;

        for i in BitIterator::new(exp) {
            if found_one { res.cyclotomic_square(); }
            else { found_one = i; }
            if i { res.mul_assign(self); }
        }
        res
    }

    // When the Fq4 element is known to be an r-th root of 
    // unity, we can use this function instead of squaring
    #[inline(always)]
    pub fn cyclotomic_square(&mut self) {
        swap(&mut self.c0, &mut self.c1);
        self.c1.add_assign(&self.c0);
        self.c1.square();
        self.c0.square();
        self.c1.sub_assign(&self.c0);
        P::mul_by_nonresidue(&mut self.c0);
        self.c1.sub_assign(&self.c0);
        self.c0.double();
        let one = Fp2::<P::Fp2P>::one();
        self.c0.add_assign(&one);
        self.c1.sub_assign(&one);
    }
}

impl<P: Fp4Extension> Field for Fp4<P> {
    
    #[inline(always)]
    fn zero() -> Self {
        Self {
            c0: Fp2::<P::Fp2P>::zero(),
            c1: Fp2::<P::Fp2P>::zero(),
        }
    }

    #[inline(always)]
    fn one() -> Self {
        Self {
            c0: Fp2::<P::Fp2P>::one(),
            c1: Fp2::<P::Fp2P>::zero(),
        }
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    #[inline(always)]
    fn double(&mut self) {
        self.c0.double();
        self.c1.double();
    }

    #[inline(always)]
    fn negate(&mut self) {
        self.c0.negate();
        self.c1.negate();
    }

    #[inline(always)]
    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }

    #[inline(always)]
    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }

    #[inline(always)]
    fn frobenius_map(&mut self, power: usize) {
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);

        self.c1.c0.mul_assign(&P::FROBENIUS_COEFFICIENTS_C1[power % 4]);
        self.c1.c1.mul_assign(&P::FROBENIUS_COEFFICIENTS_C1[power % 4]);
    }

    #[inline(always)]
    fn square(&mut self) {
        // Devegili OhEig Scott Dahab
        // --- 
        // Multiplication and Squaring on Pairing-Friendly Fields.pdf; 
        // Section 3 (Karatsuba squaring)
        let mut aa = self.c0;
        aa.square();
        let mut bb = self.c1;
        bb.square();
        self.c1.mul_assign(&self.c0);
        self.c1.double();
        self.c0 = aa;
        P::mul_by_nonresidue(&mut bb);
        self.c0.add_assign(&bb);
    }

    #[inline(always)]
    fn mul_assign(&mut self, other: &Self) {
        let mut aa = self.c0;
        aa.mul_assign(&other.c0);
        let mut bb = self.c1;
        bb.mul_assign(&other.c1);
        let mut o = other.c0;
        o.add_assign(&other.c1);
        self.c1.add_assign(&self.c0);
        self.c1.mul_assign(&o);
        self.c1.sub_assign(&aa);
        self.c1.sub_assign(&bb);
        self.c0 = aa;
        P::mul_by_nonresidue(&mut bb);
        self.c0.add_assign(&bb);
    }

    #[inline(always)]
    fn inverse(&self) -> Option<Self> {
        let mut c0s = self.c0;
        c0s.square();
        let mut c1s = self.c1;
        c1s.square();
        P::mul_by_nonresidue(&mut c1s);
        c0s.sub_assign(&c1s);

        c0s.inverse().map(|t| {
            let mut tmp = Self { c0: t, c1: t };
            tmp.c0.mul_assign(&self.c0);
            tmp.c1.mul_assign(&self.c1);
            tmp.c1.negate();

            tmp
        })
    }
}