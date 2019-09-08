use ff::{
    LegendreSymbol::{
        Zero,
        QuadraticResidue,
        QuadraticNonResidue
    },
    Field, 
    SqrtField,
};

use rand::{Rand, Rng};
use std::{
    cmp::Ordering,
    fmt::Debug,
};

pub trait Fp2Extension: 'static + Copy + Debug + Eq {
    
    type Fp: Field;

    const FROBENIUS_COEFFICIENTS: [Self::Fp; 2];
    const NON_RESIDUE: Self::Fp;

    #[inline(always)]
    fn mul_by_nonresidue(x: &mut Self::Fp) {
        x.mul_assign(&Self::NON_RESIDUE);
    }
}

/// An element of Fp2, represented by c0 + c1 * u.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Fp2<P: Fp2Extension> {
    pub c0: P::Fp,
    pub c1: P::Fp,
}

impl<P: Fp2Extension> ::std::fmt::Display for Fp2<P> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fp2({} + {} * u)", self.c0, self.c1)
    }
}

impl<P: Fp2Extension> PartialOrd for Fp2<P> where P::Fp: Ord + Field {
    #[inline(always)]
    fn partial_cmp(&self, other: &Fp2<P>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// `Fp2` elements are ordered lexicographically.
impl<P: Fp2Extension> Ord for Fp2<P> where P::Fp: Ord + Field {
    #[inline(always)]
    fn cmp(&self, other: &Fp2<P>) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
        }
    }
}

impl<P: Fp2Extension> Rand for Fp2<P> {
    fn rand<R: Rng>(rng: &mut R) -> Self {
        Self {
            c0: rng.gen(),
            c1: rng.gen(),
        }
    }
}

impl<P: Fp2Extension> Fp2<P> {

    /// Norm of Fp2 as extension field in i over P::Fp
    #[inline(always)]
    pub fn mul_assign_by_fp(&mut self, f: &P::Fp) {
        self.c0.mul_assign(f);
        self.c1.mul_assign(f);
    }

    /// Norm of Fp2 as extension field in i over P::Fp
    #[inline(always)]
    pub fn norm(&self) -> P::Fp {
        let mut t0 = self.c0;
        t0.square();
        let mut t1 = self.c1;
        t1.square();
        P::mul_by_nonresidue(&mut t1);
        t0.sub_assign(&t1);
        t0
    }
}

impl<P: Fp2Extension> Field for Fp2<P> {
    #[inline(always)]
    fn zero() -> Self {
        Self {
            c0: P::Fp::zero(),
            c1: P::Fp::zero(),
        }
    }

    #[inline(always)]
    fn one() -> Self {
        Self {
            c0: P::Fp::one(),
            c1: P::Fp::zero(),
        }
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
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
        let mut t0 = self.c0;
        t0.square();
        let mut t1 = self.c1;
        t1.square();
        P::mul_by_nonresidue(&mut t1);
        t0.sub_assign(&t1);
        t0.inverse().map(|t1| {
            let mut tmp = Self {
                c0: self.c0,
                c1: self.c1,
            };
            tmp.c0.mul_assign(&t1);
            tmp.c1.mul_assign(&t1);
            tmp.c1.negate();

            tmp
        })
    }

    #[inline(always)]
    fn frobenius_map(&mut self, power: usize) {
        self.c1.mul_assign(&P::FROBENIUS_COEFFICIENTS[power % 2]);
    }
}

impl<P: Fp2Extension> SqrtField for Fp2<P> where P::Fp: SqrtField {

    #[inline(always)]
    fn legendre(&self) -> ::ff::LegendreSymbol {
        self.norm().legendre()
    }

    #[inline(always)]
    fn sqrt(&self) -> Option<Self> {
        if self.c1.is_zero() {
            return self.c0.sqrt().map(|c0| Self{c0: c0, c1: P::Fp::zero()});
        }
        match self.legendre() {
            // Square root based on the complex method. See
            // https://eprint.iacr.org/2012/685.pdf (page 15, algorithm 8)
            Zero => Some(*self),
            QuadraticNonResidue => None,
            QuadraticResidue => {
                let mut two_inv = P::Fp::one();
                two_inv.double();
                two_inv = two_inv.inverse().expect("Two should always have an inverse");

                let alpha = self.norm().sqrt().expect("We are in the QR case, the norm should have a square root");
                //let mut delta = (alpha + &self.c0) * &two_inv;
                let mut delta = self.c0;
                delta.add_assign(&alpha);
                delta.mul_assign(&two_inv);
                if delta.legendre() == QuadraticNonResidue {
                    delta.sub_assign(&alpha);
                }

                let c0 = delta.sqrt().expect("Delta must have a square root");
                let c0_inv = c0.inverse().expect("c0 must have an inverse");

                let mut c1 = self.c1;
                c1.mul_assign(&two_inv);
                c1.mul_assign(&c0_inv);
                Some(Self{c0: c0, c1: c1})
            },
        }
    }
}