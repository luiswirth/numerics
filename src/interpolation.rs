use crate::{altsign, fx};

use std::f64::consts::PI;

/// 1D Polynomial interpolator using barycentric interpolation formula
/// for the reference interval [-1,+1]
pub struct BaryPolyInterpolator1dRef {
  nodes: na::DVector<fx>,
  weights: na::DVector<fx>,
}
impl BaryPolyInterpolator1dRef {
  pub fn new_chebychev(n: usize) -> Self {
    let mut nodes = na::DVector::zeros(n);
    let mut weights = na::DVector::zeros(n);
    for k in 0..n {
      let angle = (k as fx + 0.5) * PI / n as fx;
      nodes[k] = angle.cos();
      // 2^(q-1) removed, because of cancellation
      weights[k] = altsign(k) * angle.sin();
    }
    Self { nodes, weights }
  }

  pub fn new_generic(nodes: na::DVector<fx>) -> Self {
    let n = nodes.len();
    let mut weights = na::DVector::zeros(n);
    for i in 0..n {
      weights[i] = 1.0;
      for j in 0..n {
        if i == j {
          continue;
        }
        weights[i] *= nodes[i] - nodes[j];
      }
      weights[i] = weights[i].recip();
    }
    Self { nodes, weights }
  }

  /// polynomial interpolant evaluation using barycentric interpolation formula
  pub fn evaluate(
    &self,
    // function values at interpolation points
    interp_values: &na::DVector<fx>,
    // evaluation points
    eval_points: &na::DVector<fx>,
  ) -> na::DVector<fx> {
    const EPS: f64 = 1e-14;

    assert_eq!(interp_values.len(), self.nodes.len());
    let q = self.nodes.len();

    let mut eval_values = na::DVector::zeros(eval_points.len());
    for k in 0..eval_points.len() {
      let eval_point = eval_points[k];

      let mut numerator = 0.0;
      let mut denominator = 0.0;
      for i in 0..q {
        let interp_point = self.nodes[i];
        let interp_value = interp_values[i];
        let weight = self.weights[i];

        let diff = eval_point - interp_point;

        if diff.abs() < EPS {
          numerator = interp_value;
          denominator = 1.0;
          break;
        }

        let term = weight / diff;
        numerator += term * interp_value;
        denominator += term;
      }
      numerator /= denominator;
      eval_values[k] = numerator;
    }
    eval_values
  }
}
