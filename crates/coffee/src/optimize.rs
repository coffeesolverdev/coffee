use crate::extras::{OptimizerArgs, OptimizerError, OptimizerResults};
use crate::format::{conclude_message, process_message, start_message};
use crate::steihaug::Steihaug;
use chrono::Utc;
use core::f64;
use ndarray::{Array1, Array2, ArrayView1, Axis};
use std::error::Error;

/// Cuts off values smaller than e^(this value) due to lack of precision in f64.
const SMALLEST_EXP_VALUE: f64 = -230.0;

pub struct Optimizer {
    monomers: Array1<f64>,
    polymers: Array2<f64>,
    polymers_q: Array1<f64>,
    max_iterations: usize,
    curr_iteration: usize,
    time_us: usize,
    delta: f64,
    max_delta: f64,
    eta: f64,
    norm_ratio_threshold: f64,
    rho_thresholds: [f64; 2],
    scale_factors: [f64; 2],
    optimal_lambda: Array1<f64>,
    optimal_x: Array1<f64>,
    optimal_lagrangian: f64,
    steihaug_trust_region: Steihaug,
    use_terminal: bool,
    verbose: bool,
    log_msgs: Vec<String>,
    scalarity: bool,
    temp_celsius: f64,
}

/// Caclulates the density of water at a given temperature.
/// This is based on the IAPWS-95 formulation for the density of water.
/// We use this density as the mininum volume unit for the monomers and polymers.
///
/// # Arguments
///
/// * `temp` - The temperature in Kelvin.
///
/// # Returns
///
/// The density of water at the given temperature in g/cm^3.
fn density_water(t: f64) -> f64 {
    let a1 = -3.983035;
    let a2 = 301.797;
    let a3 = 522528.9;
    let a4 = 69.34881;
    let a5 = 999.974950;
    a5 * (1. - (t + a1) * (t + a1) * (t + a2) / a3 / (t + a4)) / 18.0152
}

/// Creates a new `Optimizer` instance with the given parameters.
///
/// # Arguments
///
/// * `monomers` - A reference to a 1-dimensional array of monomer concentrations.
/// * `polymers` - A reference to a 2-dimensional array representing the polymer matrix.
/// * `polymers_q` - A reference to a 1-dimensional array of polymer quantities.
/// * `optional_args` - An instance of `OptimizerArgs` containing optional parameters for the optimizer.
///
/// # Returns
///
/// A new instance of `Optimizer`.
///
/// # Panics
///
/// This function will panic if:
/// - `monomers` is empty.
/// - `polymers` is empty.
/// - The number of polymers is less than the number of monomers.
impl Optimizer {
    pub fn new(
        monomers: &Array1<f64>,
        polymers: &Array2<f64>,
        polymers_q_nonexp: &Array1<f64>,
        optional_args: &OptimizerArgs,
    ) -> Result<Self, Box<dyn Error>> {
        let num_monomers = monomers.len();
        let num_polymers = polymers.len_of(Axis(0));

        if num_monomers == 0 {
            return Err(Box::new(OptimizerError(
                "Monomers array is empty.".to_string(),
            )));
        }
        if num_polymers == 0 {
            return Err(Box::new(OptimizerError(
                "Polymers array is empty.".to_string(),
            )));
        }
        if num_polymers < num_monomers {
            return Err(Box::new(OptimizerError(
                "Number of polymers is less than number of monomers.".to_string(),
            )));
        }

        /* Check sizes between the arrays. */
        if num_monomers != polymers.len_of(Axis(1)) {
            return Err(Box::new(OptimizerError(
                "Monomers and polymer compositions inconsistent.".to_string(),
            )));
        }
        if num_polymers != polymers_q_nonexp.len() {
            return Err(Box::new(OptimizerError(
                "Polymers and polymer quantities have different sizes.".to_string(),
            )));
        }

        /* Scale for water molecule volume size if necessary. */
        let temp_celsius = optional_args.temp_celsius;
        let scalarity = optional_args.scalarity;
        let k_t = if scalarity {
            0.00198717 * (temp_celsius + 273.15)
        } else {
            1.0
        };
        let scaled_monomers = if scalarity {
            monomers / density_water(temp_celsius)
        } else {
            monomers.clone()
        };
        let polymers_q = polymers_q_nonexp.mapv(|x| (-x.max(SMALLEST_EXP_VALUE) / k_t).exp());

        let max_iterations = optional_args.max_iterations;
        Ok(Optimizer {
            monomers: scaled_monomers,
            polymers: polymers.clone(),
            polymers_q,
            max_iterations,
            curr_iteration: 0,
            time_us: 0,
            delta: 1.0,
            max_delta: optional_args.max_delta,
            eta: optional_args.eta,
            norm_ratio_threshold: optional_args.norm_ratio_threshold,
            rho_thresholds: optional_args.rho_thresholds,
            scale_factors: optional_args.scale_factors,
            optimal_lambda: Array1::zeros(num_monomers),
            optimal_x: Array1::zeros(num_polymers),
            optimal_lagrangian: 0.0,
            steihaug_trust_region: Steihaug::new(max_iterations, num_monomers),
            use_terminal: optional_args.use_terminal,
            verbose: optional_args.verbose,
            log_msgs: Vec::new(),
            scalarity,
            temp_celsius,
        })
    }

    /// Calculates the norm of the given vector. Replaces ndarray-linalg crate's implementation
    /// because it isn't working with WASM. Equivalent to `ndarray_linalg::Norm::norm2`, which is
    /// the implementation of pythagorean theorem for N elements.
    ///
    /// # Arguments
    ///
    /// * `v` - A 1-dimensional array representing the vector to calculate the norm of.
    ///
    /// # Returns
    ///
    /// The norm of the vector.
    fn norm(&self, v: ArrayView1<f64>) -> f64 {
        v.iter().map(|&x| x * x).sum::<f64>().sqrt()
    }

    /// Updates the optimal x values based on the current polymer lambdas and polymer quantities.
    /// This is used to calculate the optimal concentrations of the polymers.
    /// It also scales the values based on the temperature and whether scalarity is enabled.
    /// No output is needed as it is automatically updated internally.
    fn update_optimal_x(&mut self) {
        if self.scalarity {
            self.optimal_x =
                &self.polymers_q * &self.polymer_lambdas() * density_water(self.temp_celsius);
        } else {
            self.optimal_x = &self.polymers_q * &self.polymer_lambdas();
        }
    }

    /// Returns the current polymer lambdas given the current monomer lambda values.
    /// This is used to calculate the polymer quantities.
    ///
    /// # Returns
    ///
    /// A 1-dimensional array of polymer lambdas in exponential form.
    ///
    /// # Panics
    ///
    /// This function will panic if the non-exponentiated polymer concentrations are not finite.
    /// This is to ensure that the optimization is working correctly.
    fn polymer_lambdas(&self) -> Array1<f64> {
        (self.polymers.dot(&self.optimal_lambda)).exp()
    }

    /// Calculates the Lagrangian of the optimization using the current lambda and also:
    /// - The polymer quantities.
    /// - The polymer energies.
    /// - The monomer concentrations.
    ///
    /// # Returns
    ///
    /// The Lagrangian scalar value.
    ///
    /// # Panics
    ///
    /// This function will panic if the Lagrangian value is not finite.
    /// This is to ensure that the optimization is working correctly.
    fn lagrangian(&self, polymer_lambdas: &Array1<f64>) -> f64 {
        let after_energies = self.polymers_q.dot(polymer_lambdas);
        let after_initial = self.optimal_lambda.dot(&self.monomers);

        (after_energies - after_initial).ln()
    }

    fn jacobian(&self, polymer_lambdas: &Array1<f64>, lagrangian: f64) -> Array1<f64> {
        let after_energies = &self.polymers_q * polymer_lambdas;
        let jacobian = self.polymers.t().dot(&after_energies) - &self.monomers;
        jacobian / lagrangian.exp()
    }

    /// Calculates the Hessian of the optimization using the current lambda and also:
    /// - The polymer quantities.
    /// - The polymer energies.
    /// - The monomer concentrations.
    ///
    /// # Returns
    ///
    /// A 2-dimensional array representing the Hessian. Its size is M x M, where M is the number of monomers.
    fn hessian(
        &self,
        polymer_lambdas: &Array1<f64>,
        lagrangian: f64,
        jacobian: &Array1<f64>,
    ) -> Array2<f64> {
        let first_part = 1. / lagrangian.exp();

        let after_energies = &self.polymers_q * polymer_lambdas;
        /* Truncate it to M x M, and element wise multiply with the polymers matrix. */
        let polymerization = &self.polymers * &after_energies.insert_axis(Axis(1));

        let second_part = self.polymers.t().dot(&polymerization);

        let fourth_part = jacobian.view().insert_axis(Axis(1));

        let fifth_part = jacobian.view().insert_axis(Axis(0));

        first_part * second_part - fourth_part.dot(&fifth_part)
    }

    /// Optimizes the given function using the Steihaug trust region method.
    /// Requires an initial delta value to start the optimization.
    /// Initialized with the monomer concentrations, exponentiated polymer energies, and the polymer quantities.
    ///
    /// # Arguments
    ///
    /// * `initial_delta` - The initial delta value to start the optimization.
    ///
    /// # Returns
    ///
    /// A 1-dimensional array representing the optimal x values. Its size is N, where N is the number of polymers.
    ///
    /// # Panics
    ///
    /// This function will panic if the calculations are not finite.
    /// This is to ensure that the optimization is working correctly.
    pub fn optimize(&mut self, initial_delta: f64) -> Result<bool, Box<dyn Error>> {
        /* Error Check for delta value. */
        if initial_delta <= 0.0 || !initial_delta.is_finite() {
            return Err(Box::new(OptimizerError(
                "Initial delta value is not valid.".to_string(),
            )));
        }

        self.print(&start_message());

        /* Initialization and resetting from previous optimizations. */
        self.delta = initial_delta;
        let mut final_it = 0;
        self.reset();
        let start_time = Utc::now();

        /* Start of optimization. */
        for it in 0..self.max_iterations {
            /* Calculate mathematical values to generate predictions for changes. */
            let polymer_lambdas = self.polymer_lambdas();
            self.optimal_lagrangian = self.lagrangian(&polymer_lambdas);

            let function = self.optimal_lagrangian;
            let gradient = self.jacobian(&polymer_lambdas, self.optimal_lagrangian);
            let hessian = self.hessian(&polymer_lambdas, self.optimal_lagrangian, &gradient);

            /* Get miscellaneous difference values for the steihaug calculation. */
            let step = self.norm(gradient.view());
            let epsilon = step.sqrt().min(0.5f64) * step;

            /* Find predicted next step through the steihaug trust method. */
            let success = self
                .steihaug_trust_region
                .iterate(&gradient, &hessian, epsilon, self.delta);
            if !success {
                /* Conclude the optimization prematurely as it failed. */
                self.time_us = (Utc::now() - start_time)
                    .num_microseconds()
                    .unwrap_or_default() as usize;
                self.print(&conclude_message(
                    it,
                    success,
                    self.time_us,
                    self.verbose,
                    None,
                ));

                return Err(Box::new(OptimizerError(
                    "The Steihaug optimization did not succeed".to_string(),
                )));
            }
            let update_step = self.steihaug_trust_region.get_result();

            /* Pre-emptively update optimal lambdas and their math calcs to find whether reduction is accurate. */
            self.optimal_lambda = &self.optimal_lambda + &update_step;
            self.optimal_lagrangian = self.lagrangian(&self.polymer_lambdas());

            /* Find predicted and actual reductions to see how significant the optimizing change is. */
            let pred_reduction =
                -(gradient.dot(&update_step) + 0.5 * update_step.dot(&hessian.dot(&update_step)));
            let actual_reduction = function - self.optimal_lagrangian;

            /* No more optimization is needed as there is no optimizing change. */
            if actual_reduction == 0.0 {
                final_it = it;
                break;
            }

            /* Ratio calculation to determine next iteration's parameters. */
            let rho = if pred_reduction != 0.0 {
                actual_reduction / pred_reduction
            } else {
                0.0
            };

            /* Change delta based on whether reductions is too small or too high. */
            if rho < self.rho_thresholds[0] {
                /* Actual reduction is much less than predicted --> scale down delta param. */
                self.delta *= self.scale_factors[0];
            } else if rho > self.rho_thresholds[1]
                && self.norm(update_step.view()) >= self.norm_ratio_threshold * self.delta
            {
                /* Actual reduction is close to predicted --> scale up delta param up to a point. */
                self.delta = self.max_delta.min(self.scale_factors[1] * self.delta);
            }

            /* Actual reduction is scary less than predicted --> can't trust steihaug update value. */
            if rho <= self.eta {
                /* Remove the update from lambda if quadratic isn't reliable. Update relevant values for debugging. */
                self.optimal_lambda = &self.optimal_lambda - &update_step;
                self.optimal_lagrangian = self.lagrangian(&self.polymer_lambdas());
            }

            /* Calculate backtrack (error) by updating optimal_x to latest vals. */
            self.update_optimal_x();
            self.print(&process_message(it, self.optimal_lagrangian, self.error()));

            /* Update iteration. */
            final_it = it;
            self.curr_iteration += 1;
        }

        /* Find the optimal concentrations. */
        self.update_optimal_x();

        /* Calculate optimization time and print concluding results. */
        self.time_us = (Utc::now() - start_time)
            .num_microseconds()
            .unwrap_or_default() as usize;

        self.print(&conclude_message(
            final_it,
            true,
            self.time_us,
            self.verbose,
            Some(&OptimizerResults {
                optimal_x: self.optimal_x.to_vec(),
                optimal_lagrangian: self.optimal_lagrangian,
                optimal_lambda: self.optimal_lambda.to_vec(),
                concentration_error: self.error(),
                log_messages: self.log_msgs.clone(),
                elapsed_time: self.time_us,
            }),
        ));

        Ok(true)
    }

    /// Resets the optimizer to its initial state.
    /// This is useful when reusing the optimizer for multiple optimizations.
    /// It resets the lambda values and the x values.
    pub fn reset(&mut self) {
        self.curr_iteration = 0;
        self.time_us = 0;
        self.optimal_lambda.fill(0.);
        self.optimal_x.fill(0.);
        self.optimal_lagrangian = 0.0;
        self.log_msgs.clear();
    }

    /// Returns the optimal results of the optimization.
    /// This is a owned version, so it deepcopies the results.
    ///
    /// # Returns
    ///
    /// (Optimal X Values, Optimal Lagrangian, Optimal Lambda Values)
    pub fn get_results(&self) -> OptimizerResults {
        OptimizerResults {
            optimal_x: self.optimal_x.to_vec(),
            optimal_lagrangian: self.optimal_lagrangian,
            optimal_lambda: self.optimal_lambda.to_vec(),
            concentration_error: self.error(),
            log_messages: self.log_msgs.clone(),
            elapsed_time: self.time_us,
        }
    }

    fn print(&mut self, msg: &str) {
        if self.use_terminal {
            print!("{}", msg);
        } else {
            self.log_msgs.push(msg.to_string());
        }
    }

    /// Returns the time taken for the optimization in microseconds.
    /// This is useful for benchmarking the optimization.
    ///
    /// # Returns
    ///
    /// The time taken for the optimization in ms.
    pub fn benchmark(&self) -> usize {
        self.time_us
    }

    /// Finds the backtrack error of each configuration value in polymers, returns the maximum error.
    /// This is used to determine how accurate the optimization is compared to
    /// the original monomer concentrations (e.g. how well is mass conserved).
    /// Based on the current optimal x values.
    ///
    /// # Returns
    ///
    /// A f64 value representing the maximum error .
    /// This is the maximum difference between the monomer concentrations and the polymer concentrations.
    fn error(&self) -> f64 {
        let concs = self
            .polymers
            .t()
            .dot(&self.optimal_x.view().insert_axis(Axis(1)));
        let scaling = if self.scalarity {
            density_water(self.temp_celsius)
        } else {
            1.0
        };
        let backtrack = (&self.monomers * scaling).insert_axis(Axis(1)) - concs;
        backtrack
            .iter()
            .fold(f64::NEG_INFINITY, |a, &b| f64::max(a, b.abs()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_wrong_size_params() {
        /* Mismatch between polymers and monomers. */
        let monomers = array![1.0e-3, 2.0e-3];
        let polymers = array![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
            [1.0, 1.0, 1.0]
        ];
        let polymers_q = array![0.0, 0.0, -1.0e+3, -2.0e+3];
        let args = OptimizerArgs::default();
        let result = Optimizer::new(&monomers, &polymers, &polymers_q, &args);
        assert!(result.is_err());

        /* Mismatch between polymers and polymers_q. */
        let polymers = array![[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
        let polymers_q = array![0.0, 0.0, -1.0e+3, -2.0e+3];
        let result = Optimizer::new(&monomers, &polymers, &polymers_q, &args);
        assert!(result.is_err());

        /* Polymers must be greater than monomers. */
        let polymers = array![[1.0, 0.0],];
        let polymers_q = array![0.0];
        let result = Optimizer::new(&monomers, &polymers, &polymers_q, &args);
        assert!(result.is_err());
    }
}
