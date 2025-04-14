use std::error::Error;
use std::fmt;

/// Struct containing optional parameters for the optimizer.
/// These parameters can be used to customize the optimization process.
/// Has Default values that can be overridden.
#[derive(Clone)]
pub struct OptimizerArgs {
    pub max_iterations: usize,
    pub max_delta: f64,
    pub eta: f64,
    pub norm_ratio_threshold: f64,
    pub rho_thresholds: [f64; 2],
    pub scale_factors: [f64; 2],
    pub use_terminal: bool,
    pub scalarity: bool,
    pub temp_celsius: f64,
    pub verbose: bool,
}

#[derive(Clone)]
pub struct OptimizerResults {
    pub optimal_x: Vec<f64>,
    pub optimal_lagrangian: f64,
    pub optimal_lambda: Vec<f64>,
    pub concentration_error: f64,
    pub log_messages: Vec<String>,
    pub elapsed_time: usize,
}

/// Default implementation for `OptimizerArgs`.
/// Contains default values for all optional parameters.
impl Default for OptimizerArgs {
    fn default() -> Self {
        OptimizerArgs {
            max_iterations: 250,
            max_delta: 1000.0,
            eta: 0.15,
            norm_ratio_threshold: 0.95,
            rho_thresholds: [0.25, 0.75],
            scale_factors: [0.25, 2.0],
            use_terminal: true,
            scalarity: true,
            temp_celsius: 37.0,
            verbose: false,
        }
    }
}

#[derive(Debug)]
pub struct OptimizerError(pub String);

impl fmt::Display for OptimizerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Error for OptimizerError {}
