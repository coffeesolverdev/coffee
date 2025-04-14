use ndarray::{Array1, Array2, ArrayView1};

pub struct Steihaug {
    curr_iterations: usize,
    max_iterations: usize,
    vector_size: usize,
    curr_zstep: Array1<f64>,
    curr_rstep: Array1<f64>,
    curr_dstep: Array1<f64>,
}

/// Steihaug's method for solving trust region subproblems.
/// This method is used to solve the trust region subproblem:
/// min f(x) = 1/2 x^T H x + g^T x
/// s.t. ||x|| <= delta
/// where H is the Hessian matrix and g is the gradient.
/// The method is based on the paper by Steihaug (1983).
///
/// # Arguments
///
/// * `tolerance` - The tolerance for the norm of the gradient.
/// * `max_iterations` - The maximum number of iterations.
/// * `vector_size` - The size of the vectors.
impl Steihaug {
    pub fn new(max_iterations: usize, vector_size: usize) -> Self {
        Self {
            curr_iterations: 0,
            max_iterations,
            vector_size,
            curr_zstep: Array1::zeros(vector_size),
            curr_rstep: Array1::zeros(vector_size),
            curr_dstep: Array1::zeros(vector_size),
        }
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

    /// Solve the quadratic equation for the curvature.
    ///
    /// # Arguments
    ///
    /// * `delta` - The trust region radius.
    ///
    /// # Returns
    ///
    /// * `Option<f64>` - The solution to the quadratic equation.
    fn solve_curvature_quadratic(&self, delta: f64) -> Option<f64> {
        /* Operation is d * d, returning scalar. */
        let a = self.curr_dstep.dot(&self.curr_dstep);

        /* Operation is 2 * (z * d), returning scalar. */
        let b = 2.0 * self.curr_zstep.dot(&self.curr_dstep);

        /* Operation is z * z - delta^2, returning scalar. */
        let c = self.curr_zstep.dot(&self.curr_zstep) - delta * delta;

        /* Solve for real solution for quadratic equation given coefficients. */
        let t = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);

        if t.is_finite() {
            Some(t)
        } else {
            None
        }
    }

    /// Update the zstep early.
    ///
    /// # Arguments
    ///
    /// * `delta` - The trust region radius.
    ///
    /// # Returns
    ///
    /// * `bool` - Whether the update was successful.
    fn early_update_zstep(&mut self, delta: f64) -> bool {
        /* Find weight solution. */
        let tau = match self.solve_curvature_quadratic(delta) {
            Some(t) => t,
            None => {
                println!("Failed to find weight solution.");
                return false;
            }
        };

        self.curr_zstep += &(tau * &self.curr_dstep);
        true
    }

    /// Iterate the trust region subproblem.
    ///
    /// # Arguments
    ///
    /// * `gradient` - The gradient vector.
    /// * `hessian` - The Hessian matrix.
    /// * `delta` - The trust region radius.
    ///
    /// # Returns
    ///
    /// * `bool` - Whether the iteration was successful.
    ///
    /// # Panics
    ///
    /// If the calculations are not finite or they are NaN.
    /// If the vectors are not sized correctly.
    pub fn iterate(
        &mut self,
        gradient: &Array1<f64>,
        hessian: &Array2<f64>,
        eps: f64,
        delta: f64,
    ) -> bool {
        /* Limit number of iterations. */
        if self.curr_iterations >= self.max_iterations {
            println!("Exceeded maximum iterations.");
            return false;
        }

        /* Assert that vectors are sized correctly. */
        assert_eq!(gradient.dim(), self.vector_size);
        assert_eq!(hessian.dim(), (self.vector_size, self.vector_size));

        /* Reset ztep. */
        self.curr_zstep.fill(0.0);

        /* Copy over gradient into rstep and dstep (the negative). */
        self.curr_rstep = gradient.clone();
        self.curr_dstep = gradient.iter().map(|&x| -x).collect();

        /* Stop early if the magnitude of the gradient is within tolerance. */
        if self.norm(self.curr_rstep.view()) < eps {
            return true;
        }

        for _i in 0..self.vector_size {
            /* Calculate the curvature. Matrix operation is d^T @ hessian @ d, returns a scalar. */
            let curvature = self.curr_dstep.t().dot(&hessian.dot(&self.curr_dstep));

            /* Find new zstep, wait if it's needed for next iteration. */
            let alpha = (self.curr_rstep.dot(&self.curr_rstep)) / curvature;
            let new_zstep = &self.curr_zstep + alpha * &self.curr_dstep;

            if self.norm(new_zstep.view()) >= delta {
                return self.early_update_zstep(delta);
            }

            /* Find new rstep, wait if it's needed for next iteration. */
            let new_rstep = &self.curr_rstep + alpha * &hessian.dot(&self.curr_dstep);
            if self.norm(new_rstep.view()) < eps {
                self.curr_zstep = new_zstep;
                return true;
            }

            /* Find new dstep and assign it back for next iteration. */
            let beta = (new_rstep.dot(&new_rstep)) / (self.curr_rstep.dot(&self.curr_rstep));

            self.curr_dstep = beta * &self.curr_dstep - &new_rstep;

            self.curr_zstep = new_zstep;
            self.curr_rstep = new_rstep;
        }

        self.curr_iterations += 1;
        true
    }

    /// Get the latest read-only result.
    ///
    /// # Returns
    ///
    /// * `ArrayView1<f64>` - The latest result read-only.
    pub fn get_result_readonly(&self) -> ArrayView1<f64> {
        self.curr_zstep.view()
    }

    /// Get the current number of iterations.
    ///
    /// # Returns
    ///
    /// * `usize` - The current number of iterations.
    pub fn get_curr_iterations(&self) -> usize {
        self.curr_iterations
    }

    /// Get the result.
    ///     
    /// # Returns
    ///
    /// * `Array1<f64>` - The copy of the result.
    pub fn get_result(&self) -> Array1<f64> {
        self.curr_zstep.clone()
    }
}
