pub mod extras;
pub mod fileparse;
pub mod format;
pub mod optimize;
pub mod steihaug;

use std::fs::File;
use std::io::Read;

use extras::{OptimizerArgs, OptimizerResults};
use fileparse::{parse_float, read_inputs_to_dataframe};
use format::results_message;
use ndarray::{Array1, Array2};
use optimize::Optimizer;

use core::result::Result;
use std::error::Error;
use std::io::Write;

use polars::prelude::DataType;

fn run_coffee_computation(
    cfe_bytes: &[u8],
    con_bytes: &[u8],
    optimizer_args: &OptimizerArgs,
) -> Result<OptimizerResults, Box<dyn Error>> {
    // Call fileparse to read the inputs and create a dataframe
    let table = match read_inputs_to_dataframe(cfe_bytes, con_bytes) {
        Ok(table) => table,
        Err(e) => {
            return Err(format!("Error reading files: {}", e).into());
        }
    };

    let polymer_rows = table.0.height();
    let polymer_cols = table.0.width();

    let mut polymer_data = Vec::<f64>::new();
    for col in table.0.get_columns() {
        let col_f64 = col.cast(&DataType::Float64)?;
        let series = col_f64.f64()?;
        polymer_data.extend(series.into_iter().map(|v| v.unwrap_or(0.0)));
    }

    let monomer_series_f64 = table.2.cast(&DataType::Float64)?;
    let monomers_vec = monomer_series_f64
        .f64()?
        .into_iter()
        .flatten()
        .collect::<Vec<f64>>();
    let polymer_energy_vec = parse_float(&table.1)?;

    // Create the optimizer
    let monomers = Array1::from_vec(monomers_vec);
    let mut polymers = match Array2::from_shape_vec((polymer_cols, polymer_rows), polymer_data) {
        Ok(polymers) => polymers,
        Err(e) => {
            return Err(format!("Failed to create polymers array: {}", e).into());
        }
    };
    polymers.swap_axes(0, 1);
    let polymers_energies = Array1::from_vec(polymer_energy_vec);

    let initial_delta = 1.0;
    let mut optimizer =
        match Optimizer::new(&monomers, &polymers, &polymers_energies, optimizer_args) {
            Ok(opt) => opt,
            Err(e) => {
                return Err(format!("Failed to create optimizer: {}", e).into());
            }
        };

    // Call the optimizer
    if let Err(e) = optimizer.optimize(initial_delta) {
        return Err(format!("Optimization failed: {}", e).into());
    }

    Ok(optimizer.get_results())
}

pub fn run_coffee_server(cfe_bytes: &[u8], con_bytes: &[u8]) -> Result<String, Box<dyn Error>> {
    let args = OptimizerArgs {
        use_terminal: true, // print to logs for websocket version
        verbose: true,
        ..Default::default()
    };
    let optimizer_results = match run_coffee_computation(cfe_bytes, con_bytes, &args) {
        Ok(optimizer_results) => optimizer_results,
        Err(e) => {
            eprintln!("Error during optimization: {}", e);
            return Err(e);
        }
    };

    Ok(results_message(&optimizer_results))
}

pub fn run_coffee(
    file_path_cfe: &str,
    file_path_con: &str,
    file_path_log: Option<&str>,
    file_path_out: Option<&str>,
    optimizer_args: &OptimizerArgs,
) -> Result<String, Box<dyn Error>> {
    // Read the file contents
    let mut file = File::open(file_path_cfe)?;
    let mut file_content_cfe = Vec::new();
    if let Err(e) = file.read_to_end(&mut file_content_cfe) {
        return Err(format!("Error reading monomer/polymer file: {}", e).into());
    }

    file = File::open(file_path_con)?;
    let mut file_content_con = Vec::new();
    if let Err(e) = file.read_to_end(&mut file_content_con) {
        return Err(format!("Error reading concentration file: {}", e).into());
    }

    let mut log_file = None;
    if let Some(log_path) = file_path_log {
        log_file = Some(File::create(log_path)?);
    }
    let mut out_file = None;
    if let Some(out_path) = file_path_out {
        out_file = Some(File::create(out_path)?);
    }

    let optimizer_results =
        match run_coffee_computation(&file_content_cfe, &file_content_con, optimizer_args) {
            Ok(optimizer_results) => optimizer_results,
            Err(e) => {
                eprintln!("Error during optimization: {}", e);
                return Err(e);
            }
        };

    let results_string = results_message(&optimizer_results);

    if let Some(ref mut log_file) = log_file {
        for message in &optimizer_results.log_messages {
            log_file.write_all(message.as_bytes())?;
        }
        log_file.write_all(results_string.as_bytes())?;
        log_file.flush()?;
    } else {
        println!("{}", results_string);
    }

    if let Some(ref mut out_file) = out_file {
        out_file.write_all(results_string.as_bytes())?;
        out_file.flush()?;
    };

    Ok(results_string)
}
