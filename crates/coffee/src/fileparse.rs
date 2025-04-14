use std::error::Error;
use std::io::Cursor;
use std::result::Result;

use polars::prelude::{CsvReader, DataFrame, DataType, PolarsError, SerReader, Series};

type ParsedData = (DataFrame, Series, Series);

pub fn read_inputs_to_dataframe(
    file_content_cfe: &[u8],
    file_content_con: &[u8],
) -> Result<ParsedData, Box<dyn Error>> {
    let cfe_delimiter = detect_delimiter(file_content_cfe)?;
    let cfe_cursor = Cursor::new(file_content_cfe);

    let mut cfe_df = CsvReader::new(cfe_cursor)
        .has_header(false)
        .with_delimiter(cfe_delimiter)
        .finish()?;

    let num_columns = cfe_df.width();
    let num_rows = cfe_df.height();

    let mut is_nupack = true;
    let sample_size = 20.min(num_rows);

    // Check first min(rows, sample_size) entries for NUPACK formatting
    for index in 0..sample_size {
        let row = cfe_df.get(index).ok_or("Failed to get row")?;
        let value1: i64 = row.first().unwrap().extract::<i64>().unwrap_or(0);
        let value2: i64 = row.get(1).unwrap().extract::<i64>().unwrap_or(0);

        // Any failed comparison sets is_nupack to false
        is_nupack &= value1 == (index + 1) as i64;
        is_nupack &= value2 == 1;

        if !is_nupack {
            break;
        }
    }

    let float_col = cfe_df
        .select_at_idx(num_columns - 1)
        .ok_or("Failed to find float_col column")?
        .clone();
    cfe_df = cfe_df.drop(cfe_df.get_column_names()[num_columns - 1])?;

    if is_nupack {
        cfe_df = cfe_df.drop(cfe_df.get_column_names()[0])?;
        cfe_df = cfe_df.drop(cfe_df.get_column_names()[0])?;
    }

    // Parse .con file
    let con_cursor = Cursor::new(file_content_con);
    let con_df = CsvReader::new(con_cursor).has_header(false).finish()?;

    if con_df.width() != 1 {
        return Err(PolarsError::ComputeError("Invalid .con file".into()).into());
    }

    let con_vector = con_df
        .select_at_idx(0)
        .ok_or("Failed to select column")?
        .clone();

    Ok((cfe_df, float_col, con_vector))
}

pub fn parse_float(series: &Series) -> Result<Vec<f64>, Box<dyn Error>> {
    /* Normal format, what Rust can natively handle. */
    if series.dtype() == &DataType::Float64 {
        let float_series = match series.f64() {
            Ok(s) => s,
            Err(e) => {
                return Err(format!("Error converting series to f64: {}", e).into());
            }
        };
        Ok(float_series.into_iter().flatten().collect())
    }
    /* Integer format, in case all values were given as ints. */
    else if series.dtype() == &DataType::Int64 {
        let int_series = match series.i64() {
            Ok(s) => s,
            Err(e) => {
                return Err(format!("Error converting series to i64: {}", e).into());
            }
        };
        return Ok(int_series
            .into_iter()
            .map(|v| v.unwrap_or(0) as f64)
            .collect());
    }
    /* Default format as str, which means series has combination on weird number formats that Rust can't handle natively. */
    else if series.dtype() == &DataType::Utf8 {
        /* Convert to utf8, or strings. */
        let utf8_series = match series.utf8() {
            Ok(s) => s,
            Err(e) => {
                return Err(format!("Error converting series to utf8: {}", e).into());
            }
        };
        let mut float_values = Vec::new();
        for (index, value) in utf8_series.into_iter().enumerate() {
            /* Get string from wrap. */
            let value_str = match value {
                Some(v) => v,
                None => return Err(format!("Error parsing number at index: {} ", index).into()),
            };

            /* Check whether it's in exponent format with the e. */
            if value_str.contains('e') {
                let (base, exponent) = value_str.split_once('e').unwrap();
                if let Ok(base_value) = base.trim().parse::<f64>() {
                    if let Ok(exp_value) = exponent.trim().parse::<i32>() {
                        float_values.push(base_value * 10f64.powi(exp_value));
                    } else {
                        return Err(format!(
                            "Error parsing number: {} at index {}.",
                            value_str, index
                        )
                        .into());
                    }
                }
            } else {
                /* Else it's in decimal form. */
                let (base_value, decimal_value) =
                    value_str.split_once('.').unwrap_or((value_str, ""));

                if let Ok(base_value_f64) = base_value.trim().parse::<f64>() {
                    if decimal_value.is_empty() {
                        /* No decimal value, just push base. */
                        float_values.push(base_value_f64);
                    } else {
                        /* Do calculation for number = base + decimal / 10^len.  */
                        if let Ok(decimal_value_f64) = decimal_value.trim().parse::<f64>() {
                            float_values.push(
                                base_value_f64
                                    + decimal_value_f64 / 10f64.powi(decimal_value.len() as i32),
                            );
                        } else {
                            return Err(format!(
                                "Error parsing decimal part: {} at index {}.",
                                decimal_value, index
                            )
                            .into());
                        }
                    }
                } else {
                    return Err(
                        format!("Error parsing number: {} at index {}.", value_str, index).into(),
                    );
                }
            }
        }
        return Ok(float_values);
    } else {
        return Err("Unsupported number types given in files for float conversion".into());
    }
}

fn detect_delimiter(file_content: &[u8]) -> std::result::Result<u8, Box<dyn Error>> {
    file_content
        .iter()
        .copied()
        .find(|&b| !b.is_ascii_alphanumeric())
        .ok_or_else(|| "Failed to detect delimiter".into())
}
