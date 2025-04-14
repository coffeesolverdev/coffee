use clap::{Arg, Command};
use coffee::extras::OptimizerArgs;
use coffee::run_coffee;

fn command() -> Command {
    Command::new("coffee_cli")
        .version("1.0")
        .author("UT Austin Senior Design Group FH12, 2024-2025")
        .about("CLI for COFFEE optimization")
        .arg(
            Arg::new("cfe")
                .help("The file path containing the input file for compositions and free energies.")
                .required(true)
                .index(1)
                .value_parser(|file: &str| {
                    let allowed_extensions = [".cfe", ".ocx", ".txt", ".csv", ".tsv"];
                    if !allowed_extensions.iter().any(|ext| file.ends_with(ext)) {
                        return Err("File must be a .cfe, .ocx, .txt, .csv, or .tsv file".to_string());
                    }
                    Ok(file.to_string())
                }),
        )
        .arg(
            Arg::new("con")
                .help("The file path containing the input file for concentrations.")
                .required(true)
                .index(2)
                .value_parser(|file: &str| {
                    let allowed_extensions = [".con", ".txt", ".csv", ".tsv"];
                    if !allowed_extensions.iter().any(|ext| file.ends_with(ext)) {
                        return Err("File must be a .con, .txt, .csv, or .tsv file".to_string());
                    }
                    Ok(file.to_string())
                }),
        )
        .arg(
            Arg::new("log")
                .short('l')
                .long("log")
                .help("The file path to output the log, including the results. If this is not provided, log will print to stdout by default.")
                .required(false)
                .value_parser(|file: &str| {
                    let allowed_extensions = [".txt", ".log"];
                    if !allowed_extensions.iter().any(|ext| file.ends_with(ext)) {
                        return Err("File must be a .txt or .log file".to_string());
                    }
                    Ok(file.to_string())
                }),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("The file path to output only the results. If this is not provided, results will not be saved to a file and does not affect log printing.")
                .required(false)
                .value_parser(|file: &str| {
                    let allowed_extensions = [".txt", ".log"];
                    if !allowed_extensions.iter().any(|ext| file.ends_with(ext)) {
                        return Err("File must be a .txt or .log file".to_string());
                    }
                    Ok(file.to_string())
                }),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .required(false)
                .action(clap::ArgAction::SetTrue)
                .help("Enable verbose output"),
        )
}

struct CoffeeArgs {
    pub desc: clap::ArgMatches,
}

impl CoffeeArgs {
    pub fn new() -> CoffeeArgs {
        let matches = command().get_matches();

        CoffeeArgs { desc: matches }
    }

    pub fn get_file(&self, arg: &str) -> Option<String> {
        self.desc.get_one::<String>(arg).cloned()
    }

    pub fn verbose(&self) -> bool {
        self.desc.get_flag("verbose")
    }
}

fn main() {
    let args = CoffeeArgs::new();

    let cfe_path = if let Some(path) = args.get_file("cfe") {
        path
    } else {
        eprintln!("CFE file path not provided.");
        return;
    };
    let con_path = if let Some(path) = args.get_file("con") {
        path
    } else {
        eprintln!("CON file path not provided.");
        return;
    };

    let log_path = args.get_file("log");
    let out_path = args.get_file("output");
    let verbose = args.verbose();

    let optimizer_args = OptimizerArgs {
        verbose,
        use_terminal: log_path.is_none(),
        ..OptimizerArgs::default()
    };

    // Call run_coffee with the file paths and get the result
    let coffee_result = run_coffee(
        &cfe_path,
        &con_path,
        log_path.as_deref(),
        out_path.as_deref(),
        &optimizer_args,
    );

    // Pass the result to print it
    match coffee_result {
        Ok(_) => return,
        Err(e) => format!("Error: {}", e),
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_required_args() {
        /* Most simple case, no extra args. */
        let mut matches = command().try_get_matches_from(vec![
            "coffee_cli",
            "~/coffee-internal/testcases/0/input.ocx",
            "~/coffee-internal/testcases/0/input.con",
        ]);
        assert!(matches.is_ok());

        let args = CoffeeArgs {
            desc: matches.unwrap(),
        };
        assert_eq!(
            args.get_file("cfe"),
            Some("~/coffee-internal/testcases/0/input.ocx".to_string())
        );
        assert_eq!(
            args.get_file("con"),
            Some("~/coffee-internal/testcases/0/input.con".to_string())
        );
        assert_eq!(args.get_file("log"), None);
        assert_eq!(args.get_file("output"), None);
        assert!(!args.verbose());

        /* Test 0 and 1 args, which should fail. */
        matches = command().try_get_matches_from(vec!["coffee_cli"]);
        assert!(matches.is_err());

        matches = command().try_get_matches_from(vec![
            "coffee_cli",
            "~/coffee-internal/testcases/0/input.ocx",
        ]);
        assert!(matches.is_err());

        matches = command().try_get_matches_from(vec![
            "coffee_cli",
            "~/coffee-internal/testcases/0/input.con",
        ]);
        assert!(matches.is_err());
    }

    #[test]
    fn test_optional_args() {
        /* Test optional args with valid inputs, long version. */
        let mut matches = command().try_get_matches_from(vec![
            "coffee_cli",
            "~/coffee-internal/testcases/0/input.ocx",
            "~/coffee-internal/testcases/0/input.con",
            "--log",
            "~/coffee-internal/testcases/0/log.txt",
            "--output",
            "~/coffee-internal/testcases/0/output.txt",
            "--verbose",
        ]);
        assert!(matches.is_ok());

        let args = CoffeeArgs {
            desc: matches.unwrap(),
        };
        assert_eq!(
            args.get_file("log"),
            Some("~/coffee-internal/testcases/0/log.txt".to_string())
        );
        assert_eq!(
            args.get_file("output"),
            Some("~/coffee-internal/testcases/0/output.txt".to_string())
        );
        assert!(args.verbose());

        /* Test optional args with valid inputs, short version. */
        matches = command().try_get_matches_from(vec![
            "coffee_cli",
            "~/coffee-internal/testcases/0/input.ocx",
            "~/coffee-internal/testcases/0/input.con",
            "-l",
            "~/coffee-internal/testcases/0/log.txt",
            "-o",
            "~/coffee-internal/testcases/0/output.txt",
            "-v",
        ]);
        assert!(matches.is_ok());

        let args = CoffeeArgs {
            desc: matches.unwrap(),
        };
        assert_eq!(
            args.get_file("log"),
            Some("~/coffee-internal/testcases/0/log.txt".to_string())
        );
        assert_eq!(
            args.get_file("output"),
            Some("~/coffee-internal/testcases/0/output.txt".to_string())
        );
        assert!(args.verbose());

        /* Test whether optional arguments are correctly parsed */
        matches = command().try_get_matches_from(vec![
            "coffee_cli",
            "~/coffee-internal/testcases/0/input.ocx",
            "~/coffee-internal/testcases/0/input.con",
            "--log",
            "~/coffee-internal/testcases/0/log.txt",
        ]);
        assert!(matches.is_ok());

        let args = CoffeeArgs {
            desc: matches.unwrap(),
        };
        assert_eq!(
            args.get_file("log"),
            Some("~/coffee-internal/testcases/0/log.txt".to_string())
        );
        assert_eq!(args.get_file("output"), None);
        assert!(!args.verbose());

        matches = command().try_get_matches_from(vec![
            "coffee_cli",
            "~/coffee-internal/testcases/0/input.ocx",
            "~/coffee-internal/testcases/0/input.con",
            "--output",
            "~/coffee-internal/testcases/0/out.txt",
        ]);
        assert!(matches.is_ok());

        let args = CoffeeArgs {
            desc: matches.unwrap(),
        };
        assert_eq!(args.get_file("log"), None);
        assert_eq!(
            args.get_file("output"),
            Some("~/coffee-internal/testcases/0/out.txt".to_string())
        );
        assert!(!args.verbose());

        matches = command().try_get_matches_from(vec![
            "coffee_cli",
            "~/coffee-internal/testcases/0/input.ocx",
            "~/coffee-internal/testcases/0/input.con",
            "--verbose",
        ]);
        assert!(matches.is_ok());
        assert!(CoffeeArgs {
            desc: matches.unwrap()
        }
        .verbose());
    }
}
