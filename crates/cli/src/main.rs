//! `nanocad` command line entry point.
#![forbid(unsafe_code)]

use std::io::Write;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut stdout = std::io::stdout();
    let mut stderr = std::io::stderr();
    let code = nanocad_cli::run(&args, &mut stdout, &mut stderr);
    let _ = stdout.flush();
    let _ = stderr.flush();
    ExitCode::from(code)
}
