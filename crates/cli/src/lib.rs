//! The `nanocad` command line.
//!
//! This crate holds the argument handling and the command dispatch. The
//! `nanocad` binary is a thin wrapper around [`run`]. Keeping the logic in a
//! library lets the tests call it without a subprocess.
#![forbid(unsafe_code)]

use std::collections::BTreeMap;
use std::fmt;
use std::fs;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

use nanocad_format::{
    export_mmp, export_xyz, import_mmp, import_xyz, read_ncz, write_ncz, FormatError,
};
use nanocad_model::{Document, Element};
use thiserror::Error;

/// The help text.
pub const USAGE: &str = "\
nanocad <command> [arguments]

Commands:
  info <file>            Print a summary of a document.
  convert <in> <out>     Read <in> and write <out>.

Options:
  -h, --help             Print this help.
  -V, --version          Print the version.

The file format comes from the extension: .ncz, .mmp, or .xyz
";

/// Returns the crate version string.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// A file format the command line can read and write.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FileFormat {
    /// The native NCZ archive.
    Ncz,
    /// The NanoEngineer MMP text.
    Mmp,
    /// The plain XYZ text.
    Xyz,
}

impl FileFormat {
    /// Returns the display name of the format.
    pub const fn label(self) -> &'static str {
        match self {
            Self::Ncz => "ncz",
            Self::Mmp => "mmp",
            Self::Xyz => "xyz",
        }
    }

    /// Detects a format from a path extension. The match ignores case.
    pub fn from_path(path: impl AsRef<Path>) -> Option<Self> {
        let extension = path.as_ref().extension()?.to_str()?.to_ascii_lowercase();
        match extension.as_str() {
            "ncz" => Some(Self::Ncz),
            "mmp" => Some(Self::Mmp),
            "xyz" => Some(Self::Xyz),
            _ => None,
        }
    }
}

/// An error from the command line.
///
/// The binary prints one of these to stderr and exits with a nonzero code. It
/// never panics on user input.
#[derive(Debug, Error)]
pub enum CliError {
    /// The command or its arguments are not valid.
    #[error("{0}")]
    Usage(String),
    /// The file extension is not one of the supported formats.
    #[error("unrecognized file extension in {path:?}; expected .ncz, .mmp, or .xyz")]
    UnknownExtension {
        /// The offending path.
        path: String,
    },
    /// A file could not be read.
    #[error("cannot read {path}: {source}")]
    Read {
        /// The offending path.
        path: String,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
    /// A file could not be written.
    #[error("cannot write {path}: {source}")]
    Write {
        /// The offending path.
        path: String,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
    /// A file could not be parsed.
    #[error("cannot parse {path}: {source}")]
    Parse {
        /// The offending path.
        path: String,
        /// The underlying error.
        #[source]
        source: FormatError,
    },
    /// A document could not be encoded.
    #[error("cannot encode {path}: {source}")]
    Encode {
        /// The offending path.
        path: String,
        /// The underlying error.
        #[source]
        source: FormatError,
    },
    /// An error on the program's own output stream.
    #[error("cannot write output: {0}")]
    Output(#[from] std::io::Error),
}

impl CliError {
    /// Reports whether the error is a usage error.
    pub const fn is_usage(&self) -> bool {
        matches!(self, Self::Usage(_))
    }

    /// Returns the process exit code for this error.
    pub const fn exit_code(&self) -> u8 {
        if self.is_usage() {
            2
        } else {
            1
        }
    }
}

/// Runs the command line and returns the process exit code.
///
/// `args` excludes the program name. Normal output goes to `out`; errors and
/// the usage message go to `err`.
pub fn run(args: &[String], out: &mut impl Write, err: &mut impl Write) -> u8 {
    let result = match args.first().map(String::as_str) {
        None => Err(CliError::Usage("missing command".to_owned())),
        Some("-h") | Some("--help") => out.write_all(USAGE.as_bytes()).map_err(CliError::Output),
        Some("-V") | Some("--version") => {
            writeln!(out, "nanocad {}", version()).map_err(CliError::Output)
        }
        Some("info") => info(&args[1..], out),
        Some("convert") => convert(&args[1..]),
        Some(other) => Err(CliError::Usage(format!("unknown command {other:?}"))),
    };

    match result {
        Ok(()) => 0,
        Err(error) => {
            let _ = writeln!(err, "error: {error}");
            if error.is_usage() {
                let _ = write!(err, "{USAGE}");
            }
            error.exit_code()
        }
    }
}

fn info(args: &[String], out: &mut impl Write) -> Result<(), CliError> {
    let path = require_path(args, "info")?;
    let format = detect(&path)?;
    let document = read_document(&path, format)?;
    out.write_all(summarize(&document, format).as_bytes())
        .map_err(CliError::Output)
}

fn convert(args: &[String]) -> Result<(), CliError> {
    if args.len() != 2 {
        return Err(CliError::Usage(
            "convert needs an input file and an output file".to_owned(),
        ));
    }
    let input = &args[0];
    let output = &args[1];
    let input_format = detect(input)?;
    let output_format = detect(output)?;
    let document = read_document(input, input_format)?;
    write_document(output, output_format, &document)
}

fn require_path(args: &[String], command: &str) -> Result<String, CliError> {
    if args.len() != 1 {
        return Err(CliError::Usage(format!("{command} needs exactly one file")));
    }
    Ok(args[0].clone())
}

fn detect(path: &str) -> Result<FileFormat, CliError> {
    FileFormat::from_path(path).ok_or_else(|| CliError::UnknownExtension {
        path: path.to_owned(),
    })
}

fn read_document(path: &str, format: FileFormat) -> Result<Document, CliError> {
    match format {
        FileFormat::Ncz => {
            let file = fs::File::open(path).map_err(|source| CliError::Read {
                path: path.to_owned(),
                source,
            })?;
            read_ncz(BufReader::new(file)).map_err(|source| CliError::Parse {
                path: path.to_owned(),
                source,
            })
        }
        FileFormat::Mmp => {
            let text = fs::read_to_string(path).map_err(|source| CliError::Read {
                path: path.to_owned(),
                source,
            })?;
            import_mmp(&text).map_err(|source| CliError::Parse {
                path: path.to_owned(),
                source,
            })
        }
        FileFormat::Xyz => {
            let text = fs::read_to_string(path).map_err(|source| CliError::Read {
                path: path.to_owned(),
                source,
            })?;
            import_xyz(&text).map_err(|source| CliError::Parse {
                path: path.to_owned(),
                source,
            })
        }
    }
}

fn write_document(path: &str, format: FileFormat, document: &Document) -> Result<(), CliError> {
    match format {
        FileFormat::Ncz => {
            let file = fs::File::create(path).map_err(|source| CliError::Write {
                path: path.to_owned(),
                source,
            })?;
            let mut writer =
                write_ncz(BufWriter::new(file), document).map_err(|source| CliError::Encode {
                    path: path.to_owned(),
                    source,
                })?;
            writer.flush().map_err(|source| CliError::Write {
                path: path.to_owned(),
                source,
            })
        }
        FileFormat::Mmp => {
            let text = export_mmp(document).map_err(|source| CliError::Encode {
                path: path.to_owned(),
                source,
            })?;
            fs::write(path, text).map_err(|source| CliError::Write {
                path: path.to_owned(),
                source,
            })
        }
        FileFormat::Xyz => {
            let text = export_xyz(document).map_err(|source| CliError::Encode {
                path: path.to_owned(),
                source,
            })?;
            fs::write(path, text).map_err(|source| CliError::Write {
                path: path.to_owned(),
                source,
            })
        }
    }
}

/// Builds the `info` summary for a document.
pub fn summarize(document: &Document, format: FileFormat) -> String {
    let mut out = String::new();
    let name = if document.name.is_empty() {
        "(unnamed)"
    } else {
        document.name.as_str()
    };
    out.push_str(&format!("format: {}\n", format.label()));
    out.push_str(&format!("name: {name}\n"));
    out.push_str(&format!("parts: {}\n", document.part_count()));

    let atom_total: usize = document.parts.iter().map(|part| part.atom_count()).sum();
    let bond_total: usize = document.parts.iter().map(|part| part.bond_count()).sum();
    out.push_str(&format!("atoms: {atom_total}\n"));
    out.push_str(&format!("bonds: {bond_total}\n"));

    let mut histogram: BTreeMap<Element, usize> = BTreeMap::new();
    for part in &document.parts {
        for atom in part.topology.atoms() {
            *histogram.entry(atom.element).or_insert(0) += 1;
        }
    }
    out.push_str("elements:\n");
    if histogram.is_empty() {
        out.push_str("  (none)\n");
    } else {
        for (element, count) in &histogram {
            out.push_str(&format!("  {element}: {count}\n"));
        }
    }

    out.push_str("metadata:\n");
    if document.metadata.is_empty() {
        out.push_str("  (none)\n");
    } else {
        for (key, value) in &document.metadata {
            out.push_str(&format!("  {key}: {value}\n"));
        }
    }
    out
}

impl fmt::Display for FileFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.label())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_format_comes_from_a_case_insensitive_extension() {
        assert_eq!(FileFormat::from_path("part.ncz"), Some(FileFormat::Ncz));
        assert_eq!(FileFormat::from_path("PART.MMP"), Some(FileFormat::Mmp));
        assert_eq!(FileFormat::from_path("./a/b.xyz"), Some(FileFormat::Xyz));
        assert_eq!(FileFormat::from_path("part.pdb"), None);
        assert_eq!(FileFormat::from_path("noext"), None);
    }

    #[test]
    fn the_version_is_not_empty() {
        assert!(!version().is_empty());
    }

    #[test]
    fn a_usage_error_exits_with_two() {
        let mut out = Vec::new();
        let mut err = Vec::new();
        let code = run(&[], &mut out, &mut err);
        assert_eq!(code, 2);
        assert!(String::from_utf8_lossy(&err).contains("missing command"));
    }

    #[test]
    fn a_bad_extension_exits_with_one() {
        let mut out = Vec::new();
        let mut err = Vec::new();
        let args = vec!["convert".to_owned(), "a.pdb".to_owned(), "b.ncz".to_owned()];
        let code = run(&args, &mut out, &mut err);
        assert_eq!(code, 1);
        assert!(String::from_utf8_lossy(&err).contains("unrecognized file extension"));
    }

    #[test]
    fn the_summary_reports_the_counts_and_the_histogram() {
        let document = import_mmp(
            "mol (p) def\natom 1 (6) (0, 0, 0) def\natom 2 (1) (1000, 0, 0) def\nbond1 1\n",
        )
        .expect("import");
        let summary = summarize(&document, FileFormat::Mmp);
        assert!(summary.contains("format: mmp"));
        assert!(summary.contains("parts: 1"));
        assert!(summary.contains("atoms: 2"));
        assert!(summary.contains("bonds: 1"));
        assert!(summary.contains("C: 1"));
        assert!(summary.contains("H: 1"));
    }
}
