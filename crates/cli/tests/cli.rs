//! Integration tests for the `nanocad` command line.
//!
//! Most tests call [`nanocad_cli::run`] with in-memory output buffers. One test
//! runs the built binary to prove the wiring.

use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

const SAMPLE_MMP: &str = "\
mmpformat 050920 required; 051103 preferred
mol (Synthetic) def
atom 1 (6) (0, 0, 0) def
atom 2 (8) (1200, 0, 0) def
info atom atomtype = sp2
bond2 1
atom 3 (6) (2400, 1000, 0) def
bond1 2
atom 4 (1) (3600, 1000, 0) def
bond1 3
egroup (Synthetic)
end molecular machine part Synthetic
";

static COUNTER: AtomicU32 = AtomicU32::new(0);

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new() -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock after epoch")
            .as_nanos();
        let unique = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "nanocad-cli-{}-{nanos}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&path).expect("create temp dir");
        Self { path }
    }

    fn join(&self, name: &str) -> PathBuf {
        self.path.join(name)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

fn run(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|value| (*value).to_owned()).collect();
    let mut out = Vec::new();
    let mut err = Vec::new();
    let code = nanocad_cli::run(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).expect("stdout is utf-8"),
        String::from_utf8(err).expect("stderr is utf-8"),
    )
}

#[test]
fn convert_an_mmp_to_ncz_and_then_to_xyz() {
    let dir = TempDir::new();
    let mmp = dir.join("part.mmp");
    let ncz = dir.join("part.ncz");
    let xyz = dir.join("part.xyz");
    fs::write(&mmp, SAMPLE_MMP).expect("write mmp");

    let (code, _, err) = run(&[
        "convert",
        mmp.to_str().expect("path"),
        ncz.to_str().expect("path"),
    ]);
    assert_eq!(code, 0, "mmp to ncz failed: {err}");
    assert!(ncz.exists(), "ncz was not written");

    let (code, info, err) = run(&["info", ncz.to_str().expect("path")]);
    assert_eq!(code, 0, "info failed: {err}");
    assert!(info.contains("format: ncz"), "info was {info:?}");
    assert!(info.contains("atoms: 4"), "info was {info:?}");
    assert!(info.contains("bonds: 3"), "info was {info:?}");
    assert!(info.contains("C: 2"), "info was {info:?}");
    assert!(info.contains("O: 1"), "info was {info:?}");

    let (code, _, err) = run(&[
        "convert",
        ncz.to_str().expect("path"),
        xyz.to_str().expect("path"),
    ]);
    assert_eq!(code, 0, "ncz to xyz failed: {err}");
    let xyz_text = fs::read_to_string(&xyz).expect("read xyz");
    assert!(xyz_text.starts_with("4\n"), "xyz was {xyz_text:?}");
    assert!(xyz_text.contains("C "), "xyz was {xyz_text:?}");
    assert!(xyz_text.contains("O "), "xyz was {xyz_text:?}");

    let (code, info, _) = run(&["info", xyz.to_str().expect("path")]);
    assert_eq!(code, 0);
    assert!(info.contains("format: xyz"), "info was {info:?}");
    assert!(info.contains("atoms: 4"), "info was {info:?}");
}

#[test]
fn convert_an_mmp_to_xyz_directly() {
    let dir = TempDir::new();
    let mmp = dir.join("part.mmp");
    let xyz = dir.join("part.xyz");
    fs::write(&mmp, SAMPLE_MMP).expect("write mmp");

    let (code, _, err) = run(&[
        "convert",
        mmp.to_str().expect("path"),
        xyz.to_str().expect("path"),
    ]);
    assert_eq!(code, 0, "convert failed: {err}");
    let xyz_text = fs::read_to_string(&xyz).expect("read xyz");
    assert_eq!(xyz_text.lines().next(), Some("4"));
}

#[test]
fn a_bad_input_extension_is_a_clear_error() {
    let dir = TempDir::new();
    let input = dir.join("part.pdb");
    let output = dir.join("part.ncz");
    fs::write(&input, "not a real pdb").expect("write input");

    let (code, _, err) = run(&[
        "convert",
        input.to_str().expect("path"),
        output.to_str().expect("path"),
    ]);
    assert_ne!(code, 0);
    assert!(
        err.contains("unrecognized file extension"),
        "stderr was {err:?}"
    );
    assert!(!output.exists(), "no output must be written");
}

#[test]
fn a_missing_input_file_is_a_clear_error() {
    let dir = TempDir::new();
    let missing = dir.join("missing.ncz");
    let (code, _, err) = run(&["info", missing.to_str().expect("path")]);
    assert_ne!(code, 0);
    assert!(err.contains("cannot read"), "stderr was {err:?}");
}

#[test]
fn a_malformed_mmp_is_a_clear_error() {
    let dir = TempDir::new();
    let input = dir.join("bad.mmp");
    fs::write(&input, "mol (bad) def\natom 1 (6) (oops, 0, 0) def\n").expect("write");
    let (code, _, err) = run(&["info", input.to_str().expect("path")]);
    assert_ne!(code, 0);
    assert!(err.contains("cannot parse"), "stderr was {err:?}");
}

#[test]
fn a_missing_command_is_a_usage_error() {
    let (code, _, err) = run(&[]);
    assert_eq!(code, 2);
    assert!(err.contains("missing command"), "stderr was {err:?}");
    assert!(err.contains("nanocad <command>"), "stderr was {err:?}");
}

#[test]
fn help_and_version_report_success() {
    let (code, out, _) = run(&["--help"]);
    assert_eq!(code, 0);
    assert!(out.contains("nanocad <command>"), "stdout was {out:?}");

    let (code, out, _) = run(&["--version"]);
    assert_eq!(code, 0);
    assert!(out.starts_with("nanocad "), "stdout was {out:?}");
}

#[test]
fn the_built_binary_prints_the_version() {
    let output = Command::new(env!("CARGO_BIN_EXE_nanocad"))
        .arg("--version")
        .output()
        .expect("run the nanocad binary");
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf-8");
    assert!(stdout.starts_with("nanocad "), "stdout was {stdout:?}");
}
