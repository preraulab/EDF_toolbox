// build.rs -- captures the current git SHA + working-tree dirty flag at
// compile time and exposes them as env vars (`GIT_SHA`, `GIT_BUILD_DATE`)
// so main.rs can stitch them into clap's --version string.
//
// Falls back to "unknown" when git isn't available or the source isn't a
// checkout (e.g. building from a release tarball).

use std::process::Command;

fn main() {
    let sha = run_git(&["rev-parse", "--short=10", "HEAD"]).unwrap_or_else(|| "unknown".into());

    // "" if clean working tree, "-dirty" if there are uncommitted changes.
    let dirty = match Command::new("git").args(["status", "--porcelain"]).output() {
        Ok(out) if out.status.success() && !out.stdout.is_empty() => "-dirty",
        _ => "",
    };

    let date = run_git(&["log", "-1", "--format=%cs"]).unwrap_or_else(|| "unknown".into());

    println!("cargo:rustc-env=GIT_SHA={sha}{dirty}");
    println!("cargo:rustc-env=GIT_BUILD_DATE={date}");

    // Re-run if HEAD or any ref changes so cached builds don't go stale.
    println!("cargo:rerun-if-changed=../.git/HEAD");
    println!("cargo:rerun-if-changed=../.git/index");
    println!("cargo:rerun-if-changed=../.git/refs/heads");
}

fn run_git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).trim().to_string())
}
