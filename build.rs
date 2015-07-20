use std::env;
use std::process::Command;

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();

    Command::new("make").args(&["-C", "cephes-double", "lib"])
        .status().unwrap();
    Command::new("make").args(&["-C", "cephes-single", "lib"])
        .status().unwrap();

    println!("cargo:rustc-link-search=native={}", out_dir);
    println!("cargo:rustc-link-lib=static=md");
    println!("cargo:rustc-link-lib=static=mf");
}
