//! Prints the catalog of scenes that the viewer can show.
//!
//! The output is one JSON object with a `scenes` array. Each entry has an id,
//! a name, a kind (`machine` or `part`) and a category.

fn main() {
    match serde_json::to_string(&nanocad_jigs::scene_catalog()) {
        Ok(list) => println!("{{\"scenes\":{list}}}"),
        Err(error) => {
            eprintln!("the catalog was not written: {error}");
            std::process::exit(1);
        }
    }
}
