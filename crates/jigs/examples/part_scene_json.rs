//! Prints one library part as a scene, for the viewer.
//!
//! Usage: `part_scene_json <part_id>`. The output is one JSON object with a
//! `scene` key in the shared scene schema. An unknown id exits with an error.

fn main() {
    let Some(id) = std::env::args().nth(1) else {
        eprintln!("usage: part_scene_json <part_id>");
        std::process::exit(2);
    };
    let scene = match nanocad_jigs::build_part_scene(&id) {
        Ok(scene) => scene,
        Err(error) => {
            eprintln!("the part scene failed: {error}");
            std::process::exit(1);
        }
    };
    let json = match nanocad_jigs::scene_to_json(&scene) {
        Ok(json) => json,
        Err(error) => {
            eprintln!("the part scene was not written: {error}");
            std::process::exit(1);
        }
    };
    println!("{{\"scene\":{json}}}");
}
