//! A design document: the parameters, the parts and the measurements.
//!
//! The document holds a bounded history of snapshots and a cursor. A commit
//! appends the current snapshot, undo and redo walk the cursor, and a new
//! commit after an undo drops the redo tail. The history is the undo and redo
//! stack, so one structure serves both.
//!
//! A snapshot is plain data, so it serializes to JSON and needs no generator to
//! load. The parameters are a name and a number list, because the generator
//! parameter set is an input to a part, not a stored record.

use serde::{Deserialize, Serialize};

use nanocad_meter::{Fidelity, MetricValue};
use nanocad_parts::PartSchema;

use crate::scene::SceneError;

/// The most snapshots the history keeps.
pub const MAX_HISTORY: usize = 64;

/// One named parameter value.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NamedValue {
    /// The parameter name.
    pub name: String,
    /// The parameter value, in the unit of the generator parameter.
    pub value: f64,
}

/// One recorded part.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PartRecord {
    /// The part name.
    pub name: String,
    /// The generator id that produced the part.
    pub generator_id: String,
    /// The number of atoms.
    pub atom_count: usize,
    /// The number of bonds.
    pub bond_count: usize,
    /// The material label.
    pub material: String,
}

impl PartRecord {
    /// Builds a record from a generated part schema.
    pub fn from_schema(schema: &PartSchema) -> Self {
        Self {
            name: schema.name.clone(),
            generator_id: schema.generator_id.clone(),
            atom_count: schema.atom_count,
            bond_count: schema.bond_count,
            material: schema.material.clone(),
        }
    }
}

/// One recorded measurement.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Measurement {
    /// The metric name.
    pub name: String,
    /// The measured value.
    pub value: f64,
    /// The unit symbol.
    pub unit: String,
    /// The fidelity tier, as a short label.
    pub fidelity: String,
    /// A short note on the method or the limit of the value.
    pub note: String,
}

/// Returns the short label of a fidelity tier.
fn fidelity_label(fidelity: Fidelity) -> &'static str {
    match fidelity {
        Fidelity::Geometric => "geometric",
        Fidelity::QuasiStatic => "quasi_static",
        Fidelity::Harmonic => "harmonic",
        Fidelity::Dynamics => "dynamics",
    }
}

impl Measurement {
    /// Builds a measurement from its fields.
    pub fn new(
        name: impl Into<String>,
        value: f64,
        unit: impl Into<String>,
        fidelity: impl Into<String>,
        note: impl Into<String>,
    ) -> Self {
        Self {
            name: name.into(),
            value,
            unit: unit.into(),
            fidelity: fidelity.into(),
            note: note.into(),
        }
    }

    /// Builds a measurement from a metric value.
    pub fn from_metric(metric: &MetricValue) -> Self {
        Self {
            name: metric.name.clone(),
            value: metric.value,
            unit: metric.unit.clone(),
            fidelity: fidelity_label(metric.fidelity).to_string(),
            note: metric.note.clone(),
        }
    }
}

/// One state of the design.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct DesignSnapshot {
    /// The parameters, in a stable order.
    pub parameters: Vec<NamedValue>,
    /// The parts, in a stable order.
    pub parts: Vec<PartRecord>,
    /// The measurements, in a stable order.
    pub measurements: Vec<Measurement>,
}

impl DesignSnapshot {
    /// An empty snapshot.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds or replaces one parameter value.
    pub fn with_parameter(mut self, name: impl Into<String>, value: f64) -> Self {
        self.set_parameter(name, value);
        self
    }

    /// Appends one part record.
    pub fn with_part(mut self, part: PartRecord) -> Self {
        self.parts.push(part);
        self
    }

    /// Appends one measurement.
    pub fn with_measurement(mut self, measurement: Measurement) -> Self {
        self.measurements.push(measurement);
        self
    }

    /// Adds or replaces one parameter value.
    pub fn set_parameter(&mut self, name: impl Into<String>, value: f64) {
        let name = name.into();
        if let Some(entry) = self.parameters.iter_mut().find(|entry| entry.name == name) {
            entry.value = value;
            return;
        }
        self.parameters.push(NamedValue { name, value });
    }

    /// Returns a parameter value by name, or `None`.
    pub fn parameter(&self, name: &str) -> Option<f64> {
        self.parameters
            .iter()
            .find(|entry| entry.name == name)
            .map(|entry| entry.value)
    }

    /// Returns a measurement by name, or `None`.
    pub fn measurement(&self, name: &str) -> Option<&Measurement> {
        self.measurements
            .iter()
            .find(|measurement| measurement.name == name)
    }
}

/// A design document with a bounded undo and redo history.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DesignDocument {
    snapshots: Vec<DesignSnapshot>,
    cursor: usize,
}

impl DesignDocument {
    /// Creates a document with one snapshot.
    pub fn new(initial: DesignSnapshot) -> Self {
        Self {
            snapshots: vec![initial],
            cursor: 0,
        }
    }

    /// Returns the snapshot at the cursor.
    pub fn current(&self) -> &DesignSnapshot {
        &self.snapshots[self.cursor]
    }

    /// Returns the number of snapshots.
    pub fn depth(&self) -> usize {
        self.snapshots.len()
    }

    /// Returns the cursor position.
    pub fn cursor(&self) -> usize {
        self.cursor
    }

    /// Returns the snapshot history, oldest first.
    pub fn snapshots(&self) -> &[DesignSnapshot] {
        &self.snapshots
    }

    /// Reports whether an undo is available.
    pub fn can_undo(&self) -> bool {
        self.cursor > 0
    }

    /// Reports whether a redo is available.
    pub fn can_redo(&self) -> bool {
        self.cursor + 1 < self.snapshots.len()
    }

    /// Records a new snapshot.
    ///
    /// Returns `false` and does nothing when the snapshot equals the current
    /// one. Otherwise it drops the redo tail, appends the snapshot, and moves
    /// the cursor. When the history is full it drops the oldest snapshot.
    pub fn commit(&mut self, snapshot: DesignSnapshot) -> bool {
        if snapshot == *self.current() {
            return false;
        }
        self.snapshots.truncate(self.cursor + 1);
        self.snapshots.push(snapshot);
        if self.snapshots.len() > MAX_HISTORY {
            self.snapshots.remove(0);
        }
        self.cursor = self.snapshots.len() - 1;
        true
    }

    /// Moves the cursor back one snapshot. Returns `false` when it cannot.
    pub fn undo(&mut self) -> bool {
        if !self.can_undo() {
            return false;
        }
        self.cursor -= 1;
        true
    }

    /// Moves the cursor forward one snapshot. Returns `false` when it cannot.
    pub fn redo(&mut self) -> bool {
        if !self.can_redo() {
            return false;
        }
        self.cursor += 1;
        true
    }

    /// Writes the document as JSON.
    pub fn to_json(&self) -> Result<String, SceneError> {
        serde_json::to_string(self).map_err(|error| SceneError::Json(error.to_string()))
    }

    /// Writes the document as indented JSON.
    pub fn to_json_pretty(&self) -> Result<String, SceneError> {
        serde_json::to_string_pretty(self).map_err(|error| SceneError::Json(error.to_string()))
    }

    /// Reads a document from JSON.
    ///
    /// A document with no snapshot or with a cursor outside the history is
    /// refused.
    pub fn from_json(text: &str) -> Result<Self, SceneError> {
        let document: Self =
            serde_json::from_str(text).map_err(|error| SceneError::Json(error.to_string()))?;
        if document.snapshots.is_empty() {
            return Err(SceneError::Json("the document has no snapshot".to_string()));
        }
        if document.cursor >= document.snapshots.len() {
            return Err(SceneError::Json(format!(
                "the cursor {} is outside the {} snapshots",
                document.cursor,
                document.snapshots.len()
            )));
        }
        Ok(document)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn snapshot(value: f64) -> DesignSnapshot {
        DesignSnapshot::new().with_parameter("sun_teeth", value)
    }

    #[test]
    fn a_new_document_has_one_snapshot() {
        let document = DesignDocument::new(snapshot(12.0));
        assert_eq!(document.depth(), 1);
        assert_eq!(document.cursor(), 0);
        assert!(!document.can_undo());
        assert!(!document.can_redo());
        assert_eq!(document.current().parameter("sun_teeth"), Some(12.0));
    }

    #[test]
    fn a_commit_adds_a_snapshot_and_moves_the_cursor() {
        let mut document = DesignDocument::new(snapshot(12.0));
        assert!(document.commit(snapshot(13.0)));
        assert_eq!(document.depth(), 2);
        assert_eq!(document.cursor(), 1);
        assert!(document.can_undo());
    }

    #[test]
    fn an_identical_commit_is_refused() {
        let mut document = DesignDocument::new(snapshot(12.0));
        assert!(!document.commit(snapshot(12.0)));
        assert_eq!(document.depth(), 1);
    }

    #[test]
    fn undo_and_redo_walk_the_history() {
        let mut document = DesignDocument::new(snapshot(12.0));
        document.commit(snapshot(13.0));
        document.commit(snapshot(14.0));
        assert!(document.undo());
        assert_eq!(document.current().parameter("sun_teeth"), Some(13.0));
        assert!(document.undo());
        assert_eq!(document.current().parameter("sun_teeth"), Some(12.0));
        assert!(!document.undo());
        assert!(document.redo());
        assert_eq!(document.current().parameter("sun_teeth"), Some(13.0));
        assert!(document.redo());
        assert!(!document.redo());
    }

    #[test]
    fn a_new_commit_drops_the_redo_tail() {
        let mut document = DesignDocument::new(snapshot(12.0));
        document.commit(snapshot(13.0));
        document.commit(snapshot(14.0));
        document.undo();
        assert!(document.can_redo());
        assert!(document.commit(snapshot(15.0)));
        assert!(!document.can_redo());
        assert_eq!(document.depth(), 3);
        assert_eq!(document.current().parameter("sun_teeth"), Some(15.0));
    }

    #[test]
    fn the_history_is_bounded() {
        let mut document = DesignDocument::new(snapshot(0.0));
        for step in 1..=(MAX_HISTORY + 10) {
            document.commit(snapshot(step as f64));
        }
        assert_eq!(document.depth(), MAX_HISTORY);
        assert_eq!(document.cursor(), MAX_HISTORY - 1);
        assert!(!document.can_redo());
    }

    #[test]
    fn a_design_document_round_trips_through_json() {
        let mut document = DesignDocument::new(snapshot(12.0));
        document.commit(snapshot(13.0));
        let text = document.to_json().expect("write");
        let read = DesignDocument::from_json(&text).expect("read");
        assert_eq!(read, document);
    }

    #[test]
    fn a_bad_document_is_rejected() {
        assert!(DesignDocument::from_json("not json").is_err());
        let empty = "{\"snapshots\":[],\"cursor\":0}";
        assert!(DesignDocument::from_json(empty).is_err());
        let bad_cursor =
            "{\"snapshots\":[{\"parameters\":[],\"parts\":[],\"measurements\":[]}],\"cursor\":3}";
        assert!(DesignDocument::from_json(bad_cursor).is_err());
    }

    #[test]
    fn a_measurement_reads_a_metric() {
        let metric = MetricValue {
            name: "clearance".to_string(),
            value: 2.5e-10,
            unit: "m".to_string(),
            fidelity: Fidelity::Geometric,
            note: "measured".to_string(),
        };
        let measurement = Measurement::from_metric(&metric);
        assert_eq!(measurement.name, "clearance");
        assert_eq!(measurement.fidelity, "geometric");
        assert_eq!(measurement.value, 2.5e-10);
    }

    #[test]
    fn a_part_record_reads_a_schema() {
        let schema = PartSchema {
            name: "plate".to_string(),
            generator_id: "plate".to_string(),
            parameters: Default::default(),
            atom_count: 10,
            bond_count: 8,
            material: "diamondoid".to_string(),
        };
        let record = PartRecord::from_schema(&schema);
        assert_eq!(record.name, "plate");
        assert_eq!(record.atom_count, 10);
        assert_eq!(record.material, "diamondoid");
    }
}
