//! Axis-aligned and analytic solids, in SI metres.
//!
//! A [`Solid`] is a closed region. The module holds primitives and the boolean
//! operations that combine them. Every quantity carries its unit in the name.
//! A [`Profile`] has no finite z bounds, because the outline is extruded
//! through all z. The caller must set the z limits when it combines a profile.

use std::f64::consts::PI;
use std::fmt;

/// An axis-aligned bounding box in metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bounds {
    /// The low corner, component by component, in metres.
    pub min_m: [f64; 3],
    /// The high corner, component by component, in metres.
    pub max_m: [f64; 3],
}

impl Bounds {
    /// Returns a box with no points. The low corner is `+inf` and the high
    /// corner is `-inf`, so [`Bounds::include_point_m`] grows it from nothing.
    pub fn empty() -> Self {
        Bounds {
            min_m: [f64::INFINITY; 3],
            max_m: [f64::NEG_INFINITY; 3],
        }
    }

    /// Returns the smallest box that holds one point.
    pub fn from_point_m(point: [f64; 3]) -> Self {
        Bounds {
            min_m: point,
            max_m: point,
        }
    }

    /// Grows the box so that it holds one more point.
    pub fn include_point_m(&mut self, point: [f64; 3]) {
        self.min_m = std::array::from_fn(|axis| self.min_m[axis].min(point[axis]));
        self.max_m = std::array::from_fn(|axis| self.max_m[axis].max(point[axis]));
    }

    /// Returns the smallest box that holds both `self` and `other`.
    ///
    /// An empty operand has no effect. The infinities of [`Bounds::empty`] make
    /// the componentwise minimum and maximum return the other operand.
    pub fn union(&self, other: &Bounds) -> Bounds {
        Bounds {
            min_m: std::array::from_fn(|axis| self.min_m[axis].min(other.min_m[axis])),
            max_m: std::array::from_fn(|axis| self.max_m[axis].max(other.max_m[axis])),
        }
    }

    /// Returns the box grown by `margin_m` on every side. A negative margin
    /// shrinks the box. An empty box stays empty.
    pub fn expanded_m(&self, margin_m: f64) -> Bounds {
        Bounds {
            min_m: std::array::from_fn(|axis| self.min_m[axis] - margin_m),
            max_m: std::array::from_fn(|axis| self.max_m[axis] + margin_m),
        }
    }

    /// Returns the centre of the box, or `[0.0; 3]` when the box is empty.
    pub fn center_m(&self) -> [f64; 3] {
        if self.is_empty() {
            [0.0; 3]
        } else {
            std::array::from_fn(|axis| 0.5 * (self.min_m[axis] + self.max_m[axis]))
        }
    }

    /// Returns true when any low corner is above the matching high corner.
    pub fn is_empty(&self) -> bool {
        self.min_m[0] > self.max_m[0]
            || self.min_m[1] > self.max_m[1]
            || self.min_m[2] > self.max_m[2]
    }
}

impl Solid for Bounds {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        if self.is_empty() {
            return false;
        }
        self.min_m[0] <= point_m[0]
            && point_m[0] <= self.max_m[0]
            && self.min_m[1] <= point_m[1]
            && point_m[1] <= self.max_m[1]
            && self.min_m[2] <= point_m[2]
            && point_m[2] <= self.max_m[2]
    }

    fn bounds_m(&self) -> Bounds {
        *self
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(*self)
    }
}

/// A closed solid region in metres.
pub trait Solid {
    /// Returns true when `point_m` lies in the closed region.
    fn contains_m(&self, point_m: [f64; 3]) -> bool;

    /// Returns an axis-aligned box that holds the whole region.
    fn bounds_m(&self) -> Bounds;

    /// Returns a new boxed copy of this solid.
    ///
    /// A box of a trait object cannot use [`Clone`], so the trait carries this
    /// method. It lets [`Placed`] and the boolean operations clone their inner
    /// shapes.
    fn clone_box(&self) -> Box<dyn Solid>;
}

/// A box aligned with the axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Box3 {
    /// The low corner, component by component, in metres.
    pub min_m: [f64; 3],
    /// The high corner, component by component, in metres.
    pub max_m: [f64; 3],
}

impl Box3 {
    /// Returns the box with the given centre and side lengths, in metres.
    pub fn from_center_m(center_m: [f64; 3], size_m: [f64; 3]) -> Self {
        Box3 {
            min_m: std::array::from_fn(|axis| center_m[axis] - 0.5 * size_m[axis]),
            max_m: std::array::from_fn(|axis| center_m[axis] + 0.5 * size_m[axis]),
        }
    }
}

impl Solid for Box3 {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        self.min_m[0] <= point_m[0]
            && point_m[0] <= self.max_m[0]
            && self.min_m[1] <= point_m[1]
            && point_m[1] <= self.max_m[1]
            && self.min_m[2] <= point_m[2]
            && point_m[2] <= self.max_m[2]
    }

    fn bounds_m(&self) -> Bounds {
        Bounds {
            min_m: self.min_m,
            max_m: self.max_m,
        }
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(*self)
    }
}

/// A regular hexagonal prism. The axis is z.
///
/// The hexagon has a vertex at angle `rotation_rad` and six vertices at
/// `rotation_rad + k * 60°`. `radius_m` is the circumradius, so a vertex lies
/// on it. The six edge normals point at `rotation_rad + 30° + k * 60°`, and the
/// apothem is `radius_m * cos(30°)`. A point is inside when its z coordinate
/// lies within `center_m[2] ± height_m / 2` and it lies inside all six edge
/// half-planes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HexPrism {
    /// The centre of the prism, in metres.
    pub center_m: [f64; 3],
    /// The circumradius of the hexagon, in metres.
    pub radius_m: f64,
    /// The length of the prism along z, in metres.
    pub height_m: f64,
    /// The rotation of the hexagon about z, in radians.
    pub rotation_rad: f64,
}

impl Solid for HexPrism {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        let dz_m = point_m[2] - self.center_m[2];
        if dz_m.abs() > 0.5 * self.height_m {
            return false;
        }
        let dx_m = point_m[0] - self.center_m[0];
        let dy_m = point_m[1] - self.center_m[1];
        let apothem_m = self.radius_m * (PI / 6.0).cos();
        for k in 0..6 {
            let angle_rad = self.rotation_rad + PI / 6.0 + (k as f64) * (PI / 3.0);
            let (sin_rad, cos_rad) = angle_rad.sin_cos();
            if dx_m * cos_rad + dy_m * sin_rad > apothem_m {
                return false;
            }
        }
        true
    }

    fn bounds_m(&self) -> Bounds {
        let mut bounds = Bounds::empty();
        for k in 0..6 {
            let angle_rad = self.rotation_rad + (k as f64) * (PI / 3.0);
            bounds.include_point_m([
                self.center_m[0] + self.radius_m * angle_rad.cos(),
                self.center_m[1] + self.radius_m * angle_rad.sin(),
                self.center_m[2],
            ]);
        }
        bounds.min_m[2] = self.center_m[2] - 0.5 * self.height_m;
        bounds.max_m[2] = self.center_m[2] + 0.5 * self.height_m;
        bounds
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(*self)
    }
}

/// A circular cylinder. The axis is z.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cylinder {
    /// The centre of the cylinder, in metres.
    pub center_m: [f64; 3],
    /// The radius of the circular face, in metres.
    pub radius_m: f64,
    /// The length of the cylinder along z, in metres.
    pub height_m: f64,
}

impl Solid for Cylinder {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        let dz_m = point_m[2] - self.center_m[2];
        if dz_m.abs() > 0.5 * self.height_m {
            return false;
        }
        let dx_m = point_m[0] - self.center_m[0];
        let dy_m = point_m[1] - self.center_m[1];
        dx_m * dx_m + dy_m * dy_m <= self.radius_m * self.radius_m
    }

    fn bounds_m(&self) -> Bounds {
        Bounds {
            min_m: [
                self.center_m[0] - self.radius_m,
                self.center_m[1] - self.radius_m,
                self.center_m[2] - 0.5 * self.height_m,
            ],
            max_m: [
                self.center_m[0] + self.radius_m,
                self.center_m[1] + self.radius_m,
                self.center_m[2] + 0.5 * self.height_m,
            ],
        }
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(*self)
    }
}

/// A solid rotated about the z axis and then offset.
pub struct Placed {
    /// The solid before the transform.
    pub inner: Box<dyn Solid>,
    /// The rotation about the z axis, in radians.
    pub rotation_rad: f64,
    /// The translation after the rotation, in metres.
    pub offset_m: [f64; 3],
}

impl Solid for Placed {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        let dx_m = point_m[0] - self.offset_m[0];
        let dy_m = point_m[1] - self.offset_m[1];
        let dz_m = point_m[2] - self.offset_m[2];
        let (sin_rad, cos_rad) = (-self.rotation_rad).sin_cos();
        let local_m = [
            cos_rad * dx_m - sin_rad * dy_m,
            sin_rad * dx_m + cos_rad * dy_m,
            dz_m,
        ];
        self.inner.contains_m(local_m)
    }

    fn bounds_m(&self) -> Bounds {
        let inner = self.inner.bounds_m();
        if inner.is_empty() {
            return inner;
        }
        let (sin_rad, cos_rad) = self.rotation_rad.sin_cos();
        let xs_m = [inner.min_m[0], inner.max_m[0]];
        let ys_m = [inner.min_m[1], inner.max_m[1]];
        let zs_m = [inner.min_m[2], inner.max_m[2]];
        let mut bounds = Bounds::empty();
        for &x_m in &xs_m {
            for &y_m in &ys_m {
                for &z_m in &zs_m {
                    bounds.include_point_m([
                        cos_rad * x_m - sin_rad * y_m + self.offset_m[0],
                        sin_rad * x_m + cos_rad * y_m + self.offset_m[1],
                        z_m + self.offset_m[2],
                    ]);
                }
            }
        }
        bounds
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(self.clone())
    }
}

impl Clone for Placed {
    fn clone(&self) -> Self {
        Placed {
            inner: self.inner.clone_box(),
            rotation_rad: self.rotation_rad,
            offset_m: self.offset_m,
        }
    }
}

impl fmt::Debug for Placed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Placed")
    }
}

/// The set difference `outer \ inner`.
pub struct Difference {
    /// The solid that defines the region.
    pub outer: Box<dyn Solid>,
    /// The solid that is removed from `outer`.
    pub inner: Box<dyn Solid>,
}

impl Solid for Difference {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        self.outer.contains_m(point_m) && !self.inner.contains_m(point_m)
    }

    fn bounds_m(&self) -> Bounds {
        self.outer.bounds_m()
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(self.clone())
    }
}

impl Clone for Difference {
    fn clone(&self) -> Self {
        Difference {
            outer: self.outer.clone_box(),
            inner: self.inner.clone_box(),
        }
    }
}

impl fmt::Debug for Difference {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Difference")
    }
}

/// The set intersection.
pub struct Intersection {
    /// The first solid.
    pub a: Box<dyn Solid>,
    /// The second solid.
    pub b: Box<dyn Solid>,
}

impl Solid for Intersection {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        self.a.contains_m(point_m) && self.b.contains_m(point_m)
    }

    fn bounds_m(&self) -> Bounds {
        let a = self.a.bounds_m();
        let b = self.b.bounds_m();
        Bounds {
            min_m: std::array::from_fn(|axis| a.min_m[axis].max(b.min_m[axis])),
            max_m: std::array::from_fn(|axis| a.max_m[axis].min(b.max_m[axis])),
        }
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(self.clone())
    }
}

impl Clone for Intersection {
    fn clone(&self) -> Self {
        Intersection {
            a: self.a.clone_box(),
            b: self.b.clone_box(),
        }
    }
}

impl fmt::Debug for Intersection {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Intersection")
    }
}

/// The set union.
pub struct Union {
    /// The first solid.
    pub a: Box<dyn Solid>,
    /// The second solid.
    pub b: Box<dyn Solid>,
}

impl Solid for Union {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        self.a.contains_m(point_m) || self.b.contains_m(point_m)
    }

    fn bounds_m(&self) -> Bounds {
        self.a.bounds_m().union(&self.b.bounds_m())
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(self.clone())
    }
}

impl Clone for Union {
    fn clone(&self) -> Self {
        Union {
            a: self.a.clone_box(),
            b: self.b.clone_box(),
        }
    }
}

impl fmt::Debug for Union {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Union")
    }
}

/// A 2D polygon extruded through all z.
///
/// The region is unbounded in z, so [`Solid::bounds_m`] returns
/// `[NEG_INFINITY, INFINITY]` on z. The caller must set the z limits when it
/// combines a profile.
///
/// When `internal` is false, the region is the polygon interior. When
/// `internal` is true, the region is the disc of `outer_radius_m` about the
/// origin, minus the polygon interior. An internal profile with no
/// `outer_radius_m` is empty.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile {
    /// The polygon corners in the xy plane, in metres.
    pub outline_m: Vec<[f64; 2]>,
    /// True for an internal ring: the region is outside the polygon.
    pub internal: bool,
    /// The outer radius of the disc for an internal profile, in metres.
    pub outer_radius_m: Option<f64>,
}

impl Solid for Profile {
    fn contains_m(&self, point_m: [f64; 3]) -> bool {
        let inside_outline = point_in_polygon([point_m[0], point_m[1]], &self.outline_m);
        if !self.internal {
            return inside_outline;
        }
        match self.outer_radius_m {
            None => false,
            Some(radius_m) => {
                let dx_m = point_m[0];
                let dy_m = point_m[1];
                dx_m * dx_m + dy_m * dy_m <= radius_m * radius_m && !inside_outline
            }
        }
    }

    fn bounds_m(&self) -> Bounds {
        let mut bounds = Bounds::empty();
        for point_m in &self.outline_m {
            bounds.include_point_m([point_m[0], point_m[1], 0.0]);
        }
        if self.internal {
            if let Some(radius_m) = self.outer_radius_m {
                for sign_x in [-1.0, 1.0] {
                    for sign_y in [-1.0, 1.0] {
                        bounds.include_point_m([sign_x * radius_m, sign_y * radius_m, 0.0]);
                    }
                }
            }
        }
        bounds.min_m[2] = f64::NEG_INFINITY;
        bounds.max_m[2] = f64::INFINITY;
        bounds
    }

    fn clone_box(&self) -> Box<dyn Solid> {
        Box::new(self.clone())
    }
}

/// Returns true when `point_xy_m` lies inside or on `outline_m`.
///
/// The test casts a ray along +x and counts the edges it crosses. The code
/// checks the boundary first, so a point exactly on an edge is inside.
fn point_in_polygon(point_xy_m: [f64; 2], outline_m: &[[f64; 2]]) -> bool {
    if outline_m.len() < 3 {
        return false;
    }
    let mut previous_m = match outline_m.last() {
        Some(point_m) => *point_m,
        None => return false,
    };
    let mut inside = false;
    for &current_m in outline_m {
        if point_on_segment(point_xy_m, previous_m, current_m) {
            return true;
        }
        let crosses = (current_m[1] > point_xy_m[1]) != (previous_m[1] > point_xy_m[1]);
        if crosses {
            let denominator = previous_m[1] - current_m[1];
            let x_cross = (previous_m[0] - current_m[0]) * (point_xy_m[1] - current_m[1])
                / denominator
                + current_m[0];
            if point_xy_m[0] < x_cross {
                inside = !inside;
            }
        }
        previous_m = current_m;
    }
    inside
}

/// Returns true when `point_m` lies on the segment from `a_m` to `b_m`.
///
/// The cross product measures the distance from the line, and the dot product
/// keeps the point between the two ends. The tolerance scales with the segment
/// length so that the test uses the floating-point resolution of the query.
fn point_on_segment(point_m: [f64; 2], a_m: [f64; 2], b_m: [f64; 2]) -> bool {
    let ab_x = b_m[0] - a_m[0];
    let ab_y = b_m[1] - a_m[1];
    let length_sq = ab_x * ab_x + ab_y * ab_y;
    if length_sq <= 0.0 {
        return false;
    }
    let ap_x = point_m[0] - a_m[0];
    let ap_y = point_m[1] - a_m[1];
    let tolerance = f64::EPSILON * length_sq * 8.0;
    let cross = ab_x * ap_y - ab_y * ap_x;
    if cross.abs() > tolerance {
        return false;
    }
    let dot = ap_x * ab_x + ap_y * ab_y;
    dot >= -tolerance && dot <= length_sq + tolerance
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_box_contains_its_inside_and_not_its_outside() {
        let cube = Box3::from_center_m([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        assert!(cube.contains_m([0.0, 0.0, 0.0]));
        assert!(!cube.contains_m([2.0, 0.0, 0.0]));
    }

    #[test]
    fn bounds_of_a_box_are_the_box_corners() {
        let cube = Box3::from_center_m([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        assert_eq!(
            cube.bounds_m(),
            Bounds {
                min_m: [-1.0, -2.0, -3.0],
                max_m: [1.0, 2.0, 3.0],
            }
        );
    }

    #[test]
    fn two_boxes_union_and_intersect() {
        let a = Box3::from_center_m([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = Box3::from_center_m([2.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let union = Union {
            a: Box::new(a),
            b: Box::new(b),
        };
        let intersection = Intersection {
            a: Box::new(a),
            b: Box::new(b),
        };
        assert_eq!(
            union.bounds_m(),
            Bounds {
                min_m: [-1.0, -1.0, -1.0],
                max_m: [3.0, 1.0, 1.0],
            }
        );
        assert_eq!(
            intersection.bounds_m(),
            Bounds {
                min_m: [1.0, -1.0, -1.0],
                max_m: [1.0, 1.0, 1.0],
            }
        );
        assert!(union.contains_m([2.5, 0.0, 0.0]));
        assert!(union.contains_m([-0.5, 0.0, 0.0]));
        assert!(intersection.contains_m([1.0, 0.0, 0.0]));
        assert!(!intersection.contains_m([2.0, 0.0, 0.0]));
    }

    #[test]
    fn a_placed_shape_rotates_about_z() {
        let cube = Box3::from_center_m([2.0e-9, 0.0, 0.0], [1.0e-9, 1.0e-9, 1.0e-9]);
        let placed = Placed {
            inner: Box::new(cube),
            rotation_rad: std::f64::consts::FRAC_PI_2,
            offset_m: [0.0, 0.0, 0.0],
        };
        assert!(placed.contains_m([0.0, 2.0e-9, 0.0]));
        assert!(!cube.contains_m([0.0, 2.0e-9, 0.0]));
    }

    #[test]
    fn a_hex_prism_has_six_sides() {
        let radius_m = 1.0e-9;
        let hex = HexPrism {
            center_m: [0.0, 0.0, 0.0],
            radius_m,
            height_m: 1.0e-9,
            rotation_rad: 0.0,
        };
        assert!(hex.contains_m([radius_m, 0.0, 0.0]));
        let apothem_m = radius_m * (PI / 6.0).cos();
        let beyond_m = apothem_m + 1.0e-11;
        assert!(!hex.contains_m([
            beyond_m * (PI / 6.0).cos(),
            beyond_m * (PI / 6.0).sin(),
            0.0
        ]));
    }

    #[test]
    fn a_difference_removes_a_hole() {
        let outer = Box3::from_center_m([0.0; 3], [2.0; 3]);
        let inner = Box3::from_center_m([0.0; 3], [1.0; 3]);
        let difference = Difference {
            outer: Box::new(outer),
            inner: Box::new(inner),
        };
        assert!(!difference.contains_m([0.0, 0.0, 0.0]));
        assert!(difference.contains_m([0.75, 0.0, 0.0]));
    }

    #[test]
    fn an_internal_profile_is_outside_the_outline() {
        let half_m = 2.0e-10;
        let outline_m = vec![
            [-half_m, -half_m],
            [half_m, -half_m],
            [half_m, half_m],
            [-half_m, half_m],
        ];
        let profile = Profile {
            outline_m,
            internal: true,
            outer_radius_m: Some(1.0e-9),
        };
        assert!(!profile.contains_m([0.0, 0.0, 0.0]));
        assert!(profile.contains_m([0.8e-9, 0.0, 0.0]));
    }

    #[test]
    fn an_empty_bounds_unions_with_a_point() {
        let empty = Bounds::empty();
        let point = Bounds::from_point_m([1.0, 2.0, 3.0]);
        assert!(empty.is_empty());
        assert_eq!(empty.center_m(), [0.0, 0.0, 0.0]);
        assert_eq!(empty.union(&point), point);
    }
}
