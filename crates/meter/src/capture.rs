//! The sidewall capture of a bound guest.
//!
//! A sorting pocket holds its guest anywhere in the pocket void. The guest can
//! sit at the pocket axis, or it can bind against the sidewall. The ejection
//! rod has one pushing face. This module measures how much of the pocket the
//! face can reach, and it states the dead zone that the face cannot reach.
//!
//! The model is rigid and geometric. The pocket void is a cylinder of
//! `pocket_radius_m`. The guest is a sphere of `guest_radius_m`. The guest
//! centre can sit anywhere within `pocket_radius_m - guest_radius_m` of the
//! pocket axis. The pushing face is a disc of `piston_radius_m` that advances
//! along the axis. The face touches the guest when the centre offset is at
//! most `piston_radius_m + guest_radius_m`.
//!
//! A face that reaches every allowed offset captures every bound guest. A face
//! that is too narrow leaves a dead zone against the wall. The report states
//! the dead zone and the face radius that closes it.
//!
//! The model ignores the guest shape, its orientation, the wall chemistry, the
//! friction and the solvent. It is a coverage check, not a force model. See
//! ADR-0072 and ADR-0066. All lengths are SI metres.

/// The geometry of one pocket and one pushing face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaptureTarget {
    /// The radius of the pocket void, in metres.
    pub pocket_radius_m: f64,
    /// The radius of the bound guest, in metres.
    pub guest_radius_m: f64,
    /// The radius of the pushing face, in metres.
    pub piston_radius_m: f64,
}

impl Default for CaptureTarget {
    fn default() -> Self {
        Self {
            pocket_radius_m: 1.0e-9,
            guest_radius_m: 2.24e-10,
            piston_radius_m: 2.0e-10,
        }
    }
}

/// The measured sidewall capture.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaptureReport {
    /// The radius of the pocket void, in metres.
    pub pocket_radius_m: f64,
    /// The radius of the bound guest, in metres.
    pub guest_radius_m: f64,
    /// The radius of the pushing face, in metres.
    pub piston_radius_m: f64,
    /// The largest guest-centre offset from the pocket axis, in metres.
    pub allowed_offset_m: f64,
    /// The largest guest-centre offset the face can touch, in metres.
    pub reach_m: f64,
    /// The largest offset that the face both allows and touches, in metres.
    pub covered_offset_m: f64,
    /// The radial width of the wall region that the face misses, in metres.
    pub dead_zone_m: f64,
    /// The covered fraction of the allowed offsets, from zero to one.
    pub coverage: f64,
    /// True when the face reaches every allowed offset.
    pub complete: bool,
    /// The face radius that closes the dead zone, in metres.
    pub required_piston_radius_m: f64,
}

/// Measures how much of the pocket one pushing face can reach.
///
/// A non-finite or non-positive pocket radius reports a complete zero void.
pub fn capture_report(target: &CaptureTarget) -> CaptureReport {
    let pocket_radius_m = target.pocket_radius_m;
    let guest_radius_m = target.guest_radius_m.max(0.0);
    let piston_radius_m = target.piston_radius_m.max(0.0);

    if !pocket_radius_m.is_finite() || pocket_radius_m <= 0.0 {
        return CaptureReport {
            pocket_radius_m,
            guest_radius_m,
            piston_radius_m,
            allowed_offset_m: 0.0,
            reach_m: piston_radius_m + guest_radius_m,
            covered_offset_m: 0.0,
            dead_zone_m: 0.0,
            coverage: 1.0,
            complete: true,
            required_piston_radius_m: 0.0,
        };
    }

    let allowed_offset_m = (pocket_radius_m - guest_radius_m).max(0.0);
    let reach_m = piston_radius_m + guest_radius_m;
    let covered_offset_m = reach_m.min(allowed_offset_m);
    let dead_zone_m = (allowed_offset_m - covered_offset_m).max(0.0);
    let coverage = if allowed_offset_m > 0.0 {
        covered_offset_m / allowed_offset_m
    } else {
        1.0
    };
    let complete = covered_offset_m >= allowed_offset_m;

    CaptureReport {
        pocket_radius_m,
        guest_radius_m,
        piston_radius_m,
        allowed_offset_m,
        reach_m,
        covered_offset_m,
        dead_zone_m,
        coverage,
        complete,
        required_piston_radius_m: (pocket_radius_m - 2.0 * guest_radius_m).max(0.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_narrow_face_leaves_a_dead_zone() {
        let report = capture_report(&CaptureTarget::default());
        assert!(!report.complete);
        assert!(report.dead_zone_m > 0.0);
        assert!(report.coverage < 1.0);
        assert!((report.allowed_offset_m - 7.76e-10).abs() < 1.0e-15);
        assert!((report.reach_m - 4.24e-10).abs() < 1.0e-15);
    }

    #[test]
    fn the_required_face_closes_the_dead_zone() {
        let target = CaptureTarget::default();
        let report = capture_report(&target);
        let wide = CaptureTarget {
            piston_radius_m: report.required_piston_radius_m + 1.0e-15,
            ..target
        };
        let closed = capture_report(&wide);
        assert!(closed.complete);
        assert_eq!(closed.dead_zone_m, 0.0);
        assert_eq!(closed.coverage, 1.0);
    }

    #[test]
    fn a_face_as_wide_as_the_pocket_always_captures() {
        let target = CaptureTarget {
            pocket_radius_m: 8.0e-10,
            guest_radius_m: 2.0e-10,
            piston_radius_m: 8.0e-10,
        };
        let report = capture_report(&target);
        assert!(report.complete);
        assert_eq!(report.required_piston_radius_m, 4.0e-10);
    }

    #[test]
    fn a_guest_that_fills_the_pocket_is_captured() {
        let target = CaptureTarget {
            pocket_radius_m: 3.0e-10,
            guest_radius_m: 3.0e-10,
            piston_radius_m: 0.0,
        };
        let report = capture_report(&target);
        assert!(report.complete);
        assert_eq!(report.allowed_offset_m, 0.0);
    }

    #[test]
    fn a_zero_pocket_is_complete() {
        let report = capture_report(&CaptureTarget {
            pocket_radius_m: 0.0,
            ..CaptureTarget::default()
        });
        assert!(report.complete);
        assert_eq!(report.coverage, 1.0);
    }
}
