use std::f64::consts::PI;

// WGS-84 constants
const A: f64 = 6_378_137.0; // semi-major axis, meters
const E2: f64 = 6.694_379_9901413165e-3; // first eccentricity squared

fn deg_to_rad(deg: f64) -> f64 {
    deg * PI / 180.0
}

/// Convert geodetic (lat, lon, h) to ECEF (X, Y, Z)
/// lat, lon in degrees; h in meters
pub fn geodetic_to_ecef(lat_deg: f64, lon_deg: f64, h: f64) -> (f64, f64, f64) {
    let lat = deg_to_rad(lat_deg);
    let lon = deg_to_rad(lon_deg);

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let cos_lon = lon.cos();
    let sin_lon = lon.sin();

    let n = A / (1.0 - E2 * sin_lat * sin_lat).sqrt();

    let x = (n + h) * cos_lat * cos_lon;
    let y = (n + h) * cos_lat * sin_lon;
    let z = (n * (1.0 - E2) + h) * sin_lat;

    (x, y, z)
}

/// Convert ECEF -> ENU relative to reference point (lat0, lon0, h0)
/// lat0, lon0 in degrees; h0 in meters
/// Returns (east, north, up) in meters
pub fn ecef_to_enu(x: f64, y: f64, z: f64, lat0_deg: f64, lon0_deg: f64, h0: f64) -> (f64, f64, f64) {
    // ECEF of the reference point
    let (x0, y0, z0) = geodetic_to_ecef(lat0_deg, lon0_deg, h0);

    // Deltas
    let dx = x - x0;
    let dy = y - y0;
    let dz = z - z0;

    // Angles of the reference point (radians)
    let lat0 = deg_to_rad(lat0_deg);
    let lon0 = deg_to_rad(lon0_deg);

    let sin_lat0 = lat0.sin();
    let cos_lat0 = lat0.cos();
    let sin_lon0 = lon0.sin();
    let cos_lon0 = lon0.cos();

    // Rotation matrix for ENU
    let east = -sin_lon0 * dx + cos_lon0 * dy;
    let north = -sin_lat0 * cos_lon0 * dx - sin_lat0 * sin_lon0 * dy + cos_lat0 * dz;
    let up = cos_lat0 * cos_lon0 * dx + cos_lat0 * sin_lon0 * dy + sin_lat0 * dz;

    (east, north, up)
}

/// Convenience wrapper: (lat, lon, h) -> (e, n, u) relative to reference (lat0, lon0, h0)
pub fn llh_to_enu(lat_deg: f64, lon_deg: f64, h: f64, lat0_deg: f64, lon0_deg: f64, h0: f64) -> (f64, f64, f64) {
    let (x, y, z) = geodetic_to_ecef(lat_deg, lon_deg, h);
    ecef_to_enu(x, y, z, lat0_deg, lon0_deg, h0)
}

/// Compute centroid latitude/longitude for a list of points (lat, lon) in degrees.
/// Simple average; suitable for relatively small areas.
pub fn centroid_lat_lon(points: &[(f64, f64)]) -> (f64, f64) {
    let n = points.len() as f64;
    let mut sum_lat = 0.0;
    let mut sum_lon = 0.0;
    for &(lat, lon) in points {
        sum_lat += lat;
        sum_lon += lon;
    }
    (sum_lat / n, sum_lon / n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ecef_enu_roundtrip_center_zero() {
        let lat = 50.0;
        let lon = 10.0;
        let h = 200.0;
        let (x, y, z) = geodetic_to_ecef(lat, lon, h);
        let (e, n, u) = ecef_to_enu(x, y, z, lat, lon, h);
        assert!(e.abs() < 1e-9);
        assert!(n.abs() < 1e-9);
        assert!(u.abs() < 1e-9);
    }

    #[test]
    fn small_offset_check() {
        let lat0 = 42.680067;
        let lon0 = 3.034061;
        let h0 = 0.0;

        let lat1 = 42.680499;
        let lon1 = 3.035775;
        let h1 = 1.0;

        let (e, n, _u) = llh_to_enu(lat1, lon1, h1, lat0, lon0, h0);
        let approx_n = (super::deg_to_rad(lat1 - lat0)) * A;
        println!("e= {:.3}, n={:.3} approx={:.3}",e, n, approx_n);
        assert!((n - approx_n).abs() < 5.0);
        assert!(e.abs() > 0.0);
    }
}
