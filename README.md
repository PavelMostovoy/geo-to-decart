# geo-to-decart

A small Rust library to convert geodetic coordinates to ECEF and to local ENU coordinates, suitable for small areas.

## Usage

Add to Cargo.toml:

[dependencies]
geo-to-decart = { path = "." }

Example:

```rust
use geo_to_decart::{geodetic_to_ecef, llh_to_enu, centroid_lat_lon};

fn main() {
    // reference point from centroid
    let points = vec![(55.0005, 37.0007), (54.9960, 36.9950)];
    let (lat0, lon0) = centroid_lat_lon(&points);
    let h0 = 120.0;

    let (e, n, u) = llh_to_enu(55.0005, 37.0007, 120.0, lat0, lon0, h0);
    println!("E={:.3} N={:.3} U={:.3}", e, n, u);

    let (_x, _y, _z) = geodetic_to_ecef(55.0, 37.0, 120.0);
}
```
