#[macro_export]
macro_rules! assert_approx_eq {
        ($a:expr, $b:expr) => {
            assert_approx_eq!($a, $b, 1e-3);
        };
        ($a:expr, $b:expr, $eps:expr) => {{
            let a = $a;
            let b = $b;
            let eps:f64 = $eps;
            let diff = (a - b).abs();
            let threshold = eps.max(eps * a.abs().max(b.abs()));
            assert!(
                diff < threshold,
                "assert_approx_eq failed:\n  a = {a:.15}\n  b = {b:.15}\n  |a-b| = {diff:.15}\n  threshold = {threshold:.15}"
            );
        }};
    }
