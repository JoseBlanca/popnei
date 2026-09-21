//! The errors of the core crate on their way to a JavaScript `Error`.
//!
//! `impl From<popnei::Error> for JsValue` cannot be written here, because
//! neither type belongs to this crate, so every function of the crate fails
//! with [`JsPopneiError`], which does belong to it, and wasm-bindgen throws
//! it as an `Error`. `?` on a call of the core crate works everywhere, which
//! is what `.claude/skills/coding/pyo3.md` asks of the Python binding crate
//! and what section 11 of `docs/architecture.md` asks here: one place turns
//! an error of the core into what JavaScript throws.

use wasm_bindgen::{JsError, JsValue};

/// What a function of this crate fails with.
///
/// It has no case for an argument this crate refuses on its own: the
/// TypeScript package checks every argument before the call, in
/// `js/popnei/src/arguments.ts`, because a number of JavaScript reaches a
/// whole number of the core as 32 bits with no error, and it throws the
/// `Error` itself.
pub enum JsPopneiError {
    /// Something the core crate refused: an argument it takes, or what it
    /// found in the bytes it was given.
    Core(popnei::Error),
    /// Something the core read that JavaScript does not hold: a position
    /// above 2^53, which a float64 rounds.
    NotInJavaScript(String),
    /// Something that cannot happen unless this crate has a defect: a
    /// chromosome whose number is not in the table of the reader that gave
    /// it.
    Broken(String),
}

impl From<popnei::Error> for JsPopneiError {
    fn from(error: popnei::Error) -> JsPopneiError {
        JsPopneiError::Core(error)
    }
}

impl From<JsPopneiError> for JsValue {
    /// The `Error` that JavaScript catches, with the message the error has
    /// in Rust.
    ///
    /// JavaScript has one exception for everything a library refuses, so
    /// the three cases are one `Error`, where Python tells a `ValueError`
    /// from an `OSError`.
    fn from(error: JsPopneiError) -> JsValue {
        let message = match error {
            JsPopneiError::Core(error) => error.to_string(),
            JsPopneiError::NotInJavaScript(message) | JsPopneiError::Broken(message) => message,
        };
        JsError::new(&message).into()
    }
}
