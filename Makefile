.PHONY: fmt fix 

all: check

fix:
	cargo clippy --fix --allow-dirty

fmt:
	cargo fmt --all

check:
	cargo clippy

wasm-web:
	wasm-pack build --target web --mode no-install --release --out-dir=../frontend/pkg wasm_wrapper