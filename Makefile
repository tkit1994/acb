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

wasm-node:
	wasm-pack build --target nodejs --mode no-install --release --out-dir=../frontend/pkg wasm_wrapper
mex:
	cargo build -p mex_wrapper --release

build:
	cargo build --release