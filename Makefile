.PHONY: fmt fix 

all: check

fix:
	cargo clippy --fix --allow-dirty

fmt:
	cargo fmt --all

check:
	cargo clippy