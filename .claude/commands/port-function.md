# /port-function

Port a single VarDictJava function to Rust.

Usage: `/port-function <Class.method>` (e.g. `/port-function CigarModifier.stripLeadingTrailingDel`)

## Steps
1. Read the cited Java source file and line range
2. Write a one-sentence doc-comment restating what the function does
3. Cite the source: `/// ports CigarModifier.java:61-71`
4. Implement the given Rust signature
5. Add the unit test named in the task
6. Run: `cargo fmt --all && cargo clippy --all-targets -- -D warnings && cargo test`
7. If all pass, commit: `git add -A && git commit -m "Txxx: <desc>"`
