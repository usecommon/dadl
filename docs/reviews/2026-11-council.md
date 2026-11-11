# Third review: ten reviewers and a council (2026-11-11)

DADL has no outside human reviewers, so this review was run by agents.

1. **Reviewers.** Ten independent agents each reviewed one area of v0.0.2 (the design, the plan and the code), read-only. Several ran probe scripts against the code:

   | # | Area | # | Area |
   | --- | --- | --- | --- |
   | 1 | physics | 6 | guardrail red team |
   | 2 | statistics | 7 | data and provenance |
   | 3 | architecture | 8 | plan and scope |
   | 4 | Python code | 9 | security and supply chain |
   | 5 | Rust core | 10 | user experience |

2. **Council.** A council with three seats (physicist, maintainer, user advocate) read all ten reports. It re-checked the main claims numerically and adopted a change when two of three seats voted for it.

3. **Result.** v0.0.3 applies everything the council adopted for now. The [design doc](../design.md) records the deferred items in its phase table.

The review is still agent-only. A human physicist's review remains the top open item.

## Checked numerically by the council

**Series dead time.** The model is a paralyzable stage followed by a non-paralyzable one:

    m = m1 / (1 + m1 (tau_np - tau_p)),   where m1 = n e^(-n tau_p)

- This matches a simulated chain to within 0.5%.
- The formula the doc used before was off by up to 31%.

**Corrected fringe visibility (Thompson case: 145, 55, A = 45.12).** The code reported ±0.115. The correct values are:

| Accidentals A | Delta method | Monte Carlo |
| --- | --- | --- |
| Exact | ±0.1245 | ±0.128 |
| Poisson | ±0.160 | ±0.178 |

**Attribution example.**
- A baseline of S = 2.6 requires A/C ≤ 0.088.
- Reaching S = 2.3 requires A/C = 0.230, a factor of 2.6 higher.
- A 40% rise in dark counts gives at most S = 2.52, so the example's claim is wrong.

**Bugs reproduced.**
- `pair_rate_per_s: 1.0e5` loaded as a string.
- Calling `analysis="hom_visibility"` on Bell data turned 1450/550 with A = 450 into a corrected V of 0.818.
- Visibility inputs (100, 10, A = 20) gave a corrected V of 1.29.

## Adopted and applied in v0.0.3

**Correctness**
- Dead time uses the exact series model. The empirical hybrid is kept, documented as "not an oracle".
- Visibility is split into `bell_visibility` and `hom_visibility`:
  - `hom_visibility` uses the dip convention.
  - Both use delta-method errors that include the accidentals' variance.
  - Inputs are checked so V cannot exceed 1.
- Threshold labels now use intervals:
  - within 1.96σ of the threshold;
  - zero Gaussian error at a boundary.
- Input checks run in every build, in both languages:
  - sorted input;
  - integer picoseconds;
  - int64 overflow.
- Rust functions return `Result` instead of panicking.
- Shared Rust/Python test vectors live in `tests/vectors/core.json`.

**Guardrails**
- Every result records the policy fingerprint, the overrides with their reasons, and the window.
- The fingerprint covers the thresholds and is checked against a pinned constant.
- The rule tables are read-only.
- Looser rules are allowed only through `Policy.allow`, which validates the correction name.
- CHSH gets threshold labels.

**Starter kit**
- `setup.yaml` drives the quickstart. It ships in the wheel.
- New commands:
  - `dadl init`
  - `dadl simulate`
  - `dadl analyze` (generic `.npz`/`.csv` tags)
- The config loader is hardened:
  - types and units checked;
  - unknown and duplicate keys rejected;
  - size cap.
- CLI:
  - new `--version` flag;
  - `python -m dadl` works;
  - errors exit with code 1 and a message.

**Docs**
- Consistent physics, citations and numbers. 82.84% is now credited to Garg–Mermin.
- Current status is stated honestly:
  - "what works today";
  - "declared, not verified";
  - "statistical only".
- One canonical phase table.
- A schema example that CI validates.

**Packaging and CI**
- Licence files ship in the sdist and wheel.
- SPDX licence expression and a dynamic version.
- New CI jobs:
  - sdist/wheel test;
  - MSRV 1.74 check;
  - release-mode tests;
  - clippy;
  - pinned toolchain;
  - a 200-seed pull study.
- Dependabot is configured.
- `persist-credentials: false` is set on checkouts.
- The release workflow builds in one job and publishes from another.
- A release tag must match both package versions and must sit on main.

**Repository files**
- SECURITY.md is scoped to what exists and gives response targets.
- CONTRIBUTING has a review rule a single maintainer can actually follow.
- New CHANGELOG and CITATION.cff.
- More patterns in `.gitignore`.

### Where the implementation differs from the ruling

**`bell_visibility` keeps an `accidentals` argument.**
- The ruling asked for none.
- The default policy still refuses it.
- It is kept so the published guardrail regression test ("a loosened policy is never silent") keeps working through `Policy.allow`.
- The relabelling bypass is closed either way, because the analysis type can no longer be passed in.

**The example artifact uses kind `dark-afterpulse-deadtime`.**
- The ruling suggested adding a standalone `afterpulse` kind.
- The design requires parameters fitted from the same data to share one joint posterior, so a separate kind would contradict it.

## Deferred (recorded in the design doc's phase table)

**Phase 2**
- Dataset/TagSet and ChannelId types.
- DetectorModel protocol.
- Wilson/Beta intervals, and an interval type on `Estimate`.
- SeedSequence streams.
- Afterpulse closed-form test and crosstalk simulation.
- Units and transforms on posterior parameters.
- PyO3 bindings crate:
  - Rust parameter names aligned with Python;
  - full `# Errors` docs.
- `setup-1` JSON Schema.
- Adapter template and the hom-snspd example.
- Photon-HDF5 output under `/user`, after a check with phconvert.

**Phase 3**
- Streaming O(n) coincidences.
- Signed analysis plans with trial counting.
- Adjusted 95% attribution intervals.
- JCS test vectors and a canonical posterior hash.

**v0.2+ / seeking contributors**
- Signing.
- MCP and approval tokens.
- Audit log.
- HSM signer.
- QKD attestation.
- Monitor daemon.
- Arrow C Stream adapters.
- Plugin hash pinning.

## Rejected

| Proposal | Reason |
| --- | --- |
| Self-hosted fonts and a CSP on the static report | No user data on the page. |
| Wiring the Rust core in before PyO3 | Python stays the reference. |
| Signing in phase 1 | Moved later. |
| Code of conduct, issue templates and a notebook now | No outside users yet. |
| `aria-describedby` polish | Optional. |
| cargo-audit and pip-audit jobs | Dependabot covers them. |
| `#![deny(missing_docs)]` before the API freeze | Premature. |

## Repository settings only the owner can change

These were not changed by the review. Each needs the repository owner.

- [ ] **Private vulnerability reporting.** SECURITY.md relies on it, and it is currently off.
  - Repo → Settings → Security (Code security / Advanced Security) → Private vulnerability reporting → Enable.
- [ ] **Protect the `pypi` environment.**
  - Settings → Environments → pypi:
    - Required reviewers: add yourself.
    - Deployment branches and tags: Selected, with a tag rule `v*`.
    - Untick "Allow administrators to bypass".
- [ ] **Tag ruleset.**
  - Settings → Rules → Rulesets → New tag ruleset.
  - Target `v*`.
  - Restrict creations, updates and deletions to the owner (via the bypass list).
- [ ] **Branch ruleset on `main`.**
  - Require a pull request.
  - Require status checks (the CI jobs).
  - Require linear history.
  - Block force pushes and deletion.
- [ ] **Actions settings.**
  - Settings → Actions → General.
  - Workflow permissions: Read.
  - Require approval for fork pull requests.
- [ ] **crates.io.**
  - crates.io → dadl → Settings → Trusted Publishing: add repository `usecommon/dadl`, workflow `release.yml`, environment `crates`.
  - Then Account Settings → API Tokens: revoke the manual token, and run `cargo logout`.
  - After that, a `rust-lang/crates-io-auth-action` job can replace manual publishing.
- [ ] **PyPI.**
  - Confirm the Trusted Publisher names environment `pypi`.
  - Confirm no API tokens exist (Account → API tokens).
- [ ] **Dependabot.** Enable Dependabot alerts and security updates: Settings → Security.
