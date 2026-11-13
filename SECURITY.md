# Security policy

## Supported versions

DADL is pre-alpha. Only the latest 0.x release receives fixes.

## Reporting

Use GitHub's private vulnerability reporting ("Report a vulnerability" on the
repository's Security tab). Please do not open a public issue. If that button
is not available, open a public issue that says only "security contact
requested", without details, and a maintainer will reach out.

Include what you found, how to reproduce it, and which version you tested. We
aim to acknowledge reports within 7 days and to agree a disclosure date with
you.

## Scope

In scope for the code that exists today:

- Guardrail bypasses: any way to obtain a corrected or threshold-crossing
  value without the refusal or labels the policy requires, other than
  deliberately constructing or editing a result in your own process (which
  the docs state cannot be prevented).
- Parsing of setup files, tag files (npz, csv, PTU) and calibration artifacts
  (crashes, resource exhaustion, unsafe YAML loading, decompression bombs).
- Artifact verification: a file that `read_artifact` accepts although its
  bytes, hashes or self-consistency checks should fail. (Artifacts are not
  signed: anyone who can write a file can mint a self-consistent one.)
- Output paths: the CLI writing outside `--out`, through symlinks, or over
  existing files without `--force`.
- The release supply chain: workflows, published packages and their
  provenance.

Approval tokens, artifact signing, the audit log, plugin loading and the MCP
server are designed but not built; they will join this scope when they ship.

DADL has no QKD monitor today: no blinding, after-gate, time-shift, dead-time
attack or efficiency-mismatch detection, and no QKD monitoring record. Nothing
it produces belongs in a QKD security argument. When monitoring records
arrive, they will be evidence of what the monitor observed, not proof of
security, and reports that they cannot detect an unknown attack will be
expected behaviour, not vulnerabilities.
