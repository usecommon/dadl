# Security policy

## Supported versions

DADL is pre-alpha. Only the latest 0.0.x release receives fixes.

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
- Parsing of setup files and tag files (crashes, resource exhaustion, unsafe
  YAML loading).
- The release supply chain: workflows, published packages and their
  provenance.

Approval tokens, artifact signing, the audit log, plugin loading and the MCP
server are designed but not built; they will join this scope when they ship.

QKD monitoring records are evidence of what the monitor observed. They are not
proof of security, and reports that they cannot detect an unknown attack are
expected behaviour, not vulnerabilities.
