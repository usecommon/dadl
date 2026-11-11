# Security policy

DADL drives lab hardware and signs calibration evidence, so please report
vulnerabilities privately.

## Reporting

Use GitHub's private vulnerability reporting ("Report a vulnerability" on the
repository's Security tab). Please do not open a public issue.

Include what you found, how to reproduce it, and which version you tested.

## Scope

In scope: the approval-token flow, hardware limit enforcement, guardrail
bypasses, artifact signing and verification, the audit log, plugin loading,
and the MCP server.

QKD monitoring records are evidence of what the monitor observed. They are not
proof of security, and reports that they cannot detect an unknown attack are
expected behaviour, not vulnerabilities.
