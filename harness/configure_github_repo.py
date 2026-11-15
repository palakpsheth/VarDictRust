#!/usr/bin/env python3
"""
Automated GitHub Repository & Branch Protection Configuration.
Configures repository metadata, security features, and branch protection on `main`
via GitHub REST API v3 using GITHUB_TOKEN.
"""

import json
import os
import sys
import urllib.error
import urllib.request

REPO = "palakpsheth/VarDictRust"
API_BASE = f"https://api.github.com/repos/{REPO}"


def get_token(dry_run=False):
    token = os.environ.get("GITHUB_TOKEN")
    if not token:
        if dry_run:
            return "dry-run-token"
        print("ERROR: GITHUB_TOKEN environment variable not set.", file=sys.stderr)
        sys.exit(1)
    return token


def api_request(endpoint, method="GET", data=None, dry_run=False):
    if dry_run:
        print(f"  [dry-run] {method} {endpoint}")
        return {"current_user_can_bypass": "always"}
    token = get_token(dry_run=dry_run)
    url = endpoint if endpoint.startswith("http") else f"{API_BASE}{endpoint}"
    headers = {
        "Authorization": f"token {token}",
        "Accept": "application/vnd.github+json",
        "User-Agent": "VarDictRust-Configurator",
    }
    body = json.dumps(data).encode("utf-8") if data is not None else None
    if data is not None:
        headers["Content-Type"] = "application/json"

    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req) as resp:
            content = resp.read().decode("utf-8")
            return json.loads(content) if content else {}
    except urllib.error.HTTPError as e:
        err_msg = e.read().decode("utf-8")
        print(f"HTTPError {e.code} on {method} {url}: {err_msg}", file=sys.stderr)
        raise


def configure_metadata(dry_run=False):
    print("[1/6] Configuring repository metadata and settings...")
    payload = {
        "description": "Fast, memory-safe, idiomatic Rust rewrite of VarDictJava v1.8.3 for deep NGS variant calling",
        "homepage": "https://vardictrust.readthedocs.io/",
        "has_issues": True,
        "has_projects": False,
        "has_wiki": True,
        "allow_squash_merge": True,
        "allow_merge_commit": False,
        "allow_rebase_merge": True,
        "delete_branch_on_merge": True,
        "security_and_analysis": {
            "secret_scanning_non_provider_patterns": {"status": "enabled"},
            "secret_scanning_validity_checks": {"status": "enabled"},
        },
    }
    res = api_request("", method="PATCH", data=payload, dry_run=dry_run)
    print(f"  ✓ Description: {res.get('description', payload['description'])}")
    print(f"  ✓ Homepage: {res.get('homepage', payload['homepage'])}")
    print("  ✓ Merge commits disabled (allow_merge_commit=False)")


def configure_topics(dry_run=False):
    print("[2/6] Configuring repository topics...")
    topics = [
        "bioinformatics",
        "genomics",
        "variant-caller",
        "ngs",
        "somatic",
        "rust",
        "cancer-genomics",
        "rewrites-bio",
    ]
    res = api_request("/topics", method="PUT", data={"names": topics}, dry_run=dry_run)
    print(f"  ✓ Topics set: {', '.join(res.get('names', topics))}")


def enable_dependabot(dry_run=False):
    print("[3/6] Enabling vulnerability alerts and automated security fixes (Dependabot)...")
    try:
        api_request("/vulnerability-alerts", method="PUT", dry_run=dry_run)
        print("  ✓ Vulnerability alerts enabled.")
    except urllib.error.HTTPError as e:
        print(f"  ⚠ Note: vulnerability alerts status: {e.code}")
    try:
        api_request("/automated-security-fixes", method="PUT", dry_run=dry_run)
        print("  ✓ Dependabot automated security fixes enabled.")
    except urllib.error.HTTPError as e:
        print(f"  ⚠ Note: automated security fixes status: {e.code}")


def enable_private_vulnerability_reporting(dry_run=False):
    print("[4/6] Enabling Private Vulnerability Reporting...")
    try:
        api_request("/private-vulnerability-reporting", method="PUT", dry_run=dry_run)
        print("  ✓ Private vulnerability reporting enabled.")
    except urllib.error.HTTPError as e:
        print(f"  ⚠ Note: private vulnerability reporting status: {e.code}")


def configure_branch_protection(dry_run=False):
    print("[5/6] Configuring branch protection on `main`...")
    payload = {
        "required_status_checks": {
            "strict": True,
            "contexts": [
                "Build, Unit, Integration, Parity & Coverage Gate",
                "Dual-Engine Differential Parity (Tier-2)",
                "Security & Supply-Chain Scan",
                "MkDocs Strict Build Gate",
                "plan",
            ],
        },
        "enforce_admins": False,
        "required_pull_request_reviews": {
            "dismiss_stale_reviews": True,
            "require_code_owner_reviews": False,
            "required_approving_review_count": 0,
            "require_last_push_approval": False,
        },
        "restrictions": None,
        "required_linear_history": True,
        "allow_force_pushes": False,
        "allow_deletions": False,
        "block_creations": False,
        "required_conversation_resolution": True,
    }
    res = api_request("/branches/main/protection", method="PUT", data=payload, dry_run=dry_run)
    checks = res.get("required_status_checks", {}).get("contexts", payload["required_status_checks"]["contexts"])
    strict = res.get("required_status_checks", {}).get("strict", True)
    print(f"  ✓ Branch protection enabled on `main`.")
    print(f"  ✓ Strict branch up-to-date enforcement: {strict}")
    print(f"  ✓ Required status checks: {', '.join(checks)}")
    print(f"  ✓ Force pushes allowed: {res.get('allow_force_pushes', {}).get('enabled', False)}")
    print(f"  ✓ Deletions allowed: {res.get('allow_deletions', {}).get('enabled', False)}")


def configure_tag_ruleset(dry_run=False):
    RULESET_ID = 24867690
    print(f"[6/6] Configuring tag protection ruleset on `~refs/tags/v*` (ruleset {RULESET_ID})...")
    base_rules = [
        {"type": "deletion"},
        {"type": "update"},
        {"type": "non_fast_forward"},
    ]
    payload_strict = {
        "name": "tag-protection-v*",
        "target": "tag",
        "enforcement": "active",
        "conditions": {
            "ref_name": {
                "include": ["~refs/tags/v*"],
                "exclude": [],
            }
        },
        "rules": base_rules + [{"type": "creation"}],
        "bypass_actors": [
            {
                "actor_id": 5,
                "actor_type": "RepositoryRole",
                "bypass_mode": "always",
            }
        ],
    }

    if dry_run:
        api_request(f"/rulesets/{RULESET_ID}", method="PUT", data=payload_strict, dry_run=True)
        print(f"  [dry-run] Configured ruleset {RULESET_ID} for ~refs/tags/v* with deletion/update/non_fast_forward and creation lockout check.")
        return

    # Apply strict ruleset with creation block
    api_request(f"/rulesets/{RULESET_ID}", method="PUT", data=payload_strict, dry_run=False)

    # Immediately verify whether current user can bypass tag creation
    status = api_request(f"/rulesets/{RULESET_ID}", method="GET", dry_run=False)
    can_bypass = status.get("current_user_can_bypass")
    if can_bypass == "always":
        print(f"  ✓ Tag ruleset {RULESET_ID} active with creation protection (admin bypass confirmed: {can_bypass}).")
    else:
        print(f"  ⚠ Current user bypass mode is {can_bypass!r} != 'always'; rolling back creation restriction to prevent lockout...")
        payload_safe = {
            "name": "tag-protection-v*",
            "target": "tag",
            "enforcement": "active",
            "conditions": {
                "ref_name": {
                    "include": ["~refs/tags/v*"],
                    "exclude": [],
                }
            },
            "rules": base_rules,
            "bypass_actors": [],
        }
        api_request(f"/rulesets/{RULESET_ID}", method="PUT", data=payload_safe, dry_run=False)
        print(f"  ✓ Safe rollback applied: ruleset {RULESET_ID} enforces deletion/update/non_fast_forward without tag creation lockout.")


def main():
    dry_run = "--dry-run" in sys.argv
    mode_str = " (DRY-RUN)" if dry_run else ""
    print(f"=== VarDictRust: GitHub Repository Configuration ({REPO}){mode_str} ===")
    configure_metadata(dry_run=dry_run)
    configure_topics(dry_run=dry_run)
    enable_dependabot(dry_run=dry_run)
    enable_private_vulnerability_reporting(dry_run=dry_run)
    configure_branch_protection(dry_run=dry_run)
    configure_tag_ruleset(dry_run=dry_run)
    print(f"=== All GitHub repository configurations applied successfully!{mode_str} ===")


if __name__ == "__main__":
    main()
