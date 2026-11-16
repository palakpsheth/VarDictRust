#!/usr/bin/env python3
"""
Automated GitHub Repository, Branch Protection & Tag Ruleset Configuration.
Configures repository metadata, security features, `main` branch protection, and
`tag-protection-v*` (`refs/tags/v*`) ruleset via GitHub REST API using GITHUB_TOKEN.
"""

import argparse
import json
import os
import sys
import urllib.error
import urllib.request

REPO = "palakpsheth/VarDictRust"
API_BASE = f"https://api.github.com/repos/{REPO}"
TAG_RULESET_ID = 24867690
TAG_RULESET_NAME = "tag-protection-v*"
TAG_REF_PATTERN = "refs/tags/v*"


def get_token():
    token = os.environ.get("GITHUB_TOKEN")
    if not token:
        print("ERROR: GITHUB_TOKEN environment variable not set.", file=sys.stderr)
        sys.exit(1)
    return token


def api_request(endpoint, method="GET", data=None):
    token = get_token()
    url = endpoint if endpoint.startswith("http") else f"{API_BASE}{endpoint}"
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "User-Agent": "VarDictRust-Configurator",
    }
    body = json.dumps(data).encode("utf-8") if data is not None else None
    if data is not None:
        headers["Content-Type"] = "application/json"

    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
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
    if dry_run:
        print(f"  [DRY-RUN] PATCH {API_BASE} allow_merge_commit={payload['allow_merge_commit']}")
        return payload
    res = api_request("", method="PATCH", data=payload)
    print(f"  ✓ Description: {res.get('description')}")
    print(f"  ✓ Homepage: {res.get('homepage')}")
    print(f"  ✓ Allow merge commits: {res.get('allow_merge_commit')}")
    return res


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
    payload = {"names": topics}
    if dry_run:
        print(f"  [DRY-RUN] PUT {API_BASE}/topics ({len(topics)} topics)")
        return payload
    res = api_request("/topics", method="PUT", data=payload)
    print(f"  ✓ Topics set: {', '.join(res.get('names', []))}")
    return res


def enable_dependabot(dry_run=False):
    print("[3/6] Enabling vulnerability alerts and automated security fixes (Dependabot)...")
    if dry_run:
        print(f"  [DRY-RUN] PUT {API_BASE}/vulnerability-alerts")
        print(f"  [DRY-RUN] PUT {API_BASE}/automated-security-fixes")
        return
    try:
        api_request("/vulnerability-alerts", method="PUT")
        print("  ✓ Vulnerability alerts enabled.")
    except urllib.error.HTTPError as e:
        print(f"  ⚠ Note: vulnerability alerts status: {e.code}")
    try:
        api_request("/automated-security-fixes", method="PUT")
        print("  ✓ Dependabot automated security fixes enabled.")
    except urllib.error.HTTPError as e:
        print(f"  ⚠ Note: automated security fixes status: {e.code}")


def enable_private_vulnerability_reporting(dry_run=False):
    print("[4/6] Enabling Private Vulnerability Reporting...")
    if dry_run:
        print(f"  [DRY-RUN] PUT {API_BASE}/private-vulnerability-reporting")
        return
    try:
        api_request("/private-vulnerability-reporting", method="PUT")
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
    if dry_run:
        print(
            f"  [DRY-RUN] PUT {API_BASE}/branches/main/protection "
            f"checks={payload['required_status_checks']['contexts']}"
        )
        return payload
    res = api_request("/branches/main/protection", method="PUT", data=payload)
    checks = res.get("required_status_checks", {}).get("contexts", [])
    strict = res.get("required_status_checks", {}).get("strict", False)
    print("  ✓ Branch protection enabled on `main`.")
    print(f"  ✓ Strict branch up-to-date enforcement: {strict}")
    print(f"  ✓ Required status checks: {', '.join(checks)}")
    print(f"  ✓ Force pushes allowed: {res.get('allow_force_pushes', {}).get('enabled', False)}")
    print(f"  ✓ Deletions allowed: {res.get('allow_deletions', {}).get('enabled', False)}")
    return res


def configure_tag_ruleset(dry_run=False):
    print(f"[6/6] Configuring `{TAG_RULESET_NAME}` ruleset ({TAG_RULESET_ID}) for `{TAG_REF_PATTERN}`...")
    base_rules = [
        {"type": "deletion"},
        {"type": "update"},
        {"type": "non_fast_forward"},
    ]
    strict_payload = {
        "name": TAG_RULESET_NAME,
        "target": "tag",
        "enforcement": "active",
        "conditions": {
            "ref_name": {
                "include": [TAG_REF_PATTERN],
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
    fallback_payload = {**strict_payload, "rules": base_rules, "bypass_actors": []}
    endpoint = f"/rulesets/{TAG_RULESET_ID}"
    if dry_run:
        print(
            f"  [DRY-RUN] PUT {API_BASE}{endpoint} include={[TAG_REF_PATTERN]} "
            f"rules={[r['type'] for r in strict_payload['rules']]} "
            "with post-PUT GET current_user_can_bypass=='always' lockout rollback guard"
        )
        return strict_payload

    verify_res = None
    try:
        api_request(endpoint, method="PUT", data=strict_payload)
    except urllib.error.HTTPError as exc:
        if exc.code == 422:
            print(
                f"  ⚠ Warning: strict tag ruleset PUT on {endpoint} returned HTTP 422 ({exc}); "
                "falling back to ruleset without RepositoryRole bypass_actors / creation block.",
                file=sys.stderr,
            )
            can_bypass = None
        else:
            raise
    else:
        try:
            verify_res = api_request(endpoint, method="GET")
            can_bypass = verify_res.get("current_user_can_bypass") if isinstance(verify_res, dict) else None
        except Exception as exc:
            print(
                f"  ⚠ Warning: post-PUT verification GET on {endpoint} failed ({exc}); "
                "treating current_user_can_bypass as None and rolling back creation restriction.",
                file=sys.stderr,
            )
            can_bypass = None

    if can_bypass != "always":
        print(
            f"  ⚠ Lockout guard triggered (current_user_can_bypass={can_bypass!r} != 'always'); "
            "rolling back to deletion/update/non_fast_forward without creation block.",
            file=sys.stderr,
        )
        rollback_res = api_request(endpoint, method="PUT", data=fallback_payload)
        print("  ✓ Safe fallback tag ruleset applied (deletion, update, non_fast_forward).")
        return rollback_res

    print(f"  ✓ Tag ruleset `{TAG_RULESET_NAME}` active with verified admin bypass (current_user_can_bypass='always').")
    return verify_res


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Configure GitHub repository metadata, branch protection, and tag rulesets."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned GitHub REST API calls without mutating repository settings (default).",
    )
    mode.add_argument(
        "--apply",
        action="store_true",
        help="Apply live GitHub repository, branch protection, and tag ruleset settings via GITHUB_TOKEN.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    dry_run = not args.apply
    mode_label = "DRY-RUN" if dry_run else "APPLY"
    print(f"=== VarDictRust: GitHub Repository Configuration ({REPO}) [{mode_label}] ===")
    configure_metadata(dry_run=dry_run)
    configure_topics(dry_run=dry_run)
    enable_dependabot(dry_run=dry_run)
    enable_private_vulnerability_reporting(dry_run=dry_run)
    configure_branch_protection(dry_run=dry_run)
    configure_tag_ruleset(dry_run=dry_run)
    if dry_run:
        print("=== Dry-run complete (no live GitHub settings mutated; pass --apply to execute). ===")
    else:
        print("=== All GitHub repository configurations applied successfully! ===")


if __name__ == "__main__":
    main()
