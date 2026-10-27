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
        "Authorization": f"token {token}",
        "Accept": "application/vnd.github.v3+json",
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


def configure_metadata():
    print("[1/5] Configuring repository metadata and settings...")
    payload = {
        "description": "Fast, memory-safe, idiomatic Rust rewrite of VarDictJava v1.8.3 for deep NGS variant calling",
        "homepage": "https://rewrites.bio",
        "has_issues": True,
        "has_projects": False,
        "has_wiki": True,
        "allow_squash_merge": True,
        "allow_merge_commit": True,
        "allow_rebase_merge": True,
        "delete_branch_on_merge": True,
        "security_and_analysis": {
            "secret_scanning_non_provider_patterns": {"status": "enabled"},
            "secret_scanning_validity_checks": {"status": "enabled"},
        },
    }
    res = api_request("", method="PATCH", data=payload)
    print(f"  ✓ Description: {res.get('description')}")
    print(f"  ✓ Homepage: {res.get('homepage')}")


def configure_topics():
    print("[2/5] Configuring repository topics...")
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
    res = api_request("/topics", method="PUT", data={"names": topics})
    print(f"  ✓ Topics set: {', '.join(res.get('names', []))}")


def enable_dependabot():
    print("[3/5] Enabling vulnerability alerts and automated security fixes (Dependabot)...")
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


def enable_private_vulnerability_reporting():
    print("[4/5] Enabling Private Vulnerability Reporting...")
    try:
        api_request("/private-vulnerability-reporting", method="PUT")
        print("  ✓ Private vulnerability reporting enabled.")
    except urllib.error.HTTPError as e:
        print(f"  ⚠ Note: private vulnerability reporting status: {e.code}")


def configure_branch_protection():
    print("[5/5] Configuring branch protection on `main`...")
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
    res = api_request("/branches/main/protection", method="PUT", data=payload)
    checks = res.get("required_status_checks", {}).get("contexts", [])
    strict = res.get("required_status_checks", {}).get("strict", False)
    print(f"  ✓ Branch protection enabled on `main`.")
    print(f"  ✓ Strict branch up-to-date enforcement: {strict}")
    print(f"  ✓ Required status checks: {', '.join(checks)}")
    print(f"  ✓ Force pushes allowed: {res.get('allow_force_pushes', {}).get('enabled', False)}")
    print(f"  ✓ Deletions allowed: {res.get('allow_deletions', {}).get('enabled', False)}")


def main():
    print(f"=== VarDictRust: GitHub Repository Configuration ({REPO}) ===")
    configure_metadata()
    configure_topics()
    enable_dependabot()
    enable_private_vulnerability_reporting()
    configure_branch_protection()
    print("=== All GitHub repository configurations applied successfully! ===")


if __name__ == "__main__":
    main()
