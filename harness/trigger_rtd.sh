#!/usr/bin/env bash
# harness/trigger_rtd.sh — ReadTheDocs API v3 CLI helper
#
# Usage:
#   bash harness/trigger_rtd.sh status [project_slug]
#   bash harness/trigger_rtd.sh build [project_slug] [branch]
#   bash harness/trigger_rtd.sh builds [project_slug]

set -euo pipefail

CMD="${1:-status}"
PROJECT_SLUG="${2:-vardictrust}"
BRANCH="${3:-main}"

# Source RTD_TOKEN from ~/.bashrc if not present in current shell environment
if [ -z "${RTD_TOKEN:-}" ] && [ -f "$HOME/.bashrc" ]; then
    RTD_TOKEN="$(python3 -c 'import os, re; p=os.path.expanduser("~/.bashrc"); m=re.search(r"export RTD_TOKEN=[\"\x27]?([^\s\"\x27]+)", open(p).read()) if os.path.exists(p) else None; print(m.group(1) if m else "")')"
fi

if [ -z "${RTD_TOKEN:-}" ]; then
    echo "ERROR: RTD_TOKEN environment variable is not set." >&2
    echo "Export RTD_TOKEN in your environment or ~/.bashrc." >&2
    exit 1
fi

BASE_URL="https://readthedocs.org/api/v3/projects/${PROJECT_SLUG}"

case "$CMD" in
    status)
        echo "=== Querying ReadTheDocs Project Status: ${PROJECT_SLUG} ==="
        RESPONSE="$(curl -s -w "\nHTTP_STATUS:%{http_code}" \
            -H @- <<< "Authorization: Token ${RTD_TOKEN}" \
            "${BASE_URL}/")"
        STATUS_CODE="$(echo "$RESPONSE" | grep 'HTTP_STATUS:' | cut -d':' -f2)"
        BODY="$(echo "$RESPONSE" | grep -v 'HTTP_STATUS:')"
        if [ "$STATUS_CODE" -eq 200 ]; then
            echo "$BODY" | python3 -m json.tool
        else
            echo "HTTP Status: ${STATUS_CODE}" >&2
            echo "$BODY" | python3 -m json.tool || echo "$BODY"
            exit 1
        fi
        ;;
    build)
        VERSION="${3:-latest}"
        echo "=== Triggering ReadTheDocs Build: ${PROJECT_SLUG} (version: ${VERSION}) ==="
        RESPONSE="$(curl -s -w "\nHTTP_STATUS:%{http_code}" -X POST \
            -H @- <<< "Authorization: Token ${RTD_TOKEN}" \
            -H "Content-Type: application/json" \
            "${BASE_URL}/versions/${VERSION}/builds/")"
        STATUS_CODE="$(echo "$RESPONSE" | grep 'HTTP_STATUS:' | cut -d':' -f2)"
        BODY="$(echo "$RESPONSE" | grep -v 'HTTP_STATUS:')"
        if [ "$STATUS_CODE" -eq 202 ]; then
            echo "Build triggered successfully (HTTP 202 Accepted):"
            echo "$BODY" | python3 -m json.tool
        else
            echo "Build trigger returned HTTP Status: ${STATUS_CODE}" >&2
            echo "$BODY" | python3 -m json.tool || echo "$BODY"
            exit 1
        fi
        ;;
    builds)
        echo "=== Querying Recent Builds: ${PROJECT_SLUG} ==="
        RESPONSE="$(curl -s -w "\nHTTP_STATUS:%{http_code}" \
            -H @- <<< "Authorization: Token ${RTD_TOKEN}" \
            "${BASE_URL}/builds/")"
        STATUS_CODE="$(echo "$RESPONSE" | grep 'HTTP_STATUS:' | cut -d':' -f2)"
        BODY="$(echo "$RESPONSE" | grep -v 'HTTP_STATUS:')"
        if [ "$STATUS_CODE" -eq 200 ]; then
            echo "$BODY" | python3 -m json.tool
        else
            echo "HTTP Status: ${STATUS_CODE}" >&2
            echo "$BODY" | python3 -m json.tool || echo "$BODY"
            exit 1
        fi
        ;;
    *)
        echo "Usage: $0 [status|build|builds] [project_slug] [version]" >&2
        exit 1
        ;;
esac
