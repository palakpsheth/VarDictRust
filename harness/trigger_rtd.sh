#!/usr/bin/env bash
# harness/trigger_rtd.sh — ReadTheDocs API v3 CLI helper
#
# Usage:
#   bash harness/trigger_rtd.sh status [project_slug]
#   bash harness/trigger_rtd.sh sync [project_slug] [expected_stable_ref]
#   bash harness/trigger_rtd.sh build [project_slug] [branch]
#   bash harness/trigger_rtd.sh builds [project_slug]

set -euo pipefail

CMD="${1:-status}"
PROJECT_SLUG="${2:-vardictrust}"
VERSION="${3:-latest}"
SYNC_REF="${3:-}"

if [[ "$PROJECT_SLUG" == "." || "$PROJECT_SLUG" == ".." || ! "$PROJECT_SLUG" =~ ^[a-zA-Z0-9._-]+$ ]]; then
    echo "ERROR: Invalid ReadTheDocs project slug '${PROJECT_SLUG}' (allowed: [a-zA-Z0-9._-]+, excluding '.' and '..')." >&2
    exit 1
fi
if [[ "$VERSION" == "." || "$VERSION" == ".." || ! "$VERSION" =~ ^[a-zA-Z0-9._-]+$ ]]; then
    echo "ERROR: Invalid ReadTheDocs version '${VERSION}' (allowed: [a-zA-Z0-9._-]+, excluding '.' and '..')." >&2
    exit 1
fi

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
        RESPONSE="$(curl --proto '=https' --tlsv1.2 -s -w "\nHTTP_STATUS:%{http_code}" \
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
    sync)
        echo "=== Syncing ReadTheDocs Versions: ${PROJECT_SLUG} ==="
        RESPONSE="$(curl --proto '=https' --tlsv1.2 -s -w "\nHTTP_STATUS:%{http_code}" -X POST \
            -H @- <<< "Authorization: Token ${RTD_TOKEN}" \
            "${BASE_URL}/sync-versions/")"
        STATUS_CODE="$(echo "$RESPONSE" | grep 'HTTP_STATUS:' | cut -d':' -f2)"
        BODY="$(echo "$RESPONSE" | grep -v 'HTTP_STATUS:')"
        if [ "$STATUS_CODE" -eq 202 ] || [ "$STATUS_CODE" -eq 204 ]; then
            echo "Version sync triggered successfully (HTTP ${STATUS_CODE}):"
            if [ -n "$BODY" ]; then
                echo "$BODY" | python3 -m json.tool || echo "$BODY"
            fi
        else
            echo "Version sync returned HTTP Status: ${STATUS_CODE}" >&2
            echo "$BODY" | python3 -m json.tool || echo "$BODY"
            exit 1
        fi
        if [ -n "$SYNC_REF" ]; then
            MAX_ATTEMPTS="${RTD_SYNC_MAX_ATTEMPTS:-15}"
            POLL_SECS="${RTD_SYNC_POLL_SECS:-2}"
            if [[ ! "$MAX_ATTEMPTS" =~ ^[1-9][0-9]*$ ]] || [[ ! "$POLL_SECS" =~ ^[0-9]+$ ]]; then
                echo "ERROR: Invalid RTD_SYNC_MAX_ATTEMPTS='${MAX_ATTEMPTS}' or RTD_SYNC_POLL_SECS='${POLL_SECS}'." >&2
                exit 1
            fi
            attempt=1
            CUR_REF=""
            while [ "$attempt" -le "$MAX_ATTEMPTS" ]; do
                POLL_RESP="$(curl --proto '=https' --tlsv1.2 -s -w "\nHTTP_STATUS:%{http_code}" \
                    -H @- <<< "Authorization: Token ${RTD_TOKEN}" \
                    "${BASE_URL}/versions/stable/")"
                POLL_CODE="$(echo "$POLL_RESP" | grep 'HTTP_STATUS:' | cut -d':' -f2)"
                POLL_BODY="$(echo "$POLL_RESP" | grep -v 'HTTP_STATUS:')"
                if [ "$POLL_CODE" -eq 200 ]; then
                    CUR_REF="$(python3 -c 'import json, sys; d = json.load(sys.stdin); print(d.get("ref") or d.get("identifier") or "")' <<<"$POLL_BODY")"
                    if [ "$CUR_REF" = "$SYNC_REF" ]; then
                        echo "Verified ReadTheDocs 'stable' converged to '${SYNC_REF}' (attempt ${attempt}/${MAX_ATTEMPTS})."
                        break
                    fi
                fi
                if [ "$attempt" -eq "$MAX_ATTEMPTS" ]; then
                    echo "ERROR: ReadTheDocs 'stable' did not converge to '${SYNC_REF}' after ${MAX_ATTEMPTS} attempts (last ref='${CUR_REF}', HTTP ${POLL_CODE})." >&2
                    exit 1
                fi
                if [ "$POLL_SECS" -gt 0 ]; then
                    sleep "$POLL_SECS"
                fi
                attempt=$((attempt + 1))
            done
        fi
        ;;
    build)
        echo "=== Triggering ReadTheDocs Build: ${PROJECT_SLUG} (version: ${VERSION}) ==="
        RESPONSE="$(curl --proto '=https' --tlsv1.2 -s -w "\nHTTP_STATUS:%{http_code}" -X POST \
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
        RESPONSE="$(curl --proto '=https' --tlsv1.2 -s -w "\nHTTP_STATUS:%{http_code}" \
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
        echo "Usage: $0 [status|sync|build|builds] [project_slug] [version]" >&2
        exit 1
        ;;
esac
