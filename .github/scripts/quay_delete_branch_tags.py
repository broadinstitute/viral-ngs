#!/usr/bin/env python3
"""Delete a deleted branch's image tags from quay.io/broadinstitute/viral-ngs.

Run by `.github/workflows/cleanup-images.yml` when a branch is deleted. It can also be
run by hand (that workflow's workflow_dispatch, or locally) to clean up leftovers.
See issue #1125.

How tags get deleted matters. quay.io rejects the `delete` token scope that
`crane delete` asks for, so this script gets a `pull,push` token and calls
`DELETE /v2/<repo>/manifests/<tag>`. On quay that removes only the named tag. Deleting
by *digest* removes every tag on that manifest. That is how the 3.0.10-baseimage
incident happened: branch `-baseimage` builds usually reproduce main's (and the latest
release's) digest. So this script never deletes by digest.

Safety rails:
  * Refuses (exit 2) a branch whose tags would include main's tags, the bare flavor
    aliases, `latest`, or anything version-like (e.g. branch `vmain` or `core`).
  * Skips (exit 0) if a live branch still maps to the same tag prefix
    (`feature/foo` and `feature-foo`, or `foo` and `vfoo`), since the images are now
    that branch's.
  * Checks main's tags before and after deleting, and exits 1 if any vanished.
  * Dry run unless --execute. A dry run needs no credentials.

Usage:
    quay_delete_branch_tags.py --branch feature/foo              # dry run
    QUAY_USERNAME=... QUAY_TOKEN=... \\
        quay_delete_branch_tags.py --branch feature/foo --execute

Exit status: 0 ok (or skipped), 1 delete/auth errors or a main tag vanished, 2 refused.
"""

import argparse
import base64
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request

REGISTRY = 'quay.io'
DEFAULT_REPO = 'broadinstitute/viral-ngs'

# One tag per flavor. Must match the deploy-to-quay matrix in docker.yml (mega has no suffix).
FLAVOR_SUFFIXES = ('-baseimage', '-core', '-assemble', '-classify', '-phylo', '')
# Tags docker.yml points at main's images, besides main-* ("latest" is mega).
ALIAS_TAGS = ('latest', 'baseimage', 'core', 'assemble', 'classify', 'phylo')
CANARY_TAGS = tuple('main' + s for s in FLAVOR_SUFFIXES) + ALIAS_TAGS
VERSION_LIKE = re.compile(r'^\d+\.\d+')

MANIFEST_TYPES = ', '.join((
    'application/vnd.oci.image.index.v1+json',
    'application/vnd.docker.distribution.manifest.list.v2+json',
    'application/vnd.oci.image.manifest.v1+json',
    'application/vnd.docker.distribution.manifest.v2+json',
))
RETRY_STATUSES = {429, 500, 502, 503, 504}
ATTEMPTS = 4


class AuthError(Exception):
    pass


def tag_prefix(branch):
    """Same sanitization as docker.yml's "Compute Docker image tag prefix" step."""
    prefix = branch.replace('/', '-')
    return prefix[1:] if prefix.startswith('v') else prefix


def branch_tags(prefix):
    return [prefix + suffix for suffix in FLAVOR_SUFFIXES]


def is_protected(tag):
    return tag in CANARY_TAGS or VERSION_LIKE.match(tag) is not None


def refusal_reason(prefix):
    if not prefix:
        return 'empty tag prefix'
    protected = [t for t in branch_tags(prefix) if is_protected(t)]
    if protected:
        return 'its tags would include protected tags: ' + ', '.join(protected)
    return None


def parse_ls_remote(text):
    heads = []
    for line in text.splitlines():
        _, _, ref = line.partition('\t')
        if ref.startswith('refs/heads/'):
            heads.append(ref[len('refs/heads/'):])
    return heads


def live_branches(repo):
    result = subprocess.run(['git', 'ls-remote', '--heads', f'https://github.com/{repo}'],
                            check=True, capture_output=True, text=True)
    return parse_ls_remote(result.stdout)


def granted_actions(token, repo):
    payload = token.split('.')[1]
    claims = json.loads(base64.urlsafe_b64decode(payload + '=' * (-len(payload) % 4)))
    return {action for entry in claims.get('access', []) if entry.get('name') == repo
            for action in entry.get('actions', [])}


def authenticate(repo, execute, env=os.environ, opener=urllib.request.urlopen):
    """Registry token: pull,push with robot credentials to delete; anonymous pull otherwise."""
    actions = 'pull,push' if execute else 'pull'
    req = urllib.request.Request(
        f'https://{REGISTRY}/v2/auth?service={REGISTRY}&scope=repository:{repo}:{actions}')
    if execute:
        user, password = env.get('QUAY_USERNAME'), env.get('QUAY_TOKEN')
        if not user or not password:
            raise AuthError('--execute needs QUAY_USERNAME and QUAY_TOKEN in the environment')
        basic = base64.b64encode(f'{user}:{password}'.encode()).decode()
        req.add_header('Authorization', f'Basic {basic}')
    try:
        with opener(req, timeout=60) as resp:
            token = json.load(resp)['token']
    except urllib.error.HTTPError as e:
        raise AuthError(f'token request for {actions} on {repo} failed: HTTP {e.code}') from e
    if execute and 'push' not in granted_actions(token, repo):
        raise AuthError(f'quay did not grant push on {repo} to {user}; check its permissions')
    return token


class Quay:
    """Just enough of the registry v2 API to check and delete tags in one repo."""

    def __init__(self, repo, token, opener=urllib.request.urlopen, sleep=time.sleep):
        self.repo = repo
        self.token = token
        self.opener = opener
        self.sleep = sleep

    def head(self, tag):
        return self._status('HEAD', tag)

    def delete(self, tag):
        # By tag, never by digest: on quay this removes only this tag.
        return self._status('DELETE', tag)

    def _status(self, method, tag):
        req = urllib.request.Request(
            f'https://{REGISTRY}/v2/{self.repo}/manifests/{tag}', method=method,
            headers={'Authorization': f'Bearer {self.token}', 'Accept': MANIFEST_TYPES})
        status = 0  # 0 = no response
        for attempt in range(ATTEMPTS):
            try:
                with self.opener(req, timeout=60) as resp:
                    status = resp.status
            except urllib.error.HTTPError as e:
                status = e.code
                e.close()
            except urllib.error.URLError:
                status = 0
            if status and status not in RETRY_STATUSES:
                return status
            if attempt < ATTEMPTS - 1:
                self.sleep(2 ** attempt)
        return status


def annotate(out, level, message):
    if os.environ.get('GITHUB_ACTIONS') == 'true':
        print(f'::{level}::{message}', file=out)
    else:
        print(f'{level.upper()}: {message}', file=out)


def run(branch, quay, live_branches, execute, out=sys.stdout, summary=None):
    prefix = tag_prefix(branch)
    reason = refusal_reason(prefix)
    if reason:
        annotate(out, 'error', f'Refusing to delete tags for branch {branch!r}: {reason}')
        return 2

    twins = [b for b in live_branches if tag_prefix(b) == prefix]
    if twins:
        annotate(out, 'notice', f'Not deleting {prefix}* tags: live branch(es) {", ".join(twins)} '
                                f'use the same tag prefix.')
        return 0

    canary_before = [t for t in CANARY_TAGS if quay.head(t) == 200] if execute else []

    mode = 'execute' if execute else 'dry run'
    print(f'Branch {branch!r} -> tag prefix {prefix!r} ({mode})', file=out)
    rows = []
    for tag in branch_tags(prefix):
        if execute:
            status = quay.delete(tag)
            outcome = {200: 'deleted', 202: 'deleted', 204: 'deleted', 404: 'absent'}.get(status, 'ERROR')
        else:
            status = quay.head(tag)
            outcome = {200: 'present, would delete', 404: 'absent'}.get(status, 'ERROR')
        rows.append((tag, outcome, status))
        print(f'  {tag:<60} {outcome:<22} HTTP {status}', file=out)

    rc = 0
    errors = [tag for tag, outcome, _ in rows if outcome == 'ERROR']
    if errors:
        annotate(out, 'error', f'Failed on {len(errors)} tag(s): {", ".join(errors)}')
        rc = 1
    vanished = [t for t in canary_before if quay.head(t) == 404]
    if vanished:
        annotate(out, 'error',
                 f'Possible cascade: {", ".join(vanished)} resolved before the deletes and are gone now. '
                 f'Restore main-* and aliases by rerunning docker.yml deploy-to-quay on main; '
                 f'a repo admin can also restore them from quay\'s time machine.')
        rc = 1

    if summary is not None:
        summary.write(f'### quay.io branch tag cleanup: `{branch}` ({mode})\n\n'
                      '| tag | result | HTTP |\n|---|---|---|\n')
        summary.writelines(f'| `{tag}` | {outcome} | {status} |\n' for tag, outcome, status in rows)
        if vanished:
            summary.write(f'\n**Main tags vanished:** {", ".join(vanished)}\n')
    return rc


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--branch', default=os.environ.get('BRANCH'),
                        help='deleted branch name as pushed, e.g. feature/foo (default: $BRANCH)')
    parser.add_argument('--repo', default=DEFAULT_REPO, help='quay.io repository (default: %(default)s)')
    parser.add_argument('--execute', action='store_true', help='delete the tags (default: dry run)')
    args = parser.parse_args(argv)
    if args.branch is None:
        parser.error('--branch (or $BRANCH) is required')

    reason = refusal_reason(tag_prefix(args.branch))
    if reason:
        annotate(sys.stdout, 'error', f'Refusing to delete tags for branch {args.branch!r}: {reason}')
        return 2
    try:
        live = live_branches(args.repo)
        token = authenticate(args.repo, args.execute)
    except (AuthError, subprocess.CalledProcessError, urllib.error.URLError) as e:
        annotate(sys.stdout, 'error', str(e))
        return 1

    summary_path = os.environ.get('GITHUB_STEP_SUMMARY')
    if summary_path:
        with open(summary_path, 'a') as summary:
            return run(args.branch, Quay(args.repo, token), live, args.execute, summary=summary)
    return run(args.branch, Quay(args.repo, token), live, args.execute)


if __name__ == '__main__':
    sys.exit(main())
