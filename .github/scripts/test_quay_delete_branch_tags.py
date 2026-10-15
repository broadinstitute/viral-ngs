"""Unit tests for quay_delete_branch_tags.py (no network access).

Run with: python3 -m unittest discover -s .github/scripts -p 'test_*.py' -v
"""

import base64
import http.client
import io
import json
import unittest
import urllib.error
from unittest import mock

import quay_delete_branch_tags as qd

REPO = 'broadinstitute/viral-ngs'
MAIN_TAGS = ['main', 'main-baseimage', 'main-core', 'main-assemble', 'main-classify', 'main-phylo',
             'latest', 'baseimage', 'core', 'assemble', 'classify', 'phylo']


class FakeQuay:
    """Stands in for qd.Quay: tracks live tags and records every call."""

    def __init__(self, tags, delete_status=None, cascade=None, head_status=None, head_status_after_delete=None):
        self.tags = set(tags)
        self.delete_status = delete_status or {}
        self.cascade = cascade or {}
        self.head_status = head_status or {}
        self.head_status_after_delete = head_status_after_delete or {}
        self.calls = []

    def head(self, tag):
        self.calls.append(('HEAD', tag))
        if self.deleted() and tag in self.head_status_after_delete:
            return self.head_status_after_delete[tag]
        if tag in self.head_status:
            return self.head_status[tag]
        return 200 if tag in self.tags else 404

    def delete(self, tag):
        self.calls.append(('DELETE', tag))
        if tag in self.delete_status:
            return self.delete_status[tag]
        if tag not in self.tags:
            return 404
        self.tags.discard(tag)
        self.tags.difference_update(self.cascade.get(tag, ()))
        return 202

    def deleted(self):
        return [t for m, t in self.calls if m == 'DELETE']


def run(branch, quay, live_branches=(), execute=True):
    out = io.StringIO()
    rc = qd.run(branch, quay, live_branches=list(live_branches), execute=execute, out=out)
    return rc, out.getvalue()


class TestTagNames(unittest.TestCase):

    def test_prefix_matches_docker_yml(self):
        self.assertEqual(qd.tag_prefix('feature/foo'), 'feature-foo')
        self.assertEqual(qd.tag_prefix('cve-fix/2026-10-15-38710458144'), 'cve-fix-2026-10-15-38710458144')
        self.assertEqual(qd.tag_prefix('fix-quay-branch-cleanup-1125'), 'fix-quay-branch-cleanup-1125')

    def test_prefix_strips_one_leading_v(self):
        self.assertEqual(qd.tag_prefix('vfoo'), 'foo')
        self.assertEqual(qd.tag_prefix('vvfoo'), 'vfoo')
        self.assertEqual(qd.tag_prefix('dev/vfoo'), 'dev-vfoo')

    def test_branch_tags_one_per_flavor_mega_unsuffixed(self):
        self.assertEqual(
            qd.branch_tags('feature-foo'),
            ['feature-foo-baseimage', 'feature-foo-core', 'feature-foo-assemble',
             'feature-foo-classify', 'feature-foo-phylo', 'feature-foo'])


class TestRefusals(unittest.TestCase):

    def assert_refused(self, branch):
        quay = FakeQuay(MAIN_TAGS + ['3.0.27-core', '3.0-core'])
        rc, out = run(branch, quay)
        self.assertEqual(rc, 2, f'{branch!r} should be refused; output:\n{out}')
        self.assertEqual(quay.deleted(), [], f'{branch!r} must not delete anything')

    def test_refuses_main_and_aliases_of_main(self):
        for branch in ('main', 'vmain', 'main/core', 'main-core'):
            with self.subTest(branch=branch):
                self.assert_refused(branch)

    def test_refuses_version_like(self):
        for branch in ('3.0.27', 'v3.0', 'v3.0.27-rc1', '3.1'):
            with self.subTest(branch=branch):
                self.assert_refused(branch)

    def test_refuses_branches_whose_mega_tag_is_an_alias(self):
        for branch in ('core', 'latest', 'baseimage', 'phylo', 'vassemble', 'classify'):
            with self.subTest(branch=branch):
                self.assert_refused(branch)

    def test_refuses_empty(self):
        self.assert_refused('')
        self.assert_refused('v')

    def test_refuses_invalid_docker_tags(self):
        # '#' is legal in git branch names; urllib would drop it as a URL fragment and
        # send DELETE /manifests/main. A digest reference would delete by digest.
        for branch in ('main#cleanup', 'feature#1', 'foo%23bar', 'sha256:' + 'a' * 64,
                       '-leading-dash', '.leading-dot', 'x' * 120):
            with self.subTest(branch=branch):
                self.assert_refused(branch)

    def test_longest_valid_tag_not_refused(self):
        self.assertIsNone(qd.refusal_reason('x' * (128 - len('-baseimage'))))

    def test_ordinary_branch_not_refused(self):
        self.assertIsNone(qd.refusal_reason('feature-foo'))
        self.assertIsNone(qd.refusal_reason('maintenance-docs'))
        self.assertIsNone(qd.refusal_reason('core-refactor'))


class TestCollisions(unittest.TestCase):

    def test_skips_when_another_live_branch_has_same_prefix(self):
        quay = FakeQuay(qd.branch_tags('feature-foo'))
        rc, out = run('feature-foo', quay, live_branches=['main', 'feature/foo'])
        self.assertEqual(rc, 0)
        self.assertEqual(quay.deleted(), [])
        self.assertIn('feature/foo', out)

    def test_skips_v_prefixed_twin(self):
        quay = FakeQuay(qd.branch_tags('foo'))
        rc, _ = run('foo', quay, live_branches=['vfoo'])
        self.assertEqual(rc, 0)
        self.assertEqual(quay.deleted(), [])

    def test_skips_when_branch_itself_still_exists(self):
        quay = FakeQuay(qd.branch_tags('feature-foo'))
        rc, _ = run('feature/foo', quay, live_branches=['main', 'feature/foo'])
        self.assertEqual(rc, 0)
        self.assertEqual(quay.deleted(), [])

    def test_keeps_tag_that_is_a_live_branch_mega_tag(self):
        # foo's -core tag is the same name as live branch foo-core's mega tag
        quay = FakeQuay(MAIN_TAGS + qd.branch_tags('foo') + qd.branch_tags('foo-core'))
        rc, out = run('foo', quay, live_branches=['main', 'foo-core'])
        self.assertEqual(rc, 0, out)
        self.assertNotIn('foo-core', quay.deleted())
        self.assertIn('foo-core', quay.tags)
        self.assertEqual(sorted(quay.deleted()),
                         sorted(t for t in qd.branch_tags('foo') if t != 'foo-core'))

    def test_keeps_mega_tag_that_is_a_live_branch_flavor_tag(self):
        # foo-core's mega tag is the same name as live branch foo's -core tag
        quay = FakeQuay(MAIN_TAGS + qd.branch_tags('foo') + qd.branch_tags('foo-core'))
        rc, out = run('foo-core', quay, live_branches=['main', 'foo'])
        self.assertEqual(rc, 0, out)
        self.assertNotIn('foo-core', quay.deleted())
        self.assertIn('foo-core', quay.tags)
        self.assertEqual(sorted(quay.deleted()),
                         sorted(t for t in qd.branch_tags('foo-core') if t != 'foo-core'))

    def test_parse_ls_remote(self):
        text = ('abc123\trefs/heads/main\n'
                'def456\trefs/heads/feature/foo\n')
        self.assertEqual(qd.parse_ls_remote(text), ['main', 'feature/foo'])


class TestDeletion(unittest.TestCase):

    def test_deletes_all_six_tags(self):
        tags = qd.branch_tags('feature-foo')
        quay = FakeQuay(MAIN_TAGS + tags)
        rc, _ = run('feature/foo', quay, live_branches=['main'])
        self.assertEqual(rc, 0)
        self.assertEqual(quay.deleted(), tags)
        self.assertTrue(set(MAIN_TAGS) <= quay.tags)
        self.assertFalse(set(tags) & quay.tags)

    def test_absent_tags_are_ok(self):
        quay = FakeQuay(MAIN_TAGS + ['feature-foo-baseimage'])
        rc, out = run('feature/foo', quay, live_branches=['main'])
        self.assertEqual(rc, 0, out)

    def test_error_statuses_fail_after_trying_every_tag(self):
        tags = qd.branch_tags('feature-foo')
        for status in (400, 401, 403):
            with self.subTest(status=status):
                quay = FakeQuay(MAIN_TAGS + tags, delete_status={'feature-foo-core': status})
                rc, out = run('feature/foo', quay, live_branches=['main'])
                self.assertEqual(rc, 1)
                self.assertEqual(quay.deleted(), tags)
                self.assertIn(str(status), out)

    def test_canary_catches_cascade(self):
        tags = qd.branch_tags('feature-foo')
        quay = FakeQuay(MAIN_TAGS + tags,
                        cascade={'feature-foo-baseimage': ['main-baseimage', 'baseimage']})
        rc, out = run('feature/foo', quay, live_branches=['main'])
        self.assertEqual(rc, 1)
        self.assertIn('main-baseimage', out)

    def test_canary_unreadable_before_aborts_without_deleting(self):
        tags = qd.branch_tags('feature-foo')
        for status in (0, 401, 403, 503):
            with self.subTest(status=status):
                quay = FakeQuay(MAIN_TAGS + tags, head_status={'main-core': status})
                rc, out = run('feature/foo', quay, live_branches=['main'])
                self.assertEqual(rc, 1)
                self.assertEqual(quay.deleted(), [])
                self.assertIn('main-core', out)

    def test_canary_unreadable_after_fails(self):
        tags = qd.branch_tags('feature-foo')
        for status in (0, 401, 503):
            with self.subTest(status=status):
                quay = FakeQuay(MAIN_TAGS + tags, head_status_after_delete={'latest': status})
                rc, out = run('feature/foo', quay, live_branches=['main'])
                self.assertEqual(rc, 1)
                self.assertIn('latest', out)

    def test_dry_run_never_deletes(self):
        tags = qd.branch_tags('feature-foo')
        quay = FakeQuay(MAIN_TAGS + tags)
        rc, out = run('feature/foo', quay, live_branches=['main'], execute=False)
        self.assertEqual(rc, 0)
        self.assertEqual(quay.deleted(), [])
        self.assertEqual(quay.tags, set(MAIN_TAGS + tags))
        for tag in tags:
            self.assertIn(tag, out)


def _jwt(actions):
    payload = {'access': [{'type': 'repository', 'name': REPO, 'actions': actions}]}
    body = base64.urlsafe_b64encode(json.dumps(payload).encode()).decode().rstrip('=')
    return f'header.{body}.signature'


class FakeResponse(io.BytesIO):

    def __init__(self, status=200, body=b''):
        super().__init__(body)
        self.status = status


class FakeOpener:
    """Stands in for urllib.request.urlopen: replays responses, records requests."""

    def __init__(self, responses):
        self.responses = list(responses)
        self.requests = []

    def __call__(self, req, timeout=None):
        self.requests.append(req)
        r = self.responses.pop(0)
        if isinstance(r, Exception):
            raise r
        if isinstance(r, int) and r >= 400:
            raise urllib.error.HTTPError(req.full_url, r, 'error', {}, io.BytesIO(b''))
        if isinstance(r, int):
            return FakeResponse(r)
        return r


class TestQuayClient(unittest.TestCase):

    def test_retries_server_errors(self):
        opener = FakeOpener([503, 202])
        quay = qd.Quay(REPO, 'tok', opener=opener, sleep=lambda s: None)
        self.assertEqual(quay.delete('feature-foo-core'), 202)
        self.assertEqual(len(opener.requests), 2)
        self.assertEqual(opener.requests[0].get_method(), 'DELETE')
        self.assertTrue(opener.requests[0].full_url.endswith(
            '/v2/broadinstitute/viral-ngs/manifests/feature-foo-core'))
        self.assertEqual(opener.requests[0].get_header('Authorization'), 'Bearer tok')

    def test_retries_timeouts_and_dropped_connections(self):
        # urlopen raises these unwrapped (not as URLError) while waiting for the response
        for exc in (TimeoutError('timed out'), http.client.RemoteDisconnected('closed'),
                    urllib.error.URLError('unreachable')):
            with self.subTest(exc=type(exc).__name__):
                opener = FakeOpener([exc, 202])
                quay = qd.Quay(REPO, 'tok', opener=opener, sleep=lambda s: None)
                self.assertEqual(quay.delete('feature-foo-core'), 202)
                self.assertEqual(len(opener.requests), 2)

    def test_no_response_after_retries_is_status_0(self):
        opener = FakeOpener([TimeoutError('timed out')] * 10)
        quay = qd.Quay(REPO, 'tok', opener=opener, sleep=lambda s: None)
        self.assertEqual(quay.delete('feature-foo-core'), 0)
        self.assertLessEqual(len(opener.requests), 4)

    def test_timeout_on_one_tag_does_not_skip_other_tags_or_canary(self):
        def opener(req, timeout=None):
            opener.requests.append((req.get_method(), req.full_url.rsplit('/', 1)[1]))
            if req.get_method() == 'DELETE' and req.full_url.endswith('/feature-foo-core'):
                raise TimeoutError('timed out')
            return FakeResponse(202 if req.get_method() == 'DELETE' else 200)
        opener.requests = []
        quay = qd.Quay(REPO, 'tok', opener=opener, sleep=lambda s: None)
        out = io.StringIO()
        rc = qd.run('feature/foo', quay, live_branches=['main'], execute=True, out=out)
        self.assertEqual(rc, 1)
        self.assertIn('feature-foo-core', out.getvalue())
        deletes = [t for m, t in opener.requests if m == 'DELETE']
        self.assertEqual(sorted(set(deletes)), sorted(qd.branch_tags('feature-foo')))
        last_delete = max(i for i, (m, _) in enumerate(opener.requests) if m == 'DELETE')
        canary_after = {t for m, t in opener.requests[last_delete:] if m == 'HEAD'}
        self.assertEqual(canary_after, set(qd.CANARY_TAGS))

    def test_rejects_invalid_tag_without_a_request(self):
        opener = FakeOpener([])
        quay = qd.Quay(REPO, 'tok', opener=opener, sleep=lambda s: None)
        for tag in ('main#cleanup', 'sha256:' + 'a' * 64, '', 'a/b'):
            with self.subTest(tag=tag):
                with self.assertRaises(ValueError):
                    quay.delete(tag)
                with self.assertRaises(ValueError):
                    quay.head(tag)
        self.assertEqual(opener.requests, [])

    def test_does_not_retry_client_errors(self):
        opener = FakeOpener([404])
        quay = qd.Quay(REPO, 'tok', opener=opener, sleep=lambda s: None)
        self.assertEqual(quay.head('nope'), 404)
        self.assertEqual(len(opener.requests), 1)

    def test_gives_up_after_retries(self):
        opener = FakeOpener([503] * 10)
        quay = qd.Quay(REPO, 'tok', opener=opener, sleep=lambda s: None)
        self.assertEqual(quay.delete('feature-foo-core'), 503)
        self.assertLessEqual(len(opener.requests), 4)


class TestAuthenticate(unittest.TestCase):

    def token_response(self, actions):
        return FakeResponse(200, json.dumps({'token': _jwt(actions)}).encode())

    def test_execute_requires_push(self):
        opener = FakeOpener([self.token_response(['pull'])])
        env = {'QUAY_USERNAME': 'robot', 'QUAY_TOKEN': 'secret'}
        with self.assertRaises(qd.AuthError):
            qd.authenticate(REPO, execute=True, env=env, opener=opener)

    def test_execute_requests_pull_push_with_basic_auth(self):
        opener = FakeOpener([self.token_response(['pull', 'push'])])
        env = {'QUAY_USERNAME': 'robot', 'QUAY_TOKEN': 'secret'}
        token = qd.authenticate(REPO, execute=True, env=env, opener=opener)
        self.assertEqual(token, _jwt(['pull', 'push']))
        req = opener.requests[0]
        self.assertIn('scope=repository:broadinstitute/viral-ngs:pull,push', req.full_url)
        self.assertNotIn('delete', req.full_url)
        expected = 'Basic ' + base64.b64encode(b'robot:secret').decode()
        self.assertEqual(req.get_header('Authorization'), expected)

    def test_execute_without_credentials_fails(self):
        with self.assertRaises(qd.AuthError):
            qd.authenticate(REPO, execute=True, env={}, opener=FakeOpener([]))

    def test_main_reports_token_timeout_cleanly(self):
        for exc in (TimeoutError('timed out'), http.client.RemoteDisconnected('closed')):
            with self.subTest(exc=type(exc).__name__):
                with mock.patch.object(qd, 'live_branches', return_value=['main']), \
                     mock.patch.object(qd, 'authenticate', side_effect=exc), \
                     mock.patch('sys.stdout', new_callable=io.StringIO) as out:
                    self.assertEqual(qd.main(['--branch', 'feature/foo', '--execute']), 1)
                self.assertIn(str(exc), out.getvalue())

    def test_dry_run_is_anonymous_pull(self):
        opener = FakeOpener([self.token_response(['pull'])])
        env = {'QUAY_USERNAME': 'robot', 'QUAY_TOKEN': 'secret'}
        qd.authenticate(REPO, execute=False, env=env, opener=opener)
        req = opener.requests[0]
        self.assertIn('scope=repository:broadinstitute/viral-ngs:pull', req.full_url)
        self.assertNotIn('push', req.full_url)
        self.assertIsNone(req.get_header('Authorization'))


if __name__ == '__main__':
    unittest.main()
