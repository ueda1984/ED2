# Read the Docs publishing workflow

This note describes the GitHub/Read the Docs separation checked in
November 2026. It does not claim that development-branch code has
already been released as upstream ED2.

## Source locations

| Location | Role |
| --- | --- |
| [QS-Cube/ED2](https://github.com/QS-Cube/ED2) | Upstream code and repository URL configured in the Read the Docs project |
| [ueda1984/ED2](https://github.com/ueda1984/ED2) | Development fork for Phase 1 and later engineering changes |
| [ed2.readthedocs.io](https://ed2.readthedocs.io/en/latest/) | Public documentation site |

MkDocs source lives in `docs/`, with `mkdocs.yml` for navigation
and `.readthedocs.yaml` for deployment settings. `latest` was
the only enabled public version when inspected.

A push to a branch in `ueda1984/ED2` **does not directly update**
the public site configured for `QS-Cube/ED2`.

## Suggested author-review-publish loop

1. Update code and the matching user documentation in a feature branch
   of the development fork. Keep documentation for unmerged features
   explicitly labeled *development only*.
2. Preview the site locally from the source tree:
   ```bash
   python3 -m pip install -r docs/requirements.txt
   python3 -m mkdocs build --strict
   python3 -m mkdocs serve
   ```
3. Review scientific signs, file formats, numerical claims and
   internal page links against the **same source SHA**.
4. After explicit approval, integrate the changes into the intended
   upstream branch. Do not update or merge `main` implicitly.
5. Trigger or verify a new build of the correct Read the Docs version;
   compare deployed content and displayed source revision with the
   source committed to GitHub.

Development-only material should not be presented as available in the
upstream executable until the corresponding code is integrated.

## Webhook and permissions

The Read the Docs Community project has an existing GitHub webhook
integration but its activity page displayed **no recent exchanges**
during the review. An empty activity log alone does not prove failure.

Read the Docs project maintenance and GitHub repository administration
are **separate permissions**. An upstream repository administrator
must verify GitHub Webhooks and recent delivery results. When a push
hook is not usable, a Read the Docs project maintainer can request a
manual build of the approved upstream revision using the Builds page.

Never commit, paste or screenshot webhook secrets or access tokens.
If a credential has been exposed, rotate it through the service
settings and update the corresponding GitHub integration.

## Publication gate

Before calling the current docs a stable release, check:

- Source SHA and docs site build revision agree.
- Installation and example commands run from a clean source archive.
- Output file extensions and column order match the actual source.
- Hamiltonian signs, symmetry claims and limitations match code.
- New solver controls are not advertised for older executables.
- A fresh-site build completes without missing MkDocs navigation
  targets or broken internal links.

`make distcheck` and end-to-end release validation are handled in
the separate Phase 2-A workstream, not established by this
documentation-only branch.
