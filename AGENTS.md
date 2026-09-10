# Palisade HDA II fork

This is a fork of [sunnypilot](https://github.com/sunnypilot/sunnypilot) that adds Hyundai
Palisade / Kia Telluride 2023-24 HDA II (CAN/CAN FD blended) support, which sunnypilot does not
carry. `opendbc_repo` points at [mudchobo/opendbc](https://github.com/mudchobo/opendbc), a fork of
sunnypilot/opendbc kept in the same shape.

| | repo | branch | upstream |
|---|---|---|---|
| openpilot | `mudchobo/openpilot` | `pal23-hda2` | `sunnypilot/sunnypilot` master |
| opendbc | `mudchobo/opendbc` | `pal23-hda2` | `sunnypilot/opendbc` master |

## Branch shape: upstream history, plus exactly one commit

**Both branches are always `upstream/master` plus a single squashed fork commit.**

```
451256b2 Palisade HDA II fork          <- ours, the only commit we own
6135084c modeld_v2: ... (#1993)        <- upstream/master, untouched
047ae41c modeld_v2: ... (#1990)
...
```

Why: `git diff upstream/master` is then the entire fork, upstream commits keep their original
hashes and GitHub signatures, and there are no merge commits of ours cluttering the log.

### Syncing

Do **not** `git merge upstream/master`. That leaves our merge commits in the history and grows
the log every sync. Instead, rebuild the single commit on top of the new upstream:

```bash
# 1. opendbc first, since openpilot's submodule pointer depends on its hash
cd ../opendbc
git fetch upstream master
git rebase upstream/master        # one commit, so conflicts are resolved once
# resolve, then run the checks below

# 2. openpilot
cd ../openpilot
git fetch upstream master
git rebase upstream/master
git -C opendbc_repo fetch ../opendbc pal23-hda2 && git -C opendbc_repo checkout <new opendbc sha>
git add opendbc_repo && git rebase --continue     # if the rebase stopped on the submodule
git submodule update --init                       # upstream often moves tinygrad/panda
```

If a rebase gets messy, the tree-preserving alternative is exact and cannot lose changes:

```bash
NEW=$(git commit-tree HEAD^{tree} -p upstream/master -F msg.txt)
git reset --hard $NEW
# verify: the only diff against the old tip should be opendbc_repo
git diff-tree -r --name-only $NEW <old tip>
```

This means **every sync force-pushes** (`git push --force-with-lease`). A device that already
pulled the branch needs `git fetch && git reset --hard origin/pal23-hda2 && git submodule update
--init --recursive`, or a fresh install.

Keep the squashed commit message as the fork's manifest: one section per feature, updated when a
feature is added or removed. It is the only record of why anything here differs from upstream.

## Gotchas that will bite

- **LFS is read-only.** `.lfsconfig` points at `gitlab.com/sunnypilot/public/sunnypilot-new-lfs.git`
  (~2 GB) over ssh, which this fork cannot write to. Pushing anything LFS-tracked (`*.otf`, `*.png`,
  `*.ttf`, `*.wav`, `*.svg`, `*.onnx`) fails with `git@gitlab.com: Permission denied (publickey)`.
  - For **upstream** LFS files pulled in by a sync: the object already sits in sunnypilot's public
    store, so push with `git push --no-verify` to skip the LFS pre-push hook. Clones fetch it from
    the read-only https endpoint.
  - For **our own** binaries: add a `!filter !diff !merge` exception in `.gitattributes` so they are
    stored as plain git objects. That is what the Kakao fonts do.
- **Never `git filter-branch` over a range containing upstream commits.** It strips their GitHub
  `gpgsig` and changes their hashes, forking them away from upstream.
- **Commits are authored as `mudchobo`**, set per-repo in `.git/config`
  (`1011214+mudchobo@users.noreply.github.com`). The global gitconfig is a work identity — do not
  let it leak into this repo.
- **Submodules do not follow a merge or rebase.** Run `git submodule update --init` afterwards;
  `openpilot/sunnypilot/models/tests` fails loudly when `tinygrad_repo` is stale.
- Some tests read blobs via `git show HEAD:...`, so they fail while a sync is uncommitted. Commit
  first, then judge test results.
- `opendbc/sunnypilot/car/car_list.json` is generated — run `opendbc/sunnypilot/car/platform_list.py`.
- `docs/CARS.md` in opendbc is stale relative to master; regenerating it adds unrelated noise, so
  leave it to sunnypilot's automation.

## Checks before pushing

```bash
# opendbc
cd ../opendbc
uv run --with ruff ruff check . && uv run ty check
opendbc/safety/tests/misra/test_misra.sh          # the hyundai.h safety changes trip MISRA easily
uv run unittest-parallel -j4

# openpilot
cd ../openpilot
uv run ruff check openpilot/
uv run scons -j8 openpilot/common/                # needed after params_keys.h changes
uv run python -m pytest openpilot/sunnypilot openpilot/selfdrive/ui/tests -q \
  --ignore=openpilot/sunnypilot/selfdrive/car/tests/test_custom_cruise.py \
  --ignore=openpilot/sunnypilot/selfdrive/controls/lib/dec/tests/test_dec_planner_gate.py
```

The two ignored tests need an acados build that isn't set up locally. `uv run` rewrites `uv.lock`
as a side effect — `git checkout uv.lock` before committing.

## Deliberate decisions

- **Torque Control Tune Version stays v0.0.** v0 (`LatControlTorqueV0`) keeps a derivative term
  that upstream's v1 dropped, sunnypilot itself defaults away from v1 ("FIXME-SP: tuning issues"),
  and the Friction Reduction tuning is baselined on v0.
- **The driver attention warning (DAW) is not ported.** It came from the original bryangerlach
  branch as a test UI. `ALERTS_364`/`DAW_Status` stays in the Palisade DBC as bus documentation,
  but nothing reads it and there is no `CarState.dawStatus`. Do not re-add it.
- **`radarUnavailable` is `True` on this platform** unless ESCC hardware is fitted, because the
  Palisade DBC has no `Bus.radar` entry. That disables the Hyundai jerk-limited integrator, pins
  the braking jerk limit at 5.0, and leaves lead tracking vision-only. Relevant to any
  longitudinal-behaviour complaint.
