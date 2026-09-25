# The board

This directory holds the messages that the sessions working on popnei at
the same time leave for each other. Each session works on its own branch
in its own worktree, so a file it commits is invisible to the others until
its branch is merged into `main`. The board is outside the branches for
that reason: this README is the only file of the directory in git, and
the messages are ignored by git and reachable from every worktree at the
same path, `/Users/jose/devel/popnei/.claude/board/`.

A message is one markdown file, named by the time it was posted, in UTC,
and the branch that posted it, so that `ls` lists the messages in the
order they were written:

    2026-09-22T1405-plan-dists-kosman.md

The file is under 15 lines. It opens with a sentence that says what it is
about, then says which branch posted it and what a reader on another
branch would do with it. A message is one of:

- **A start.** Which plan or spec the session works on, its branch, and
  the files that other sessions may also change: `docs/glossary.md`,
  `docs/architecture.md`, the `Variants` struct, a public function whose
  signature changes.
- **A finding.** Something learned that another session would otherwise
  find again: a behaviour of pyNei, a value that surprised, a trap of the
  build or of a tool.
- **A change to a shared file**, when it happens, with what changed.

A message is information for the session that reads it and not an order.
It weighs the message as it weighs what a reviewer reports: it decides
what to do with it, and when two branches will change one file it says
so to the owner, who orders the merges.

When a session's branch has been merged or dropped, its messages come
out. The session that merges removes its own, and a session that finds a
message whose branch no longer exists removes it and says so to the
owner. What was in a message and still matters goes where the next
reader looks for it: the report of the plan, the spec, or
`docs/rust_core.md`.

The rule for when a session reads and posts is in `CLAUDE.md`, for every
session that writes in the repository, and in the `following-plans` skill,
for the sessions that carry out a plan.
