# popnei: how the assistant works here

What popnei is and what it is for is in `docs/objectives.md`, the data flow
and the modules in `docs/architecture.md`, and the measurements behind the
decisions in `docs/rust_core.md`.

The skills are under `.claude/skills/` and the subagents under
`.claude/agents/`. The `writing` skill is read before anything a person
will read is written, and a document or a GitHub issue goes to the
`first-reader` subagent before it is handed over.

## Replies in chat

The owner reads a reply to decide something or to learn something they
need for their own work. They did not see the session. The principles of
the writing skill, `.claude/skills/writing/SKILL.md`, hold in chat, and
these are the ones that fail most often there:

- A reply opens with a sentence that says what it is about, also when it
  follows a question: what was asked or what was being worked on, with the
  options or the things named by what they are. The owner may come back to
  it after hours of other work. The answer comes next, in a full sentence.
  When several things are true, the one that changes what the owner does
  goes first.
- A longer reply is sketched before it is written, as the writing skill
  says, so that no term is used before the sentence that explains it.
- A reply holds the decisions that are needed from the owner and what they
  need to know. How the work went stays out, unless it changes what they
  do next.
- A reply carries the numbers the answer turns on and leaves the rest of
  what was measured for a document, with a line that says where it is or
  that it can be written. A reply that needs headings is usually a
  document.
- A request for a decision has the options, what each one costs, and a
  recommendation. What can be decided without the owner is decided, and is
  mentioned only when they need to know of it.
- No name before its explanation. Labels from the session, and the names
  of files, types and scripts, are not words the owner has, unless the
  owner used them first.
- A number where an adjective would go, with what it was measured on. A
  comparison has both sides in the same units and says which is better.
- Nothing about the reply itself, and nothing about how interesting or
  important a thing is.
- What failed, or was not done, is said as plainly as what worked.
