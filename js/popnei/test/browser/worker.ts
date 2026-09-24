/**
 * The web worker of the browser tests: it loads the WebAssembly of popnei
 * once and runs the case the page asks it for.
 *
 * A case is a module of `cases/` that exports `run`, which reads what it
 * needs over http and throws an `Error` when a value is not the one the
 * node tests assert. The worker answers the page with the message of that
 * error, or with `null` when the case ran through. A test file of its own
 * is added to `cases/` and asked for by name, and nothing of the harness
 * changes.
 *
 * `FileReaderSync`, which reads a range of a file, exists only inside a
 * web worker, so every browser test of popnei runs in here.
 *
 * It is not a test file: node's test runner runs the files whose name ends
 * in `.test.ts`, and Playwright the ones whose name ends in `.browser.ts`.
 */

import { init } from "../../dist/web.js";

/** What this worker is, of the global object a module worker runs in. */
interface WorkerScope {
  postMessage(message: unknown): void;
  addEventListener(
    kind: "message",
    listener: (event: { data: unknown }) => void,
  ): void;
}

/** What the page asks for: the case to run, under the number of the ask. */
interface Ask {
  id: number;
  case: string;
}

const worker = globalThis as unknown as WorkerScope;

/** The WebAssembly, fetched once for every case this worker runs. */
const loading = init();

/**
 * Runs the case `name` of `cases/`.
 *
 * @throws {Error} When the name is not one a file of `cases/` could be
 * called, when that file exports no `run`, or when the case itself throws,
 * which is what an assertion of a case that failed does.
 */
async function runTheCase(name: string): Promise<void> {
  if (!/^[a-z][a-z0-9_]*$/.test(name)) {
    throw new Error(`\`${name}\` is no case of \`cases/\``);
  }
  await loading;
  const module: unknown = await import(`./cases/${name}.ts`);
  const run = (module as { run?: unknown }).run;
  if (typeof run !== "function") {
    throw new Error(`\`cases/${name}.ts\` exports no \`run\``);
  }
  await (run as () => Promise<void>)();
}

worker.addEventListener("message", (event) => {
  const ask = event.data as Ask;
  runTheCase(ask.case).then(
    () => {
      worker.postMessage({ id: ask.id, failure: null });
    },
    (error: unknown) => {
      const failure = error instanceof Error ? error.message : String(error);
      worker.postMessage({ id: ask.id, failure });
    },
  );
});
