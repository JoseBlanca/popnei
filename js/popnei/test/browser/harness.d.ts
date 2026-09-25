/**
 * What `harness.html`, the page of the browser tests, gives the test that
 * opened it: one function that runs a case of `cases/` inside the web
 * worker of `worker.ts` and whose promise is rejected with the message of
 * the assertion that failed there.
 */

declare global {
  interface Window {
    runInTheWorker(name: string): Promise<void>;
  }
}

export {};
