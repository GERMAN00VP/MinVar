import argparse
import glob
import logging
import multiprocessing
from multiprocessing import Pool, Manager, Queue
from .main import run_minvar

# ------------------- Logger centralizado -------------------
def listener_configurer():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] [%(processName)s] %(message)s",
    )

def listener_process(queue):
    listener_configurer()
    while True:
        try:
            record = queue.get()
            if record is None:
                break
            logger = logging.getLogger()
            logger.handle(record)
        except Exception:
            import sys, traceback
            print("Logging listener failed:", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)

def worker_configurer(queue):
    h = logging.handlers.QueueHandler(queue)
    root = logging.getLogger()
    root.addHandler(h)
    root.setLevel(logging.INFO)

# ------------------- Wrapper para multiprocessing -------------------
def run_folder(args):
    folder, min_freq, gff, min_depth, min_base_quality = args
    return run_minvar(folder, min_freq, gff, min_depth, min_base_quality)

# ------------------- Main CLI -------------------
def main():

    parser = argparse.ArgumentParser(prog="MinVar")

    parser.add_argument("--input", required=True,
                        help="Path or pattern to IRMA folders (str:'path/pattern')")
    parser.add_argument("--min_freq", type=float, default=0.1,
                        help="Minimum frequency of reads for minority variant")
    parser.add_argument("--min_depth", type=int, default=20,
                        help="Minimum supporting reads for minority variant")
    parser.add_argument("--min_base_quality", type=float, default=25,
                        help="Minimum mean base quality for variant support")
    parser.add_argument("--gff", default=None,
                        help="Path to annotations file")
    parser.add_argument("--threads", type=int, default=1)

    args = parser.parse_args()

    folders = glob.glob(args.input)
    if not folders:
        raise FileNotFoundError(f"No folders matched the pattern: {args.input}")

    # ------------------- Setup logging con Queue -------------------
    manager = Manager()
    log_queue = manager.Queue()

    listener = multiprocessing.Process(target=listener_process, args=(log_queue,))
    listener.start()

    worker_args = [
        (f, args.min_freq, args.gff, args.min_depth, args.min_base_quality)
        for f in folders
    ]

    # ------------------- Ejecutar multiprocessing -------------------
    with Pool(args.threads, initializer=worker_configurer, initargs=(log_queue,)) as pool:
        results = pool.map(run_folder, worker_args)

    # ------------------- Cerrar listener -------------------
    log_queue.put(None)
    listener.join()

    # ------------------- Resumen final -------------------
    total_folders = len(results)
    success_count = sum(1 for r in results if r["status"] == "OK")
    warning_count = sum(1 for r in results if r["status"] == "WARNING")
    failed_count = total_folders - success_count - warning_count

    print("\n=== MinVar Summary ===")
    print(f"Total folders processed: {total_folders}")
    print(f"Successful: {success_count}")
    print(f"Warnings: {warning_count}")
    print(f"Failed: {failed_count}\n")

    for r in results:
        folder = r["folder"]
        status = r["status"]
        time_elapsed = r.get("time", 0)
        warnings = len(r.get("warnings", []))
        print(f"{folder}: {status} | Time: {time_elapsed:.2f}s | Warnings: {warnings}")
        if warnings:
            for w in r["warnings"]:
                print(f"  - {w}")

    print("\nPipeline completed.")

if __name__ == "__main__":
    main()
