#!/usr/bin/env python3
"""Run a command with process-tree and system-memory safety limits.

The guard is intentionally dependency-free so it can protect environment setup,
native builds, and numerical runs before the project Python environment exists.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import platform
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def _system_free_percent() -> float | None:
    """Return an approximate available-memory percentage when supported."""
    if platform.system() == "Darwin":
        try:
            result = subprocess.run(
                ["memory_pressure", "-Q"],
                check=False,
                capture_output=True,
                text=True,
                timeout=3,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        match = re.search(r"free percentage:\s*([0-9.]+)%", result.stdout)
        return float(match.group(1)) if match else None

    meminfo = Path("/proc/meminfo")
    if meminfo.exists():
        values: dict[str, int] = {}
        for line in meminfo.read_text(encoding="utf-8").splitlines():
            key, _, raw_value = line.partition(":")
            match = re.search(r"\d+", raw_value)
            if match:
                values[key] = int(match.group())
        total = values.get("MemTotal")
        available = values.get("MemAvailable")
        if total and available is not None:
            return 100.0 * available / total
    return None


def _process_tree_rss_kib(root_pid: int) -> tuple[int, set[int]]:
    """Return aggregate RSS and PIDs for root_pid and all descendants."""
    try:
        result = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="],
            check=True,
            capture_output=True,
            text=True,
            timeout=3,
        )
    except (OSError, subprocess.SubprocessError):
        return 0, {root_pid}

    rows: list[tuple[int, int, int]] = []
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            rows.append(tuple(map(int, fields)))
        except ValueError:
            continue

    pids = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, ppid, _ in rows:
            if ppid in pids and pid not in pids:
                pids.add(pid)
                changed = True
    rss_kib = sum(rss for pid, _, rss in rows if pid in pids)
    return rss_kib, pids


def _timestamp() -> str:
    return dt.datetime.now(dt.UTC).isoformat(timespec="seconds")


def _emit(message: str, log_handle: object | None) -> None:
    line = f"[memory-guard {_timestamp()}] {message}"
    print(line, file=sys.stderr, flush=True)
    if log_handle is not None:
        print(line, file=log_handle, flush=True)


def _terminate_tree(
    process: subprocess.Popen[bytes], pids: set[int], grace_seconds: float
) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        pass
    for pid in sorted(pids, reverse=True):
        try:
            os.kill(pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            pass

    try:
        process.wait(timeout=grace_seconds)
        return
    except subprocess.TimeoutExpired:
        pass

    try:
        os.killpg(process.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass
    for pid in sorted(pids, reverse=True):
        try:
            os.kill(pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
    process.wait()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a command under process-tree and system-memory limits."
    )
    parser.add_argument("--max-rss-gb", type=float, default=12.0)
    parser.add_argument("--min-free-percent", type=float, default=12.0)
    parser.add_argument("--interval", type=float, default=2.0)
    parser.add_argument("--breach-samples", type=int, default=3)
    parser.add_argument("--grace-seconds", type=float, default=10.0)
    parser.add_argument("--status-every", type=float, default=30.0)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--log", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command[:1] == ["--"]:
        args.command = args.command[1:]
    if not args.command:
        parser.error("a command is required after --")
    if args.max_rss_gb <= 0 or args.interval <= 0 or args.breach_samples <= 0:
        parser.error("memory, interval, and breach limits must be positive")
    if not 0 <= args.min_free_percent <= 100:
        parser.error("--min-free-percent must be between 0 and 100")
    if args.threads <= 0:
        parser.error("--threads must be positive")
    return args


def main() -> int:
    args = _parse_args()
    log_handle = None
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        log_handle = args.log.open("a", encoding="utf-8")

    environment = os.environ.copy()
    environment["WAVELET_QC_MEMORY_GUARD"] = "1"
    for name in THREAD_ENV_VARS:
        environment[name] = str(args.threads)

    max_rss_kib = int(args.max_rss_gb * 1024 * 1024)
    command_display = " ".join(args.command)
    _emit(
        f"start: {command_display}; RSS cap={args.max_rss_gb:.2f} GiB; "
        f"minimum system free={args.min_free_percent:.1f}%; threads={args.threads}",
        log_handle,
    )

    process = subprocess.Popen(
        args.command,
        env=environment,
        start_new_session=True,
    )
    peak_rss_kib = 0
    consecutive_breaches = 0
    next_status = time.monotonic()

    try:
        while process.poll() is None:
            rss_kib, pids = _process_tree_rss_kib(process.pid)
            peak_rss_kib = max(peak_rss_kib, rss_kib)
            free_percent = _system_free_percent()
            over_rss = rss_kib > max_rss_kib
            under_free = (
                free_percent is not None
                and free_percent < args.min_free_percent
            )
            consecutive_breaches = (
                consecutive_breaches + 1 if over_rss or under_free else 0
            )

            now = time.monotonic()
            if now >= next_status:
                free_text = "unknown" if free_percent is None else f"{free_percent:.1f}%"
                _emit(
                    f"running: tree RSS={rss_kib / 1024 / 1024:.2f} GiB; "
                    f"system free={free_text}; processes={len(pids)}",
                    log_handle,
                )
                next_status = now + args.status_every

            if consecutive_breaches >= args.breach_samples:
                reason = (
                    f"tree RSS {rss_kib / 1024 / 1024:.2f} GiB exceeded cap"
                    if over_rss
                    else f"system free memory {free_percent:.1f}% fell below floor"
                )
                _emit(f"terminating: {reason}", log_handle)
                _terminate_tree(process, pids, args.grace_seconds)
                _emit(
                    f"stopped; peak tree RSS={peak_rss_kib / 1024 / 1024:.2f} GiB",
                    log_handle,
                )
                return 137

            time.sleep(args.interval)
    except KeyboardInterrupt:
        _, pids = _process_tree_rss_kib(process.pid)
        _emit("interrupted; terminating command tree", log_handle)
        _terminate_tree(process, pids, args.grace_seconds)
        return 130
    else:
        return_code = process.wait()
        _emit(
            f"exit={return_code}; "
            f"peak tree RSS={peak_rss_kib / 1024 / 1024:.2f} GiB",
            log_handle,
        )
        return return_code
    finally:
        if log_handle is not None:
            log_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
