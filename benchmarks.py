#!/usr/bin/env python3
from mupq import mupq
from interface import parse_arguments, get_platform


BENCH_TARGETS = [
        '--stack', '--speed', '--hashing', '--size', '--speed', '--subspeed', '--all'
        ]

class SubSpeedBenchmark(mupq.StackBenchmark):
    test_type = 'speed_sub'

if __name__ == "__main__":
    args, rest = parse_arguments()
    platform, settings = get_platform(args)
    with platform:
        schemes = [s for s in rest if not s.startwith('--')]
        if not set(rest).intersection(set(BENCH_TARGETS)):
            print("No benchmark target selected.")
            print(f"Choose some of {', '.join(BENCH_TARGETS)}.")
        else
            if "--stack" in rest or '--all' in rest:
                test = mupq.StackBenchmark(settings, platform)
                test.test_all(schemes)
            if "--speed" in rest or '--all' in rest:
                test = mupq.SpeedBenchmark(settings, platform)
                test.test_all(schemes)
            if "--hashing" in rest or '--all' in rest:
                test = mupq.HashingBenchmark(settings, platform)
                test.test_all(schemes)
            if "--size" in rest or '--all' in rest:
                test = mupq.SizeBenchmark(settings, platform)
                test.test_all(schemes)
            if "--subspeed" in rest or '--all' in rest:
                test = SubSpeedBenchmark(settings, platform)
                test.test_all(schemes)
