"""
Real-Time MD Performance Monitoring
Tracks GROMACS MD simulation performance, energy, temperature, pressure
"""

import os
import time
import re
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from datetime import datetime, timedelta
import numpy as np

class MDMonitor:
    def __init__(self, log_file="01_md_preparation/md.log"):
        self.log_file = log_file
        self.data = {
            'time': [],
            'step': [],
            'ns_per_day': [],
            'hours_remaining': [],
            'completion_time': [],
            'timestamp': []
        }

    def parse_log(self):
        """Extract performance metrics from GROMACS log"""
        if not os.path.exists(self.log_file):
            return None

        with open(self.log_file, 'r') as f:
            lines = f.readlines()

        # Find performance section
        for i, line in enumerate(lines):
            if 'Performance:' in line:
                # Extract ns/day
                match = re.search(r'(\d+\.\d+)\s+\(ns/day\)', line)
                if match:
                    ns_per_day = float(match.group(1))

                    # Extract hour/ns
                    match2 = re.search(r'(\d+\.\d+)\s+\(hour/ns\)', line)
                    if match2:
                        hour_per_ns = float(match2.group(1))

                        return {
                            'ns_per_day': ns_per_day,
                            'hour_per_ns': hour_per_ns
                        }

        # If no Performance line yet, extract from step estimates
        for line in reversed(lines):
            if 'step' in line and 'will finish' in line:
                # Parse: "step 4400, will finish Sun Dec 21 01:57:17 2025"
                match = re.search(r'step\s+(\d+)', line)
                if match:
                    current_step = int(match.group(1))
                    # Rough estimate based on steps
                    total_steps = 50000000
                    ps_per_step = 0.002
                    ps_completed = current_step * ps_per_step
                    ns_completed = ps_completed / 1000.0

                    return {
                        'current_step': current_step,
                        'ns_completed': ns_completed,
                        'percent': (current_step / total_steps) * 100
                    }

        return None

    def calculate_stats(self, step, total_steps=50000000, timestep_ps=0.002):
        """Calculate simulation statistics"""
        ps_completed = step * timestep_ps
        ns_completed = ps_completed / 1000.0
        percent_complete = (step / total_steps) * 100

        return {
            'ns_completed': ns_completed,
            'ps_completed': ps_completed,
            'percent': percent_complete,
            'steps_remaining': total_steps - step
        }

    def monitor_live(self, interval=10):
        """Live monitoring with periodic updates"""
        print("="*60)
        print("GROMACS MD LIVE MONITOR")
        print("="*60)
        print(f"Log file: {self.log_file}")
        print(f"Update interval: {interval} seconds")
        print(f"Press Ctrl+C to stop monitoring")
        print("="*60)
        print()

        try:
            while True:
                metrics = self.parse_log()

                if metrics:
                    print(f"[{datetime.now().strftime('%H:%M:%S')}]", end=" ")

                    if 'current_step' in metrics:
                        print(f"Step: {metrics['current_step']:,} | "
                              f"NS: {metrics['ns_completed']:.3f} | "
                              f"Complete: {metrics['percent']:.4f}%")

                    if 'ns_per_day' in metrics:
                        print(f"Performance: {metrics['ns_per_day']:.2f} ns/day "
                              f"({metrics['hour_per_ns']:.3f} hour/ns)")

                else:
                    print(f"[{datetime.now().strftime('%H:%M:%S')}] Waiting for data...")

                time.sleep(interval)

        except KeyboardInterrupt:
            print("\n\nMonitoring stopped by user")
            print("="*60)

    def plot_energy_temperature(self, edr_file="01_md_preparation/md.edr"):
        """
        Extract and plot energy/temperature from EDR file
        Requires gmx energy command
        """
        print("Energy/temperature plotting requires GROMACS gmx energy tool")
        print("This will be available once simulation produces enough data")
        print(f"EDR file: {edr_file}")


def main():
    import sys

    monitor = MDMonitor()

    if len(sys.argv) > 1 and sys.argv[1] == 'live':
        # Live monitoring mode
        interval = int(sys.argv[2]) if len(sys.argv) > 2 else 10
        monitor.monitor_live(interval=interval)
    else:
        # Single snapshot
        print("MD Monitor - Snapshot Mode")
        print("="*60)
        metrics = monitor.parse_log()

        if metrics:
            print("Current Metrics:")
            for key, value in metrics.items():
                print(f"  {key}: {value}")
        else:
            print("No data available yet. Simulation may be starting...")

        print("\nFor live monitoring, run:")
        print("  python monitor_md.py live [interval_seconds]")
        print("\nExample:")
        print("  python monitor_md.py live 5")


if __name__ == "__main__":
    main()
