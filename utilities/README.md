# Utility Scripts

**Purpose:** Helper scripts for monitoring and managing MD simulations

---

## Scripts

### monitor_md.py
**Use Case:** Real-time monitoring of running MD simulation

**What it does:**
- Connects to running Docker container
- Extracts performance metrics
- Displays progress, speed, ETA
- Checks for errors or crashes

**How to use:**

```bash
# For free peptide MD (if it were still running)
python utilities/monitor_md.py --container md_production_peptide

# For fusion protein MD (currently running)
python utilities/monitor_md.py --container fusion_md_production
```

**Output:**
```
=== MD Simulation Monitor ===
Container: fusion_md_production
Status: Running
Progress: 1.2 ns / 100 ns (1.2%)
Performance: 6.3 ns/day
ETA: 15.7 days
Temperature: 300.1 K
Last update: 2026-01-02 12:15:30
```

**Manual monitoring (without script):**

```bash
# Check if container is running
docker ps | grep fusion_md_production

# View latest logs
docker logs --tail 50 fusion_md_production

# Follow logs in real-time
docker logs -f fusion_md_production
```

---

## Monitoring Commands Reference

### Docker Container Status

```bash
# List all running containers
docker ps

# List all containers (including stopped)
docker ps -a

# Check specific container
docker ps -a --filter "name=fusion_md_production"
```

### View Logs

```bash
# Last 50 lines
docker logs --tail 50 fusion_md_production

# Last 100 lines
docker logs --tail 100 fusion_md_production

# Follow in real-time (Ctrl+C to stop)
docker logs -f fusion_md_production

# Search logs for specific info
docker logs fusion_md_production 2>&1 | grep "Performance"
docker logs fusion_md_production 2>&1 | grep "ns/day"
```

### Container Management

```bash
# Stop container
docker stop fusion_md_production

# Start stopped container
docker start fusion_md_production

# Restart container
docker restart fusion_md_production

# Remove container (WARNING: deletes progress!)
docker rm fusion_md_production
```

### Check Resource Usage

```bash
# Container CPU/memory usage
docker stats fusion_md_production

# All containers
docker stats
```

---

## Interpreting MD Progress

### Performance Metrics

**ns/day (nanoseconds per day):**
- Typical range: 5-20 ns/day depending on system size
- Free peptide: ~13.6 ns/day (small system)
- Fusion protein: ~6.3 ns/day (larger system)

**ETA Calculation:**
```
Target: 100 ns
Current speed: 6.3 ns/day
Time needed: 100 / 6.3 = 15.9 days
```

### Warning Signs

**Temperature fluctuations:**
- Normal: ± 3 K around target (300 K)
- Warning: > ± 10 K
- Action: Check logs for "LINCS warning" or "step too large"

**Crashes:**
- Look for: "Segmentation fault" or container stopped unexpectedly
- Check: `docker ps -a` to see if exited
- Recovery: Often need to restart from last checkpoint

**Slow performance:**
- Expected: Larger systems are slower
- Check: Docker resource allocation
- Optimize: Increase CPU cores if available

---

## Checkpoints and Recovery

### Understanding Checkpoints

GROMACS saves checkpoints periodically:
- File: `md.cpt` (checkpoint file)
- Contains: Full system state at specific time
- Purpose: Resume if simulation crashes

### Resuming from Checkpoint

If simulation crashes or is stopped:

```bash
# Restart from last checkpoint
docker run --rm \
  -v "path/to/md_folder:/work" \
  -w /work \
  gromacs/gromacs:latest \
  gmx mdrun -v -deffnm md -cpi md.cpt
```

The `-cpi md.cpt` flag tells GROMACS to continue from checkpoint.

---

## Expected Output Files During MD

As simulation runs, these files grow:

```
md.xtc       # Trajectory (coordinates over time)
md.edr       # Energy file
md.log       # Log file
md.cpt       # Checkpoint file
```

**File sizes:**
- Free peptide (23 ns): ~150 MB trajectory
- Fusion protein (100 ns, est.): ~5 GB trajectory

---

## Quick Health Check

**Is my simulation healthy?**

```bash
# Check container is running
docker ps | grep fusion_md_production

# Check last log entries
docker logs --tail 20 fusion_md_production

# Look for:
# - "Step XXXXX" (progress)
# - Performance metrics (ns/day)
# - No error messages
```

**Good signs:**
- Step number increasing
- Regular performance output
- Temperature near 300 K
- No warnings

**Bad signs:**
- Container stopped
- "LINCS WARNING" repeated
- "Step too large" errors
- No output for >1 hour

---

## Troubleshooting

### Container not found

```bash
# List all containers
docker ps -a

# Check if it exited
docker ps -a | grep fusion
```

### Can't access logs

```bash
# Make sure using correct container name
docker ps

# Try without filtering
docker logs <container_id>
```

### Simulation seems stuck

```bash
# Check if making progress
docker logs --tail 5 fusion_md_production

# Wait 5 minutes, check again
docker logs --tail 5 fusion_md_production

# Compare step numbers - should increase
```

### Want to pause simulation

```bash
# Stop container (can resume later)
docker stop fusion_md_production

# Resume later
docker start fusion_md_production
```

---

## Future Utility Scripts (Coming Soon)

**cluster_trajectory.py**
- Cluster conformations using RMSD
- Extract representative structures
- Prepare for RFDiffusion3

**extract_frames.py**
- Extract specific frames from trajectory
- Convert to PDB format
- For visualization or analysis

**compare_trajectories.py**
- Statistical comparison of free vs fusion
- Generate comparison plots
- Decision support for ensemble selection

---

**Last Updated:** January 2, 2026
