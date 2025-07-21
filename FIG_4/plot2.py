import glob, re
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from collections import defaultdict

# 1. Build a list of all your data files
pattern = re.compile(r"c_(\d+)_(\d+)_([+-]?\d*\.?\d+)_([0-9eE\+\-\.]+)\.txt")
records = []

for fn in glob.glob("c_*.txt"):
    m = pattern.match(fn)
    if not m: 
        continue
    inp, dk, ph, coup = m.groups()
    inp, dk = int(inp), int(dk)
    ph, coup = float(ph), float(coup)
    sig = np.loadtxt(fn)[:,0]
    records.append({
        'input': inp,
        'deltak': dk,
        'phase': ph,
        'coupling': coup,
        'signal': sig
    })

# 2. Group by (deltak, phase, coupling)
groups = defaultdict(dict)
for rec in records:
    key = (rec['deltak'], rec['phase'], rec['coupling'])
    groups[key][rec['input']] = rec['signal']

# 3. Compute metrics for each group that has both input=1 and input=2
results = []
nskip = 1
dtN   = 25
dtE   = dtN/25
conv  = dtE/41.34137   # your energy→time factor

for (dk, ph, coup), sigs in groups.items():
    if 1 in sigs and 2 in sigs:
        A = sigs[1].copy();  B = sigs[2].copy()
        # normalize
        A = (A - A.mean())/A.std()
        B = (B - B.mean())/B.std()
        xcorr = correlate(A, B)
        nsamp = A.size
        xcorr /= nsamp

        # time‐lag axis
        lags = np.arange(1 - nsamp, nsamp) * nskip * conv
        i0   = xcorr.argmax()

        # metrics
        dt_rec    = lags[i0]
        phi_rec   = np.angle(xcorr[i0]*nsamp)   # multiply back to get complex value
        results.append({
            'Δk':          dk,
            'phase_in':    ph,
            'coupling':    coup,
            'time_shift':  dt_rec,
            'phase_shift': phi_rec
        })

# 4. Split into “no phonon” vs. “with phonon” and compare
no_phonon   = [r for r in results if r['coupling'] == 0]
with_phonon = [r for r in results if r['coupling'] != 0]

print("Comparison of recovered shifts (with phonon minus no phonon):")
print(" Δk │  Δtime (ps) │ Δphase (rad)")
print("────┼─────────────┼─────────────")
for dk in sorted({r['Δk'] for r in results}):
    r0 = next((r for r in no_phonon   if r['Δk']==dk), None)
    r1 = next((r for r in with_phonon if r['Δk']==dk), None)
    if r0 and r1:
        ddt  = r1['time_shift']  - r0['time_shift']
        dphi = r1['phase_shift'] - r0['phase_shift']
        print(f" {dk:>2} │ {ddt:>+9.4e} │ {dphi:>+9.4e}")

# 5. (Optional) Plot your results
# e.g. time_shift vs. coupling for each Δk
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
for dk in sorted({r['Δk'] for r in results}):
    pts = sorted([r for r in results if r['Δk']==dk], key=lambda x: x['coupling'])
    cs  = [r['coupling']   for r in pts]
    ts  = [r['time_shift'] for r in pts]
    ax.plot(cs, ts, marker='o', label=f"Δk={dk}")
ax.set(xlabel="phonon coupling", ylabel="recovered time shift (ps)")
ax.legend()
plt.show()