#!/usr/bin/env python3
"""
render_brain_video.py

Rotating 3D video (or four-view still) of a fluorescent signal in a cleared
brain imaged by light sheet microscopy, with the brain outline taken from an
autofluorescence channel and, optionally, the largest bright region ("core")
highlighted in a second color.

HOW TO USE
  1. Install packages:   pip install -r requirements.txt
  2. Edit the USER SETTINGS below.
  3. Set STILL = True and run:   python render_brain_video.py
     Adjust the display settings until the still looks right.
  4. Set STILL = False and run again to make the video.

Each channel is ONE of:
  * a folder of single-plane 2D TIFFs (one file per z plane), or
  * one 3D level of a Zarr array (z, y, x). If the Zarr metadata file is
    missing, 128 x 128 x 128 uint16 Blosc chunks in z/y/x folders are ASSUMED.

Required axis order: z = dorsoventral, y = anteroposterior with anterior at
low y, x = mediolateral. Input files are only read, never changed.

INTENSITY UNITS
  All intensity settings are in image counts (arbitrary units, a.u.): the
  16-bit pixel values (0-65535) in the files provided. They are not calibrated
  and are affected by laser power, exposure, depth, stitching, destriping and
  downsampling. Do not compare them between brains without matched acquisition
  and processing. CORE_CUTOFF is applied after downsampling and Gaussian
  smoothing (sigma = 25 um), so it is not a cutoff on unsmoothed pixel values.

  The rendering is for visualization. The core is an intensity-defined region
  and does not by itself establish an injection-site boundary.
"""

# =============================================================================
# USER SETTINGS - edit this section
# =============================================================================

# --- Files -------------------------------------------------------------------
SIGNAL_PATH  = "/path/to/signal_channel"           # e.g. tdTomato, 561 nm
OUTLINE_PATH = "/path/to/autofluorescence_channel" # e.g. 488 nm
OUTPUT_PATH  = "brain_still.png"                   # .png for a still, .mp4 for a video
                                                   # (an existing file is overwritten)

# --- Image geometry ----------------------------------------------------------
VOXEL_UM     = None     # (z, y, x) um of the input files, e.g. (16, 14.4, 14.4)
TARGET_UM    = 15       # approximate rendered voxel size (um); larger = less memory
DORSAL_FIRST = False    # True if the first z plane is the dorsal side.
                        # Rotates 180 deg about y; does NOT mirror left/right.
                        # Check the hemisphere against your injection records.

# --- Output ------------------------------------------------------------------
STILL        = True     # True = four-view PNG; False = rotating video
FRAMES       = 360      # frames per full rotation
FPS          = 30
SIZE         = (1920, 1080)   # pixels (width, height), multiples of 8
ZOOM         = 1.0      # 1.25 was used for videos; 1.0 for stills
ELEVATION    = 20       # camera angle above horizontal (degrees)
SCALEBAR_UM  = 2000     # 0 hides it

# --- Signal display ----------------------------------------------------------
SIGNAL_RANGE   = (250, 10000)   # counts (a.u.): (invisible below, saturated above)
GAMMA          = 0.6            # <1 brightens faint signal
SIGNAL_OPACITY = 0.85
BLEND          = "maximum"      # "maximum" (crisp fibers) or "composite"

# --- Highlighted core --------------------------------------------------------
CORE_CUTOFF  = None     # counts (a.u.); None or 0 turns it off.
                        # Largest connected region above this value is drawn,
                        # always on top, even when behind other signal.
CORE_COLOR   = "#ffd633"
CORE_OPACITY = 0.85

# --- Brain outline -----------------------------------------------------------
OUTLINE_CUTOFF  = 5     # counts (a.u.) in the autofluorescence channel
OUTLINE_OPACITY = 0.2
OUTLINE_OPEN    = 5     # trims spikes (iterations at ~30 um)
OUTLINE_CLOSE   = 5     # fills gaps (iterations at ~30 um, maximum 13)

# =============================================================================
# End of user settings
# =============================================================================

import glob, os, re, sys
import imageio.v2 as imageio
import matplotlib
import numpy as np
import pyvista as pv
import tifffile
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi


def files(d):
    f = [p for p in glob.glob(os.path.join(d, "*")) if p.lower().endswith((".tif", ".tiff"))
         and not os.path.basename(p).startswith("._")]
    if not f: sys.exit(f"No TIFFs in {d}")
    return sorted(f, key=lambda p: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", os.path.basename(p))])

def load(fs, f, rotate=False):
    """Read TIFF planes and block-average by f = (z, y, x)."""
    if rotate: fs = fs[::-1]
    fz, fy, fx = f; out = []; shape = None
    for k in range(0, len(fs) // fz * fz, fz):
        acc = 0
        for p in fs[k:k + fz]:
            a = tifffile.imread(p)
            if shape is None: shape = a.shape
            if a.ndim != 2 or a.shape != shape: sys.exit(f"Plane {p} is {a.shape}, expected 2D {shape}")
            if rotate: a = a[:, ::-1]
            h, w = a.shape[0] // fy * fy, a.shape[1] // fx * fx
            acc = acc + a[:h, :w].astype(np.float32).reshape(h // fy, fy, w // fx, fx).mean((1, 3))
        out.append(acc / fz)
    return np.stack(out)

def is_zarr(d):
    return any(os.path.exists(os.path.join(d, n)) for n in (".zarray", "zarr.json")) or \
        (os.path.isdir(d) and bool(glob.glob(os.path.join(d, "0", "0", "0"))))

def load_zarr(d):
    """One 3D zarr level (z,y,x). Without metadata, assumes 128^3 uint16 Blosc chunks
    in z/y/x folders and trims all-zero trailing padding."""
    if any(os.path.exists(os.path.join(d, n)) for n in (".zarray", "zarr.json")):
        import zarr
        v = np.asarray(zarr.open(d, mode="r"))
    else:
        from numcodecs import Blosc
        print("  no zarr metadata in", d, "- assuming 128^3 uint16 chunks", flush=True)
        ps = glob.glob(os.path.join(d, "[0-9]*", "[0-9]*", "[0-9]*"))
        ix = [tuple(map(int, p.split(os.sep)[-3:])) for p in ps]
        if len(ix) != np.prod(np.max(ix, 0) + 1): print("  WARNING: some chunks are missing (shown as empty)", flush=True)
        v = np.zeros([(m + 1) * 128 for m in np.max(ix, 0)], np.uint16)
        for p, (i, j, k) in zip(ps, ix):
            raw = Blosc().decode(open(p, "rb").read())
            if len(raw) != 128 ** 3 * 2: sys.exit(f"Unexpected chunk size in {p}")
            v[i*128:(i+1)*128, j*128:(j+1)*128, k*128:(k+1)*128] = np.frombuffer(raw, "<u2").reshape(128, 128, 128)
        e = [np.flatnonzero(v.any(tuple(q for q in range(3) if q != ax)))[-1] + 1 for ax in range(3)]
        v = v[:e[0], :e[1], :e[2]]
    v = np.squeeze(v)
    if v.ndim != 3: sys.exit(f"{d} has shape {v.shape}; point at one 3D level")
    return v

def block(v, f):
    if tuple(f) == (1, 1, 1): return v
    n = [s // k for s, k in zip(v.shape, f)]
    v = v[:n[0]*f[0], :n[1]*f[1], :n[2]*f[2]]
    return v.reshape(n[0], f[0], n[1], f[1], n[2], f[2]).mean((1, 3, 5), dtype=np.float32)

def grid(v, sp, origin=(0, 0, 0)):
    g = pv.ImageData(dimensions=v.shape[::-1], spacing=sp[::-1], origin=tuple(origin)[::-1])
    g.point_data["v"] = v.ravel()      # z,y,x C-order == x-fastest
    g.set_active_scalars("v")
    return g

def main():
    for label, d in (("SIGNAL_PATH", SIGNAL_PATH), ("OUTLINE_PATH", OUTLINE_PATH)):
        if not os.path.isdir(d): sys.exit(f"{label} is not a folder: {d}  (edit USER SETTINGS)")
    if not VOXEL_UM or len(VOXEL_UM) != 3: sys.exit("Set VOXEL_UM = (z, y, x) in USER SETTINGS")
    if not OUTPUT_PATH.lower().endswith(".png" if STILL else ".mp4"):
        sys.exit("OUTPUT_PATH must end in .png if STILL is True, or .mp4 if it is False")
    if not 0 <= OUTLINE_CLOSE <= 13: sys.exit("OUTLINE_CLOSE must be 0-13")
    if SIGNAL_RANGE[1] <= SIGNAL_RANGE[0]: sys.exit("SIGNAL_RANGE must be (low, high)")
    if is_zarr(SIGNAL_PATH) != is_zarr(OUTLINE_PATH): sys.exit("Both channels must be the same format")
    f = tuple(max(1, round(TARGET_UM / v)) for v in VOXEL_UM); sp = [v * k for v, k in zip(VOXEL_UM, f)]
    fm = tuple(max(1, round(30 / v)) for v in VOXEL_UM); sm = [v * k for v, k in zip(VOXEL_UM, fm)]
    print("factors", f, "voxel um", sp, flush=True)
    if is_zarr(SIGNAL_PATH):
        inj, af = load_zarr(SIGNAL_PATH), load_zarr(OUTLINE_PATH)
        if inj.shape != af.shape: sys.exit(f"Channel shapes differ: {inj.shape} vs {af.shape}")
        if DORSAL_FIRST: inj, af = inj[::-1, :, ::-1], af[::-1, :, ::-1]   # 180 deg about y, not a mirror
        inj, af = block(inj, f), block(af, fm)
    else:
        fi, fa = files(SIGNAL_PATH), files(OUTLINE_PATH)
        if len(fi) != len(fa): sys.exit(f"Plane counts differ: {len(fi)} vs {len(fa)}")
        if tifffile.imread(fi[0]).shape != tifffile.imread(fa[0]).shape: sys.exit("Channel plane sizes differ")
        inj, af = load(fi, f, DORSAL_FIRST), load(fa, fm, DORSAL_FIRST)   # outline is always built at ~30 um
    print("shape", inj.shape, flush=True)
    m = ndi.gaussian_filter(af.astype(np.float32), 2) > OUTLINE_CUTOFF
    m = ndi.binary_opening(m, iterations=OUTLINE_OPEN) if OUTLINE_OPEN else m
    lab, n = ndi.label(m)
    if not n: sys.exit("Empty brain mask; lower OUTLINE_CUTOFF")
    m = lab == np.bincount(lab.ravel())[1:].argmax() + 1
    m = np.pad(m, 14)
    if OUTLINE_CLOSE: m = ndi.binary_closing(m, iterations=OUTLINE_CLOSE)
    m = ndi.binary_fill_holes(m)[14:-14, 14:-14, 14:-14]
    vol = m.sum() * np.prod(sm) / 1e9
    print(f"outline volume {vol:.0f} mm3", flush=True)
    idx = np.argwhere(m); lo, hi = idx.min(0), idx.max(0) + 1
    ms = ndi.gaussian_filter(np.pad(m, 4).astype(np.float32), 1.5)
    mg = grid(ms, sm)
    # voxel-center offset so both channels stay registered at different resolutions
    mg.origin = tuple(-4 * b + (b - s) / 2 for b, s in zip(sm[::-1], sp[::-1]))
    surf = mg.contour([0.5], scalars="v").smooth_taubin(n_iter=40, pass_band=0.05)
    surf = surf.compute_normals(auto_orient_normals=True, consistent_normals=True, split_vertices=False)
    core = None
    if CORE_CUTOFF:
        fc = tuple(max(1, round(20 / q)) for q in sp); sc = [q * k for q, k in zip(sp, fc)]
        g = ndi.gaussian_filter(block(inj, fc).astype(np.float32), 25.0 / np.array(sc))
        cl, cn = ndi.label(g > CORE_CUTOFF)
        if not cn: sys.exit("Nothing above CORE_CUTOFF")
        k = np.bincount(cl.ravel())[1:].argmax() + 1
        print(f"core >{CORE_CUTOFF:g}: {(cl == k).sum() * np.prod(sc) / 1e9:.2f} mm3, centroid z,y,x voxel "
              f"{[int(q) for q in ndi.center_of_mass(cl == k)]}", flush=True)
        near = ndi.binary_dilation(cl == k, iterations=2)
        core = grid(np.where(near, g, 0).astype(np.float32), sc, [(q - w) / 2 for q, w in zip(sc, sp)]).contour([CORE_CUTOFF], scalars="v").smooth_taubin(n_iter=20, pass_band=0.1)
        del g, cl, near
    # crop the signal to the brain (+0.3 mm) and convert plane by plane to save memory
    c0 = [max(0, int((l * b - 300) // q)) for l, b, q in zip(lo, sm, sp)]
    c1 = [min(n, int((h * b + 300) // q) + 1) for h, b, q, n in zip(hi, sm, sp, inj.shape)]
    w0, w1 = SIGNAL_RANGE
    iv = np.stack([(np.clip((inj[z, c0[1]:c1[1], c0[2]:c1[2]].astype(np.float32) - w0) / (w1 - w0), 0, 1) ** GAMMA * 255
                    ).astype(np.uint8) for z in range(c0[0], c1[0])])
    del inj, af
    print("rendered signal block", iv.shape, flush=True)
    pl = pv.Plotter(off_screen=True, window_size=SIZE); pl.set_background("black")
    asurf = pl.add_mesh(surf, color="lightgray", opacity=OUTLINE_OPACITY, smooth_shading=True, specular=0.2,
                culling="back")
    if core is not None:
        acore = pl.add_mesh(core, color=CORE_COLOR, smooth_shading=True, ambient=0.35, specular=0.2)
    cm = LinearSegmentedColormap.from_list("t", ["#8a0f0a", "#c4161c", "#f0202a", "#ff4038"] if CORE_CUTOFF else ["#8a0f0a", "#ff2d20", "#ff7a50", "#ffd9b0"])
    avol = pl.add_volume(grid(iv, sp, [c * q for c, q in zip(c0, sp)]), scalars="v", blending=BLEND, cmap=cm, clim=(0, 255), opacity=list(np.linspace(0, SIGNAL_OPACITY, 64)),
                  opacity_unit_distance=64, show_scalar_bar=False, shade=False)
    c = (lo + hi) / 2 * sm; ext = (hi - lo) * sm; R = 4e4
    pl.enable_parallel_projection()
    W, H = SIZE; aspect = W / H
    e = np.radians(ELEVATION)
    half_h = (ext[0] * np.cos(e) + np.hypot(ext[1], ext[2]) * np.sin(e)) / 2
    pl.camera.focal_point = c[::-1]
    ps = max(half_h, np.hypot(ext[1], ext[2]) / 2 / aspect) * 1.1 / ZOOM
    font = ImageFont.truetype(os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans.ttf"), H // 40)
    def shot(az, el=ELEVATION):
        t, e = np.radians(az), np.radians(el)
        d = np.array([np.sin(t) * np.cos(e), -np.cos(t) * np.cos(e), np.sin(e)])   # az 0 = from low y
        pl.camera.position = c[::-1] + R * d
        pl.camera.up = (0, 0, 1) if abs(el) < 89 else (0, -1, 0)
        pl.camera.parallel_scale = ps; pl.camera.clipping_range = (1e3, 8e4)
        def grab(others, c):
            asurf.SetVisibility(others); avol.SetVisibility(others)
            if core is not None: acore.SetVisibility(c)
            pl.render(); return pl.screenshot(return_img=True)[..., :3].astype(np.float32)
        img = grab(True, False)
        if core is not None:   # core is overlaid so it is never hidden or tinted by signal in front
            k = grab(False, True); al = (k.max(2, keepdims=True) > 8) * CORE_OPACITY
            img = img * (1 - al) + k * al
        im = Image.fromarray(img.astype(np.uint8)).resize((W, H))
        if SCALEBAR_UM:
            L = SCALEBAR_UM * H / (2 * ps); dr = ImageDraw.Draw(im); x, y = W * .05, H * .94
            dr.rectangle([x, y, x + L, y + H // 150], fill="white")
            dr.text((x + L / 2, y - 6), f"{SCALEBAR_UM / 1000:g} mm", font=font, fill="white", anchor="ms")
        return np.asarray(im)
    pl.show(auto_close=False)
    if STILL:
        s = [shot(0), shot(50), shot(90), shot(0, 89.9)]
        imageio.imwrite(OUTPUT_PATH, np.vstack([np.hstack(s[:2]), np.hstack(s[2:])])); return
    w = imageio.get_writer(OUTPUT_PATH, fps=FPS, codec="libx264", quality=8, macro_block_size=8)
    for i in range(FRAMES):
        w.append_data(shot(360 * i / FRAMES))
        if i % 60 == 0: print("frame", i, flush=True)
    w.close(); print("Done")


if __name__ == "__main__":
    main()
