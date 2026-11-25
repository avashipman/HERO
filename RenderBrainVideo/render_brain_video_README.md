# Light sheet brain video

Makes a rotating 3D video or four-view still of fluorescent signal in a cleared
brain, with an autofluorescence outline and an optional highlighted bright region.

## Run

    pip install -r requirements.txt
    python render_brain_video.py

First edit **USER SETTINGS** at the top of `render_brain_video.py` (paths, voxel
size, orientation, output). Use `STILL = True` to tune, then `STILL = False`
and an `.mp4` output name.
A working OpenGL display is needed. Allow roughly 16-32 GB RAM for ~10 um volumes.

## Settings used

| Setting | BNST | Insula |
|---|---|---|
| Input | TIFF planes | Zarr, level 2 |
| VOXEL_UM (z, y, x) | (16, 14.4, 14.4) | (10, 10, 10) |
| TARGET_UM | 15 | 10 |
| DORSAL_FIRST | False | True |
| SIGNAL_RANGE | (250, 6000) | (250, 12000) |
| CORE_CUTOFF | 4000 | 12000 |
| Resulting core volume | 0.81 mm3 | 1.21 mm3 |

All other settings: defaults in the script. Videos used `ZOOM = 1.25`,
`SIZE = (1920, 1080)`; stills used `ZOOM = 1.0`, `SIZE = (1200, 900)` per view.
Raw voxel size was 4 x 1.8 x 1.8 um (z, y, x). The Insula Zarr metadata were not
available; its 10 um spacing was inferred from array dimensions and the raw voxel size.

## Units and interpretation

- Intensities are image counts (arbitrary units) from 16-bit files; uncalibrated.
- CORE_CUTOFF is applied after downsampling and 25 um (sigma) Gaussian smoothing.
  It was chosen by eye for each brain. The core is the largest connected region
  above it, drawn on top of other signal.
- Brains were processed differently and displayed over different ranges. Do not
  compare intensities or core volumes between them.
- The core is not, by itself, a verified injection-site boundary.
- Scale bars are not corrected for tissue expansion or shrinkage.
- Outline smoothing and gap filling can create surfaces that are not anatomy.

## AI assistance

The visualization code was written with Claude (Anthropic) under the authors'
direction. The authors reviewed the outputs, confirmed the injected hemisphere and
raw voxel size, and take responsibility for the code and its use.

## License

MIT; see `LICENSE`.

## Citation

TO BE COMPLETED BY THE AUTHORS.
