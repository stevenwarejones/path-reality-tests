# Verified Dryad originals and download instructions

The authoritative source is [Dryad DOI 10.5061/dryad.x0k6djj14](https://doi.org/10.5061/dryad.x0k6djj14),
**version 3, version ID 450489**, published June 30, 2026. The
[exact version metadata](https://datadryad.org/api/v2/versions/450489) identifies
that revision. The deposited files are released under
[CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).

This repository includes byte-for-byte copies of the 19 verified originals in
`dataset/` to make automated checks independent of download availability. The
repository's Apache license does not replace their CC0 status. The deposit's own
`dataset/README.md` is preserved unchanged. This instruction file is separate.
Paper and supplement PDFs are not part of this redistribution.

## Verify the included files

From `studies/wen-2026-propagator/`, run:

```sh
python fetch_data.py --verify-only
```

The command must report **all 19 files as verified** and exit successfully. A
missing file, changed byte count, changed SHA-256 or unexpected dataset entry
fails verification. Stop and report any mismatch; do not modify the manifest or
silently replace an existing file to make a check pass.

## Obtain independent copies from Dryad

1. Open the [dataset page](https://doi.org/10.5061/dryad.x0k6djj14) in a browser
   and select the June 30, 2026 revision corresponding to version 3 / 450489.
2. Click **Download full dataset**, or download each listed file separately.
3. Unzip into `raw/dataset/` in a fresh checkout or a separate temporary folder,
   not directly into `raw/`. Keep filenames and bytes unchanged. Do not put the ZIP
   itself in `dataset/`; that folder must contain only the 19 expected originals.
4. Run the verification command above. For a separate folder, use
   `python fetch_data.py --verify-only --directory /path/to/downloaded/files`.
5. If any file fails, stop and report its name and the error. The committed
   manifest remains fixed to the audited version.

Scripted downloads may fail with 401/403 because of Dryad's access controls. A
normal browser download supplied the verified archive used here. Use the browser
route if necessary; do not try to bypass access controls. Automatic download is a
best-effort fallback, not a CI dependency. Downloading a newer revision does not
make it interchangeable with the pinned files.

## Expected files

Sizes are exact byte counts from the unchanged SHA-256 manifest.

| Filename | Bytes |
|---|---:|
| `Fig3.xlsx` | 50,582 |
| `Fig4A.xlsx` | 10,698 |
| `Fig4B.xlsx` | 11,615 |
| `Fig4C.xlsx` | 15,061 |
| `FigS1B-L.xlsx` | 142,424 |
| `FigS1B-Minus.xlsx` | 142,786 |
| `FigS1B-Plus.xlsx` | 142,855 |
| `FigS1B-PSI.xlsx` | 144,703 |
| `FigS1B-R.xlsx` | 142,739 |
| `FigS1C.xlsx` | 20,020 |
| `FigS2A.xlsx` | 16,297 |
| `FigS2B.xlsx` | 14,005 |
| `FigS2C.xlsx` | 15,937 |
| `FigS2D.xlsx` | 13,783 |
| `FigS3A.xlsx` | 11,676 |
| `FigS3B.xlsx` | 16,790 |
| `FigS4A.xlsx` | 13,495 |
| `FigS4B.xlsx` | 25,197 |
| `README.md` | 11,480 |

## Data credit

Wen, Yong-Li; Tian, Li-Man; Wang, Yunfei; Zhang, Shanchao; Li, Chang; Li, Jianfeng;
Wang, Enke; Yan, Hui; Zhu, Shi-Liang (2026). *Direct experimental test of Feynman's
path integral postulates with single photons* [Dataset]. Dryad.
[https://doi.org/10.5061/dryad.x0k6djj14](https://doi.org/10.5061/dryad.x0k6djj14).
Use of these data does not imply endorsement of this independent study.
