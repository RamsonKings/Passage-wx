#!/usr/bin/env python3
"""
Passage weather relay.

Runs on GitHub Actions every 6 hours (and whenever route.txt changes). Downloads the latest
free forecasts, samples them along the route in route.txt and writes one compact file,
_site/relay.json, that the Passage Weather Planner reads in a single request.

Sources (all free):
  Currents  Copernicus Marine GLOBAL_ANALYSISFORECAST_PHY_001_024  merged-uv hourly (incl. tides)   ~10 days
  Waves     Copernicus Marine GLOBAL_ANALYSISFORECAST_WAV_001_027  MFWAM 3-hourly, wind sea + swell ~10 days
            ECMWF open data, IFS wave model (swh, mwd, mwp, pp1d)                               15 days
            NOAA GFS-Wave 0.25 (AWS open data)                                                  ~10 days used
  Wind/MSLP ECMWF open data IFS 0.25 (10u, 10v, gust, msl)                                       15 days
            ECMWF open data AIFS single (AI model; 10u, 10v, msl)                               15 days
            NOAA GFS 0.25 (AWS open data; 10 m u/v, gust, PRMSL)                                ~10 days used

Every source is optional: when one fails the others are still written, and the log says why.
"""
import datetime as dt, json, math, os, re, sys, time, traceback
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import requests

OUT_DIR = os.environ.get("RELAY_OUT", "_site")
WORK = os.environ.get("RELAY_WORK", "work")
NOW = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0)
LOG = []


def log(msg):
    line = f"{dt.datetime.now(dt.timezone.utc):%H:%M:%S} {msg}"
    print(line, flush=True)
    LOG.append(line)


# ---------------------------------------------------------------- route & points
def read_route(path="route.txt"):
    wps = []
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = [x.strip() for x in line.split(";")]
        try:
            wps.append((p[0], float(p[1]), float(p[2])))
        except (IndexError, ValueError):
            log(f"route: skipped line '{line}'")
    if len(wps) < 2:
        raise SystemExit("route.txt needs at least two waypoints: name; lat; lon")
    return wps


def gc_dist(a1, o1, a2, o2):
    p1, p2 = math.radians(a1), math.radians(a2)
    dp, dl = p2 - p1, math.radians(o2 - o1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * 3440.065 * math.asin(math.sqrt(min(1, a)))


def dlon(o1, o2):
    d = o2 - o1
    return d - 360 if d > 180 else d + 360 if d < -180 else d


def norm_lon(o):
    return ((o + 540) % 360) - 180


def points_along(wps, short_nm, long_nm):
    """same idea as the planner: dense in straits / short legs, sparser in open sea"""
    pts = [(wps[0][1], wps[0][2])]
    for (n1, a1, o1), (n2, a2, o2) in zip(wps, wps[1:]):
        d = gc_dist(a1, o1, a2, o2)
        n = max(1, math.ceil(d / (short_nm if d < 150 else long_nm)))
        for k in range(1, n + 1):
            f = k / n
            pts.append((round(a1 + (a2 - a1) * f, 3), round(norm_lon(o1 + dlon(o1, o2) * f), 3)))
    return pts


# ---------------------------------------------------------------- grid sampling
class Grid:
    """regular lat/lon grid; bilinear sampling at points (NaN-aware)"""

    def __init__(self, lats, lons, values):
        lats, lons, values = np.asarray(lats, float), np.asarray(lons, float), np.asarray(values, float)
        if lats[0] > lats[-1]:
            lats, values = lats[::-1], values[::-1, :]
        self.lat0, self.dlat, self.nlat = lats[0], (lats[-1] - lats[0]) / max(1, len(lats) - 1), len(lats)
        self.lon0, self.dlon, self.nlon = lons[0], (lons[-1] - lons[0]) / max(1, len(lons) - 1), len(lons)
        self.glob = abs(self.dlon * self.nlon - 360) < 1e-6
        self.v = values

    def at(self, lat, lon):
        y = (lat - self.lat0) / self.dlat
        x = ((lon - self.lon0) % 360) / self.dlon if self.glob else (lon - self.lon0) / self.dlon
        if not self.glob and x < 0:
            x = (lon + 360 - self.lon0) / self.dlon
        if y < 0 or y > self.nlat - 1 or x < 0 or (not self.glob and x > self.nlon - 1):
            return None
        y0, x0 = min(int(y), self.nlat - 2), int(x)
        fy, fx = y - y0, x - x0
        acc = wsum = 0.0
        for dy, dx, w in ((0, 0, (1 - fy) * (1 - fx)), (0, 1, (1 - fy) * fx), (1, 0, fy * (1 - fx)), (1, 1, fy * fx)):
            xx = (x0 + dx) % self.nlon if self.glob else min(x0 + dx, self.nlon - 1)
            val = self.v[y0 + dy, xx]
            if w > 0 and np.isfinite(val):
                acc += val * w
                wsum += w
        return acc / wsum if wsum >= 0.3 else None


def grib_messages(path):
    """yield (shortName, step_hours, Grid) for every GRIB message in a file (regular_ll only)"""
    import eccodes
    with open(path, "rb") as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                if eccodes.codes_get(gid, "gridType") != "regular_ll":
                    continue
                name = eccodes.codes_get(gid, "shortName")
                try:
                    step = int(eccodes.codes_get(gid, "endStep"))
                except Exception:
                    step = int(eccodes.codes_get(gid, "step"))
                ni, nj = eccodes.codes_get(gid, "Ni"), eccodes.codes_get(gid, "Nj")
                la1, la2 = eccodes.codes_get(gid, "latitudeOfFirstGridPointInDegrees"), eccodes.codes_get(gid, "latitudeOfLastGridPointInDegrees")
                lo1, lo2 = eccodes.codes_get(gid, "longitudeOfFirstGridPointInDegrees"), eccodes.codes_get(gid, "longitudeOfLastGridPointInDegrees")
                vals = eccodes.codes_get_values(gid).astype(float)
                if eccodes.codes_get(gid, "bitmapPresent"):
                    vals[vals == eccodes.codes_get(gid, "missingValue")] = np.nan
                if lo2 < lo1:
                    lo2 += 360
                yield name, step, Grid(np.linspace(la1, la2, nj), np.linspace(lo1, lo2, ni), vals.reshape(nj, ni))
            finally:
                eccodes.codes_release(gid)


def series_from_grids(grids, pts, scale):
    """grids: {step: Grid}; returns (steps, [[int per step] per point])"""
    steps = sorted(grids)
    out = []
    for la, lo in pts:
        row = []
        for s in steps:
            v = grids[s].at(la, lo)
            row.append(None if v is None or not math.isfinite(v) else int(round(v * scale)))
        out.append(row)
    return steps, out


# ---------------------------------------------------------------- ECMWF open data
ECMWF_SOURCES = ["ecmwf", "google", "azure", "aws"]   # data.ecmwf.int first; cloud mirrors as backup


def ecmwf_block(model, stream, params, steps, tag, pts, scales, rename):
    """latest 00z / 12z run (those go to 360 h; 06z / 18z stop at 144 h)"""
    from ecmwf.opendata import Client
    target = os.path.join(WORK, f"{tag}.grib2")
    base, last = None, None
    cand = NOW.replace(hour=12 if NOW.hour >= 12 else 0)
    for _ in range(4):
        # several mirrors: AWS S3 throttles ("503 Slow Down"); small chunks of steps, few retries, then the next mirror
        ok = False
        for source in ECMWF_SOURCES:
            try:
                if os.path.exists(target):
                    os.remove(target)
                for k in range(0, len(steps), 12):
                    part = target + ".part"
                    kw = dict(type="fc", step=steps[k:k + 12], param=params, target=part, date=cand.strftime("%Y%m%d"), time=cand.hour)
                    if stream:
                        kw["stream"] = stream
                    Client(source=source, model=model, maximum_retries=3, retry_after=15).retrieve(**kw)
                    with open(target, "ab") as out_f, open(part, "rb") as in_f:
                        out_f.write(in_f.read())
                    os.remove(part)
                ok = True
                log(f"{tag}: downloaded run {cand:%Y-%m-%d %HZ} from {source}")
                break
            except Exception as e:
                last = f"{source}: {e}"
                log(f"{tag}: {source} failed for run {cand:%Y-%m-%d %HZ} ({str(e)[:120]})")
        if ok:
            base = cand
            break
        cand -= dt.timedelta(hours=12)
    if base is None:
        raise RuntimeError(f"no {tag} run found ({last})")
    grids = {}
    for name, step, g in grib_messages(target):
        grids.setdefault(rename.get(name, name), {})[step] = g
    out = {"base": base.isoformat(), "vars": {}}
    for var, by_step in grids.items():
        if var not in scales:
            continue
        st, rows = series_from_grids(by_step, pts, scales[var])
        out["vars"][var] = {"steps": st, "scale": scales[var], "v": rows}
    os.remove(target)
    log(f"{tag}: run {out['base']}, vars {sorted(out['vars'])}, {len(next(iter(out['vars'].values()))['steps']) if out['vars'] else 0} steps")
    return out


def ifs_steps():
    return list(range(0, 145, 3)) + list(range(150, 361, 6))


def get_ecmwf(pts):
    res = {}
    sc = {"u10": 10, "v10": 10, "gust": 10, "msl": 0.1, "swh": 100, "mwd": 1, "mwp": 10, "pp1d": 10}
    rn = {"10u": "u10", "10v": "v10", "10fg": "gust", "fg10": "gust", "i10fg": "gust", "msl": "msl", "swh": "swh", "mwd": "mwd", "mwp": "mwp", "pp1d": "pp1d"}
    for tag, model, stream, params, steps in (
        ("ifs", "ifs", None, ["10u", "10v", "msl", "10fg"], ifs_steps()),
        ("ifs_wave", "ifs", "wave", ["swh", "mwd", "mwp", "pp1d"], ifs_steps()),
        ("aifs", "aifs-single", None, ["10u", "10v", "msl"], list(range(0, 361, 6))),
    ):
        try:
            res[tag] = ecmwf_block(model, stream, params, steps, tag, pts, sc, rn)
        except Exception as e:
            if tag == "ifs":     # gust can be missing at step 0 – retry without it
                try:
                    res[tag] = ecmwf_block(model, stream, params[:3], steps, tag, pts, sc, rn)
                    continue
                except Exception as e2:
                    e = e2
            log(f"{tag}: FAILED {e}")
    return res


# ---------------------------------------------------------------- NOAA GFS / GFS-Wave (AWS open data, byte ranges)
GFS_BASE = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
GFS_ATM = {"u10": ("UGRD", "10 m above ground"), "v10": ("VGRD", "10 m above ground"), "gust": ("GUST", "surface"), "msl": ("PRMSL", "mean sea level")}
GFS_WAV = {"swh": ("HTSGW", "surface"), "mwd": ("DIRPW", "surface"), "pp1d": ("PERPW", "surface"),
           "shww": ("WVHGT", "surface"), "mdww": ("WVDIR", "surface"), "mpww": ("WVPER", "surface"),
           "shts": ("SWELL", "1 in sequence"), "mdts": ("SWDIR", "1 in sequence"), "mpts": ("SWPER", "1 in sequence")}
GFS_SCALE = {"u10": 10, "v10": 10, "gust": 10, "msl": 0.1, "swh": 100, "mwd": 1, "pp1d": 10, "shww": 100, "mdww": 1, "mpww": 10, "shts": 100, "mdts": 1, "mpts": 10}
GFS_LAST = 240          # hours used from GFS / GFS-Wave (ECMWF covers to 360 h)
S = requests.Session()


def gfs_urls(cyc, fh, kind):
    d, h = cyc.strftime("%Y%m%d"), cyc.strftime("%H")
    if kind == "atm":
        return f"{GFS_BASE}/gfs.{d}/{h}/atmos/gfs.t{h}z.pgrb2.0p25.f{fh:03d}"
    return f"{GFS_BASE}/gfs.{d}/{h}/wave/gridded/gfswave.t{h}z.global.0p25.f{fh:03d}.grib2"


def latest_gfs_cycle(last_fh):
    cyc = NOW.replace(hour=NOW.hour // 6 * 6)
    for _ in range(6):
        r = S.head(gfs_urls(cyc, last_fh, "wave") + ".idx", timeout=30)
        if r.status_code == 200:
            return cyc
        cyc -= dt.timedelta(hours=6)
    raise RuntimeError("no complete GFS cycle found in the last 36 h")


def idx_ranges(url, wanted):
    lines = S.get(url + ".idx", timeout=60).text.strip().splitlines()
    rows = [l.split(":") for l in lines]
    found = {}
    for i, r in enumerate(rows):
        for key, (var, lev) in wanted.items():
            if key not in found and r[3] == var and r[4] == lev:
                start = int(r[1])
                end = int(rows[i + 1][1]) - 1 if i + 1 < len(rows) else ""
                found[key] = (start, end)
    return found


def gfs_field(url, rng):
    start, end = rng
    data = S.get(url, headers={"Range": f"bytes={start}-{end}"}, timeout=120).content
    path = os.path.join(WORK, f"g{abs(hash((url, start)))}.grib2")
    with open(path, "wb") as f:
        f.write(data)
    try:
        for _, _, g in grib_messages(path):
            return g
    finally:
        os.remove(path)


def get_gfs(pts_coarse):
    hours = list(range(0, GFS_LAST + 1, 3))
    cyc = latest_gfs_cycle(hours[-1])
    log(f"gfs: cycle {cyc:%Y-%m-%d %HZ}")
    out = {}
    for kind, wanted, tag in (("atm", GFS_ATM, "gfs"), ("wave", GFS_WAV, "gfs_wave")):
        grids = {k: {} for k in wanted}

        def job(fh):
            url = gfs_urls(cyc, fh, kind)
            try:
                for key, rng in idx_ranges(url, wanted).items():
                    grids[key][fh] = gfs_field(url, rng)
            except Exception as e:
                log(f"{tag} f{fh:03d}: {e}")

        with ThreadPoolExecutor(12) as ex:
            list(ex.map(job, hours))
        vars_ = {}
        for key, by in grids.items():
            if not by:
                continue
            st, rows = series_from_grids(by, pts_coarse, GFS_SCALE[key])
            vars_[key] = {"steps": st, "scale": GFS_SCALE[key], "v": rows}
        out[tag] = {"base": cyc.isoformat(), "vars": vars_}
        log(f"{tag}: vars {sorted(vars_)}")
    return out


# ---------------------------------------------------------------- Copernicus Marine (free account)
def boxes_for(pts, max_area=30.0, pad=0.15):
    """split the route into lat/lon boxes small enough to download quickly"""
    boxes, cur = [], []
    for p in pts:
        trial = cur + [p]
        la = [q[0] for q in trial]
        lo = [q[1] for q in trial]
        if cur and (max(la) - min(la) + 2 * pad) * (max(lo) - min(lo) + 2 * pad) > max_area:
            boxes.append(cur)
            cur = [cur[-1], p]
        else:
            cur = trial
    boxes.append(cur)
    out = []
    for b in boxes:
        la = [q[0] for q in b]
        lo = [q[1] for q in b]
        out.append((min(la) - pad, max(la) + pad, min(lo) - pad, max(lo) + pad))
    return out


def cmems_subset(dataset, variables, box, start, end, tag):
    import copernicusmarine as cm
    import xarray as xr
    s, n, w, e = box
    fn = f"{tag}_{abs(hash(box))}.nc"
    cm.subset(dataset_id=dataset, variables=variables, minimum_latitude=s, maximum_latitude=n,
              minimum_longitude=w, maximum_longitude=e, start_datetime=start.strftime("%Y-%m-%dT%H:%M:%S"),
              end_datetime=end.strftime("%Y-%m-%dT%H:%M:%S"), output_directory=WORK, output_filename=fn,
              overwrite=True, disable_progress_bar=True)
    ds = xr.open_dataset(os.path.join(WORK, fn))
    return ds, os.path.join(WORK, fn)


def sample_ds(ds, var, la, lo):
    da = ds[var]
    for dim in ("depth", "elevation"):
        if dim in da.dims:
            da = da.isel({dim: 0})
    lat_name = "latitude" if "latitude" in da.dims else "lat"
    lon_name = "longitude" if "longitude" in da.dims else "lon"
    try:
        v = da.interp({lat_name: la, lon_name: lo})
        arr = v.values
        if np.all(~np.isfinite(arr)):                 # coastal point: nearest sea cell
            arr = da.sel({lat_name: la, lon_name: lo}, method="nearest").values
        return arr
    except Exception:
        return None


def get_cmems(pts_fine, pts_coarse):
    if not os.environ.get("COPERNICUSMARINE_SERVICE_USERNAME"):
        log("copernicus: no account set (secrets CMEMS_USER / CMEMS_PASS) – skipped")
        return {}
    out = {}
    start = NOW - dt.timedelta(hours=3)
    end = NOW + dt.timedelta(days=10)
    jobs = (
        ("cur", "cmems_mod_glo_phy_anfc_merged-uv_PT1H-i", ["utotal", "vtotal"], pts_fine, {"utotal": 100, "vtotal": 100}),
        ("mfwam", "cmems_mod_glo_wav_anfc_0.083deg_PT3H-i",
         ["VHM0", "VMDR", "VHM0_WW", "VMDR_WW", "VTM01_WW", "VHM0_SW1", "VMDR_SW1", "VTM01_SW1"], pts_coarse,
         {"VHM0": 100, "VMDR": 1, "VHM0_WW": 100, "VMDR_WW": 1, "VTM01_WW": 10, "VHM0_SW1": 100, "VMDR_SW1": 1, "VTM01_SW1": 10}),
    )
    for tag, dataset, variables, pts, scales in jobs:
        try:
            rows = {v: [None] * len(pts) for v in variables}
            times = None
            for box in boxes_for(pts):
                ds, path = cmems_subset(dataset, variables, box, start, end, tag)
                tvals = ds["time"].values
                if times is None:
                    times = [str(np.datetime_as_string(t, unit="m")) + "Z" for t in tvals]
                s, n, w, e = box
                for i, (la, lo) in enumerate(pts):
                    if rows[variables[0]][i] is not None or not (s <= la <= n and w <= lo <= e):
                        continue
                    for v in variables:
                        arr = sample_ds(ds, v, la, lo)
                        if arr is None:
                            continue
                        arr = list(arr)[: len(times)]
                        rows[v][i] = [None if not np.isfinite(x) else int(round(float(x) * scales[v])) for x in arr]
                ds.close()
                os.remove(path)
            out[tag] = {"dataset": dataset, "t0": times[0] if times else None, "times": times,
                        "vars": {v: {"scale": scales[v], "v": rows[v]} for v in variables}}
            got = sum(1 for r in rows[variables[0]] if r)
            log(f"copernicus {tag}: {got}/{len(pts)} points, {len(times or [])} times")
        except Exception as e:
            log(f"copernicus {tag}: FAILED {e}")
            traceback.print_exc()
    return out


# ---------------------------------------------------------------- main
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    wps = read_route()
    fine = points_along(wps, 8, 20)        # currents (ETA)
    coarse = points_along(wps, 40, 40)     # wind / waves
    log(f"route: {len(wps)} waypoints, {len(fine)} current points, {len(coarse)} weather points")
    data = {"version": 1, "generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
            "route": [{"name": n, "lat": a, "lon": o} for n, a, o in wps],
            "fine": fine, "coarse": coarse, "models": {}}
    t0 = time.time()
    data["models"].update(get_cmems(fine, coarse))
    try:
        data["models"].update(get_ecmwf(coarse))
    except Exception as e:
        log(f"ecmwf: FAILED {e}")
    try:
        data["models"].update(get_gfs(coarse))
    except Exception as e:
        log(f"gfs: FAILED {e}")
    data["log"] = LOG
    data["seconds"] = round(time.time() - t0)
    path = os.path.join(OUT_DIR, "relay.json")
    with open(path, "w") as f:
        json.dump(data, f, separators=(",", ":"))
    log(f"wrote {path}: {os.path.getsize(path) / 1024:.0f} kB, models {sorted(data['models'])}")
    if not data["models"]:
        sys.exit("no source worked – see the log above")


if __name__ == "__main__":
    main()
