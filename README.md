# Passage weather relay

A free helper that runs on GitHub every 6 hours. It downloads the latest forecasts from the
data centres, samples them along the route in `route.txt`, and publishes one file,
`relay.json`, next to the Passage Weather Planner on your own GitHub Pages address.

On the ship the planner then downloads that one file (a few hundred kB) instead of making
hundreds of requests to Open-Meteo and Copernicus: no Starlink shared-IP limits, no time-outs.

| Data | Source | Range |
|---|---|---|
| Currents (incl. tides) | Copernicus Marine, merged-uv hourly | ~10 days |
| Waves, wind sea, swell | Copernicus MFWAM · ECMWF WAM · NOAA GFS-Wave | 10 / 15 / 10 days |
| Wind, gusts, MSLP | ECMWF IFS · ECMWF AIFS (AI) · NOAA GFS | 15 / 15 / 10 days |

The planner shows the median of the models and the low–high range ("Models" columns).

## One-time setup (about 20 minutes)

1. **GitHub account** – sign up free at github.com.
2. **New repository** – top right **+ → New repository**. Name: `passage-wx`. Choose **Public**
   (public repositories get free GitHub Actions and Pages). Create.
3. **Upload these files** – in the new repository: **Add file → Upload files**, drag in
   `relay.py`, `route.txt`, `index.html`, `README.md` and the folder `.github`
   (it contains `workflows/relay.yml`). **Commit changes**.
   If the `.github` folder does not upload by drag-and-drop: **Add file → Create new file**,
   type the name `.github/workflows/relay.yml`, paste the contents of that file, commit.
4. **Copernicus Marine account** – sign up free at marine.copernicus.eu (top right, Register).
5. **Store the Copernicus login as secrets** – in the repository: **Settings → Secrets and variables →
   Actions → New repository secret**. Add `CMEMS_USER` (your Copernicus username) and
   `CMEMS_PASS` (your Copernicus password). Nobody can read them back, not even in the logs.
6. **Turn on Pages** – **Settings → Pages → Build and deployment → Source: GitHub Actions**.
7. **First run** – **Actions** tab → **Weather relay** → **Run workflow**. It takes 15–40 minutes.
   A green tick means done. Your planner is then at
   `https://<your-github-name>.github.io/passage-wx/` and the data at `…/passage-wx/relay.json`.

## Each voyage

1. In the planner: **Route waypoints → Copy route for relay**.
2. On GitHub open `route.txt` → pencil icon (Edit) → replace the waypoint lines → **Commit changes**.
3. The relay rebuilds by itself (about 20–30 minutes). After that it refreshes every 6 hours.

## Checking it

* Planner status line after Run forecast: `relay 1.4 h old · …` means the relay was used.
* **Copy diagnostics** in the planner includes the relay's own log (which sources worked).
* On GitHub: **Actions** tab shows each run; a red cross opens the log with the reason.
* If the relay is older than 13 h, or the route differs, the planner falls back to the direct
  sources automatically and says so in the alerts.

Note: GitHub pauses scheduled runs in a repository with no activity for 60 days. Editing
`route.txt` for a new voyage counts as activity; otherwise press **Run workflow** once.

Data licences: Copernicus Marine (free, attribution), ECMWF open data (CC-BY-4.0),
NOAA (public domain). For personal planning; not a substitute for official warnings.
