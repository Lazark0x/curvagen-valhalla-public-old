// See valhalla/thor/road_twin_index.h.
#include "thor/road_twin_index.h"

#include "baldr/graphconstants.h"
#include "midgard/constants.h"
#include "midgard/logging.h"
#include "midgard/pointll.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <unordered_map>

using namespace valhalla::baldr;
using namespace valhalla::midgard;

namespace valhalla {
namespace thor {

namespace {

// Shape sampling pitch.  30 m is the twin radius, so a 25 m pitch keeps a sample inside
// every twin-radius disc along the segment; the cap bounds the memory of a long rural
// edge (a 5 km edge would otherwise carry 200 samples).
constexpr double kSampleSpacingM = 25.0;
constexpr uint32_t kMaxSamplesPerEdge = 12;
// "over most of the shorter edge": the share of one segment's samples that must sit
// inside the radius of the other segment's polyline.
constexpr double kCoverFraction = 0.70;
// The share of the SHORTER segment's length that the matched run must span along the
// other segment.  Without it a continuation reads as a twin (see run() below).
constexpr double kSpanFraction = 0.60;
// Collinearity tolerance on the chord bearings (the ticket's "~180 deg +/- tolerance").
constexpr double kBearingTolDeg = 30.0;
// proto/v4-p1.1 SWITCHBACK TEST (curvagen-valhalla#12 item 1).  Planimetric geometry
// cannot tell a mountain hairpin pair from a dual carriageway (P1 report §7.6: block C
// paid for it with is_lollipop c0.8 1.30 -> 3.91 % and vlasina-50 km 48.3 -> 62.7 %).
// Two physical facts separate them:
//   * a dual carriageway is ALWAYS two different OSM ways (each carriageway is its own
//     one-way way); the two arms either side of a hairpin are usually one way;
//   * a carriageway pair holds a roughly CONSTANT lateral offset along the matched run,
//     while a hairpin's arms converge to ~0 at the apex and diverge away from it.
// A pair qualifies only if it is a different way OR its offset is constant.
// "Constant" = the covered samples' offsets span at most kOffsetSpreadM, or the
// smallest is at least kOffsetConstFrac of the largest (which admits a wide-median
// motorway whose separation swings 20-60 m but never collapses).
constexpr double kOffsetSpreadM = 12.0;
constexpr double kOffsetConstFrac = 0.40;
// Spatial bin pitch; must be >= the widest radius so a 1-ring neighbourhood is complete.
constexpr double kCellM = 120.0;
// Serbia-scale constant for the lon/metre conversion of the BIN GRID only (distances
// themselves use the pair's own latitude).  cos(46 deg) keeps every cell >= kCellM in
// the north of the tileset, which is what the 1-ring completeness argument needs.
constexpr double kGridCosLat = 0.695;

struct Pt {
  int32_t lat_e6;
  int32_t lon_e6;
};

inline double norm_deg(double d) {
  while (d < 0.0)
    d += 360.0;
  while (d >= 360.0)
    d -= 360.0;
  return d;
}

// Angle between two undirected chord bearings, 0..90 after folding at 180 (collinear
// covers both the 0 deg and the 180 deg presentations — see the header's note).
inline double collinear_angle(double b1, double b2) {
  double d = std::fabs(norm_deg(b1) - norm_deg(b2));
  if (d > 180.0)
    d = 360.0 - d;
  return d > 90.0 ? 180.0 - d : d;
}

// Squared distance, metres, from p to the segment ab in a local flat frame, plus the
// clamped projection parameter (needed for the along-run span test below).
inline double seg_dist2(double px, double py, double ax, double ay, double bx, double by,
                        double& t_out) {
  const double vx = bx - ax, vy = by - ay;
  const double wx = px - ax, wy = py - ay;
  const double vv = vx * vx + vy * vy;
  double t = vv > 0.0 ? (wx * vx + wy * vy) / vv : 0.0;
  t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
  t_out = t;
  const double dx = wx - t * vx, dy = wy - t * vy;
  return dx * dx + dy * dy;
}

std::mutex& registry_mutex() {
  static std::mutex m;
  return m;
}

// ---- persisted form ----------------------------------------------------------------
// One header, then the arrays back to back in this order, each exactly its length:
//   twin_keys u64, twin_vals u64, par_keys u64, par_vals u64, filter u64,
//   twin_off u32 (twin_keys + 1), par_off u32 (par_keys + 1).
// Every u64 array starts 8-aligned because the header is 8-aligned and the u32 arrays
// come last.  Host byte order; `endian` refuses a file from the other order.
constexpr char kMagic[8] = {'C', 'V', 'G', 'R', 'T', 'W', 'I', 'N'};
constexpr uint32_t kFormatVersion = 1;
constexpr uint32_t kEndianProbe = 0x01020304u;
constexpr uint32_t kFlagParallels = 1u;
constexpr uint32_t kFlagSwitchback = 2u;
// Bump when Build's logic changes in a way the constants below do not show.
constexpr uint64_t kBuildRulesVersion = 1;

struct FileHeader {
  char magic[8];
  uint32_t format;
  uint32_t endian;
  uint64_t rules;   // build_rules_fingerprint()
  uint64_t tileset; // RoadTwinIndex::tileset_fingerprint()
  uint64_t tiles;
  double twin_radius_m;
  double parallel_radius_m; // 0 when parallels are not built
  uint32_t flags;
  uint32_t reserved;
  uint64_t twin_keys, twin_vals, par_keys, par_vals, filter_words;
};
static_assert(sizeof(FileHeader) % 8 == 0, "u64 arrays follow the header");

inline uint64_t mix64(uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ull;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebull;
  return x ^ (x >> 31);
}
inline uint64_t fold(uint64_t h, uint64_t v) {
  return mix64(h ^ mix64(v + 0x9e3779b97f4a7c15ull));
}
inline uint64_t bits_of(double d) {
  uint64_t u;
  std::memcpy(&u, &d, sizeof u);
  return u;
}

// Every rule constant Build reads, so a persisted index from other rules is refused.
uint64_t build_rules_fingerprint(uint64_t filter_bits) {
  uint64_t h = fold(0, kBuildRulesVersion);
  for (double d : {kSampleSpacingM, kCoverFraction, kSpanFraction, kBearingTolDeg, kOffsetSpreadM,
                   kOffsetConstFrac, kCellM, kGridCosLat})
    h = fold(h, bits_of(d));
  h = fold(h, kMaxSamplesPerEdge);
  return fold(h, filter_bits);
}

uint32_t flags_of(bool build_parallels, bool switchback_test) {
  return (build_parallels ? kFlagParallels : 0u) | (switchback_test ? kFlagSwitchback : 0u);
}

uint64_t payload_bytes(const FileHeader& h) {
  return 8 * (h.twin_keys + h.twin_vals + h.par_keys + h.par_vals + h.filter_words) +
         4 * (h.twin_keys + 1 + h.par_keys + 1);
}

} // namespace

struct RoadTwinIndex::Mapping {
  void* addr = nullptr;
  size_t len = 0;
  ~Mapping() {
    if (addr)
      munmap(addr, len);
  }
};

void RoadTwinIndex::point_views_at_vectors() {
  twin_filter_ = {own_filter_.data(), own_filter_.size()};
  twin_keys_ = {own_twin_keys_.data(), own_twin_keys_.size()};
  twin_vals_ = {own_twin_vals_.data(), own_twin_vals_.size()};
  twin_off_ = {own_twin_off_.data(), own_twin_off_.size()};
  par_keys_ = {own_par_keys_.data(), own_par_keys_.size()};
  par_vals_ = {own_par_vals_.data(), own_par_vals_.size()};
  par_off_ = {own_par_off_.data(), own_par_off_.size()};
}

uint64_t RoadTwinIndex::tileset_fingerprint(GraphReader& reader) {
  const auto set = reader.GetTileSet();
  std::vector<GraphId> ids(set.begin(), set.end());
  std::sort(ids.begin(), ids.end(),
            [](const GraphId& a, const GraphId& b) { return a.value < b.value; });
  uint64_t h = fold(0, ids.size());
  for (const auto& id : ids) {
    h = fold(h, id.value);
    graph_tile_ptr tile = reader.GetGraphTile(id);
    if (!tile)
      continue;
    h = fold(h, tile->header()->tile_checksum());
    h = fold(h, tile->header()->dataset_id());
    h = fold(h, tile->header()->end_offset());
  }
  return h;
}

std::unique_ptr<RoadTwinIndex> RoadTwinIndex::Scan(GraphReader& reader,
                                                   double twin_radius_m,
                                                   double parallel_radius_m,
                                                   bool build_parallels,
                                                   bool switchback_test) {
  auto idx = std::unique_ptr<RoadTwinIndex>(new RoadTwinIndex());
  idx->Build(reader, twin_radius_m, parallel_radius_m, build_parallels, switchback_test);
  return idx;
}

bool RoadTwinIndex::Save(const std::string& path, GraphReader& reader) const {
  FileHeader h{};
  std::memcpy(h.magic, kMagic, sizeof kMagic);
  h.format = kFormatVersion;
  h.endian = kEndianProbe;
  h.rules = build_rules_fingerprint(kFilterBits);
  h.tileset = tileset_fingerprint(reader);
  h.tiles = reader.GetTileSet().size();
  h.twin_radius_m = twin_radius_m_;
  h.parallel_radius_m = parallel_radius_m_;
  h.flags = flags_of(parallel_radius_m_ > 0.0, switchback_test_);
  h.twin_keys = twin_keys_.n;
  h.twin_vals = twin_vals_.n;
  h.par_keys = par_keys_.n;
  h.par_vals = par_vals_.n;
  h.filter_words = twin_filter_.n;

  // Written beside the target and renamed over it, so a reader never maps half a file.
  const std::string tmp = path + ".tmp";
  FILE* f = std::fopen(tmp.c_str(), "wb");
  if (!f)
    return false;
  bool ok = std::fwrite(&h, sizeof h, 1, f) == 1;
  auto put = [&](const auto& v) {
    if (ok && v.n)
      ok = std::fwrite(v.p, sizeof(*v.p), v.n, f) == v.n;
  };
  put(twin_keys_);
  put(twin_vals_);
  put(par_keys_);
  put(par_vals_);
  put(twin_filter_);
  put(twin_off_);
  put(par_off_);
  ok = (std::fclose(f) == 0) && ok;
  if (ok)
    ok = std::rename(tmp.c_str(), path.c_str()) == 0;
  if (!ok)
    std::remove(tmp.c_str());
  return ok;
}

std::unique_ptr<RoadTwinIndex> RoadTwinIndex::Load(const std::string& path,
                                                   GraphReader& reader,
                                                   double twin_radius_m,
                                                   double parallel_radius_m,
                                                   bool build_parallels,
                                                   bool switchback_test,
                                                   std::string* why) {
  auto refuse = [&](const std::string& reason) {
    if (why)
      *why = reason;
    return std::unique_ptr<RoadTwinIndex>();
  };
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    return refuse("cannot open " + path);
  struct stat st{};
  if (::fstat(fd, &st) != 0 || static_cast<uint64_t>(st.st_size) < sizeof(FileHeader)) {
    ::close(fd);
    return refuse("shorter than its header");
  }
  auto map = std::make_shared<Mapping>();
  map->len = static_cast<size_t>(st.st_size);
  void* addr = ::mmap(nullptr, map->len, PROT_READ, MAP_PRIVATE, fd, 0);
  ::close(fd);
  if (addr == MAP_FAILED)
    return refuse("mmap failed");
  map->addr = addr;

  FileHeader h;
  std::memcpy(&h, addr, sizeof h);
  if (std::memcmp(h.magic, kMagic, sizeof kMagic) != 0)
    return refuse("not a road identity file");
  if (h.endian != kEndianProbe)
    return refuse("other byte order");
  if (h.format != kFormatVersion)
    return refuse("format " + std::to_string(h.format) + ", this engine reads " +
                  std::to_string(kFormatVersion));
  if (h.rules != build_rules_fingerprint(kFilterBits))
    return refuse("built under other build rules");
  const double par = build_parallels ? parallel_radius_m : 0.0;
  if (h.twin_radius_m != twin_radius_m || h.parallel_radius_m != par ||
      h.flags != flags_of(build_parallels, switchback_test))
    return refuse("built with other radii or switchback setting");
  if (h.filter_words != kFilterBits / 64)
    return refuse("prefilter of the wrong size");
  if (sizeof h + payload_bytes(h) != map->len)
    return refuse("size does not match its header");
  if (h.tiles != reader.GetTileSet().size() || h.tileset != tileset_fingerprint(reader))
    return refuse("built from other tiles");

  auto idx = std::unique_ptr<RoadTwinIndex>(new RoadTwinIndex());
  const char* p = static_cast<const char*>(addr) + sizeof h;
  auto take64 = [&](View<uint64_t>& v, uint64_t n) {
    v = {reinterpret_cast<const uint64_t*>(p), static_cast<size_t>(n)};
    p += 8 * n;
  };
  auto take32 = [&](View<uint32_t>& v, uint64_t n) {
    v = {reinterpret_cast<const uint32_t*>(p), static_cast<size_t>(n)};
    p += 4 * n;
  };
  take64(idx->twin_keys_, h.twin_keys);
  take64(idx->twin_vals_, h.twin_vals);
  take64(idx->par_keys_, h.par_keys);
  take64(idx->par_vals_, h.par_vals);
  take64(idx->twin_filter_, h.filter_words);
  take32(idx->twin_off_, h.twin_keys + 1);
  take32(idx->par_off_, h.par_keys + 1);
  if (idx->twin_off_.p[h.twin_keys] != h.twin_vals || idx->par_off_.p[h.par_keys] != h.par_vals)
    return refuse("offsets do not end at their values");
  idx->mapping_ = std::move(map);
  idx->twin_radius_m_ = h.twin_radius_m;
  idx->parallel_radius_m_ = h.parallel_radius_m;
  idx->switchback_test_ = switchback_test;
  idx->stats_.tiles = h.tiles;
  idx->stats_.twin_keys = h.twin_keys;
  idx->stats_.parallel_keys = h.par_keys;
  idx->stats_.twin_radius_m = h.twin_radius_m;
  idx->stats_.parallel_radius_m = h.parallel_radius_m;
  return idx;
}

bool RoadTwinIndex::same_contents(const RoadTwinIndex& o) const {
  auto eq = [](const auto& a, const auto& b) {
    return a.n == b.n && std::equal(a.begin(), a.end(), b.begin());
  };
  return eq(twin_filter_, o.twin_filter_) && eq(twin_keys_, o.twin_keys_) &&
         eq(twin_vals_, o.twin_vals_) && eq(twin_off_, o.twin_off_) && eq(par_keys_, o.par_keys_) &&
         eq(par_vals_, o.par_vals_) && eq(par_off_, o.par_off_);
}

const RoadTwinIndex& RoadTwinIndex::get(GraphReader& reader,
                                        double twin_radius_m,
                                        double parallel_radius_m,
                                        bool build_parallels,
                                        bool switchback_test,
                                        const std::string& persisted_path) {
  // Keyed by tileset location + radii + file: gurka builds many tilesets in one process.
  static std::map<std::string, std::unique_ptr<RoadTwinIndex>> registry;
  const std::string key = reader.GetTileSetLocation() + "|" + std::to_string(twin_radius_m) + "|" +
                          std::to_string(build_parallels ? parallel_radius_m : 0.0) + "|" +
                          (switchback_test ? "sb" : "nosb") + "|" + persisted_path;
  std::lock_guard<std::mutex> lock(registry_mutex());
  auto it = registry.find(key);
  if (it != registry.end())
    return *it->second;
  std::unique_ptr<RoadTwinIndex> idx;
  if (!persisted_path.empty()) {
    const auto t0 = std::chrono::steady_clock::now();
    std::string why;
    idx = Load(persisted_path, reader, twin_radius_m, parallel_radius_m, build_parallels,
               switchback_test, &why);
    if (idx) {
      idx->stats_.build_ms =
          std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
      LOG_INFO("roundtrip road identity sidecar: mapped " + persisted_path + " in " +
               std::to_string(static_cast<int>(idx->stats_.build_ms)) + " ms — " +
               std::to_string(idx->stats_.tiles) + " tiles, twins " +
               std::to_string(idx->stats_.twin_keys) + " keyed segments, parallels " +
               std::to_string(idx->stats_.parallel_keys) + " keyed segments");
    } else {
      LOG_WARN("roundtrip road identity sidecar: " + persisted_path + " refused (" + why +
               "); scanning the tiles instead");
    }
  }
  if (!idx)
    idx = Scan(reader, twin_radius_m, parallel_radius_m, build_parallels, switchback_test);
  auto* raw = idx.get();
  registry.emplace(key, std::move(idx));
  return *raw;
}

void RoadTwinIndex::Build(GraphReader& reader,
                          double twin_radius_m,
                          double parallel_radius_m,
                          bool build_parallels,
                          bool switchback_test) {
  const auto t0 = std::chrono::steady_clock::now();
  twin_radius_m_ = twin_radius_m;
  parallel_radius_m_ = build_parallels ? parallel_radius_m : 0.0;
  switchback_test_ = switchback_test;
  stats_.twin_radius_m = twin_radius_m;
  stats_.parallel_radius_m = build_parallels ? parallel_radius_m : 0.0;
  const double max_radius = build_parallels ? std::max(twin_radius_m, parallel_radius_m)
                                            : twin_radius_m;

  // ---- pass 1: sample every canonical (forward) drivable, non-shortcut edge --------
  std::vector<uint64_t> edge_ids;
  std::vector<uint32_t> samp_off;
  std::vector<uint16_t> samp_cnt;
  std::vector<float> bearing;
  std::vector<float> len_m; // edge length (metres) — the span test's yardstick
  std::vector<uint64_t> way_id; // proto/v4-p1.1: OSM way id — the switchback test
  std::vector<Pt> samples;
  std::vector<std::pair<uint64_t, uint32_t>> cell_entries; // (cell key, edge index)

  const int64_t cell_lat_e6 =
      static_cast<int64_t>(kCellM / kMetersPerDegreeLat * 1e6);
  const int64_t cell_lon_e6 =
      static_cast<int64_t>(kCellM / (kMetersPerDegreeLat * kGridCosLat) * 1e6);
  auto cell_key = [&](const Pt& p) {
    const int64_t ci = static_cast<int64_t>(std::floor(double(p.lat_e6) / cell_lat_e6));
    const int64_t cj = static_cast<int64_t>(std::floor(double(p.lon_e6) / cell_lon_e6));
    return (static_cast<uint64_t>(ci + (1 << 22)) << 32) |
           (static_cast<uint64_t>(cj + (1 << 22)) & 0xffffffffull);
  };
  auto cell_key_ij = [&](int64_t ci, int64_t cj) {
    return (static_cast<uint64_t>(ci + (1 << 22)) << 32) |
           (static_cast<uint64_t>(cj + (1 << 22)) & 0xffffffffull);
  };

  const auto tiles = reader.GetTileSet();
  stats_.tiles = tiles.size();
  std::vector<uint64_t> my_cells;
  for (const auto& tid : tiles) {
    graph_tile_ptr tile = reader.GetGraphTile(tid);
    if (!tile)
      continue;
    const uint32_t n = tile->header()->directededgecount();
    for (uint32_t i = 0; i < n; ++i) {
      const GraphId eid(tid.tileid(), tid.level(), i);
      const DirectedEdge* de = tile->directededge(eid);
      if (de->is_shortcut() || de->IsTransitLine() || de->use() == Use::kTransitConnection)
        continue;
      if (!de->forward()) // canonical = the edge carrying the stored shape direction
        continue;
      if (!((de->forwardaccess() | de->reverseaccess()) & kAutoAccess))
        continue;
      const auto ei = tile->edgeinfo(de);
      const auto shape = ei.shape();
      if (shape.size() < 2)
        continue;
      // proto/v4-p1.1: the canonical key is now min(edge, opposing edge) — see the
      // header.  One opposing lookup per segment at BUILD time buys the F01 forest pass
      // a tile-free key at request time.
      graph_tile_ptr opp_tile = tile;
      const GraphId opp = reader.GetOpposingEdgeId(eid, opp_tile);
      const uint64_t canon = opp.is_valid() ? std::min(eid.value, opp.value) : eid.value;

      // Sample at kSampleSpacingM, endpoints always in, capped.
      const double len = static_cast<double>(de->length());
      double pitch = std::max(kSampleSpacingM, len / (kMaxSamplesPerEdge - 1));
      const uint32_t idx = static_cast<uint32_t>(edge_ids.size());
      edge_ids.push_back(canon);
      way_id.push_back(ei.wayid());
      samp_off.push_back(static_cast<uint32_t>(samples.size()));
      auto push = [&](const PointLL& p) {
        samples.push_back({static_cast<int32_t>(std::llround(p.lat() * 1e6)),
                           static_cast<int32_t>(std::llround(p.lng() * 1e6))});
      };
      push(shape.front());
      double since = 0.0;
      for (size_t s = 1; s < shape.size(); ++s) {
        const double d = shape[s - 1].Distance(shape[s]);
        since += d;
        if (since >= pitch && s + 1 < shape.size()) {
          push(shape[s]);
          since = 0.0;
        }
      }
      push(shape.back());
      samp_cnt.push_back(static_cast<uint16_t>(samples.size() - samp_off.back()));
      bearing.push_back(static_cast<float>(shape.front().Heading(shape.back())));
      len_m.push_back(static_cast<float>(len));

      my_cells.clear();
      for (uint32_t s = samp_off.back(); s < samples.size(); ++s)
        my_cells.push_back(cell_key(samples[s]));
      std::sort(my_cells.begin(), my_cells.end());
      my_cells.erase(std::unique(my_cells.begin(), my_cells.end()), my_cells.end());
      for (uint64_t c : my_cells)
        cell_entries.emplace_back(c, idx);
    }
    // Prototype: the tileset is small enough to keep resident; no eviction dance.
  }
  const uint32_t nedges = static_cast<uint32_t>(edge_ids.size());
  samp_off.push_back(static_cast<uint32_t>(samples.size()));
  stats_.edges_indexed = nedges;
  stats_.samples = samples.size();

  // ---- CSR over the bin grid -------------------------------------------------------
  std::sort(cell_entries.begin(), cell_entries.end());
  std::vector<uint64_t> cell_keys;
  std::vector<uint32_t> cell_off, cell_edges;
  cell_edges.reserve(cell_entries.size());
  for (size_t i = 0; i < cell_entries.size();) {
    const uint64_t k = cell_entries[i].first;
    cell_keys.push_back(k);
    cell_off.push_back(static_cast<uint32_t>(cell_edges.size()));
    while (i < cell_entries.size() && cell_entries[i].first == k)
      cell_edges.push_back(cell_entries[i++].second);
  }
  cell_off.push_back(static_cast<uint32_t>(cell_edges.size()));
  stats_.peak_bytes = cell_entries.size() * sizeof(cell_entries[0]) +
                      samples.size() * sizeof(Pt) + cell_edges.size() * 4 +
                      cell_keys.size() * 12 + nedges * 18;
  std::vector<std::pair<uint64_t, uint32_t>>().swap(cell_entries);

  // ---- pass 2: pair every edge against its 1-ring bin neighbourhood ----------------
  const double twin_r2 = twin_radius_m * twin_radius_m;
  const double max_r2 = max_radius * max_radius;
  std::vector<uint32_t> stamp(nedges, 0);
  std::vector<uint32_t> cand;
  std::vector<std::pair<uint64_t, uint64_t>> twin_pairs, par_pairs;

  // proto/v4-p1.1: where the switchback test bites, on a 0.5 deg grid — the terrain
  // breakdown the ticket asks for (block C mountain origins vs the Belgrade cells),
  // reported without hard-coding any corpus geography into the engine.
  std::map<std::pair<int, int>, uint32_t> drop_cells;
  auto note_drop = [&](const Pt& p) {
    ++drop_cells[{static_cast<int>(std::floor(p.lat_e6 * 1e-6 * 2.0)),
                  static_cast<int>(std::floor(p.lon_e6 * 1e-6 * 2.0))}];
  };

  for (uint32_t u = 0; u < nedges; ++u) {
    const uint32_t uo = samp_off[u], uc = samp_cnt[u];
    // candidate gather
    cand.clear();
    my_cells.clear();
    for (uint32_t s = uo; s < uo + uc; ++s) {
      const int64_t ci = static_cast<int64_t>(std::floor(double(samples[s].lat_e6) / cell_lat_e6));
      const int64_t cj = static_cast<int64_t>(std::floor(double(samples[s].lon_e6) / cell_lon_e6));
      for (int64_t di = -1; di <= 1; ++di)
        for (int64_t dj = -1; dj <= 1; ++dj)
          my_cells.push_back(cell_key_ij(ci + di, cj + dj));
    }
    std::sort(my_cells.begin(), my_cells.end());
    my_cells.erase(std::unique(my_cells.begin(), my_cells.end()), my_cells.end());
    for (uint64_t k : my_cells) {
      auto it = std::lower_bound(cell_keys.begin(), cell_keys.end(), k);
      if (it == cell_keys.end() || *it != k)
        continue;
      const size_t ci = static_cast<size_t>(it - cell_keys.begin());
      for (uint32_t j = cell_off[ci]; j < cell_off[ci + 1]; ++j) {
        const uint32_t v = cell_edges[j];
        if (v <= u || stamp[v] == u + 1)
          continue;
        stamp[v] = u + 1;
        cand.push_back(v);
      }
    }
    if (cand.empty())
      continue;

    // local flat frame at u's first sample
    const double lat0 = samples[uo].lat_e6 * 1e-6;
    const double m_lat = kMetersPerDegreeLat * 1e-6;
    const double m_lon = kMetersPerDegreeLat * std::cos(lat0 * kRadPerDeg) * 1e-6;
    auto X = [&](const Pt& p) { return p.lon_e6 * m_lon; };
    auto Y = [&](const Pt& p) { return p.lat_e6 * m_lat; };

    for (uint32_t v : cand) {
      if (collinear_angle(bearing[u], bearing[v]) > kBearingTolDeg)
        continue;
      ++stats_.pair_tests;
      const uint32_t vo = samp_off[v], vc = samp_cnt[v];

      // A pair qualifies as "the same physical road, two carriageways" on TWO tests,
      // both needed:
      //   COVER — the share of a's samples whose distance to b's polyline is <= r;
      //   SPAN  — the along-b extent those samples project onto, as a share of the
      //           shorter polyline.
      // The span test is what separates a parallel run from a CONTINUATION.  Two short
      // consecutive edges (or two edges crossing at a junction) are collinear and every
      // sample sits within 30 m of the shared node, so cover alone reads 100 % — but
      // every sample projects onto the same point of the neighbour, so the span is ~0.
      // Without it the Serbia index called 40 % of all segments twins.
      auto run = [&](uint32_t ao, uint32_t ac, uint32_t bo, uint32_t bc, double r2,
                     double& cov, double& span, double& off_min, double& off_max) {
        double blen[kMaxSamplesPerEdge + 2];
        double acc = 0.0;
        blen[0] = 0.0;
        for (uint32_t t = 0; t + 1 < bc; ++t) {
          const double dx = X(samples[bo + t + 1]) - X(samples[bo + t]);
          const double dy = Y(samples[bo + t + 1]) - Y(samples[bo + t]);
          acc += std::sqrt(dx * dx + dy * dy);
          blen[t + 1] = acc;
        }
        uint32_t n = 0;
        double smin = std::numeric_limits<double>::max(), smax = -1.0;
        off_min = std::numeric_limits<double>::max();
        off_max = -1.0;
        for (uint32_t s = ao; s < ao + ac; ++s) {
          const double px = X(samples[s]), py = Y(samples[s]);
          double best = std::numeric_limits<double>::max(), best_s = 0.0, tt = 0.0;
          for (uint32_t t = 0; t + 1 < bc; ++t) {
            const double d2 = seg_dist2(px, py, X(samples[bo + t]), Y(samples[bo + t]),
                                        X(samples[bo + t + 1]), Y(samples[bo + t + 1]), tt);
            if (d2 < best) {
              best = d2;
              best_s = blen[t] + tt * (blen[t + 1] - blen[t]);
            }
          }
          if (best <= r2) {
            ++n;
            smin = std::min(smin, best_s);
            smax = std::max(smax, best_s);
            // proto/v4-p1.1: the covered samples' lateral offsets — the switchback test.
            const double off = std::sqrt(best);
            off_min = std::min(off_min, off);
            off_max = std::max(off_max, off);
          }
        }
        cov = ac ? double(n) / ac : 0.0;
        span = (smax >= smin) ? (smax - smin) : 0.0;
      };
      // proto/v4-p1.1: `constant` is the switchback verdict for whichever direction
      // qualified; `same_way` is the cheap OSM-identity half of it.
      const bool same_way = way_id[u] == way_id[v];
      auto qualifies = [&](double r2, bool& constant) {
        double c1 = 0, s1 = 0, c2 = 0, s2 = 0;
        double o1lo = 0, o1hi = 0, o2lo = 0, o2hi = 0;
        run(uo, uc, vo, vc, r2, c1, s1, o1lo, o1hi);
        run(vo, vc, uo, uc, r2, c2, s2, o2lo, o2hi);
        const double shorter = std::min(len_m[u], len_m[v]);
        const bool q1 = c1 >= kCoverFraction && s1 >= kSpanFraction * shorter;
        const bool q2 = c2 >= kCoverFraction && s2 >= kSpanFraction * shorter;
        auto flat = [](double lo, double hi) {
          return hi >= 0.0 && ((hi - lo) <= kOffsetSpreadM || lo >= kOffsetConstFrac * hi);
        };
        constant = (q1 && flat(o1lo, o1hi)) || (q2 && flat(o2lo, o2hi));
        return q1 || q2;
      };
      // A pair is the same physical road only if the geometry says so AND the
      // switchback test does not veto it: different OSM way, or a constant offset.
      bool constant = false;
      if (qualifies(twin_r2, constant)) {
        ++stats_.geom_twin_pairs;
        if (!switchback_test || !same_way || constant) {
          twin_pairs.emplace_back(edge_ids[u], edge_ids[v]);
        } else {
          ++stats_.switchback_dropped;
          ++stats_.sb_same_way;
          if (!constant)
            ++stats_.sb_divergent;
          note_drop(samples[uo]);
        }
      } else if (build_parallels && qualifies(max_r2, constant)) {
        ++stats_.geom_par_pairs;
        if (!switchback_test || !same_way || constant)
          par_pairs.emplace_back(edge_ids[u], edge_ids[v]);
        else
          ++stats_.switchback_dropped_par;
      }
    }
  }
  stats_.twin_pairs = twin_pairs.size();
  stats_.parallel_pairs = par_pairs.size();

  // ---- emit CSR maps (both directions of each pair) ---------------------------------
  auto emit = [](std::vector<std::pair<uint64_t, uint64_t>>& pairs, std::vector<uint64_t>& keys,
                 std::vector<uint32_t>& off, std::vector<uint64_t>& vals) {
    std::vector<std::pair<uint64_t, uint64_t>> both;
    both.reserve(pairs.size() * 2);
    for (const auto& p : pairs) {
      both.emplace_back(p.first, p.second);
      both.emplace_back(p.second, p.first);
    }
    std::vector<std::pair<uint64_t, uint64_t>>().swap(pairs);
    std::sort(both.begin(), both.end());
    both.erase(std::unique(both.begin(), both.end()), both.end());
    vals.reserve(both.size());
    for (size_t i = 0; i < both.size();) {
      const uint64_t k = both[i].first;
      keys.push_back(k);
      off.push_back(static_cast<uint32_t>(vals.size()));
      while (i < both.size() && both[i].first == k)
        vals.push_back(both[i++].second);
    }
    off.push_back(static_cast<uint32_t>(vals.size()));
  };
  emit(twin_pairs, own_twin_keys_, own_twin_off_, own_twin_vals_);
  emit(par_pairs, own_par_keys_, own_par_off_, own_par_vals_);
  // proto/v4-p1.1: the Bloom prefilter over the twin key set (see the header).
  own_filter_.assign(kFilterBits / 64, 0ull);
  for (uint64_t k : own_twin_keys_) {
    const uint64_t h = mix(k) & kFilterMask;
    own_filter_[h >> 6] |= 1ull << (h & 63);
  }
  point_views_at_vectors();
  stats_.twin_keys = twin_keys_.n;
  stats_.parallel_keys = par_keys_.n;
  stats_.retained_bytes = twin_keys_.n * 8 + twin_off_.n * 4 + twin_vals_.n * 8 + par_keys_.n * 8 +
                          par_off_.n * 4 + par_vals_.n * 8;
  stats_.build_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();

  LOG_INFO("roundtrip road identity sidecar: built in " +
           std::to_string(static_cast<int>(stats_.build_ms)) + " ms — " +
           std::to_string(stats_.tiles) + " tiles, " + std::to_string(stats_.edges_indexed) +
           " segments, " + std::to_string(stats_.samples) + " shape samples, " +
           std::to_string(stats_.pair_tests) + " pair tests; twins " +
           std::to_string(stats_.twin_pairs) + " pairs / " + std::to_string(stats_.twin_keys) +
           " keyed segments (r<=" + std::to_string(static_cast<int>(twin_radius_m)) +
           " m), parallels " + std::to_string(stats_.parallel_pairs) + " pairs / " +
           std::to_string(stats_.parallel_keys) + " keyed segments (r<=" +
           std::to_string(static_cast<int>(stats_.parallel_radius_m)) + " m); retained " +
           std::to_string(stats_.retained_bytes / 1024 / 1024) + " MiB, peak ~" +
           std::to_string(stats_.peak_bytes / 1024 / 1024) + " MiB");

  // proto/v4-p1.1 switchback ledger: what the OSM-way / constant-offset test removed,
  // and the 0.5 deg cells it removed it from (top 12, descending).
  if (switchback_test) {
    std::vector<std::pair<uint32_t, std::pair<int, int>>> top;
    for (const auto& kv : drop_cells)
      top.emplace_back(kv.second, kv.first);
    std::sort(top.begin(), top.end(), std::greater<>());
    std::string cells;
    for (size_t i = 0; i < top.size() && i < 12; ++i)
      cells += " (" + std::to_string(top[i].second.first / 2.0).substr(0, 4) + "," +
               std::to_string(top[i].second.second / 2.0).substr(0, 5) + ")=" +
               std::to_string(top[i].first);
    LOG_INFO("roundtrip switchback test: dropped " + std::to_string(stats_.switchback_dropped) +
             " of " + std::to_string(stats_.geom_twin_pairs) + " geometric twin pairs (" +
             std::to_string(stats_.geom_twin_pairs
                                ? 100 * stats_.switchback_dropped / stats_.geom_twin_pairs
                                : 0) +
             " %) — same OSM way " + std::to_string(stats_.sb_same_way) +
             ", divergent offset " + std::to_string(stats_.sb_divergent) + "; parallels dropped " +
             std::to_string(stats_.switchback_dropped_par) + " of " +
             std::to_string(stats_.geom_par_pairs) + "; drop cells (0.5 deg):" + cells);
  }
}

} // namespace thor
} // namespace valhalla
