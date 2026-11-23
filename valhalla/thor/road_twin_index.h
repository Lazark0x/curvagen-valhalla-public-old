#ifndef VALHALLA_THOR_ROAD_TWIN_INDEX_H_
#define VALHALLA_THOR_ROAD_TWIN_INDEX_H_

// The round-trip v4 road-identity sidecar (ADR-0042), served by default.
//
// The census (docs/curvagen/research/2026-09-06-defect-atlas-v2.md §3) named one
// mechanism as the whole rider-visible residual: F02, "the return rides the other
// carriageway".  v3 keys exclusion, leash, rejoin and xcand on DIRECTED EDGE IDS, and
// the opposite carriageway of a divided road is a different OSM way with different
// nodes — so nothing in the engine knows the two are the same physical road.
//
// This sidecar is the missing notion of road identity, computed ONCE per process from
// the loaded tiles:
//   * TWIN     — another road segment whose shape runs within `twin_radius_m` (30 m)
//                of this one over most of the shorter of the two, and whose chord is
//                collinear (angle within tolerance of 0 deg or 180 deg).  That is the
//                opposite carriageway / one-way pair.
//   * PARALLEL — the same test at `parallel_radius_m` (80 m) but not at the twin
//                radius: the street one block over.
//
// Orientation note.  The ticket phrases the twin test as "bearing difference ~180 deg".
// Physically that is right — the return rides the twin ANTI-parallel to the forward
// leg — but OSM's storage direction for a carriageway is arbitrary, so at index time
// the same physical pair can present as 0 deg or 180 deg.  The index therefore
// qualifies on COLLINEARITY and stores segment identity; the call sites exclude both
// directions of a twin (exactly as route_leg already does for the corridor edge and
// its opposing edge), which is what makes the anti-parallel ride impossible.
//
// Persisted form.  The index depends only on the tiles and the four parameters, so a
// region pack can carry it: valhalla_build_road_identity writes it next to the pack's
// tar, and thor.roundtrip_road_identity_file names it.  get() then maps the file
// read-only instead of scanning every tile; its pages are file-backed and clean, so
// they do not count against the process footprint.  A file whose format, parameters,
// build rules or tileset fingerprint differ from what this process would build is
// refused with a warning and the index is scanned as before.
//
// Keys and values are CANONICAL directed-edge ids: the directed edge with
// DirectedEdge::forward() == true, i.e. the one whose travel direction matches the
// stored EdgeInfo shape.  Every undirected road segment has exactly one.  Call sites
// canonicalise a query id with canonical_id() below (they already hold the tile and
// the opposing id, so it is free).

#include "baldr/graphreader.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace valhalla {
namespace thor {

class RoadTwinIndex {
public:
  struct Stats {
    uint64_t tiles = 0;
    uint64_t edges_indexed = 0; // canonical (forward) directed edges sampled
    uint64_t samples = 0;
    uint64_t pair_tests = 0;
    uint64_t twin_pairs = 0;
    uint64_t parallel_pairs = 0;
    // proto/v4-p1.1 switchback test (curvagen-valhalla#12 item 1)
    uint64_t geom_twin_pairs = 0;      // pairs the P1 geometry alone called twins
    uint64_t switchback_dropped = 0;   // ...of which the switchback test rejected
    uint64_t sb_same_way = 0;          // dropped: same OSM way id
    uint64_t sb_divergent = 0;         // dropped: offset not constant along the run
    uint64_t geom_par_pairs = 0;
    uint64_t switchback_dropped_par = 0;
    uint64_t twin_keys = 0;
    uint64_t parallel_keys = 0;
    double build_ms = 0.0;
    uint64_t retained_bytes = 0;
    uint64_t peak_bytes = 0;
    double twin_radius_m = 0.0;
    double parallel_radius_m = 0.0;
  };

  // Process-wide, built once per (tile location, radii) triple.  Keyed rather than a
  // plain singleton because gurka builds many tilesets in one process.  Thread-safe.
  // A non-empty `persisted_path` is tried first (see Load); on refusal the index is
  // scanned from the tiles.
  static const RoadTwinIndex& get(baldr::GraphReader& reader,
                                  double twin_radius_m,
                                  double parallel_radius_m,
                                  bool build_parallels,
                                  bool switchback_test,
                                  const std::string& persisted_path = "");

  // An index scanned from every tile `reader` holds, outside the process registry.
  static std::unique_ptr<RoadTwinIndex> Scan(baldr::GraphReader& reader,
                                             double twin_radius_m,
                                             double parallel_radius_m,
                                             bool build_parallels,
                                             bool switchback_test);

  // Maps a file written by Save.  nullptr, with the reason in `why`, when the file is
  // missing, malformed, or was built from other tiles, other parameters or other rules.
  static std::unique_ptr<RoadTwinIndex> Load(const std::string& path,
                                             baldr::GraphReader& reader,
                                             double twin_radius_m,
                                             double parallel_radius_m,
                                             bool build_parallels,
                                             bool switchback_test,
                                             std::string* why = nullptr);

  // Writes this index for the tileset `reader` holds; false on an IO error.
  bool Save(const std::string& path, baldr::GraphReader& reader) const;

  // Identical twin and parallel maps and prefilter.
  bool same_contents(const RoadTwinIndex& other) const;

  // A fingerprint of every tile's id, size and data hash: what a persisted index is
  // checked against before it is used.
  static uint64_t tileset_fingerprint(baldr::GraphReader& reader);

  // The canonical id of a directed edge = min(edge, opposing edge).
  //
  // proto/v4-p1.1: this WAS `de->forward() ? eid : opp`.  min() is the same
  // undirected identity but needs no DirectedEdge dereference, which is what lets the
  // F01 forest pass (roundtrip_expansion.cc) key every settled label without a tile
  // lookup — see that file's canonical_from_label().  `de` is kept in the signature
  // (unused) so the call sites read the same.
  static uint64_t canonical_id(const baldr::DirectedEdge*,
                               const baldr::GraphId& eid,
                               const baldr::GraphId& opp) {
    return opp.is_valid() ? std::min(eid.value, opp.value) : eid.value;
  }
  static uint64_t canonical_id(uint64_t eid, uint64_t opp) {
    return std::min(eid, opp);
  }

  void append_twins(uint64_t canonical, std::vector<uint64_t>& out) const {
    append(twin_keys_, twin_off_, twin_vals_, canonical, out);
  }
  void append_parallels(uint64_t canonical, std::vector<uint64_t>& out) const {
    append(par_keys_, par_off_, par_vals_, canonical, out);
  }
  bool has_twins(uint64_t canonical) const {
    return lookup(twin_keys_, canonical) != twin_keys_.n;
  }

  // proto/v4-p1.1: one-probe Bloom prefilter over the twin key set.  A clear bit is a
  // PROOF of "no twins"; a set bit falls through to the CSR binary search.  The F01
  // forest pass asks this question once per settled label (hundreds of thousands per
  // request) and >98 % of the answers are "no" — a 512 KiB bit array turns those into
  // one cache probe instead of a ~17-step lower_bound over 115 k keys.
  bool maybe_has_twins(uint64_t canonical) const {
    if (twin_filter_.n == 0)
      return true;
    const uint64_t h = mix(canonical) & kFilterMask;
    return (twin_filter_.p[h >> 6] >> (h & 63)) & 1ull;
  }

  const Stats& stats() const {
    return stats_;
  }
  // True when the arrays live in a mapped file rather than in this process's heap.
  bool mapped() const {
    return mapping_ != nullptr;
  }

private:
  // A read-only array: into the vectors below when scanned, into the mapping when loaded.
  template <typename T> struct View {
    const T* p = nullptr;
    size_t n = 0;
    const T* begin() const {
      return p;
    }
    const T* end() const {
      return p + n;
    }
  };
  struct Mapping;

  RoadTwinIndex() = default;
  RoadTwinIndex(const RoadTwinIndex&) = delete; // the views point into this object
  RoadTwinIndex& operator=(const RoadTwinIndex&) = delete;
  void Build(baldr::GraphReader& reader,
             double twin_radius_m,
             double parallel_radius_m,
             bool build_parallels,
             bool switchback_test);
  void point_views_at_vectors();

  static size_t lookup(const View<uint64_t>& keys, uint64_t k) {
    auto it = std::lower_bound(keys.begin(), keys.end(), k);
    return (it == keys.end() || *it != k) ? keys.n : static_cast<size_t>(it - keys.begin());
  }
  static void append(const View<uint64_t>& keys,
                     const View<uint32_t>& off,
                     const View<uint64_t>& vals,
                     uint64_t k,
                     std::vector<uint64_t>& out) {
    const size_t i = lookup(keys, k);
    if (i == keys.n)
      return;
    out.insert(out.end(), vals.p + off.p[i], vals.p + off.p[i + 1]);
  }

  static constexpr uint64_t kFilterBits = 1ull << 22; // 4 Mi bits = 512 KiB
  static constexpr uint64_t kFilterMask = kFilterBits - 1;
  static uint64_t mix(uint64_t x) { // splitmix64 finalizer
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ull;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebull;
    return x ^ (x >> 31);
  }

  // CSR: sorted keys, [off[i], off[i+1]) into vals.
  View<uint64_t> twin_filter_, twin_keys_, twin_vals_, par_keys_, par_vals_;
  View<uint32_t> twin_off_, par_off_;

  // Storage behind the views: the scanned vectors, or the mapped file.
  std::vector<uint64_t> own_filter_, own_twin_keys_, own_twin_vals_, own_par_keys_, own_par_vals_;
  std::vector<uint32_t> own_twin_off_, own_par_off_;
  std::shared_ptr<const Mapping> mapping_;

  double twin_radius_m_ = 0.0;
  double parallel_radius_m_ = 0.0; // 0 when parallels are not built
  bool switchback_test_ = false;
  Stats stats_;
};

} // namespace thor
} // namespace valhalla

#endif // VALHALLA_THOR_ROAD_TWIN_INDEX_H_
